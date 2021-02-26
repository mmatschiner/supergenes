# m_matschiner Tue Mar 10 10:51:36 CET 2020

import re

import newick

def parse_species_tree(tree, Ne, branch_length_units="gen", generation_time=None):
    """
    This function parses a species tree in
    `Newick <https://en.wikipedia.org/wiki/Newick_format>`_ format and defines a
    simulation model according to the species tree. The tree is assumed to be
    rooted and ultrametric and branch lengths must be included and correspond to
    time, either in units of millions of years, years, or generations. Leaves must
    be named.

    After reading the input tree, this function defines a
    :class:`.PopulationConfiguration` instance for each terminal node in the tree,
    corresponding to extant species. These population configurations store the
    species' name and population size. The specified Ne is used as the size of all
    populations. Additionally, one or more :class:`.MassMigration` instances are
    defined for each internal node, with the time of the mass migration set
    according to the age of the node in the species tree. :math:`n - 1` mass
    migration events are defined for internal nodes with `n` descendants, meaning
    that a single event is defined for bifurcating nodes. For each internal node,
    the left-most of the descendant populations is arbitrarily selected as the
    destination in all mass migrations defined for that node.

    :param str tree: The tree string in Newick format, with named leaves and branch
        lengths.
    :param float Ne: The effective population size.
    :param str branch_length_units: The units of time in which the species tree's
        branch lengths are measured. Allowed branch length units are millions of
        years, years, and generations; these should be specified with the strings
        ``"myr"``, ``"yr"``, or ``"gen"``, respectively. This defaults to
        ``"gen"``.
    :param float generation_time: The number of years per generation. If and only
        if the branch lengths are not in units of generations, the generation time
        must be specified. This defaults to `None`.
    :type generation_time: float or None
    :return: A tuple of two lists of which the first contains
        :class:`.PopulationConfiguration` instances and the second contains
        :class:`.MassMigration` instances. The population configurations specify
        the size of each population and the species name corresponding to each
        population. Species names are stored as metadata in each
        :class:`.PopulationConfiguration` instance, with the metadata tag
        "species_name". Sampling configurations and growth rates are not specified
        in the population configurations. The list of population configurations is
        ordered according to the order of the corresponding extant species in a
        `post-order tree traversal
        <https://en.wikipedia.org/wiki/Tree_traversal#Post-order_(LRN)>`_. The list
        of mass migration events is ordered by the time of the events, from young
        to old events.
    :rtype: (list, list)
    :warning: This function does not modify migration matrices. When the population
        configurations and mass migration events returned by this function are used
        to simulate with the :func:`.simulate` function, it should be ensured that
        migration rates to source populations of mass migration events are zero
        after the mass migration (viewed backwards in time).
    """

    # Make sure that branch length units are either "myr", "yr", or "gen".
    allowed_branch_lenth_units = ["myr", "yr", "gen"]
    if branch_length_units not in allowed_branch_lenth_units:
        err = "The specified units for branch lengths ("
        err += f'"{branch_length_units}") are not accepted. '
        err += 'Accepted units are "myr" (millions of years), "yr" (years), '
        err += 'and "gen" (generations).'
        raise ValueError(err)

    try:
        Ne = float(Ne)
    except ValueError:
        raise ValueError("Population size Ne must be numeric.")
    if Ne <= 0:
        raise ValueError("Population size Ne must be > 0.")

    # Make sure that the generation time is either None or positive.
    if generation_time is not None:
        generation_time = check_generation_time(generation_time)

    # Make sure that the generation time is specified if and only if
    # branch lengths are not in units of generations.
    if branch_length_units == "gen":
        if generation_time is not None:
            err = 'With branch lengths in units of generations ("gen"), '
            err += "a generation time should not be specified additionally."
            raise ValueError(err)
    else:
        if generation_time is None:
            err = "With branch lengths in units of "
            err += f'"{branch_length_units}", a generation time must be '
            err += "specified additionally."
            raise ValueError(err)

    # Get the number of generations per branch length unit.
    generations_per_branch_length_unit = get_generations_per_branch_length_unit(
        branch_length_units, generation_time
    )

    # Define populations and demographic events according to the
    # specified population size and the divergence times in the species tree.
    # Per divergence event (node in the tree), a mass migration with a proportion
    # of 1 of the population is used. The destination is the left-most leaf for
    # each node. Because we are using the n leaf populations, we map each node back
    # to the leaf population that it corresponds to.
    populations = []
    demographic_events = []
    leaf_map = {}
    for node in root.walk("postorder"):
        if len(node.descendants) == 0:
            # Per extant species (= leaf node) in the tree, add a population with
            # size Ne. Species names are stored as metadata with the "species_name"
            # tag.
            populations.append(
                msprime.PopulationConfiguration(initial_size=Ne, metadata=node.name.strip())
            )
            leaf_map[node] = len(populations) - 1
        else:
            # Per internal node, add one (if the node is bifurcating) or multiple
            # (if the node is multi-furcating) MassMigrations. The parent species
            # maps implicitly to the left-most child species, so we don't generate
            # any MassMigrations for that.
            leaf_map[node] = leaf_map[node.descendants[0]]
            # For each child species after the left-most one, we create
            # a MassMigration into the left-most species.
            for child in node.descendants[1:]:
                demographic_events.append(
                    msprime.MassMigration(
                        source=leaf_map[child], dest=leaf_map[node], time=node.time
                    )
                )

    # Sort demographic events by time.
    demographic_events.sort(key=lambda de: de.time)

    return populations, demographic_events

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False


def check_generation_time(generation_time):
    try:
        generation_time = float(generation_time)
    except ValueError:
        raise ValueError("Generation time must be numeric.")
    if generation_time <= 0:
        raise ValueError("Generation time must be > 0.")
    return generation_time


def get_generations_per_branch_length_unit(branch_length_units, generation_time):
    """
    Method to calculate the number of generations per branch length
    unit, given the branch length unit and a generation time.
    """
    if branch_length_units == "gen":
        generations_per_branch_length_unit = 1
    elif branch_length_units == "myr":
        generations_per_branch_length_unit = 10 ** 6 / generation_time
    else:
        generations_per_branch_length_unit = 1 / generation_time
    return generations_per_branch_length_unit


def parse_newick(tree, branch_length_multiplier):
    """
    Parses the newick tree and annotates the resulting nodes with their
    time values, appropriately scaled.
    """
    # Parse the newick tree string.
    parsed = newick.loads(tree)
    if len(parsed) == 0:
        raise ValueError(f"Not a valid newick tree: '{tree}'")
    root = parsed[0]

    # Set node depths (distances from root).
    stack = [(root, 0)]
    num_nodes = 0
    max_depth = 0
    while len(stack) > 0:
        node, depth = stack.pop()
        if depth > max_depth:
            max_depth = depth
        num_nodes += 1
        node.depth = depth
        for child in node.descendants:
            stack.append((child, depth + child.length))
    if num_nodes < 3:
        raise ValueError("Newick tree must have at least three nodes")

    # Set node times (distances from present).
    for node in root.walk():
        node.time = (max_depth - node.depth) * branch_length_multiplier

    return root


# Import libraries.
import msprime
import sys
import newick
import random
from random import randint
import datetime
import re
import os

# Get the command line arguments.
tree_file_name = sys.argv[1]
rec_map_file_name = sys.argv[2]
mut_rate = float(sys.argv[3])
generation_time = float(sys.argv[4])
pop_size = int(sys.argv[5])
vcf_outfile_name = sys.argv[6]
chr_length = int(sys.argv[7])

# Read the recombination map.
recombination_map = msprime.RecombinationMap.read_hapmap(rec_map_file_name)

# Set simulation parameters.
bottleneck_pop_size = 0.5
bottleneck_end_time = 880800
bottleneck_start_time = 880810
sample_size = 4

# Read the species tree string.
with open(tree_file_name, "r") as f:
    species_tree_string = f.read()

# Parse the species tree.
root = parse_newick(species_tree_string, 10**6 / generation_time)

# Parse the species tree with msprime and generate population configurations and demographic evens.
parsed_tuple = parse_species_tree(
    tree=species_tree_string,
    Ne=pop_size,
    branch_length_units="myr",
    generation_time=generation_time
    )
population_configurations = parsed_tuple[0]
demographic_events = parsed_tuple[1]
for n in population_configurations:
    n.sample_size = sample_size

# Remove the last demographic event, and keep to reinsert later.
last_demographic_event = demographic_events[-1]
del demographic_events[-1]

# Add two events for the bottleneck.
demographic_events.append(
    msprime.PopulationParametersChange(
        time=bottleneck_end_time/generation_time,
        initial_size=bottleneck_pop_size,
        growth_rate=0,
        population_id=0),
    )
demographic_events.append(
    msprime.PopulationParametersChange(
        time=bottleneck_start_time/generation_time,
        initial_size=pop_size,
        growth_rate=0,
        population_id=0),
    )

# Reinsert the last demographic event.
demographic_events.append(last_demographic_event)

# Write the vcf file.
print('Simulating with msprime...', end='', flush=True)
new_tree_sequence_obj = msprime.simulate(
        population_configurations=population_configurations,
        Ne=30000,
        demographic_events=demographic_events,
        mutation_rate=mut_rate*generation_time,
        recombination_map=recombination_map,
        random_seed=random.randint(1, 10000000)
        )
print(" done.")

# Prepare the vcf header.
print('Preparing the vcf...', end='', flush=True)
vcf_string = '##fileformat=VCFv4.2\n'
now = datetime.datetime.now()
vcf_string += '##fileDate={}\n'.format(now.strftime("%Y%m%d"))
vcf_string += '##source=simulate_lg02_with_recmap.py\n'
vcf_string += '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
vcf_string += '##contig=<ID=1,length={}>\n'.format(chr_length)
vcf_string += '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT'
for pop_config in population_configurations:
    vcf_string += '\t{}_1\t{}_2'.format(pop_config.metadata, pop_config.metadata)
vcf_string += '\n'

# Prepare vcf file.
DNA_alphabet = ['A','C','G','T']
n_variants = 0
in_gene_position = 0
for variant in new_tree_sequence_obj.variants():
    vp = int(variant.site.position)
    if vp > in_gene_position:  # exclude multiple mutations at the same site (i.e. at the moment we generate only biallelic sites) TO DO: make multiallelics in these cases
        in_gene_position = vp
        n_variants = n_variants + 1
        ra = randint(0, 3)
        ancestralBase = DNA_alphabet[ra]
        derivedPossibilities = DNA_alphabet[:ra] +  DNA_alphabet[ra+1:]
        rd = randint(0, 2)
        derivedBase = derivedPossibilities[rd]
        vcf_string += '1\t{}\t.\t{}\t{}\t.\t.\t.\tGT'.format(vp,ancestralBase,derivedBase)
        # print(variant.genotypes)
        # print()
        # sys.exit(0)
        for sp in range(0,len(population_configurations)):
            vcf_string += '\t{}|{}'.format(variant.genotypes[sample_size*sp],variant.genotypes[(sample_size*sp)+1])
            if sample_size == 4:
                vcf_string += '\t{}|{}'.format(variant.genotypes[(sample_size*sp)+2],variant.genotypes[(sample_size*sp)+3])
        vcf_string += '\n'

# Write the vcf output file.
vcf_outfile = open(vcf_outfile_name, 'w')
vcf_outfile.write(vcf_string)
print(" done.")
