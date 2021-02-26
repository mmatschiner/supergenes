# m_matschiner Wed Jul 8 12:15:44 CEST 2020

# Import libraries.
import msprime
import sys
import random
from random import randint
import datetime
import re
import os

# Get the command line arguments.
recombination_rate = float(sys.argv[1])
mutation_rate = float(sys.argv[2])
generation_time = float(sys.argv[3])
# bottleneck_start_time = float(sys.argv[4])
# bottleneck_end_time = bottleneck_start_time - 20*generation_time
bottleneck_end_time = float(sys.argv[4])
population_size = int(sys.argv[5]) # The population size is understood by msprime as the number of diploid individuals.
bottleneck_population_size = 0.5
chr_length = int(sys.argv[6])

# Define the populations.
# population_configurations = [
#         msprime.PopulationConfiguration(sample_size=2,
#                 initial_size=population_size),
#         msprime.PopulationConfiguration(sample_size=2,
#                 initial_size=population_size)
# ]
population_configurations = [
        msprime.PopulationConfiguration(sample_size=2,
                initial_size=population_size)
]

# Define the bottleneck.
# demographic_events = [
#         msprime.PopulationParametersChange(
#                 time=bottleneck_end_time/generation_time,
#                 initial_size=bottleneck_population_size,
#                 growth_rate=0,
#                 population_id=0),
#         msprime.PopulationParametersChange(
#                 time=bottleneck_start_time/generation_time,
#                 initial_size=population_size,
#                 growth_rate=0,
#                 population_id=0),
#         msprime.MassMigration(
#                 time=bottleneck_start_time/generation_time, source=1, destination=0, proportion=1.0)
# ]
demographic_events = [
        msprime.PopulationParametersChange(
                time=bottleneck_end_time/generation_time,
                initial_size=bottleneck_population_size,
                growth_rate=0,
                population_id=0)
]

# Simulate.
# new_tree_sequence_obj = msprime.simulate(
#         population_configurations=population_configurations,
#         demographic_events=demographic_events,
#         mutation_rate=mutation_rate*generation_time,
#         recombination_rate=recombination_rate,
#         length=chr_length,
#         random_seed=random.randint(1, 10000000)
#         )
new_tree_sequence_obj = msprime.simulate(
        population_configurations=population_configurations,
        demographic_events=demographic_events,
        recombination_rate=recombination_rate,
        length=chr_length,
        random_seed=random.randint(1, 10000000)
        )

# sum_branch_lengths = 0
# sum_weighted_branch_lengths = 0
# tree_count = 0
# for tree in new_tree_sequence_obj.trees():
#     print(tree.draw(format="unicode"))
#     sys.exit(0)
#     tree_count += 1
#     sum_branch_lengths += tree.branch_length(0)
#     sum_weighted_branch_lengths += (tree.interval[1] - tree.interval[0]) * tree.branch_length(0)
# mean_branch_length = sum_branch_lengths/float(tree_count)
# mean_weighted_branch_length = sum_weighted_branch_lengths/chr_length
# print("{}\t{}\t{}\t{}".format(population_size,tree_count,mean_branch_length,mean_weighted_branch_length))


bp_coalescing_in_bottleneck = 0
sum_branch_lengths = 0
tree_count = 0
for tree in new_tree_sequence_obj.trees():
    tree_count += 1
    if tree.branch_length(0) > bottleneck_end_time/generation_time:
        bp_coalescing_in_bottleneck += tree.interval[1] - tree.interval[0]
    sum_branch_lengths += tree.branch_length(0) * generation_time
genome_proportion_coalescing_in_bottleneck = bp_coalescing_in_bottleneck/float(chr_length)
mean_branch_length = sum_branch_lengths/float(tree_count)
print("{}\t{}\t{}\t{}\t{}".format(bottleneck_end_time,population_size,tree_count,mean_branch_length,genome_proportion_coalescing_in_bottleneck))

# # Calculate genetic diversity in the two populations.
# variable_in_pop1 = 0
# variable_in_pop2 = 0
# for variant in new_tree_sequence_obj.variants():
#     if variant.genotypes[0] != variant.genotypes[1]:
#         variable_in_pop1 += 1
#     if variant.genotypes[2] != variant.genotypes[3]:
#         variable_in_pop2 += 1
# pi1 = variable_in_pop1/chr_length
# pi2 = variable_in_pop2/chr_length
# print("{}\t{}\t{}\t{}\t{}\t{}".format(population_size,bottleneck_start_time,variable_in_pop1,variable_in_pop2,pi1,pi2))
