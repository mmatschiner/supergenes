# m_matschiner Fri Sep 10 14:14:12 CEST 2021
# michaelmatschiner@mac.com

# This script performs sliding window analyses over a vcf, and
# reports the number of SNPs and the proportion of heterozyguous
# SNPs, Weir and Cockerham's F_st, d_xy, the proportion of fixed
# difference d_f and the level of nucleotde variation measured
# as pi.

# Import libraries and make sure we're on python 3.
import sys
if sys.version_info[0] < 3:
    print('Python 3 is needed to run this script!')
    sys.exit(0)
import argparse, textwrap, vcf, math, scipy
from scipy import special
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio import AlignIO

# Expand the class MultipleSeqAlignment in order to add more statistics.
class XMultipleSeqAlignment(MultipleSeqAlignment):

    def get_number_of_variable_sites(self, pops=[]):
        # Get all sequences of this population.
        seqs = []
        for y in range(0,len(self)):
            if pops == []:
                seqs.append(str(self[y].seq))
            else:
                for pop in pops:
                    if pop in self[y].id:
                        seqs.append(str(self[y].seq))
                        break
        number_of_variable_sites = 0
        for x in range(0,self.get_alignment_length()):
            bases = []
            for y in range(0,len(seqs)):
                if seqs[y][x] == 'A':
                    bases.append("A")
                    bases.append("A")
                elif seqs[y][x] == 'M':
                    bases.append("A")
                    bases.append("C")
                elif seqs[y][x] == 'R':
                    bases.append("A")
                    bases.append("G")
                elif seqs[y][x] == 'W':
                    bases.append("A")
                    bases.append("T")
                elif seqs[y][x] == 'C':
                    bases.append("C")
                    bases.append("C")
                elif seqs[y][x] == 'S':
                    bases.append("C")
                    bases.append("G")
                elif seqs[y][x] == 'Y':
                    bases.append("C")
                    bases.append("T")
                elif seqs[y][x] == 'G':
                    bases.append("G")
                    bases.append("G")
                elif seqs[y][x] == 'K':
                    bases.append("G")
                    bases.append("T")
                elif seqs[y][x] == 'T':
                    bases.append("T")
                    bases.append("T")
            if len(set(bases)) > 1:
                number_of_variable_sites += 1
        return number_of_variable_sites

    def get_heterozygosity(self, pops=[]):
        # Get all sequences of this population.
        seqs = []
        for y in range(0,len(self)):
            if pops == []:
                seqs.append(str(self[y].seq))
            else:
                for pop in pops:
                    if pop in self[y].id:
                        seqs.append(str(self[y].seq))
                        break
        number_of_heterozyguous_positions = 0
        number_of_homozyguous_positions = 0
        for x in range(0,self.get_alignment_length()):
            bases = []
            for y in range(0,len(seqs)):
                if seqs[y][x] == 'A':
                    number_of_homozyguous_positions += 1
                elif seqs[y][x] == 'M':
                    number_of_heterozyguous_positions += 1
                elif seqs[y][x] == 'R':
                    number_of_heterozyguous_positions += 1
                elif seqs[y][x] == 'W':
                    number_of_heterozyguous_positions += 1
                elif seqs[y][x] == 'C':
                    number_of_homozyguous_positions += 1
                elif seqs[y][x] == 'S':
                    number_of_heterozyguous_positions += 1
                elif seqs[y][x] == 'Y':
                    number_of_heterozyguous_positions += 1
                elif seqs[y][x] == 'G':
                    number_of_homozyguous_positions += 1
                elif seqs[y][x] == 'K':
                    number_of_heterozyguous_positions += 1
                elif seqs[y][x] == 'T':
                    number_of_homozyguous_positions += 1
        if number_of_heterozyguous_positions + number_of_homozyguous_positions == 0:
            heterozygosity = None
        else:
            heterozygosity = number_of_heterozyguous_positions/(number_of_heterozyguous_positions+number_of_homozyguous_positions)
        return heterozygosity
        

    def get_gc_content(self, pops=[]):
        # Get all sequences of this population.
        seqs = []
        for y in range(0,len(self)):
            if pops == []:
                seqs.append(str(self[y].seq))
            else:
                for pop in pops:
                    if pop in self[y].id:
                        seqs.append(str(self[y].seq))
                        break
        number_of_gc_positions = 0
        number_of_at_positions = 0
        for x in range(0,self.get_alignment_length()):
            bases = []
            for y in range(0,len(seqs)):
                if seqs[y][x] == 'A':
                    number_of_at_positions += 2
                elif seqs[y][x] == 'M':
                    number_of_gc_positions += 1
                    number_of_at_positions += 1
                elif seqs[y][x] == 'R':
                    number_of_gc_positions += 1
                    number_of_at_positions += 1
                elif seqs[y][x] == 'W':
                    number_of_at_positions += 1
                elif seqs[y][x] == 'C':
                    number_of_gc_positions += 2
                elif seqs[y][x] == 'S':
                    number_of_gc_positions += 1
                elif seqs[y][x] == 'Y':
                    number_of_gc_positions += 1
                    number_of_at_positions += 1
                elif seqs[y][x] == 'G':
                    number_of_gc_positions += 2
                elif seqs[y][x] == 'K':
                    number_of_gc_positions += 1
                    number_of_at_positions += 1
                elif seqs[y][x] == 'T':
                    number_of_at_positions += 2
        if (number_of_gc_positions+number_of_at_positions) == 0:
            gc_content = None
        else:
            gc_content = number_of_gc_positions/(number_of_gc_positions+number_of_at_positions)
        return gc_content

    def get_number_of_unique_genotypes(self, pop=None):
        seqs = []
        for record in self:
            if pop == None or pop in record.id:
                seqs.append(str(record.seq))
        return len(set(seqs))

    def get_proportion_of_variable_sites(self, pops=[]):
        return self.get_number_of_variable_sites(pops)/self.get_alignment_length()

    def get_number_of_invariable_sites(self, pops=[]):
        return self.get_alignment_length()-self.get_number_of_variable_sites(pops)

    def get_allele_frequencies(self, x, pops=[]):
        if isinstance(pops, list) == False:
            print("ERROR: Populations are not given as a list (get_allele_frequencies):")
            print(pops)
            sys.exit(1)
        else:
            allele_frequencies = [0, 0, 0, 0]
            for y in range(0,len(self)):
                for pop in pops:
                    if pop in self[y].id:
                        if self[y].seq[x] == 'A':
                            allele_frequencies[0] += 2
                        elif self[y].seq[x] == 'M':
                            allele_frequencies[0] += 1
                            allele_frequencies[1] += 1
                        elif self[y].seq[x] == 'R':
                            allele_frequencies[0] += 1
                            allele_frequencies[2] += 1
                        elif self[y].seq[x] == 'W':
                            allele_frequencies[0] += 1
                            allele_frequencies[3] += 1
                        elif self[y].seq[x] == 'C':
                            allele_frequencies[1] += 2
                        elif self[y].seq[x] == 'S':
                            allele_frequencies[1] += 1
                            allele_frequencies[2] += 1
                        elif self[y].seq[x] == 'Y':
                            allele_frequencies[1] += 1
                            allele_frequencies[3] += 1
                        elif self[y].seq[x] == 'G':
                            allele_frequencies[2] += 2
                        elif self[y].seq[x] == 'K':
                            allele_frequencies[2] += 1
                            allele_frequencies[3] += 1
                        elif self[y].seq[x] == 'T':
                            allele_frequencies[3] += 2
            return allele_frequencies

    def get_is_biallelic_per_site(self, x, pops=[]):
        if isinstance(pops, list) == False:
            print("ERROR: Populations are not given as a list (get_is_biallelic_per_site):")
            print(pops)
            sys.exit(1)
        else:
            allele_frequencies = [0, 0, 0, 0]
            for pop in pops:
                pop_allele_frequencies = self.get_allele_frequencies(x, [pop])
                allele_frequencies[0] += pop_allele_frequencies[0]
                allele_frequencies[1] += pop_allele_frequencies[1]
                allele_frequencies[2] += pop_allele_frequencies[2]
                allele_frequencies[3] += pop_allele_frequencies[3]
            allele_frequency_zeros = 0
            if allele_frequencies[0] == 0:
                allele_frequency_zeros += 1
            if allele_frequencies[1] == 0:
                allele_frequency_zeros += 1
            if allele_frequencies[2] == 0:
                allele_frequency_zeros += 1
            if allele_frequencies[3] == 0:
                allele_frequency_zeros += 1
            if allele_frequency_zeros == 2:
                return True
            else:
                return False

    def get_is_variable_per_site(self, x, pops=[]):
        if isinstance(pops, list) == False:
            print("ERROR: Populations are not given as a list (get_is_variable_per_site):")
            print(pops)
            sys.exit(1)
        else:
            allele_frequencies = [0, 0, 0, 0]
            for pop in pops:
                pop_allele_frequencies = self.get_allele_frequencies(x, [pop])
                allele_frequencies[0] += pop_allele_frequencies[0]
                allele_frequencies[1] += pop_allele_frequencies[1]
                allele_frequencies[2] += pop_allele_frequencies[2]
                allele_frequencies[3] += pop_allele_frequencies[3]
            allele_frequency_zeros = 0
            if allele_frequencies[0] == 0:
                allele_frequency_zeros += 1
            if allele_frequencies[1] == 0:
                allele_frequency_zeros += 1
            if allele_frequencies[2] == 0:
                allele_frequency_zeros += 1
            if allele_frequencies[3] == 0:
                allele_frequency_zeros += 1
            if allele_frequency_zeros == 3:
                return False
            else:
                return True

    def get_alleles_per_site(self, x, pops=[]):
        if isinstance(pops, list) == False:
            print("ERROR: Populations are not given as a list (get_alleles_per_site):")
            print(pops)
            sys.exit(1)
        else:
            alleles = []
            for y in range(0,len(self)):
                if pops == []:
                    if self[y].seq[x] == 'A':
                        alleles.append("A")
                        alleles.append("A")
                    elif self[y].seq[x] == 'M':
                        alleles.append("A")
                        alleles.append("C")
                    elif self[y].seq[x] == 'R':
                        alleles.append("A")
                        alleles.append("G")
                    elif self[y].seq[x] == 'W':
                        alleles.append("A")
                        alleles.append("T")
                    elif self[y].seq[x] == 'C':
                        alleles.append("C")
                        alleles.append("C")
                    elif self[y].seq[x] == 'S':
                        alleles.append("C")
                        alleles.append("G")
                    elif self[y].seq[x] == 'Y':
                        alleles.append("C")
                        alleles.append("T")
                    elif self[y].seq[x] == 'G':
                        alleles.append("G")
                        alleles.append("G")
                    elif self[y].seq[x] == 'K':
                        alleles.append("G")
                        alleles.append("T")
                    elif self[y].seq[x] == 'T':
                        alleles.append("T")
                        alleles.append("T")
                else:
                    for pop in pops:
                        if pop in self[y].id:
                            if self[y].seq[x] == 'A':
                                alleles.append("A")
                                alleles.append("A")
                            elif self[y].seq[x] == 'M':
                                alleles.append("A")
                                alleles.append("C")
                            elif self[y].seq[x] == 'R':
                                alleles.append("A")
                                alleles.append("G")
                            elif self[y].seq[x] == 'W':
                                alleles.append("A")
                                alleles.append("T")
                            elif self[y].seq[x] == 'C':
                                alleles.append("C")
                                alleles.append("C")
                            elif self[y].seq[x] == 'S':
                                alleles.append("C")
                                alleles.append("G")
                            elif self[y].seq[x] == 'Y':
                                alleles.append("C")
                                alleles.append("T")
                            elif self[y].seq[x] == 'G':
                                alleles.append("G")
                                alleles.append("G")
                            elif self[y].seq[x] == 'K':
                                alleles.append("G")
                                alleles.append("T")
                            elif self[y].seq[x] == 'T':
                                alleles.append("T")
                                alleles.append("T")
            return alleles

    def get_pi_per_site(self, x, pops=[]):
        # pi is the probability that two randomly chosen sequences from the
        # sample have different alleles at a site x.
        # Following Ruegg et al. (2014, Mol Ecol, A role for migration-linked genes and genomic
        # islands in divergence of a songbird), only biallelic (or in this case monomorphic)
        # SNPs are allowed.
        if isinstance(pops, list) == False:
            print("ERROR: Populations are not given as a list (get_pi_per_site):")
            print(pops)
            sys.exit(1)
        else:
            all_allele_frequencies = self.get_allele_frequencies(x, pops)
            two_allele_frequencies = []
            for all_allele_frequency in all_allele_frequencies:
                if all_allele_frequency != 0:
                    two_allele_frequencies.append(all_allele_frequency)
            if len(two_allele_frequencies) == 0:
                return 0
            elif len(two_allele_frequencies) == 1:
                return 0
            elif len(two_allele_frequencies) == 2:
                # Use r and a as in Ruegg et al. (2014).
                r = two_allele_frequencies[0]
                a = two_allele_frequencies[1]
                numerator = r * a
                denominator = scipy.special.binom((r+a),2)
                pi = numerator/denominator
                return pi
            elif len(two_allele_frequencies) > 2:
                return None

    def get_pi(self, pops=[]):
        pi_per_site_values = []
        l_k = self.get_alignment_length()
        for x in range(self.get_alignment_length()):
            pi_per_site = self.get_pi_per_site(x, pops)
            if pi_per_site is not None:
                pi_per_site_values.append(pi_per_site)
        if l_k == 0:
            return None
        else:
            return sum(pi_per_site_values)/l_k

    def get_F_st(self, pops=[]):
        if len(pops) != 2:
            print("ERROR: Exactly two populations must be specified to calculate pairwise F_st!")
            sys.exit(1)
        else:
            # Initiate six lists that will be needed for Fst calculation.
            # ninds_dup_1 and ninds_dup_2 and the numbers of individuals for pop1 and pop2, per locus.
            ninds_dup_1 = []
            ninds_dup_2 = []
            # p1 and p2 are the frequencies of allele A for pop1 and pop2, per locus.
            p1 = []
            p2 = []
            # oh1 and oh2 are the frequencies of heterozygotes for pop1 and pop2, per locus.
            oh1 = []
            oh2 = []

            # The number of populations, r is set to two.
            r = 2

            # Consider only biallelic loci, and only individuals without missing data.
            for x in range(self.get_alignment_length()):
                if self.get_is_biallelic_per_site(x, pops):

                    # Get all alleles at this site.
                    alleles_pop1 = self.get_alleles_per_site(x, [pops[0]])
                    alleles_pop2 = self.get_alleles_per_site(x, [pops[1]])
                    alleles_both_pops = alleles_pop1 + alleles_pop2
                    unique_alleles_both_pops = sorted(list(set(alleles_both_pops)))
                    good_unique_alleles_both_pops = []
                    for unique_allele_both_pops in unique_alleles_both_pops:
                        if unique_allele_both_pops in ["A", "C", "G", "T"]:
                            good_unique_alleles_both_pops.append(unique_allele_both_pops)

                    # Get the good (= no missing data) pairs of alleles for both pops, for this site.
                    good_pairs_of_alleles_pop1 = []
                    good_pairs_of_alleles_pop2 = []
                    for z in range(int(len(alleles_pop1)/2)):
                        if alleles_pop1[z*2] in good_unique_alleles_both_pops and alleles_pop1[(z*2)+1] in good_unique_alleles_both_pops:
                            good_pairs_of_alleles_pop1.append([alleles_pop1[z*2],alleles_pop1[(z*2)+1]])
                    for z in range(int(len(alleles_pop2)/2)):
                        if alleles_pop2[z*2] in good_unique_alleles_both_pops and alleles_pop2[(z*2)+1] in good_unique_alleles_both_pops:
                            good_pairs_of_alleles_pop2.append([alleles_pop2[z*2],alleles_pop2[(z*2)+1]])

                    # Get the number of good individuals in both pops based on the number of good pairs of alleles, for this site.
                    ninds_dup_1_this_site = len(good_pairs_of_alleles_pop1)
                    ninds_dup_2_this_site = len(good_pairs_of_alleles_pop2)

                    # Continue only if good individuals are found in both populations, for this site.
                    if ninds_dup_1_this_site > 0 and ninds_dup_2_this_site > 0:

                        # Get the frequencies of the two alleles, for this site.
                        allele_A_occurrences_pop1 = 0
                        allele_A_occurrences_pop2 = 0
                        for good_pair_of_alleles_pop1 in good_pairs_of_alleles_pop1:
                            for good_allele_pop1 in good_pair_of_alleles_pop1:
                                if good_allele_pop1 == good_unique_alleles_both_pops[0]:
                                    allele_A_occurrences_pop1 += 1
                                elif good_allele_pop1 != good_unique_alleles_both_pops[1]:
                                    print("ERROR: Unexpected allele found at site " + str(x) + ": " + good_allele_pop1 + "!")
                                    sys.exit(1)
                        for good_pair_of_alleles_pop2 in good_pairs_of_alleles_pop2:
                            for good_allele_pop2 in good_pair_of_alleles_pop2:
                                if good_allele_pop2 == good_unique_alleles_both_pops[0]:
                                    allele_A_occurrences_pop2 += 1
                                elif good_allele_pop2 != good_unique_alleles_both_pops[1]:
                                    print("ERROR: Unexpected allele found at site " + str(x) + ": " + good_allele_pop2 + "!")
                                    sys.exit(1)
                        p1_this_site = allele_A_occurrences_pop1/(len(good_pairs_of_alleles_pop1)*2)
                        p2_this_site = allele_A_occurrences_pop2/(len(good_pairs_of_alleles_pop2)*2)

                        # Get the frequencies of heterozygotes, for this site.
                        heterozygote_occurrences_pop1 = 0
                        heterozygote_occurrences_pop2 = 0
                        for good_pair_of_alleles_pop1 in good_pairs_of_alleles_pop1:
                            if good_pair_of_alleles_pop1[0] != good_pair_of_alleles_pop1[1]:
                                heterozygote_occurrences_pop1 += 1
                        for good_pair_of_alleles_pop2 in good_pairs_of_alleles_pop2:
                            if good_pair_of_alleles_pop2[0] != good_pair_of_alleles_pop2[1]:
                                heterozygote_occurrences_pop2 += 1
                        oh1_this_site = heterozygote_occurrences_pop1/len(good_pairs_of_alleles_pop1)
                        oh2_this_site = heterozygote_occurrences_pop2/len(good_pairs_of_alleles_pop2)

                        ninds_dup_1.append(ninds_dup_1_this_site)
                        ninds_dup_2.append(ninds_dup_2_this_site)
                        p1.append(p1_this_site)
                        p2.append(p2_this_site)
                        oh1.append(oh1_this_site)
                        oh2.append(oh2_this_site)

            # Calculation of Fst according to Weir & Cockerham (1984), as in method stamppFst of R package StAMPP.
            # n_bar is the average number of individuals in a population, for each locus.
            n_bar = []
            for x in range(len(ninds_dup_1)):
                n_bar.append((ninds_dup_1[x] + ninds_dup_2[x])/r)
            nc = []
            for x in range(len(ninds_dup_1)):
                nc.append((r * n_bar[x]) - (((ninds_dup_1[x]**2) + (ninds_dup_2[x]**2))/(r * n_bar[x])))
            # p_bar is the average sample frequency of allele A in a population, for each locus.
            p_bar = []
            for x in range(len(ninds_dup_1)):
                p_bar.append(((ninds_dup_1[x] * p1[x])/(r * n_bar[x])) + ((ninds_dup_2[x] * p2[x])/(r * n_bar[x])))
            # s_square is the sample variance of allele A frequencies over populations.
            s_square = []
            for x in range(len(ninds_dup_1)):
                s_square.append(((ninds_dup_1[x] * ((p1[x] - p_bar[x])**2))/n_bar[x]) + ((ninds_dup_2[x] * ((p2[x] - p_bar[x])**2))/n_bar[x]))
            # h_bar is the average heterozygote frequency for allele A.
            h_bar = []
            for x in range(len(ninds_dup_1)):
                h_bar.append(((ninds_dup_1[x] * oh1[x])/(r * n_bar[x])) + ((ninds_dup_2[x] * oh2[x])/(r * n_bar[x])))
            # Equation 2 in WC84.
            a = []
            for x in range(len(ninds_dup_1)):
                if n_bar[x] > 1:
                    a_cand = (n_bar[x]/nc[x]) * (s_square[x] - (1/(n_bar[x] - 1)) * ((p_bar[x] * (1 - p_bar[x])) - (((r - 1)/r) * s_square[x]) - ((1/4) * h_bar[x])))
                    if not math.isnan(a_cand):
                        a.append(a_cand)
            # Equation 3 in WC84.
            b = []
            for x in range(len(ninds_dup_1)):
                if n_bar[x] > 1:
                    b_cand = (n_bar[x]/(n_bar[x] - 1)) * ((p_bar[x] * (1 - p_bar[x])) - (((r - 1)/r) * s_square[x]) - (((2 * n_bar[x] - 1)/(4 * n_bar[x])) * h_bar[x]))
                    if not math.isnan(b_cand):
                        b.append(b_cand)
            # Equation 4 in WC84.
            c = []
            for x in range(len(ninds_dup_1)):
                if n_bar[x] > 1:
                    c_cand = (1/2) * h_bar[x]
                    if not math.isnan(c_cand):
                        c.append(c_cand)
            if (sum(a) + sum(b) + sum(c)) == 0:
                return None
            else:
                fst = sum(a)/(sum(a) + sum(b) + sum(c))
                return fst

    def get_d_xy(self, pops=[], l_k=-1):
        if len(pops) != 2:
            print("ERROR: Exactly two populations must be specified to calculate pairwise d_xy!")
            sys.exit(1)
        else:
            if l_k == -1:
                l_k = self.get_alignment_length()
            if l_k == 0:
                d_xy = None
            else:
                sum_of_quotients = 0
                for x in range(self.get_alignment_length()):
                    if self.get_is_biallelic_per_site(x, pops):
                        all_allele_frequencies0 = self.get_allele_frequencies(x, [pops[0]])
                        all_allele_frequencies1 = self.get_allele_frequencies(x, [pops[1]])
                        if sum(all_allele_frequencies0) > 0 and sum(all_allele_frequencies1) > 0:
                            two_allele_frequencies0 = []
                            two_allele_frequencies1 = []
                            if all_allele_frequencies0[0] != 0 or all_allele_frequencies1[0] != 0:
                                two_allele_frequencies0.append(all_allele_frequencies0[0])
                                two_allele_frequencies1.append(all_allele_frequencies1[0])
                            if all_allele_frequencies0[1] != 0 or all_allele_frequencies1[1] != 0:
                                two_allele_frequencies0.append(all_allele_frequencies0[1])
                                two_allele_frequencies1.append(all_allele_frequencies1[1])
                            if all_allele_frequencies0[2] != 0 or all_allele_frequencies1[2] != 0:
                                two_allele_frequencies0.append(all_allele_frequencies0[2])
                                two_allele_frequencies1.append(all_allele_frequencies1[2])
                            if all_allele_frequencies0[3] != 0 or all_allele_frequencies1[3] != 0:
                                two_allele_frequencies0.append(all_allele_frequencies0[3])
                                two_allele_frequencies1.append(all_allele_frequencies1[3])
                            if len(two_allele_frequencies0) != 2:
                                print("ERROR: Wrong number of allele frequencies!")
                                print(all_allele_frequencies0)
                                print(all_allele_frequencies1)
                                sys.exit(1)
                            r0 = two_allele_frequencies0[0]
                            a0 = two_allele_frequencies0[1]
                            r1 = two_allele_frequencies1[0]
                            a1 = two_allele_frequencies1[1]
                            numerator = r0*a1 + r1*a0
                            denominator = (r0+a0) * (r1+a1)
                            sum_of_quotients += numerator/denominator
                d_xy = sum_of_quotients/l_k
            return d_xy

    def get_d_f(self, pops=[]):
        if len(pops) != 2:
            print("ERROR: Exactly two populations must be specified to calculate pairwise d_f!")
            sys.exit(1)
        else:
            l_k = self.get_alignment_length()
            if l_k == 0:
                d_f = None
            else:
                number_of_fixed_snps = 0
                for x in range(self.get_alignment_length()):
                    if self.get_is_biallelic_per_site(x, pops):
                        pop0_variable = self.get_is_variable_per_site(x, [pops[0]])
                        pop1_variable = self.get_is_variable_per_site(x, [pops[1]])
                        if pop0_variable == False and pop1_variable == False:
                            number_of_fixed_snps += 1
                d_f = number_of_fixed_snps/l_k
            return d_f

    def get_afd(self, pops=[]):
        if len(pops) != 2:
            print("ERROR: Exactly two populations must be specified to calculate pairwise afd!")
            sys.exit(1)
        else:
            afds = []
            for x in range(self.get_alignment_length()):
                if self.get_is_biallelic_per_site(x, pops):
                    all_allele_frequencies0 = self.get_allele_frequencies(x, [pops[0]])
                    all_allele_frequencies1 = self.get_allele_frequencies(x, [pops[1]])
                    if sum(all_allele_frequencies0) > 0 and sum(all_allele_frequencies1) > 0:
                        two_allele_frequencies0 = []
                        two_allele_frequencies1 = []
                        if all_allele_frequencies0[0] != 0 or all_allele_frequencies1[0] != 0:
                            two_allele_frequencies0.append(all_allele_frequencies0[0])
                            two_allele_frequencies1.append(all_allele_frequencies1[0])
                        if all_allele_frequencies0[1] != 0 or all_allele_frequencies1[1] != 0:
                            two_allele_frequencies0.append(all_allele_frequencies0[1])
                            two_allele_frequencies1.append(all_allele_frequencies1[1])
                        if all_allele_frequencies0[2] != 0 or all_allele_frequencies1[2] != 0:
                            two_allele_frequencies0.append(all_allele_frequencies0[2])
                            two_allele_frequencies1.append(all_allele_frequencies1[2])
                        if all_allele_frequencies0[3] != 0 or all_allele_frequencies1[3] != 0:
                            two_allele_frequencies0.append(all_allele_frequencies0[3])
                            two_allele_frequencies1.append(all_allele_frequencies1[3])
                        if len(two_allele_frequencies0) != 2:
                            print("ERROR: Wrong number of allele frequencies!")
                            print(all_allele_frequencies0)
                            print(all_allele_frequencies1)
                            sys.exit(1)
                        af0 = two_allele_frequencies0[0]/(two_allele_frequencies0[0]+two_allele_frequencies0[1])
                        af1 = two_allele_frequencies1[0]/(two_allele_frequencies1[0]+two_allele_frequencies1[1])
                        afds.append(abs(af0-af1))
            if len(afds) == 0:
                return None
            else:
                return sum(afds)/len(afds)

# Parse the command line arguments.
parser = argparse.ArgumentParser(
    formatter_class=argparse.RawDescriptionHelpFormatter,
    description=textwrap.dedent('''\
      %(prog)s
    -----------------------------------------
      Calculate measures of genetic variation and population
      divergence over sliding windows or fixed regions in a vcf file.
    '''))
parser.add_argument('-v', '--version', action='version', version='%(prog)s 1.01')
parser.add_argument(
    '-p', '--populations',
    nargs='*',
    type=str,
    help="Exactly two population identifiers are required. The vcf file may contain samples from other populations, these will be ignored.")
parser.add_argument(
    '-c', '--chromosome',
    nargs=1,
    type=str,
    help='The name of the chromosome/linkage group (required).')
parser.add_argument(
    '-w', '--window-size',
    nargs=1,
    type=int,
    default=[1000],
    dest='window_size',
    help="Size of the sliding window in bp (default: 1000)."
    )
parser.add_argument(
    '-f', '--from',
    nargs=1,
    type=int,
    default=[1],
    dest='start',
    help="First position of vcf to be considered (count starts at 1) (default: 1)."
    )
parser.add_argument(
    '-t', '--to',
    nargs=1,
    type=int,
    default=[-1],
    dest='end',
    help="Last position of vcf to be considered (count starts at 1) (default: not used)."
    )
parser.add_argument('infile', nargs='?', type=argparse.FileType('r'), default='-', help='The input file name.')
parser.add_argument('outfile', nargs='?', type=argparse.FileType('w'), default=sys.stdout, help='The output file name.')
args = parser.parse_args()
infile = args.infile
outfile = args.outfile
pops = args.populations
window_size = args.window_size[0]
first_pos = args.start[0]
last_pos = args.end[-1]
if args.chromosome == None:
    print("ERROR: A chromosome id must be specified with argument '-c'!")
    sys.exit(1)
chromosome = args.chromosome[0]
if infile.isatty():
    print('No input file specified, and no input piped through stdin!')
    sys.exit(0)
if len(pops) != 2:
    print("ERROR: " + str(len(pops)) + " population identifiers were specified, but exactly two population identifiers were expected!")
    sys.exit(1)

# Read the infile.
vcf_data = vcf.Reader(infile)
all_sample_ids = vcf_data.samples
sample_ids = []
number_of_sample_ids_for_pop1 = 0
number_of_sample_ids_for_pop2 = 0
for sample_id in all_sample_ids:
    if pops[0] in sample_id:
        if pops[1] in sample_id:
            print("ERROR: Sample id " + sample_id + " contains both population identifiers!")
            sys.exit(1)
        else:
            number_of_sample_ids_for_pop1 += 1
            sample_ids.append(sample_id)
    elif pops[1] in sample_id:
        number_of_sample_ids_for_pop2 += 1
        sample_ids.append(sample_id)

if number_of_sample_ids_for_pop1 == 0:
    print("ERROR: No samples found for pop1!")
    print("pop1 is: " + pops[0])
    print("sample ids of pop1 are: ")
    for sample_id in all_sample_ids:
        if pops[0] in sample_id:
            print(sample_id)
    sys.exit(1)
if number_of_sample_ids_for_pop2 == 0:
    print("ERROR: No samples found for pop2!")
    print("pop2 is: " + pops[1])
    print("sample ids of pop2 are: ")
    for sample_id in all_sample_ids:
        if pops[1] in sample_id:
            print(sample_id)
    sys.exit(1)

records = []
record_chromosomes = []
record_positions = []
first_record = None
for record in vcf_data:
    # Skip this record if it is multi-allelic.
    if len(record.REF) > 1:
        continue
    if type(first_record) == type(None):
        first_record = record
    record_positions.append(record.POS)
    record_chromosomes.append(record.CHROM)
    record_genotypes = []
    if len(record.samples) != len(first_record.samples):
        print("Record lengths differ!", file=sys.stderr)
        sys.exit(1)
    for sample in record:
        if pops[0] in sample.sample or pops[1] in sample.sample:
            if sample.gt_bases == "A/A" or sample.gt_bases == "A|A":
                record_genotypes.append("A")
            elif sample.gt_bases == "A/C" or sample.gt_bases == "A|C":
                record_genotypes.append("M")
            elif sample.gt_bases == "A/G" or sample.gt_bases == "A|G":
                record_genotypes.append("R")
            elif sample.gt_bases == "A/T" or sample.gt_bases == "A|T":
                record_genotypes.append("W")
            elif sample.gt_bases == "C/A" or sample.gt_bases == "C|A":
                record_genotypes.append("M")
            elif sample.gt_bases == "C/C" or sample.gt_bases == "C|C":
                record_genotypes.append("C")
            elif sample.gt_bases == "C/G" or sample.gt_bases == "C|G":
                record_genotypes.append("S")
            elif sample.gt_bases == "C/T" or sample.gt_bases == "C|T":
                record_genotypes.append("Y")
            elif sample.gt_bases == "G/A" or sample.gt_bases == "G|A":
                record_genotypes.append("R")
            elif sample.gt_bases == "G/C" or sample.gt_bases == "G|C":
                record_genotypes.append("S")
            elif sample.gt_bases == "G/G" or sample.gt_bases == "G|G":
                record_genotypes.append("G")
            elif sample.gt_bases == "G/T" or sample.gt_bases == "G|T":
                record_genotypes.append("K")
            elif sample.gt_bases == "T/A" or sample.gt_bases == "T|A":
                record_genotypes.append("W")
            elif sample.gt_bases == "T/C" or sample.gt_bases == "T|C":
                record_genotypes.append("Y")
            elif sample.gt_bases == "T/G" or sample.gt_bases == "T|G":
                record_genotypes.append("K")
            elif sample.gt_bases == "T/T" or sample.gt_bases == "T|T":
                record_genotypes.append("T")
            elif sample.gt_bases == None:
                record_genotypes.append("N")
            else:
                print("ERROR: Unexpected genotype of sample" + str(sample.sample) + " at site " + str(record.POS) + ": " + str(sample.gt_bases), file=sys.stderr)
                sys.exit(1)
    records.append(record_genotypes)

# Test whether all records are of the same length.
for y in range(1,len(records)):
    if len(records[y]) != len(records[0]):
        print("The length of record " + str(y) + " (" + str(len(records[y])) + ") differs from that of record 0 (" + str(len(records[0])) + ")!" , file=sys.stderr)
        sys.exit(1)

# Find the maximum position on the specified chromosome.
max_position_on_chromosome = 0
for y in range(0,len(records)):
    if record_chromosomes[y] == chromosome:
        if record_positions[y] > max_position_on_chromosome:
            max_position_on_chromosome = record_positions[y]

# Run a loop over sliding windows until the window end position exceeds the maximum position.
window_start = first_pos
if last_pos != -1:
    window_end = last_pos
    window_size = last_pos - first_pos + 1
else:
    window_end = window_start + window_size - 1
outstring = "window_center     num_snps   gc_content    prop_var1    prop_var2  prop_var1+2    prop_het1    prop_het2  prop_het1+2           pi         f_st        d_xy*       d_xy**          d_f          afd\n"
outfile.write(outstring)
while window_end < max_position_on_chromosome:
    # Produce a multiple sequence alignment object with only SNPs within this window.
    seqs = []
    for x in range(0,len(records[0])):
        seq = []
        for y in range(0,len(records)):
            if record_chromosomes[y] == chromosome:
                if record_positions[y] >= window_start and record_positions[y] <= window_end:
                    this_record = records[y]
                    this_position_of_this_record = this_record[x]
                    seq.append(this_position_of_this_record)
        # seqs.append(SeqRecord(Seq("".join(seq),generic_dna),id=sample_ids[x].replace('-','_')))
        # seqs.append(SeqRecord(Seq("".join(seq),generic_dna),id=sample_ids[x]))
        seqs.append(SeqRecord(Seq("".join(seq)),id=sample_ids[x]))
    align = XMultipleSeqAlignment(seqs)
    num_snps = align.get_alignment_length()
    gc_content = align.get_gc_content()
    prop_var1 = align.get_number_of_variable_sites([pops[0]])/window_size
    prop_var2 = align.get_number_of_variable_sites([pops[1]])/window_size
    prop_var12 = align.get_number_of_variable_sites([pops[0], pops[1]])/window_size
    prop_het1 = align.get_heterozygosity([pops[0]])
    prop_het2 = align.get_heterozygosity([pops[1]])
    prop_het12 = align.get_heterozygosity([pops[0], pops[1]])
    pi = align.get_pi([pops[0], pops[1]])
    f_st = align.get_F_st([pops[0], pops[1]])
    d_xy_asterisk = align.get_d_xy([pops[0], pops[1]], window_size)
    d_xy_doubleasterisk = align.get_d_xy([pops[0], pops[1]], align.get_alignment_length())
    d_f = align.get_d_f([pops[0], pops[1]])
    afd = align.get_afd([pops[0], pops[1]])
    outstring = "{0:.1f}".format(((window_start+window_end-1)/2)).rjust(13)
    outstring += str(num_snps).rjust(13)
    if gc_content == None:
        outstring += "NA".rjust(13)
    else:
        outstring += "{0:.4f}".format(gc_content).rjust(13)
    if prop_var1 == None:
        outstring += "NA".rjust(13)
    else:
        outstring += "{0:.6f}".format(prop_var1).rjust(13)
    if prop_var2 == None:
        outstring += "NA".rjust(13)
    else:
        outstring += "{0:.6f}".format(prop_var2).rjust(13)
    if prop_var12 == None:
        outstring += "NA".rjust(13)
    else:
        outstring += "{0:.6f}".format(prop_var12).rjust(13)
    if prop_het1 == None:
        outstring += "NA".rjust(13)
    else:
        outstring += "{0:.4f}".format(prop_het1).rjust(13)
    if prop_het2 == None:
        outstring += "NA".rjust(13)
    else:
        outstring += "{0:.4f}".format(prop_het2).rjust(13)
    if prop_het12 == None:
        outstring += "NA".rjust(13)
    else:
        outstring += "{0:.4f}".format(prop_het12).rjust(13)
    if pi == None:
        outstring += "NA".rjust(13)
    else:
        outstring += "{0:.4f}".format(pi).rjust(13)
    if f_st == None:
        outstring += "NA".rjust(13)
    else:
        outstring += "{0:.4f}".format(f_st).rjust(13)
    if d_xy_asterisk == None:
        outstring += "NA".rjust(13)
    else:
        outstring += "{0:.4f}".format(d_xy_asterisk).rjust(13)
    if d_xy_doubleasterisk == None:
        outstring += "NA".rjust(13)
    else:
        outstring += "{0:.4f}".format(d_xy_doubleasterisk).rjust(13)
    if d_f == None:
        outstring += "NA".rjust(13)
    else:
        outstring += "{0:.4f}".format(d_f).rjust(13)
    if afd == None:
        outstring += "NA".rjust(13)
    else:
        outstring += "{0:.4f}".format(afd).rjust(13)
    outstring += "\n"
    outfile.write(outstring)
    if last_pos != -1:
        break
    else:
        window_start += window_size
        window_end += window_size

# footnote_string = "\n"
# footnote_string += "*  Following Ruegg et al. (2014), using the window size as 'contig length' Lk.\n"
# footnote_string += "   Thus this dxy takes into account the number of invariant sites, and will vary depending on filtering options.\n"
# footnote_string += "** As above, but using the number of bi-allelic SNPs as 'contig length' Lk.\n"
# footnote_string += "   This dxy ignores all positions that are not bi-allelic SNPs. It should be more robust to vcf filtering and\n"
# footnote_string += "   is equivalent to the mean contribution of a SNP to the first dxy measure (above).\n"
# outfile.write(footnote_string)
