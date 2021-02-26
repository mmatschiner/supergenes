# m_matschiner Mon Oct 8 20:32:56 CEST 2018

# Get the command-line arguments.
vcf_file_name = ARGV[0]
region_string = ARGV[1]
min_likelihood = ARGV[2].to_i
alignment_file_name = ARGV[3]

# Read the vcf file.
vcf_file = File.open(vcf_file_name)
vcf_lines = vcf_file.readlines

# Get the sample ids.
ids = []
vcf_lines.each do |l|
	if l[0..5] == "#CHROM"
		ids = l.split[9..-1]
		break
	end
end

# Make sure that sample IDs were found.
if ids == []
	puts "ERROR: No sample IDs could be found!"
	exit(1)
end

# Get the window coordinates from the region string.
lg = region_string.split(":")[0]
from = region_string.split(":")[1].split("-")[0].to_i
to = region_string.split(":")[1].split("-")[1].to_i

# Initiate the sequences as composed only of missing data.
seqs = []
ids.size.times do |x|
	seq = ""
	from.upto(to) { seq << "N" }
	seqs << seq
end

# Read sequence data from the vcf
vcf_lines.each do |l|
	if l[0] != "#"
		line_ary = l.split
		pos = line_ary[1].to_i
		if line_ary[0] == lg and pos >= from and pos <= to
			ref = line_ary[3]
			fmt = line_ary[8]
			if line_ary[4] == "." # If this site is invariant.
				unless fmt == "GT:AD:DP"
					rgq_index = fmt.split(":").index("RGQ")
					if rgq_index == nil
						puts "ERROR: Format string (#{fmt}) could not be read!"
						puts "  Record is:"
						puts l
						exit(1)
					end
					ids.size.times do |x|
						ary_index = x+9
						if line_ary[ary_index].split(":")[rgq_index].to_i >= min_likelihood
							seqs[x][pos-from] = ref
						end
					end
				end
			else # If this site is variant.
				alts = line_ary[4].split(",")
				all_alleles = [ref]
				alts.each { |a| all_alleles << a }
				pl_index = fmt.split(":").index("PL")
				if pl_index == nil
					puts "ERROR: Format string could not be read!"
					exit(1)
				end
				if all_alleles.size == 2
					ids.size.times do |x|
						ary_index = x+9
						pls = line_ary[ary_index].split(":")[pl_index].split(",")
						gt = "N"
						unless pls == ["."]
							if pls.size == 3
								if pls[0] == "0"
									gt = all_alleles[0] if pls[2].to_i >= min_likelihood
								elsif pls[2] == "0"
									gt = all_alleles[1] if pls[0].to_i >= min_likelihood
								elsif pls[1] == "0"
									if pls[0].to_i >= min_likelihood and pls[2].to_i >= min_likelihood
										if pls[0].to_i == pls[2].to_i
											gt = all_alleles.sample
										elsif pls[0].to_i < pls[2].to_i
											gt = all_alleles[0]
										elsif pls[0].to_i > pls[2].to_i
											gt = all_alleles[1]
										else
											puts "ERROR: Unexpected case!"
											exit(1)
										end
									elsif pls[0].to_i >= min_likelihood and pls[2].to_i < min_likelihood
										gt = all_alleles[1]
									elsif pls[0].to_i < min_likelihood and pls[2].to_i >= min_likelihood
										gt = all_alleles[0]
									end
								else
									puts "ERROR: No PL value is 0 (PL: #{pls.join(", ")})!"
									exit(1)
								end
							else
								puts "ERROR: Unexpected combination of allele number (#{all_alleles.size}) and number of likelihoods (#{pls.size})!"
								exit(1)
							end
						end
						seqs[x][pos-from] = gt if gt != "N"
					end
				elsif all_alleles.size == 3
					ids.size.times do |x|
						ary_index = x+9
						pls = line_ary[ary_index].split(":")[pl_index].split(",")
						gt = "N"
						unless pls == ["."]
							if pls.size == 6
								if pls[0] == "0"
									gt = all_alleles[0] if pls[2].to_i >= min_likelihood and pls[4].to_i >= min_likelihood and pls[5].to_i >= min_likelihood
								elsif pls[2] == "0"
									gt = all_alleles[1] if pls[0].to_i >= min_likelihood and pls[3].to_i >= min_likelihood and pls[5].to_i >= min_likelihood
								elsif pls[5] == "0"
									gt = all_alleles[2] if pls[0].to_i >= min_likelihood and pls[1].to_i >= min_likelihood and pls[2].to_i >= min_likelihood
								elsif pls[1] == "0"
									if pls[0].to_i == pls[2].to_i
										gt = [all_alleles[0],all_alleles[1]].sample
									elsif pls[0].to_i < pls[2].to_i
										gt = all_alleles[0]
									elsif pls[0].to_i > pls[2].to_i
										gt = all_alleles[1]
									else
										puts "ERROR: Unexpected case!"
										exit(1)
									end
								elsif pls[3] == "0"
									if pls[0].to_i == pls[5].to_i
										gt = [all_alleles[0],all_alleles[2]].sample
									elsif pls[0].to_i < pls[5].to_i
										gt = all_alleles[0]
									elsif pls[0].to_i > pls[5].to_i
										gt = all_alleles[2]
									else
										puts "ERROR: Unexpected case!"
										exit(1)
									end
								elsif pls[4] == "0"
									if pls[2].to_i == pls[5].to_i
										gt = [all_alleles[1],all_alleles[2]].sample
									elsif pls[2].to_i < pls[5].to_i
										gt = all_alleles[1]
									elsif pls[2].to_i > pls[5].to_i
										gt = all_alleles[2]
									else
										puts "ERROR: Unexpected case!"
										exit(1)
									end
								else
									puts "ERROR: No PL value is 0 (PL: #{pls.join(", ")})!"
									exit(1)
								end
							else
								puts "ERROR: Unexpected combination of allele number (#{all_alleles.size}) and number of likelihoods (#{pls.size}) at pos #{pos}!"
								puts "  Record is:"
								puts l
								exit(1)
							end
						end
						seqs[x][pos-from] = gt if gt != "N"
					end
				elsif all_alleles.size != 4
					puts "ERROR: Unexpected allele number (#{all_alleles.size}) at position #{pos}!"
					exit(1)
				end
			end
		end
	end
end

# Get the longest id length.
max_id_length = 0
ids.each { |i| max_id_length = i.size if i.size > max_id_length }

# Prepare the alignment string.
alignment_string = "#{ids.size} #{to-from+1}\n"
ids.size.times do |x|
	alignment_string << "#{ids[x].ljust(max_id_length)}  #{seqs[x]}\n"
end

# Write the alignment file.
alignment_file = File.open(alignment_file_name, "w")
alignment_file.write(alignment_string)
