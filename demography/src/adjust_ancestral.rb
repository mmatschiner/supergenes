# m_matschiner Mon Aug 17 13:32:00 CEST 2020

# Get the command-line arguments.
vcf_file_name = ARGV[0]
adjusted_vcf_file_name = ARGV[1]
reference_file_name = ARGV[2]
ancestor_file_name = ARGV[3]

# Read the reference.
ref_ids = []
ref_seqs = []
reference_file = File.open(reference_file_name)
reference_file.each do |l|
	if l[0] == ">"
		ref_ids << l.strip[1..-1]
		ref_seqs << ""
	elsif l.strip != ""
		ref_seqs.last << l.strip
	end
end

# Read the vcf.
vcf_file = File.open(vcf_file_name)
adjusted_vcf_file = File.open(adjusted_vcf_file_name, "w")
n_multi_nucleotide = 0
n_no_blast_hits = 0
n_weak_blast_hits = 0
n_good_blast_hits = 0
n_gts_rotated = 0
n_ancestral_gt_is_ref_gt = 0
n_ancestral_gt_is_third_gt = 0
previous_lg = nil
vcf_file.each do |l|
	if l[0] == "#"
		adjusted_vcf_file.write(l)
	else
		line_ary = l.split
		vcf_lg = line_ary[0]
		vcf_pos = line_ary[1].to_i
		# Feedback
		print "#{vcf_lg}\t#{vcf_pos}\t"
		vcf_ref_allele = line_ary[3]
		vcf_alt_allele = line_ary[4]
		vcf_gts = line_ary[9..-1]
		ref_seq = ref_seqs[ref_ids.index(vcf_lg)]
		ref_gt = ref_seq[vcf_pos-1]
		if vcf_ref_allele.size > 1 or vcf_alt_allele.size > 1
			n_multi_nucleotide += 1
			puts "multi-nucleotide"
		else
			# Make sure that the reference genotype is correctly read.
			unless ref_gt == vcf_ref_allele
				puts "ERROR: Reference genotype (#{ref_gt}) differs from reference allele in vcf (#{vcf_ref_allele}) at position #{vcf_pos}!"
				exit
			end
			# Extract the query sequence.
			query_start_pos = [0,vcf_pos-1-200].max
			query_end_pos = [ref_seq.size-1,vcf_pos-1+200].min
			query_seq = ref_seq[query_start_pos..query_end_pos]
			vcf_pos_in_query_pos = vcf_pos-1 - query_start_pos
			# Make sure that the site position in the query is correctly interpreted.
			unless ref_gt == query_seq[vcf_pos_in_query_pos]
				puts "ERROR: Reference genotype (#{ref_gt}) differs from site in query (#{query_seq[vcf_pos_in_query_pos]}) at position #{vcf_pos}!"
				exit 1
			end
			# Remove gaps from the query sequence.
			if query_seq.include?("-")
				n_gaps_between_start_and_site_in_query_seq = query_seq[0..vcf_pos_in_query_pos-1].count("-")
				query_seq.delete!("-")
				vcf_pos_in_query_pos -= n_gaps_between_start_and_site_in_query_seq
				# Make sure that the site position in the query is correctly interpreted.
				unless ref_gt == query_seq[vcf_pos_in_query_pos]
					puts "ERROR: Reference genotype (#{ref_gt}) differs from site in query (#{query_seq[vcf_pos_in_query_pos]}) at position #{vcf_pos}!"
					exit 1
				end
			end
			# Prepare a query string for blast searches.
			query_string = ">query\n"
			query_string << "#{query_seq}\n"
			query_file = File.open("tmp.query.fasta","w")		
			query_file.write(query_string)
			query_file.close
			# Run blast.
			blast_res = `blastn -query tmp.query.fasta -db #{ancestor_file_name} -culling_limit 1 -max_target_seqs 1 -outfmt '6 evalue bitscore pident qstart qseq sseq' 2>/dev/null`
			blast_res_lines = blast_res.split("\n")
			if blast_res_lines.size > 0
				blast_res_line = blast_res_lines[0]
				blast_res_ary = blast_res_line.split
				blast_evalue = blast_res_ary[0].to_f
				blast_bitscore = blast_res_ary[1].to_f
				blast_pident = blast_res_ary[2].to_f
				blast_qstart = blast_res_ary[3].to_i
				blast_qseq = blast_res_ary[4]
				blast_sseq = blast_res_ary[5]
				if blast_evalue < 10**(-5) and blast_bitscore > 50 and blast_pident > 0.95 and blast_qseq.delete("-").size > 200
					vcf_pos_in_blasted_query_pos_without_gaps = vcf_pos_in_query_pos+1 - blast_qstart
					# Get the corresponding position in the blasted query with gaps.
					vcf_pos_in_blasted_query_pos_with_gaps = 0
					count_non_missing = 0
					until count_non_missing == vcf_pos_in_blasted_query_pos_without_gaps
						vcf_pos_in_blasted_query_pos_with_gaps += 1
						count_non_missing += 1 if blast_qseq[vcf_pos_in_blasted_query_pos_with_gaps] != "-"
					end
					# Make sure that the site position in the blasted query is correctly interpreted.
					unless ref_gt == blast_qseq[vcf_pos_in_blasted_query_pos_with_gaps]
						puts "ERROR: Reference genotype (#{ref_gt}) differs from site in blasted query (#{blast_qseq[vcf_pos_in_blasted_query_pos_with_gaps]}) at position #{vcf_pos}!"
						exit 1
					end
					ancestral_genotype = blast_sseq[vcf_pos_in_blasted_query_pos_with_gaps]
					n_good_blast_hits += 1
					if ancestral_genotype == ref_gt
						n_ancestral_gt_is_ref_gt += 1
						adjusted_vcf_file.write(l)
						puts "ancestral gt is ref"
					elsif ancestral_genotype == vcf_alt_allele
						n_gts_rotated += 1
						# Rotate all genotypes.
						rotated_line = "#{vcf_lg}\t#{vcf_pos}\t.\t#{ancestral_genotype}\t#{ref_gt}\t.\tPASS\t.\tGT\t"
						vcf_gts.each do |gt|
							if gt == "0|0"
								rotated_line << "1|1\t"
							elsif gt == "0|1"
								rotated_line << "1|0\t"
							elsif gt == "1|0"
								rotated_line << "0|1\t"
							elsif gt == "1|1"
								rotated_line << "0|0\t"
							else
								puts "ERROR: Unexpected genotype: #{gt}!"
								exit 1
							end
						end
						rotated_line.chomp!
						rotated_line << "\n"
						adjusted_vcf_file.write(rotated_line)
						puts "rotated"
					else
						n_ancestral_gt_is_third_gt += 1
						adjusted_vcf_file.write(l)
						puts "ancestral is third"
					end
				else
					adjusted_vcf_file.write(l)
					n_weak_blast_hits += 1
					puts "weak blast hit"
				end
			else
				adjusted_vcf_file.write(l)
				n_no_blast_hits += 1
				puts "no blast hit"
			end
		end
	end
end

# Feedback.
puts "n_no_blast_hits: #{n_no_blast_hits}"
puts "n_weak_blast_hits: #{n_weak_blast_hits}"
puts "n_good_blast_hits: #{n_good_blast_hits}"
puts "n_gts_rotated: #{n_gts_rotated}"
puts "n_ancestral_gt_is_ref_gt: #{n_ancestral_gt_is_ref_gt}"
puts "n_ancestral_gt_is_third_gt: #{n_ancestral_gt_is_third_gt}"
