# m_matschiner Mon Aug 31 16:15:44 CEST 2020

# Get the command-line arguments.
vcf_file_name = ARGV[0]
ref_file_name = ARGV[1]
sample_name = ARGV[2]
fasta_file_name = ARGV[3]

# Read the reference.
ref_file = File.open(ref_file_name)
ref_lines = ref_file.readlines
ref_ids = []
ref_seqs = []
ref_lines.each do |l|
	if l[0] == ">"
		ref_ids << l.strip[1..-1]
		ref_seqs << ""
	elsif l.strip != ""
		ref_seqs.last << l.strip
	end
end

# Read the vcf file line by line for each reference sequence id, and only focus on those sites from that sequence.
fasta_string = ""
ref_ids.size.times do |x|
	print "Preparing the fasta sequence for id #{ref_ids[x]}..."
	fasta_string << ">#{ref_ids[x]}\n"
	fasta_seq = ""
	sample_index = nil
	vcf_file = File.open(vcf_file_name)
	vcf_file.each do |l|
		if l[0..1] == "##"
			next
		elsif l[0] == "#"
			header_ary = l.split
			if header_ary.include?(sample_name)
				sample_index = header_ary.index(sample_name)
			else
				puts "ERROR: Sample #{sample_name} not found in VCF header!"
				exit 1
			end
		else
			line_ary = l.split
			lg = line_ary[0]
			if lg == ref_ids[x]
				pos = line_ary[1].to_i
				fasta_seq << ref_seqs[x][fasta_seq.size..pos-2] if pos > 1
				ref_allele = line_ary[3]
				alt_allele = line_ary[4]
				sample_gt = line_ary[sample_index].split(":")[0]
				sample_gt.sub!("|","/")
				ref_lg_index = nil
				if ref_ids.include?(lg)
					ref_lg_index = ref_ids.index(lg)
					if sample_gt == ("0/0")
						fasta_seq << ref_allele
					elsif sample_gt == ("1/1")
						fasta_seq << alt_allele
					else
						puts "ERROR: Unexpected sample genotype (#{sample_gt})!"
						exit 1
					end
				else
					puts "ERROR! LG #{lg} not found in reference!"
					exit 1
				end
			end
		end
	end
	fasta_seq << ref_seqs[x][fasta_seq.size..-1]
	fasta_string << "#{fasta_seq}\n"
	vcf_file.close
	puts " done."
end

# Write the output file.
fasta_file = File.open(fasta_file_name, "w")
fasta_file.write(fasta_string)
puts "Wrote file #{fasta_file_name}."