# m_matschiner Sun Aug 30 14:18:29 CEST 2020

# Get the command-line arguments.
in_vcf_file_name = ARGV[0]
anc_id = ARGV[1]
out_vcf_file_name = ARGV[2]
bed_file_name = ARGV[3]

# Open the output file.
out_vcf_file = File.open(out_vcf_file_name, "w")

# Open the bed file.
bed_file = File.open(bed_file_name, "w")
bed_file.write("chrom\tchromStart\tchromEnd\n")

# Read the input vcf file.
in_vcf_file = File.open(in_vcf_file_name)
anc_index = nil
n_sites_excluded = 0
in_vcf_file.each do |l|
	if l[0..1] == "##"
		out_vcf_file.write(l)
	elsif l[0] == "#"
		header_ary = l.split
		if header_ary.include?(anc_id)
			anc_index = header_ary.index(anc_id)
		else
			puts "ERROR: Ancestor id #{anc_id} not included in header!"
			exit 1
		end
		out_vcf_file.write(l)
	else
		line_ary = l.split
		lg = line_ary[0]
		pos = line_ary[1].to_i
		if line_ary[anc_index] == "0|0" or line_ary[anc_index] == "1|1"
			out_vcf_file.write(l)
		else
			n_sites_excluded += 1
			bed_file.write("#{lg}\t#{pos-1}\t#{pos}\n")
		end
	end
end

# Report.
puts "Excluded #{n_sites_excluded} sites heterozygous in ancestor #{anc_id}."