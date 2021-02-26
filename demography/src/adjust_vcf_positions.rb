# m_matschiner Wed Jun 10 14:22:13 CEST 2020

# Get the command-line arguments.
vcf_file_name = ARGV[0]
shift = ARGV[1].to_i

# Read the vcf.
vcf_file = File.open(vcf_file_name)
vcf_lines = vcf_file.readlines
vcf_file.close
vcf_prep_string = ""
vcf_lines.each do |l|
	if l[0] == "#"
		vcf_prep_string << "#{l}"
	else
		line_ary = l.split
		line_ary[1] = (line_ary[1].to_i + shift).to_s
		vcf_prep_string << line_ary.join("\t")
		vcf_prep_string << "\n"
	end
end

# Write the adjusted vcf file.
vcf_file = File.open(vcf_file_name, "w")
vcf_file.write(vcf_prep_string)
