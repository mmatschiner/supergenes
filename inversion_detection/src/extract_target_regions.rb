# m_matschiner Mon Jan 16 16:44:00 CET 2017

# Get the command line options.
assembly_file_name = ARGV[0]
target_regions_table_file_name = ARGV[1]

# Read the target regions table.
target_regions_table_file = File.open(target_regions_table_file_name)
target_regions_table_lines = target_regions_table_file.readlines
target_lgs = []
target_froms = []
target_tos = []
target_regions_table_lines[1..-1].each do |l|
	line_ary = l.split
	target_lgs << line_ary[0]
	target_froms << line_ary[1].to_i
	target_tos << line_ary[2].to_i
end

# Read the assembly file.
assembly_file = File.open(assembly_file_name)
assembly_lines = assembly_file.readlines
fasta_ids = []
fasta_seqs = []
assembly_lines.each do |l|
	if l[0] == ">"
		fasta_ids << l[1..-1].strip
		fasta_seqs << ""
	elsif l.strip != ""
		fasta_seqs.last << l.strip
	end
end

# Extract the target regions.
target_lgs.size.times do |x|
	extract_id = ""
	extract_seq = ""
	fasta_ids.size.times do |y|
		if target_lgs[x].downcase == fasta_ids[y].downcase
			extract_id = "#{target_lgs[x]}_#{target_froms[x]}_#{target_tos[x]}"
			extract_from = target_froms[x]*1000000
			extract_to = (target_tos[x]*1000000)-1
			extract_seq = fasta_seqs[y][extract_from..extract_to]
			while extract_seq.size < (extract_to+1) - extract_from
				extract_seq << "N"
			end
		end
	end
	# Write the extracted target region to a fasta file.
	extract_fasta_string = ">#{extract_id}\n"
	extract_fasta_string << "#{extract_seq}\n"
	# while extract_seq.size > 0
	# 	extract_fasta_string << "#{extract_seq.slice!(0..79)}\n"
	# end
	extract_fasta_file_name = "#{assembly_file_name.chomp(".fasta")}_#{extract_id}.fasta"
	extract_fasta_file = File.open(extract_fasta_file_name,"w")
	extract_fasta_file.write(extract_fasta_string)
	puts "Wrote file #{extract_fasta_file_name}."
end
