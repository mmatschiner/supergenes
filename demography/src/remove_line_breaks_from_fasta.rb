# m_matschiner Mon Aug 31 23:44:25 CEST 2020

# Get the command-line arguments.
in_fasta_file_name = ARGV[0]
out_fasta_file_name = ARGV[1]

# Read the input file.
ids = []
seqs = []
in_fasta_file = File.open(in_fasta_file_name)
in_fasta_file.each do |l|
	if l[0] == ">"
		ids << l.strip[1..-1]
		seqs << ""
	elsif l.strip != ""
		seqs.last << l.strip
	end
end

# Write the output file.
out_fasta_file = File.open(out_fasta_file_name, "w")
ids.size.times do |x|
	out_fasta_file.write(">#{ids[x]}\n")
	out_fasta_file.write("#{seqs[x]}\n")
end