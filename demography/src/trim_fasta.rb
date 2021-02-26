# m_matschiner Mon Aug 31 11:38:51 CEST 2020

# Get the command-line arguments.
in_fasta_file_name = ARGV[0]
from_pos = ARGV[1].to_i - 1
to_pos = ARGV[2].to_i - 1
out_fasta_file_name = ARGV[3]

# Read the input fasta file.
ids = []
seqs = []
in_fasta_file = File.open(in_fasta_file_name)
in_fasta_lines = in_fasta_file.readlines
in_fasta_lines.each do |l|
	if l[0] == ">"
		ids << l[1..-1].strip
		seqs << ""
	elsif l.strip != ""
		seqs.last << l.strip
	end
end

# Make sure only a single sequence is in the fasta file.
if ids.size > 1 or seqs.size > 1
	puts "ERROR: Fasta file includes more than one sequence!"
	exit 1
end

# Write a trimmed sequence to the output file.
out_fasta_string = ">#{ids[0]}\n"
out_fasta_string << seqs[0][from_pos..to_pos]
out_fasta_file = File.open(out_fasta_file_name, "w")
out_fasta_file.write(out_fasta_string)