# m_matschiner Wed Jan 18 12:31:30 CET 2017

# Get the command line arguments.
fasta_input_file_name = ARGV[0]
features_output_file_name = ARGV[1]
window_size = ARGV[2].to_i
missing_threshold = ARGV[3].to_f
repetitiveness_threshold = ARGV[4].to_f

# Read the fasta input file.
fasta_input_file = File.open(fasta_input_file_name)
fasta_input_lines = fasta_input_file.readlines
seq = fasta_input_lines[1].strip.upcase

# Prepare an output string with data on missing data and low-complexity regions.
features_output_string = ""
from = 0
to = window_size-1
while to <= seq.size-1
	seq_window = seq[from..to]
	# Classify missingness.
	if seq_window.count("N") > missing_threshold * window_size
		missing = true
	else
		missing = false
	end
	# Classify repetitiveness.
	upper_end = window_size-1
	n_diagonals = upper_end
	1.upto(upper_end) do |x|
		(x+1).upto(upper_end) do |y|
			n_diagonals += 2 if seq_window[x] == seq_window[y] and seq_window[x-1] == seq_window[y-1]
		end
	end
	repetitiveness = n_diagonals/((upper_end)**2).to_f
	if repetitiveness > repetitiveness_threshold
		repetitive = true
	else
		repetitive = false
	end
	features_output_string << "#{from}\t#{to}\t#{missing}\t#{repetitive}\n"
	from += window_size
	to += window_size
end

# Write the feature output file.
features_output_file = File.open(features_output_file_name,"w")
features_output_file.write(features_output_string)