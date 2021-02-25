# m_matschiner Tue Feb 27 20:28:14 CET 2018

# Get the name of the input alignment file.
input_file_name = ARGV[0]
output_file_name = ARGV[1]

# Read the input alignment file.
input_file = File.open(input_file_name)
lines = input_file.readlines

# Initiate two sequence strings.
query_seq = ""
cod_seq = ""

# Read the two sequences.
lines[3..-1].each do |l|
	if l[0..5] == "Query:"
		query_seq << l.split[2]
	elsif l[0..5] == "Sbjct:"
		cod_seq << l.split[2]
	end
end

# Get the start and end positions of the first window.
window_size = 50000
from = 0
to = from
n_nongap = 0
until n_nongap == window_size-1
	to += 1
	n_nongap += 1 unless cod_seq[to] == "-"
end
from_without_gaps = 0
to_without_gaps = window_size - 1

# Use a shifting window to calculate gap numbers and difference.
outstring = "from\tto\tfrom_without_gaps\tto_without_gaps\tcod_gaps\tquery_gaps\tdistance\n"
until to > cod_seq.size
	query_fragment = query_seq[from..to]
	cod_fragment = cod_seq[from..to]
	fragment_difference = 0
	fragment_cod_gaps = 0
	fragment_query_gaps = 0
	query_fragment.size.times do |x|
		if query_fragment[x] == "-"
			fragment_query_gaps += 1
		elsif cod_fragment[x] == "-"
			fragment_cod_gaps += 1
		elsif cod_fragment[x] != query_fragment[x]
			fragment_difference += 1 if ["A","C","G","T"].include?(cod_fragment[x].upcase) and ["A","C","G","T"].include?(query_fragment[x].upcase)
		end
	end
	outstring << "#{from}\t#{to}\t#{from_without_gaps}\t#{to_without_gaps}\t#{fragment_cod_gaps}\t#{fragment_query_gaps}\t#{fragment_difference}\n"

	# Get the start and end positions of the next window.
	from = to + 1
	to = from
	n_nongap = 0
	until n_nongap == window_size-1
		to += 1
		n_nongap += 1 unless cod_seq[to] == "-"
	end
	from_without_gaps = to_without_gaps + 1
	to_without_gaps = from_without_gaps + window_size - 1
end

# Write the output file.
output_file = File.open(output_file_name,"w")
output_file.write(outstring)
