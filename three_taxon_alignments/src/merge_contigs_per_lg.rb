# m_matschiner Thu Jun 15 12:28:19 CEST 2017

# Get the command-line arguments.
input_file_name = ARGV[0]
output_file_name = ARGV[1]
add_ns_str = ARGV[2]

# Interpret the command-line arguments.
lg = input_file_name.split("/").last.chomp(".fasta")
add_ns = false
add_ns = true if add_ns_str == "true"

# Read the input file.
input_lines = File.open(input_file_name).readlines

# Generate a concatenated sequence.
seq = ""
input_lines.each do |l|
  if l[0] == ">"
    seq << "NNNNNNNNNN" if add_ns and seq != ""
  else
    seq << l.strip
  end
end

# Write the output file.
outfile = File.open(output_file_name,"w")
outfile.write(">#{lg}\n")
x = 0
while x < seq.size do
  outfile.write("#{seq[x..x+79].upcase}\n")
  x += 80
end
