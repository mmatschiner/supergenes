# m_matschiner Wed Jun 14 13:56:06 CEST 2017

input_file = File.open(ARGV[0])
lines = input_file.readlines
id = lines[0].strip[1..-1]
seq = lines[1].strip
rev_comp_seq = ""
(seq.size-1).downto(0) do |x|
	if seq[x] == "A"
		rev_comp_seq << "T"
	elsif seq[x] == "C"
		rev_comp_seq << "G"
	elsif seq[x] == "G"
		rev_comp_seq << "C"
	elsif seq[x] == "T"
		rev_comp_seq << "A"
	elsif seq[x] == "N"
		rev_comp_seq << "N"
	elsif seq[x] == "a"
		rev_comp_seq << "t"
	elsif seq[x] == "c"
		rev_comp_seq << "g"
	elsif seq[x] == "g"
		rev_comp_seq << "c"
	elsif seq[x] == "t"
		rev_comp_seq << "a"
	elsif seq[x] == "n"
		rev_comp_seq << "n"
	else
		raise "ERROR: Found unexpected nucleotide: #{seq[x]}"
	end
end

puts ">#{id}"
puts "#{rev_comp_seq}"