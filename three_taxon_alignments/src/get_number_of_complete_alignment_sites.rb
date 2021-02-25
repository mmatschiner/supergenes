# m_matschiner Tue Jul 3 00:39:26 CEST 2018

# Get the command-line arguments.
alignment_file_name = ARGV[0]

# Read the alignment file.
alignment_file = File.open(alignment_file_name)
alignment_lines = alignment_file.readlines
ids = []
seqs = []
alignment_lines.each do |l|
	if l[0] == ">"
		ids << l.strip[1..-1]
		seqs << ""
	elsif l.strip != ""
		seqs.last << l.strip
	end
end

# Get the number of complete alignment sites.
n_complete_sites = 0
seqs[0].size.times do |pos|
	site_complete = true
	seqs.each do |s|
		site_complete = false if s[pos] == "-"
	end
	n_complete_sites += 1 if site_complete
end

# Report the number of complete alignment sites.
puts "Complete alignment sites: #{n_complete_sites}/#{seqs[0].size}"