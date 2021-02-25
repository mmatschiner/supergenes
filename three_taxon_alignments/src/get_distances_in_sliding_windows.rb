# m_matschiner Tue Jul 3 23:48:48 CEST 2018

# Get the command-line arguments.
alignment_file_name = ARGV[0]
window_size = ARGV[1].to_i

# Read the alignment.
alignment_file = File.open(alignment_file_name)
alignment_lines = alignment_file.readlines
ids = []
seqs = []
alignment_lines.each do |l|
	strip_line = l.strip
	if strip_line[0] == ">"
		ids << strip_line[1..-1]
		seqs << ""
	elsif strip_line != ""
		seqs.last << strip_line
	end
end

# Prepare the header of the output table.
header_line = "window_center\twindow_start\twindow_end"
1.upto(seqs.size-1) do |x|
	header_line << "\tdistance_between_#{ids[x]}_and_first"
end
puts header_line
STDOUT.flush

# Go over the alignment in sliding windows and calculate distances.
window_start = -1
window_end = -1
window_seqs = []
seqs.size.times { window_seqs << [] }
pos = 0
while pos < seqs[0].size
	nucleotides_at_this_pos = []
	seqs.each do |s|
		nucleotides_at_this_pos << s[pos].downcase
	end
	unless nucleotides_at_this_pos.include?("-") or nucleotides_at_this_pos.include?("n")
		window_start = pos if window_start == -1
		seqs.size.times { |x| window_seqs[x] << nucleotides_at_this_pos[x] }
	end
	if window_seqs[0].size == window_size
		# Set the window positions.
		window_end = pos
		window_center = (window_start+window_end)/2.0
		# Calculate the pairwise distances for this window.
		pairwise_dists = []
		(seqs.size-1).times { pairwise_dists << 0 }
		window_seqs[0].size.times do |window_pos|
			(seqs.size-1).times do |x|
				if window_seqs[x+1][window_pos] != window_seqs[0][window_pos]
					if ["a", "c", "g", "t"].include?(window_seqs[x+1][window_pos]) and ["a", "c", "g", "t"].include?(window_seqs[0][window_pos])
						pairwise_dists[x] += 1
					else
						puts "ERROR: Expected a, c, g, or t but found the following nucleotides: #{window_seqs[x+1][window_pos]} and #{window_seqs[0][window_pos]}"
						exit 1
					end
				end
			end
		end
		# Write output for this window.
		line = "#{'%.1f' % (window_center + 1)}\t#{window_start + 1}\t#{window_end + 1}"
		pairwise_dists.each do |d|
			line << "\t#{d}"
		end
		puts line
		STDOUT.flush
		# Reset the window sequences and start position.
		window_seqs = []
		seqs.size.times { window_seqs << [] }
		window_start = -1
	end
	pos += 1
end
