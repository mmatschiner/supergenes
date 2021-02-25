# m_matschiner Tue Jul 3 23:48:48 CEST 2018

# Get the command-line arguments.
alignment_file_name = ARGV[0]

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

# Convert sequences to arrays.
seqs.size.times do |x|
	seqs[x] = seqs[x].chars
end

# Go over the alignment in sliding windows and calculate distances.
puts "#{alignment_file_name}"
puts "species1\tspecies2\tlength\tn_missing\tn_different\tdistance"
0.upto(ids.size-2) do |x|
	(x+1).upto(ids.size-1) do |y|
		pos = 0
		differences_for_pair = 0
		missing_for_pair = 0
		while pos < seqs[0].size
			if ["-","n","N"].include?(seqs[x][pos]) or ["-","n","N"].include?(seqs[y][pos])
				missing_for_pair += 1
			elsif ["a","c","g","t"].include?(seqs[x][pos]) and ["a","c","g","t"].include?(seqs[y][pos])
				differences_for_pair += 1 if seqs[x][pos] != seqs[y][pos]
			else
				puts "ERROR: Expected a, c, g, or t but found the following nucleotides: #{seqs[x][pos]} and #{seqs[y][pos]}!"
				exit 1
			end
			pos += 1
		end
		puts "#{ids[x]}\t#{ids[y]}\t#{seqs[0].size}\t#{missing_for_pair}\t#{differences_for_pair}\t#{differences_for_pair/(seqs[0].size-missing_for_pair).to_f}"
	end
end

# Do the same for the inversion regions if the linkage group is 01, 02, 07, or 12.
start_pos = nil
end_pos = nil
if alignment_file_name.downcase.include?("lg01")
	start_pos = 9114741
	end_pos = 26192489
elsif alignment_file_name.downcase.include?("lg02")
	start_pos = 18489307
	end_pos = 24050282
elsif alignment_file_name.downcase.include?("lg07")
	start_pos = 13606502
	end_pos = 23016726
elsif alignment_file_name.downcase.include?("lg12")
	start_pos = 589105
	end_pos = 13631347
end
if start_pos != nil
	puts
	puts "#{alignment_file_name} - inversion"
	0.upto(ids.size-2) do |x|
		(x+1).upto(ids.size-1) do |y|
			pos = start_pos-1
			differences_for_pair = 0
			missing_for_pair = 0
			while pos < end_pos
				if ["-","n","N"].include?(seqs[x][pos]) or ["-","n","N"].include?(seqs[y][pos])
					missing_for_pair += 1
				elsif ["a","c","g","t"].include?(seqs[x][pos]) and ["a","c","g","t"].include?(seqs[y][pos])
					differences_for_pair += 1 if seqs[x][pos] != seqs[y][pos]
				else
					puts "ERROR: Expected a, c, g, or t but found the following nucleotides: #{seqs[x][pos]} and #{seqs[y][pos]}!"
					exit 1
				end
				pos += 1
			end
			puts "#{ids[x]}\t#{ids[y]}\t#{seqs[0].size}\t#{missing_for_pair}\t#{differences_for_pair}\t#{differences_for_pair/(end_pos-start_pos+1-missing_for_pair).to_f}"
		end
	end
end
