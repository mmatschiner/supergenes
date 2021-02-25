# m_matschiner Wed Jul 18 09:47:33 CEST 2018

# Get the command-line arguments.
alignment_file_name = ARGV[0]
chromosome_id = ARGV[1]
mask_file_name = ARGV[2]

# Read the alignment file.
alignment_file = File.open(alignment_file_name)
alignment_lines = alignment_file.readlines
ids = []
seqs = []
alignment_lines.each do |l|
	stripped_line = l.strip
	if stripped_line[0] == ">"
		ids << stripped_line[1..-1]
		seqs << []
	elsif stripped_line != ""
		stripped_line.size.times do |pos|
			seqs.last << stripped_line[pos]
		end
	end
end

# Get the unreliable sites (those with missing data).
unreliable_sites = []
seqs[0].size.times do |pos|
	pos_reliable = true
	seqs.each do |s|
		if s[pos] == "-" or s[pos] == "n" or s[pos] == "N"
			pos_reliable = false
			break
		end
	end
	unreliable_sites << pos unless pos_reliable
end

# Prepare the mask string.
mask_string = ""
unreliable_sites.size.times do |x|
	if x == 0
		mask_string << "#{chromosome_id}\t#{unreliable_sites[x] + 1}\t"
	elsif unreliable_sites[x] != unreliable_sites[x-1] + 1
		mask_string << "#{unreliable_sites[x-1] + 1}\n"
		mask_string << "#{chromosome_id}\t#{unreliable_sites[x] + 1}\t"
	end
	if x == unreliable_sites.size-1
		mask_string << "#{unreliable_sites[x] + 1}\n"
	end
end

# Write the mask to file.
mask_file = File.open(mask_file_name, "w")
mask_file.write(mask_string)
puts "Wrote file #{mask_file_name}."