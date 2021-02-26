# m_matschiner Tue Sep 1 00:36:14 CEST 2020

# Get the command-line arguments.
bed_file_name = ARGV[0]
lg = ARGV[1]
region_from = ARGV[2].to_i-1
region_to = ARGV[3].to_i-1

# Read the bed file.
length = region_to - region_from + 1
n_masked = 0
bed_file = File.open(bed_file_name)
bed_file.each do |l|
	lineary = l.split
	if lineary[0] == lg
		fragment_from = lineary[1].to_i
		fragment_to = lineary[2].to_i-2
		length_masked = [fragment_to,region_to].min - [fragment_from,region_from].max
		n_masked += length_masked if length_masked > 0
	end
end

puts "#{lg}\t#{region_from+1}\t#{region_to+1}\t#{length}\t#{n_masked}\t#{n_masked/length.to_f}"