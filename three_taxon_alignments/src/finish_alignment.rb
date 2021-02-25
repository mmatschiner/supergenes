# m_matschiner Mon Jul 2 16:39:55 CEST 2018

# Get the command-line arguments.
refined_alignment_file_name = ARGV[0]
mapping_alignment_file_name = ARGV[1]
finished_alignment_file_name = ARGV[2]

# Read the refined input alignment file.
print "Reading the alignment file #{refined_alignment_file_name}..."
STDOUT.flush
alignment_file = File.open(refined_alignment_file_name)
alignment_lines = alignment_file.readlines
alignment_ids = []
alignment_seqs = []
alignment_lines.each do |l|
	stripped_line = l.strip
	if stripped_line[0] == ">"
		alignment_ids << stripped_line[1..-1]
		alignment_seqs << []
	elsif stripped_line != ""
		stripped_line.size.times do |pos|
			alignment_seqs.last << stripped_line[pos]
		end
	end
end
puts " done."
STDOUT.flush

# Remove alignment sites that are gaps in the first sequence (which is assumed to be the reference).
print "Removing alignment sites that are gaps in the reference..."
STDOUT.flush
tmp_alignment_seqs = []
alignment_seqs.size.times do |x|
	tmp_alignment_seqs << []
end
alignment_seqs[0].size.times do |pos|
	unless alignment_seqs[0][pos] == "-"
		alignment_seqs.size.times do |x|
			tmp_alignment_seqs[x] << alignment_seqs[x][pos]
		end
	end
end
puts " done."
STDOUT.flush

# Get the unique species IDs (the first 7 characters of all but the first sequence IDs)
species_ids = []
alignment_ids[1..-1].each do |i|
	species_ids << i[0..6]
end
uniq_species_ids = species_ids.uniq

# Set the window size.
window_size = 100

# Per unique species, make a strict consensus sequence in which all regions are set to missing if less than 100 consecutive positions are identical among all sequences of that species.
uniq_species_seqs = []
uniq_species_ids.each do |s|
        print "Generating strict consensus sequence for #{s}..."
        STDOUT.flush
	indices_for_this_uniq_species = []
	alignment_ids.size.times do |x|
		indices_for_this_uniq_species << x if alignment_ids[x][0..6] == s
	end
	uniq_species_seq = []
	tmp_alignment_seqs[0].size.times do |pos|
		uniq_species_seq << "-"
	end
 	pos = 0
	while pos < tmp_alignment_seqs[0].size-1 do
		if pos >= (window_size/2) and pos <= tmp_alignment_seqs[0].size - (window_size/2) - 1 
			windows_centered_at_this_pos_in_this_species = []
			indices_for_this_uniq_species.each do |x|
				windows_centered_at_this_pos_in_this_species << tmp_alignment_seqs[x][pos-(window_size/2)..pos+(window_size/2)]
			end
			if windows_centered_at_this_pos_in_this_species.uniq.size == 1
				(-(window_size/2)).upto(window_size/2) do |window_pos|
					uniq_species_seq[pos+window_pos] = windows_centered_at_this_pos_in_this_species[0][window_pos+(window_size/2)]
				end
			else
				uniq_species_seq[pos] = "-"
			end
		else
			uniq_species_seq[pos] = "-"
		end
                pos += 1
	end
	uniq_species_seqs << uniq_species_seq
        puts " done."
        STDOUT.flush
end

# Read the mapping-based input alignment file.
print "Reading the alignment file #{mapping_alignment_file_name}..."
STDOUT.flush
mapping_alignment_file = File.open(mapping_alignment_file_name)
mapping_alignment_lines = mapping_alignment_file.readlines
mapping_alignment_ids = []
mapping_alignment_seqs = []
mapping_alignment_lines.each do |l|
	stripped_line = l.strip
	if stripped_line[0] == ">"
		mapping_alignment_ids << stripped_line[1..-1]
		mapping_alignment_seqs << []
	elsif stripped_line != ""
		stripped_line.size.times do |pos|
			mapping_alignment_seqs.last << stripped_line[pos]
		end
	end
end
puts " done."
STDOUT.flush

# Compare the refined and mapping-based reference sequences and set any regions to missing if the two sequences differ.
unless mapping_alignment_ids.include?("gadMor2_mapping")
	puts "ERROR: Sequence ID \"gadMor2_mapping\" could not be found among the sequences in the mapping alignment!"
	exit 1
end
compare_seq = mapping_alignment_seqs[mapping_alignment_ids.index("gadMor2_mapping")]
pos_with_differences = []
tmp_alignment_seqs[0].size.times do |pos|
	compare_seq[pos] = "N" if compare_seq[pos] == nil
	if tmp_alignment_seqs[0][pos].upcase == "A" and ["C", "G", "T", "S", "Y", "K"].include?(compare_seq[pos].upcase)
		pos_with_differences << pos
	elsif tmp_alignment_seqs[0][pos].upcase == "C" and ["A", "G", "T", "R", "W", "K"].include?(compare_seq[pos].upcase)
		pos_with_differences << pos
	elsif tmp_alignment_seqs[0][pos].upcase == "G" and ["A", "C", "T", "M", "W", "Y"].include?(compare_seq[pos].upcase)
		pos_with_differences << pos
	elsif tmp_alignment_seqs[0][pos].upcase == "T" and ["A", "C", "G", "M", "R", "S"].include?(compare_seq[pos].upcase)
		pos_with_differences << pos
	end
end
pos_with_differences.each do |pos|
	-10.upto(10) {|pos_addition| tmp_alignment_seqs[0][pos+pos_addition] = "N" unless tmp_alignment_seqs[0][pos+pos_addition] == "-"}
end

# Compare the refined and mapping-based versions of the other sequences and set any regions to missing if the two sequences differ.
uniq_species_seqs.size.times do |x|
	unless mapping_alignment_ids.include?("#{uniq_species_ids[x]}_mapping")
		puts "ERROR: Sequence ID \"#{uniq_species_ids[x]}_mapping\" could not be found among the sequences in the mapping alignment!"
		exit 1
	end
	compare_seq = mapping_alignment_seqs[mapping_alignment_ids.index("#{uniq_species_ids[x]}_mapping")]
	pos_with_differences = []
	uniq_species_seqs[x].size.times do |pos|
		compare_seq[pos] = "N" if compare_seq[pos] == nil
		if uniq_species_seqs[x][pos].upcase == "A" and ["C", "G", "T", "S", "Y", "K"].include?(compare_seq[pos].upcase)
			pos_with_differences << pos
		elsif uniq_species_seqs[x][pos].upcase == "C" and ["A", "G", "T", "R", "W", "K"].include?(compare_seq[pos].upcase)
			pos_with_differences << pos
		elsif uniq_species_seqs[x][pos].upcase == "G" and ["A", "C", "T", "M", "W", "Y"].include?(compare_seq[pos].upcase)
			pos_with_differences << pos
		elsif uniq_species_seqs[x][pos].upcase == "T" and ["A", "C", "G", "M", "R", "S"].include?(compare_seq[pos].upcase)
			pos_with_differences << pos
		end
	end
	pos_with_differences.each do |pos|
		-10.upto(10) {|pos_addition| uniq_species_seqs[x][pos+pos_addition] = "N" unless uniq_species_seqs[x][pos+pos_addition] == "-"}
	end
end

# Prepare a new output alignment.
print "Preparing output alignment..."
STDOUT.flush
finished_alignment_string = ""
finished_alignment_string << ">reference\n"
finished_alignment_seq = ""
tmp_alignment_seqs[0].size.times do |pos|
	finished_alignment_seq << tmp_alignment_seqs[0][pos]
end
finished_alignment_string << "#{finished_alignment_seq}\n"
uniq_species_seqs.size.times do |x|
	finished_alignment_string << ">#{uniq_species_ids[x]}\n"
	finished_alignment_seq = ""
	uniq_species_seqs[x].size.times do |pos|
		finished_alignment_seq << uniq_species_seqs[x][pos]
	end
	finished_alignment_string << "#{finished_alignment_seq}\n"
end
puts " done."
STDOUT.flush

# Write the new output alignment.
finished_alignment_file = File.open(finished_alignment_file_name, "w")
finished_alignment_file.write(finished_alignment_string)

# Feedback.
puts "Wrote file #{finished_alignment_file_name}."
