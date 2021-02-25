# m_matschiner Tue Jun 26 17:21:34 CEST 2018

# Get the command-line argument.
input_alignment_file_name = ARGV[0]
output_alignment_file_name = ARGV[1]
reference_file_name = ARGV[2]
chromosome_id = ARGV[3]
output_sequence_id = ARGV[4]

# Read the input file.
input_alignment_file = File.open(input_alignment_file_name)
input_alignment_lines = input_alignment_file.readlines

# Read the reference file.
reference_file = File.open(reference_file_name)
reference_lines = reference_file.readlines
if reference_lines[0][0] != ">"
	puts "ERROR: Expected reference sequence in fasta format, but could not read format!"
	exit 1
end
reference_seq = reference_lines[1].strip

# Open the output file.
output_alignment_file = File.open(output_alignment_file_name,"w")

# Initialize sequence strings for the subject and the query.
subject_seq = ""
query_seq = ""

# Determine the input-file format.
if input_alignment_lines[0][0..4] == "Query"

	# Analyse the input file.
	1.upto(input_alignment_lines.size-1) do |x|
		query_line = input_alignment_lines[x]
		if query_line[0..4] == "Query"

			# Read the current query and subject lines.
			query_line_ary = query_line.strip.split
			fragment_length = query_line_ary[2].size
			already_found_short_fragment = false
			already_found_other_end_coordinate = false
			if query_line_ary.size != 4
				puts "ERROR: Expected line to contain four elements, but found this line: #{query_line.strip}!"
				exit 1
			elsif fragment_length != 60
				if already_found_short_fragment # One short fragment is allowed, this could be the very last fragment in the alignment.
					puts "ERROR: Expected sequence fragment to have 60 bp, but found this sequence fragment: #{query_line_ary[2]}!"
					exit 1
				end
				already_found_short_fragment = true
			elsif input_alignment_lines[x+2][0..4] != "Sbjct"
				puts "ERROR: Expected line with \'Query\' keyword to be followed by a line with the \'Sbjct\' keyword (after another line in between), but found this line instead: #{input_alignment_lines[x+2]}!"
				exit 1
			end
			query_fragment = query_line_ary[2]
			subject_line = input_alignment_lines[x+2]
			subject_line_ary = subject_line.strip.split
			if subject_line_ary.size != 4
				puts "ERROR: Expected line to contain four elements, but found this line: #{subject_line.strip}!"
				exit 1
			elsif subject_line_ary[2].size != fragment_length
				puts "ERROR: Expected sequence fragment to have #{fragment_length} bp, but found this sequence fragment: #{subject_line_ary[2]}!"
				exit 1
			end
			subject_fragment = subject_line_ary[2]

			# Add the subject and query sequence fragments.
			subject_seq << subject_fragment
			query_seq << query_fragment

		end
	end

elsif input_alignment_lines[0][0..4] == "##maf"

	# Initialize an array for the sequence and the sequence length.
	gaps_added_in_subject = 0
	subject_seq_as_array = []
	query_seq_as_array = []
	old_start_coordinate_for_fragment = -1
	old_end_coordinate_for_fragment = -1

	# Analyse the input file.
	1.upto(input_alignment_lines.size-1) do |x|
		line = input_alignment_lines[x]
		if line[0..1] == "a "
			subject_line = input_alignment_lines[x+1].strip
			query_line = input_alignment_lines[x+2].strip
			if query_line[0..1] != "s "
				puts "ERROR: Expected alignment line to begin with \'s \' but found #{query_line}"
				exit 1
			elsif subject_line[0..1] != "s "
				puts "ERROR: Expected alignment line to begin with \'s \' but found #{subject_line}"
				exit 1
			end
			subject_line_ary = subject_line.split
			query_line_ary = query_line.split
			if subject_line_ary.size != 7
				puts "ERROR: Expected line to contain seven elements but found this line: #{subject_line_ary}!"
				exit 1
			elsif query_line_ary.size != 7
				puts "ERROR: Expected line to contain seven elements but found this line: #{subject_line_ary}!"
				exit 1
			end
			# Read this alignment if the chromosome id matches the one specified for this analysis.
			if chromosome_id.downcase == subject_line_ary[1][1..-1].downcase

				# Extract info from the query and subject lines.
				start_coordinate_for_fragment = subject_line_ary[2].to_i
				if start_coordinate_for_fragment < old_start_coordinate_for_fragment
					puts "ERROR: Expected increasing start coordinates of fragments but found #{start_coordinate_for_fragment} following #{old_start_coordinate_for_fragment}!"
					exit 1
				elsif start_coordinate_for_fragment < old_end_coordinate_for_fragment
					puts "ERROR: Expected non-overlapping fragments but found a fragment starting at #{start_coordinate_for_fragment} following another ending at #{old_end_coordinate_for_fragment}!"
					exit 1
				end
				end_coordinate_for_fragment = start_coordinate_for_fragment + subject_line_ary[3].to_i - 1
				subject_fragment = subject_line_ary[6]
				query_fragment = query_line_ary[6]

				# Make sure the subject and the query have the same length.
				if subject_fragment.size != query_fragment.size
					puts "ERROR: Expected subject and query sequences to have the same length but found different lengths!"
					exit 1
				end
				
				# Add to the chromosome-length alignment.
				(old_end_coordinate_for_fragment + 1).upto(start_coordinate_for_fragment - 1) do |pos|
					subject_seq_as_array << reference_seq[pos]
					query_seq_as_array << "-"
				end
				subject_fragment.size.times do |pos|
					subject_seq_as_array << subject_fragment[pos]
					query_seq_as_array << query_fragment[pos]
				end

				# Memorize the start and end coordinates.
				old_start_coordinate_for_fragment = start_coordinate_for_fragment
				old_end_coordinate_for_fragment = end_coordinate_for_fragment
			end
		end
	end

	# Make sure some sequences were found for this chromosome.
	if query_seq_as_array == []
		puts "ERROR: No sequences found!"
		exit 1
	end

	# Convert the sequence array into a string.
	subject_seq_as_array.each {|i| subject_seq << i}
	query_seq_as_array.each {|i| query_seq << i}

else

	# Quit if the file format could not be read.	
	puts "ERROR: Unexpected file format!"
	exit 1

end

# Prepare the output sequence string.
output_alignment_string = ">reference\n"
output_alignment_string << "#{subject_seq}\n"
output_alignment_string << ">#{output_sequence_id}\n"
output_alignment_string << "#{query_seq}\n"
output_alignment_file.write(output_alignment_string)

# Feedback.
puts "Wrote sequence to file #{output_alignment_file_name}"
