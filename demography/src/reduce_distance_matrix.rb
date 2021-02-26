# m_matschiner Fri Jun 29 16:01:54 CEST 2018

# Add methods to enumerables.
module Enumerable
	def sum
		self.inject(0){|accum, i| accum + i }
	end
	def mean
		if self.length == 0
			nil
		else
			self.sum/self.length.to_f
		end
	end
	def sample_variance
		if self.length == 0
			nil
		else
			m = self.mean
			sum = self.inject(0){|accum, i| accum +(i-m)**2 }
			sum/(self.length - 1).to_f
		end
	end
	def standard_deviation
		if self.length == 0
			nil
		else
			return Math.sqrt(self.sample_variance)
		end
	end
	def median
		if self.length == 0
			nil
		else
			sorted_array = self.sort
			if self.size.modulo(2) == 1 
				sorted_array[self.size/2]
			else
				(sorted_array[(self.size/2)-1]+sorted_array[self.size/2])/2.0
			end
		end
	end
end

# Get the command-line arguments.
distance_matrix_file_name = ARGV[0]
sample_table_file_name = ARGV[1]
reduced_distance_matrix_file_name = ARGV[2]
exclude_sample_id_string = ARGV[3]
unless ARGV[3] == nil
	exclude_sample_ids = exclude_sample_id_string.split(",")
end

# Read the distance matrix.
distance_matrix_file = File.open(distance_matrix_file_name)
distance_matrix_lines = distance_matrix_file.readlines
distance_matrix_lines.size.times do |x|
	if distance_matrix_lines[x][0] == "#"
		distance_matrix_lines[x] = nil
	end
end
distance_matrix_lines.compact!
distance_matrix_sample_ids = distance_matrix_lines[0].split
distance_matrix = []
row_ary = []
distance_matrix_lines[1].split.size.times do |x|
	row_ary << 0.0
end
distance_matrix << row_ary
distance_matrix_lines[1..-1].each do |l|
	line_ary = l.split
	row_ary = []
	line_ary.each do |i|
		row_ary << i.to_f
	end
	distance_matrix << row_ary
end

# Make the distance matrix symmetric.
distance_matrix.size.times do |x|
	distance_matrix[x].size.times do |y|
		distance_matrix[x][y] = distance_matrix[y][x] if distance_matrix[x][y] == 0 and distance_matrix[y][x] > 0
		distance_matrix[y][x] = distance_matrix[x][y] if distance_matrix[y][x] == 0 and distance_matrix[x][y] > 0
	end
end

# Read the sample table.
sample_table_file = File.open(sample_table_file_name)
sample_table_lines = sample_table_file.readlines
sample_table_sample_ids = []
sample_table_species_ids = []
sample_table_lines.each do |l|
	line_ary = l.split
	if distance_matrix_sample_ids.include?(line_ary[0])
		sample_table_sample_ids << line_ary[0]
		sample_table_species_ids << line_ary[1]
	end
end

# Remove samples from the table that are to be excluded.
unless exclude_sample_ids == nil
  sample_table_sample_ids.size.times do |x|
    if exclude_sample_ids.include?(sample_table_sample_ids[x])
      sample_table_sample_ids[x] = nil
      sample_table_species_ids[x] = nil
    end
  end
  sample_table_sample_ids.compact!
  sample_table_species_ids.compact!
end

# Get a list of unique species ids from the table.
uniq_sample_table_species_ids = sample_table_species_ids.uniq.sort

# Prepare a new matrix with species comparisons.
reduced_distance_matrix = []
uniq_sample_table_species_ids.size.times do |x|
	reduced_distance_matrix << []
	uniq_sample_table_species_ids.size.times do |y|
		reduced_distance_matrix[x] << 0
	end
end

# Fill the matrix with species comparisons.
relative_standard_deviations = []
uniq_sample_table_species_ids.size.times do |x|
	uniq_sample_table_species_ids.size.times do |y|
		unless x == y
			species1_sample_ids = []
			species2_sample_ids = []
			sample_table_sample_ids.size.times do |z|
				if sample_table_species_ids[z] == uniq_sample_table_species_ids[x]
					species1_sample_ids << sample_table_sample_ids[z]
				elsif sample_table_species_ids[z] == uniq_sample_table_species_ids[y]
					species2_sample_ids << sample_table_sample_ids[z]
				end
			end
			distances_between_species1_and_species2 = []
			species1_sample_ids.each do |s1|
				species2_sample_ids.each do |s2|
					if distance_matrix_sample_ids.index(s1) == nil
						puts "ERROR: Species 1 sample #{s1} cannot be found!"
						exit 1
					elsif distance_matrix_sample_ids.index(s2) == nil
						puts "ERROR: Species 1 sample #{s2} cannot be found!"
						exit 1
					else
						if distance_matrix[distance_matrix_sample_ids.index(s1)][distance_matrix_sample_ids.index(s2)].to_s.strip == ""
							puts "ERROR: No distance found between samples #{s1} and #{s2}!"
							puts distance_matrix_sample_ids.index(s1)
							puts distance_matrix_sample_ids.index(s2)
							exit 1
						end
						distances_between_species1_and_species2 << distance_matrix[distance_matrix_sample_ids.index(s1)][distance_matrix_sample_ids.index(s2)]
					end
				end
			end

			distance = distances_between_species1_and_species2.mean
			relative_standard_deviation = distances_between_species1_and_species2.standard_deviation/distance.to_f
			relative_standard_deviations << relative_standard_deviation unless relative_standard_deviation.nan? # If single samples are compared, standard deviation is NaN.
			reduced_distance_matrix[x][y] = distance
			reduced_distance_matrix[y][x] = distance
		end
	end
end

# Prepare the output.
outstring = "# Median relative standard deviation of distances measured among samples of the same species: #{relative_standard_deviations.median}\n"
uniq_sample_table_species_ids.each { |i| outstring << "\t#{i}" }
outstring << "\n"
reduced_distance_matrix.size.times do |x|
	outstring << "#{uniq_sample_table_species_ids[x]}"
	reduced_distance_matrix.size.times do |y|
		outstring << "\t#{'%.2f' % reduced_distance_matrix[x][y]}"
	end
	outstring << "\n"
end

# Write the output file.
reduced_distance_matrix_file = File.open(reduced_distance_matrix_file_name,"w")
reduced_distance_matrix_file.write(outstring)

# Feedback.
puts "Wrote combined distance matrix to file #{reduced_distance_matrix_file_name} (median relative standard deviation #{relative_standard_deviations.median})."
