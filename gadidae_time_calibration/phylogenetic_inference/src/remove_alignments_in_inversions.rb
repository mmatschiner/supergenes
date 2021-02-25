# m_matschiner Tue Aug 7 15:25:27 CEST 2018

# Get the command-line arguments.
alignment_dir = ARGV[0]
exon_table_file_name = ARGV[1]
region_table_file_name = ARGV[2]
inversion_limits_table_file_name = ARGV[3]
location_table_file_name = ARGV[4]

# Get a list of alignment files.
alignment_file_names = []
Dir.entries(alignment_dir).each { |e| alignment_file_names << e if e[-4..-1] == ".nex" }

# Read the exon table.
exon_table_file = File.open(exon_table_file_name)
exon_table_lines = exon_table_file.readlines
exon_table_exon_ids = []
exon_table_gene_ids = []
exon_table_lines[1..-1].each do |l|
	line_ary = l.split
	exon_table_exon_ids << line_ary[0]
	exon_table_gene_ids << line_ary[1]
end

# Read the regions table.
regions_table_file = File.open(region_table_file_name)
regions_table_lines = regions_table_file.readlines
regions_table_exon_ids = []
regions_table_chromosome_ids = []
regions_table_positions = []
regions_table_lines.each do |l|
	line_ary = l.split
	regions_table_exon_ids << line_ary[0]
	regions_table_chromosome_ids << line_ary[1]
	regions_table_positions << (line_ary[3].to_i + line_ary[4].to_i)/2
end

# Read the inversion limit table.
inversion_limits_table_file = File.open(inversion_limits_table_file_name)
inversion_limits_table_lines = inversion_limits_table_file.readlines
inversion_limits_table_chromosome_ids = []
inversion_limits_table_froms = []
inversion_limits_table_tos = []
inversion_limits_table_lines.each do |l|
	unless l[0] == "#"
		line_ary = l.split
		inversion_limits_table_chromosome_ids << line_ary[0]
		inversion_limits_table_froms << line_ary[1].to_i
		inversion_limits_table_tos << line_ary[2].to_i
	end
end

# Get the chromosome and position for each remaining gene.
outstring = ""
alignment_file_names.each do |f|

	# Get the gene id.
	alignment_gene_id = f.chomp(".nex")
	
	# Make sure that the alignment gene id is found in the exon table.
	unless exon_table_gene_ids.include?(alignment_gene_id)
		puts "ERROR: Could not find gene id #{alignment_gene_id}!"
		exit(1)
	end
	
	# Get the chromosome id and the position of the gene.
	gene_chromosome_ids = []
	gene_positions = []
	exon_table_exon_ids.size.times do |x|
		if exon_table_gene_ids[x] == alignment_gene_id
			exon_id = exon_table_exon_ids[x]
			if regions_table_exon_ids.include?(exon_id)
				regions_table_exon_ids.size.times do |y|
					if regions_table_exon_ids[y] == exon_id
						gene_chromosome_ids << regions_table_chromosome_ids[y]
						gene_positions << regions_table_positions[y]
					end
				end
			end
		end
	end
	
	# Make sure that all exons are on the same chromosome.
	unless gene_chromosome_ids.uniq.size == 1
		puts "ERROR: Exons for gene #{alignment_gene_id} were found on more than one chromosome!"
		exit(1)
	end
	gene_chromosome_id = gene_chromosome_ids[0]
	
	# Calculate the mean gene position.
	sum_gene_positions = 0
	gene_positions.each { |pos| sum_gene_positions += pos }
	mean_gene_position = (sum_gene_positions/(gene_positions.size.to_f)).to_i

	# Find out if the gene lies within one of the inverted regions.
	gene_within_inversion = false
	if inversion_limits_table_chromosome_ids.include?(gene_chromosome_id)
		if mean_gene_position > inversion_limits_table_froms[inversion_limits_table_chromosome_ids.index(gene_chromosome_id)]
			if mean_gene_position < inversion_limits_table_tos[inversion_limits_table_chromosome_ids.index(gene_chromosome_id)]
				gene_within_inversion = true
			end
		end
	end
	outstring << "#{alignment_gene_id}\t#{gene_chromosome_id}\t#{mean_gene_position}\t#{gene_within_inversion}\n"

end

# Write the output file.
location_table_file = File.open(location_table_file_name, "w")
location_table_file.write(outstring)
