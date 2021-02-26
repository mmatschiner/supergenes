# m_matschiner Fri Mar 20 11:20:20 CET 2020

# Define a class for lines of the SVG graph.
class Line
	def initialize(x_start,x_end,y_start,y_end,color,stroke,opacity,stroke_dasharray)
		@x_start = x_start
		@x_end = x_end
		@y_start = y_start
		@y_end = y_end
		@color = color
		@stroke = stroke
		@opacity = opacity
		@stroke_dasharray = stroke_dasharray
	end
	def to_svg
		if @stroke_dasharray == "none"
			svg = "<line x1=\"#{@x_start.round(3)}\" y1=\"#{@y_start.round(3)}\" x2=\"#{@x_end.round(3)}\" y2=\"#{@y_end.round(3)}\" stroke=\"#{@color}\" stroke-width=\"#{@stroke}\" stroke-opacity=\"#{@opacity}\" stroke-linecap=\"round\" />"
		else
			svg = "<line x1=\"#{@x_start.round(3)}\" y1=\"#{@y_start.round(3)}\" x2=\"#{@x_end.round(3)}\" y2=\"#{@y_end.round(3)}\" stroke=\"#{@color}\" stroke-width=\"#{@stroke}\" stroke-opacity=\"#{@opacity}\" stroke-linecap=\"round\" stroke-dasharray=\"#{@stroke_dasharray}\" />"
		end
		svg
	end
end

# Define class for paths of the SVG graph.
class Path
	def initialize(x,y,fill,color,stroke,opacity,dashed=false)
		@x = [x]
		@y = [y]
		@fill = fill
		@color = color
		@stroke = stroke
		@opacity = opacity
		@dashed = dashed
	end
	def add_point(x,y)
		@x << x
		@y << y
	end
	def to_svg
		svg = "<path d=\"M #{@x[0]} #{@y[0]} "
		if @x.size > 1
			1.upto(@x.size-1) do |z|
				svg << "L #{@x[z]} #{@y[z]} "
			end
		end
		if @dashed
			svg << "\" fill=\"#{@fill}\" stroke=\"#{@color}\" stroke-width=\"#{@stroke}\" stroke-opacity=\"#{@opacity}\" stroke-linejoin=\"round\" stroke-linecap=\"round\" stroke-dasharray=\"1 2\"/>"
		else
			svg << "\" fill=\"#{@fill}\" stroke=\"#{@color}\" stroke-width=\"#{@stroke}\" stroke-opacity=\"#{@opacity}\" stroke-linejoin=\"round\" stroke-linecap=\"round\"/>"
		end
		svg
	end
end

# Define class for rectangles of the SVG graph.
class Rectangle
	def initialize(x,y,width,height,fill,color,stroke,opacity)
		@x = x
		@y = y
		@width = width
		@height = height
		@fill = fill
		@color = color
		@stroke = stroke
		@opacity = opacity
	end
	def to_svg
		svg = "<rect x=\"#{@x}\" y=\"#{@y}\" width=\"#{@width}\" height=\"#{@height}\" fill=\"#{@fill}\" stroke=\"#{@color}\" stroke-width=\"#{@stroke}\" fill-opacity=\"#{@opacity}\" stroke-linejoin=\"round\" />"
		svg
	end
end

# Define a class for text of the SVG graph.
class Text
	def initialize(x,y,font,font_size,weight,color,string,anchor,transform)
		@x = x
		@y = y
		@font = font
		@font_size = font_size
		@weight = weight
		@color = color
		@string = string
		@anchor = anchor
		@transform = transform
	end
	def to_svg
		if @transform == "none"
			if @anchor == "none"
				svg = "<text font-family=\"#{@font}\" font-size=\"#{@font_size}\" fill=\"#{@color}\" font-weight=\"#{@weight}\" x=\"#{@x}\" y=\"#{@y}\">#{@string}</text>"
			else
				svg = "<text text-anchor=\"#{@anchor}\" font-family=\"#{@font}\" font-size=\"#{@font_size}\" fill=\"#{@color}\" font-weight=\"#{@weight}\" x=\"#{@x}\" y=\"#{@y}\">#{@string}</text>"
			end
		else
			if @anchor == "none"
				svg = "<text font-family=\"#{@font}\" font-size=\"#{@font_size}\" fill=\"#{@color}\" font-weight=\"#{@weight}\" x=\"#{@x}\" y=\"#{@y}\" transform=\"#{@transform}\">#{@string}</text>"
			else
				svg = "<text text-anchor=\"#{@anchor}\" font-family=\"#{@font}\" font-size=\"#{@font_size}\" fill=\"#{@color}\" font-weight=\"#{@weight}\" x=\"#{@x}\" y=\"#{@y}\" transform=\"#{@transform}\">#{@string}</text>"
			end				
		end
		svg
	end
end


# Get the command line arguments.
tree_file_dir = ARGV[0]
lg = ARGV[1]
lg_length = ARGV[2].to_i
color_code_file_name = ARGV[3]
plot_file_name = ARGV[4]

# Get the tree file names.
tree_file_names = []
tree_file_dir_entries = Dir.entries(tree_file_dir)
tree_file_dir_entries.each do |e|
	if e[0..1] == "LG" and e[-4..-1] == ".tre"
		tree_file_names << e
	end
end
tree_file_names.sort!

# Read the color code file.
color_code_file = File.open(color_code_file_name)
color_code_lines = color_code_file.readlines
color_code_species = []
color_code_color = []
color_code_dashed = []
color_code_lines.each do |l|
	line_ary = l.split
	color_code_species << line_ary[0]
	color_code_color << line_ary[1]
	if line_ary[2] == "dashed"
		color_code_dashed << true
	else
		color_code_dashed << false
	end
end

# Read the tree files.
tree_window_lgs = []
tree_window_starts = []
tree_window_ends = []
tree_strings = []
trees_strings_ary = []
tree_file_names.each do |n|
	tree_window_lg = n.split("/").last.split("_")[0]
	if tree_window_lg == lg
		tree_window_starts << n.split("/").last.split("_")[1].to_f
		tree_window_ends << n.split("/").last.split("_")[2].split(".")[0].to_f
		tree_file = File.open("#{tree_file_dir}/#{n}")
		trees_file = File.open("#{tree_file_dir}/#{n}es")
		tree_strings << tree_file.read
		tree_file.close
		trees_strings_ary << trees_file.readlines
		trees_file.close
	end
end
puts "Read #{tree_strings.size} MCC tree files."
puts "Read #{trees_strings_ary.size} posterior distribution tree files."
trees_strings_ary.each do |a|
	unless a.size == trees_strings_ary[0].size
		puts "ERROR: Not all posterior distribution tree files have #{trees_strings_ary[0].size} trees!"
		exit 1
	end
end
puts "All posterior distribution tree files have #{trees_strings_ary[0].size} trees."

# Get divergence times between neighboring species in the trees.
tree_window_centers = []
tree_species_arys = []
trees_species_arys_arys = []
tree_cumulative_distance_arys = []
trees_cumulative_distance_arys_arys = []
tree_strings.size.times do |x|
	# Collect distances between neighboring species in the mcc tree.
	tree_string = tree_strings[x].gsub("(","").gsub(")","")
	tree_species = []
	tree_distances = []
	tree_ary = tree_string.split(",")
	tree_ary.each do |i|
		tip_ary = i.split(":")
		tree_species << tip_ary[0]
		tip_times = tip_ary[1..-1].map { |n| n.to_f }
		tree_distances << tip_times.sum * 2
	end
	tree_window_centers << (tree_window_starts[x]+tree_window_ends[x])/2.0
	tree_species_arys << tree_species
	tree_cumulative_distance_ary = [0]
	(tree_distances.size-1).times do |y|
		tree_cumulative_distance_ary << tree_distances[0..y].sum
	end
	tree_cumulative_distance_arys << tree_cumulative_distance_ary

	# Do the same for each set of posterior trees per window.
	trees_species_arys = []
	trees_cumulative_distance_arys = []
	trees_strings_ary[x].each do |ts|
		trees_string = ts.gsub("(","").gsub(")","")
		trees_species = []
		trees_distances = []
		trees_ary = trees_string.split(",")
		trees_ary.each do |i|
			tip_ary = i.split(":")
			trees_species << tip_ary[0]
			tip_times = tip_ary[1..-1].map { |n| n.to_f }
			trees_distances << tip_times.sum * 2
		end
		trees_species_arys << trees_species
		trees_cumulative_distance_ary = [0]
		(trees_distances.size-1).times do |y|
			trees_cumulative_distance_ary << trees_distances[0..y].sum
		end
		trees_cumulative_distance_arys << trees_cumulative_distance_ary
	end
	trees_species_arys_arys << trees_species_arys
	trees_cumulative_distance_arys_arys << trees_cumulative_distance_arys
end

# Identify the maximum cumulative tree distance.
max_cumulative_distance = 0
tree_cumulative_distance_arys.each do |a|
	max_cumulative_distance = a.last if a.last > max_cumulative_distance
	if max_cumulative_distance > 4
		puts "WARNING: The maximum cumulative distance is greater than 4 (#{max_cumulative_distance})!"
	end
end
max_cumulative_distance = 4.0

# Generate the svg paths.
svg_width = 180
svg_height = 100
left_margin = 5
right_margin = 5
top_margin = 5
bottom_margin = 11
dim_x = svg_width - left_margin - right_margin
dim_y = svg_height - top_margin - bottom_margin
paths = []
color_code_species.each do |s|
	color = color_code_color[color_code_species.index(s)]
	dashed = color_code_dashed[color_code_species.index(s)]

	# Write the path for the mcc tree.
	path = nil
	tree_window_centers.size.times do |x|
		tree_species_ary = tree_species_arys[x]
		tree_cumulative_distance_ary = tree_cumulative_distance_arys[x]
		# Identify the largest difference between cumulative distances (= the largest distance).
		largest_difference = tree_cumulative_distance_ary[0]
		root_index = 0
		(tree_cumulative_distance_ary.size-1).times do |z|
			difference = tree_cumulative_distance_ary[z+1]-tree_cumulative_distance_ary[z]
			if difference > largest_difference
				largest_difference = difference
				root_index = z
			end
		end
		cumulative_root_distance = tree_cumulative_distance_ary[root_index] + largest_difference / 2.0
		arys_species_index = tree_species_ary.index(s)
		cumulative_distance = tree_cumulative_distance_ary[arys_species_index]
		x_coord = left_margin + dim_x * tree_window_centers[x]/lg_length.to_f
		y_coord = top_margin + 0.5 * dim_y + dim_y * (cumulative_distance-cumulative_root_distance) / max_cumulative_distance
		if path == nil
			path = Path.new(x_coord,y_coord,"none",color,1,1,dashed)
		else
			path.add_point(x_coord,y_coord)
		end
	end
	paths << path

	# Do the same for each set of posterior trees per window.
	trees_species_arys_arys[0].size.times do |q| # 119, should be 100
		path = nil
		tree_window_centers.size.times do |x| # 119
			tree_species_ary = trees_species_arys_arys[x][q]
			if trees_cumulative_distance_arys_arys[x][q] == nil
				puts "ERROR: trees_cumulative_distance_arys_arys[x][q] is nil!"
				exit
			end
			tree_cumulative_distance_ary = trees_cumulative_distance_arys_arys[x][q]
			# Identify the largest difference between cumulative distances (= the largest distance).
			largest_difference = tree_cumulative_distance_ary[0]
			root_index = 0
			(tree_cumulative_distance_ary.size-1).times do |z|
				difference = tree_cumulative_distance_ary[z+1]-tree_cumulative_distance_ary[z]
				if difference > largest_difference
					largest_difference = difference
					root_index = z
				end
			end
			cumulative_root_distance = tree_cumulative_distance_ary[root_index] + largest_difference / 2.0
			arys_species_index = tree_species_ary.index(s)
			cumulative_distance = tree_cumulative_distance_ary[arys_species_index]
			x_coord = left_margin + dim_x * tree_window_centers[x]/lg_length.to_f
			y_coord = top_margin + 0.5 * dim_y + dim_y * (cumulative_distance-cumulative_root_distance) / max_cumulative_distance
			if path == nil
				path = Path.new(x_coord,y_coord,"none",color,0.1,0.1,dashed)
			else
				path.add_point(x_coord,y_coord)
			end
		end
		paths << path
	end
end

# Add the scale bar.
lines = []
texts = []
font = "Helvetica"
font_size6 = 2.1169354839 # 6pt
font_size8 = 2.8225806452 # 8pt
bar_time = 0.5
bar_length = dim_y * bar_time / max_cumulative_distance
x_coord = left_margin
y1_coord = top_margin + 5
y2_coord = top_margin + 5 + bar_length
lines << Line.new(x_coord,x_coord,y1_coord,y2_coord,"black",1,1,"none")
x_coord = left_margin + 2.5
y_coord = top_margin + 5 + 0.5 * bar_length
texts << Text.new(x_coord,y_coord,font,font_size6,"normal","black","#{bar_time} myr","start","none")

# Add the root line.
bg_lines = []
x1_coord = left_margin
x2_coord = left_margin + dim_x
y_coord = top_margin + dim_y / 2.0
bg_lines << Line.new(x1_coord,x2_coord,y_coord,y_coord,"#839496",1,1,"0 2")
[-1,-0.5,0.5,1].each do |x|
	y_coord = top_margin.to_f + (dim_y.to_f / 2.0) - (dim_y.to_f * x.to_f / max_cumulative_distance.to_f)
	bg_lines << Line.new(x1_coord,x2_coord,y_coord,y_coord,"#839496",0.5,0.5,"none")
end

# Add the axis line.
x1_coord = left_margin
x2_coord = left_margin + dim_x
y_coord = top_margin + dim_y
bg_lines << Line.new(x1_coord,x2_coord,y_coord,y_coord,"black",1,1,"none")
x_coord = left_margin + dim_x / 2.0
y_coord = top_margin + dim_y + 9
texts << Text.new(x_coord,y_coord,font,font_size8,"normal","black","#{lg} position (Mb)","middle","none")

# Add ticks on axis.
increment = 1000000
pos = 0
while pos < lg_length
	x_coord = left_margin + dim_x * pos / lg_length
	y1_coord = top_margin + dim_y
	y2_coord = top_margin + dim_y + 1
	bg_lines << Line.new(x_coord,x_coord,y1_coord,y2_coord,"black",1,1,"none")
	if (pos/2000000)*2000000 == pos or pos == 0
		y_coord = y2_coord + 3
		texts << Text.new(x_coord,y_coord,font,font_size6,"normal","black","#{pos/1000000}","middle","none")
	end
	pos += increment
end

# Make svg frame
frame = Rectangle.new(left_margin,top_margin,dim_x,dim_y,"none","black",1,1)

# Prepare the svg string.
svg_string = ""
svg_string << "<?xml version=\"1.0\" standalone=\"no\"?>\n"
svg_string << "<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\" \"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">\n"
svg_string << "<svg width=\"#{svg_width}mm\" height=\"#{svg_height}mm\" viewBox=\"0 0 #{svg_width} #{svg_height}\" xmlns=\"http://www.w3.org/2000/svg\" version=\"1.1\">\n"
bg_lines.each {|l| svg_string << "    #{l.to_svg}\n"}
paths.each {|p| svg_string << "    #{p.to_svg}\n"}
texts.each {|t| svg_string << "    #{t.to_svg}\n"}
lines.each {|l| svg_string << "    #{l.to_svg}\n"}
svg_string << "</svg>\n"

# Write the svg string to file.
svg_file = File.new(plot_file_name,"w")
svg_file.write(svg_string)

# Feedback.
puts "Wrote file #{plot_file_name}."
