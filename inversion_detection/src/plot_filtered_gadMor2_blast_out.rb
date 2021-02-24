# m_matschiner Wed Jan 18 12:10:31 CET 2017

# Define a class for SVG graphs.
class SVG
	attr_reader :width, :height
	def initialize(width, height)
		@width = width
		@height = height
		@elements = []
	end
	def add_element(element)
		@elements << element
	end
	def to_s
		svg_string = ""
		svg_string << "<?xml version=\"1.0\" standalone=\"no\"?>\n"
		svg_string << "<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\" \"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">\n"
		svg_string << "<svg width=\"#{@width}mm\" height=\"#{@height}mm\" viewBox=\"0 0 #{@width} #{@height}\" xmlns=\"http://www.w3.org/2000/svg\" version=\"1.1\">\n"
		@elements.each {|e| svg_string << "    #{e.to_s}\n"}
		svg_string << "</svg>\n"
		svg_string
	end
end

# Define a class for lines of the SVG graph.
class Line
	attr_reader :x_start, :y_start, :x_end, :y_end, :stroke_color, :stroke_width, :opacity
	def initialize(x_start,y_start,x_end,y_end,stroke_color,stroke_width,opacity)
		@x_start = x_start
		@y_start = y_start
		@x_end = x_end
		@y_end = y_end
		if stroke_color == nil
			@stroke_color = "black"
		else
			@stroke_color = stroke_color
		end
		if stroke_width == nil
			@stroke_width = 1.0
		else
			@stroke_width = stroke_width
		end
		if opacity == nil
			@opacity = 1.0
		else
			@opacity = opacity
		end
	end
	def to_s
		svg = "<line x1=\"#{@x_start.round(3)}\" y1=\"#{@y_start.round(3)}\" x2=\"#{@x_end.round(3)}\" y2=\"#{@y_end.round(3)}\" stroke=\"#{@stroke_color}\" stroke-width=\"#{@stroke_width}\" stroke-opacity=\"#{@opacity}\" />"
		# svg = "<line x1=\"#{@x_start.round(3)}\" y1=\"#{@y_start.round(3)}\" x2=\"#{@x_end.round(3)}\" y2=\"#{@y_end.round(3)}\" stroke=\"#{@stroke_color}\" stroke-width=\"#{@stroke_width}\" />"
		svg
	end
end

# Define class for circles of the SVG graph.
class Circle
	attr_reader :x, :y, :r, :fill_color, :stroke_color, :stroke_width, :opacity
	def initialize(x,y,r,fill_color,stroke_color,stroke_width,opacity)
		@x = x
		@y = y
		@r = r
		@fill_color = fill_color
		@stroke_color = stroke_color
		@stroke_width = stroke_width
		@opacity = opacity
	end
	def to_s
		svg = "<circle cx=\"#{@x.round(3)}\" cy=\"#{y.round(3)}\" r=\"#{r.round(3)}\" fill=\"#{@fill_color}\" stroke=\"#{@stroke_color}\" stroke-width=\"#{@stroke_width}\" fill-opacity=\"#{@opacity}\" stroke-opacity=\"#{@opacity}\" />"
		# svg = "<circle cx=\"#{@x.round(3)}\" cy=\"#{y.round(3)}\" r=\"#{r.round(3)}\" fill=\"#{@fill_color}\" stroke=\"#{@stroke_color}\" stroke-width=\"#{@stroke_width}\" />"
		svg
	end
end

# Define class for rectangles of the SVG graph.
class Rectangle
	attr_reader :x, :y, :width, :height, :fill_color, :stroke_color, :stroke_width, :opacity
	def initialize(x,y,width,height,fill_color,stroke_color,stroke_width,opacity)
		@x = x
		@y = y
		@width = width
		@height = height
		if fill_color == nil
			@fill_color = "none"
		else
			@fill_color = fill_color
		end
		if stroke_color == nil
			@stroke_color = "black"
		else
			@stroke_color = stroke_color
		end
		if stroke_width == nil
			@stroke_width = 1.0
		else
			@stroke_width = stroke_width
		end
		if opacity == nil
			@opacity = 1.0
		else
			@opacity = opacity
		end
	end
	def to_s
		svg = "<rect x=\"#{@x}\" y=\"#{@y}\" width=\"#{@width}\" height=\"#{@height}\" fill=\"#{@fill_color}\" stroke=\"#{@stroke_color}\" stroke-width=\"#{@stroke_width}\" fill-opacity=\"#{@opacity}\" stroke-opacity=\"#{@opacity}\" />"
		# svg = "<rect x=\"#{@x}\" y=\"#{@y}\" width=\"#{@width}\" height=\"#{@height}\" fill=\"#{@fill_color}\" stroke=\"#{@stroke_color}\" stroke-width=\"#{@stroke_width}\" />"
		svg
	end
end

# Define class for paths of the SVG graph.
class Path
	attr_reader :x, :y, :fill_color, :stroke_color, :stroke_width, :opacity
	def initialize(x,y,fill_color,stroke_color,stroke_width,opacity)
		@x = [x]
		@y = [y]
		@fill_color = fill_color
		@stroke_color = stroke_color
		@stroke_width = stroke_width
		@opacity = opacity
	end
	def add_point(x,y)
		@x << x
		@y << y
	end
	def to_s
		svg = "<path d=\"M #{@x[0]} #{@y[0]} "
		if @x.size > 1
			1.upto(@x.size-1) do |z|
				svg << "L #{@x[z]} #{@y[z]} "
			end
		end
		svg << "z\" fill=\"#{@fill_color}\" stroke=\"#{@stroke_color}\" stroke-width=\"#{@stroke_width}\" fill-opacity=\"#{@opacity}\" stroke-opacity=\"#{@opacity}\"/>"
		# svg << "z\" fill=\"#{@fill_color}\" stroke=\"#{@stroke_color}\" stroke-width=\"#{@stroke_width}\" />"
		svg
	end
end

# Define class for texts of the SVG graph.
class Text
	attr_reader :x, :y
	def initialize(x,y,font_size,string,anchor,baseline)
		@x = x
		@y = y
		@font_size = font_size
		@string = string
		if anchor == nil
			@anchor = "middle"
		else
			@anchor = anchor
		end
		if baseline == nil
			@baseline = "auto"
		else
			@baseline = baseline
		end
	end
	def to_s
		svg = "<text text-anchor=\"#{@anchor}\" dominant-baseline=\"#{@baseline}\" x=\"#{@x}\" y=\"#{@y}\" font-family=\"Helvetica\" font-size=\"#{@font_size}pt\">#{@string}</text>"
		svg
	end
end

# Get the command line arguments.
filtered_blast_output_file_name = ARGV[0]
ld_file_name = ARGV[1]
target_region_features_file_name = ARGV[2]
anomalous_links_file_name = ARGV[3]
plot_output_file_name = ARGV[4]
bitscore_threshold = ARGV[5].to_i 
n_scaffolds_for_plot = 15

# Read the BLAST output file.
filtered_blast_output_file = File.open(filtered_blast_output_file_name)
filtered_blast_output_lines = filtered_blast_output_file.readlines
sseqids = []
pidents = []
qstarts = []
qends = []
bitscores = []
dirs = []
filtered_blast_output_lines.each do |l|
	line_ary = l.split
	sseqids << line_ary[1]
	pidents << line_ary[2].to_f
	qstarts << line_ary[6].to_i
	qends << line_ary[7].to_i
	bitscores << line_ary[11].to_i
	dirs << line_ary[12]
end

# Remove all hits below the bitscore threshold.
sseqids.size.times do |x|
	if bitscores[x] < bitscore_threshold
		sseqids[x] = nil
		pidents[x] = nil
		qstarts[x] = nil
		qends[x] = nil
		bitscores[x] = nil
		dirs[x] = nil
	end
end
sseqids.compact!
pidents.compact!
qstarts.compact!
qends.compact!
bitscores.compact!
dirs.compact!

# Find the sseqids with the longest overall hit length.
uniq_sseqids = sseqids.uniq
n_positions_covered_per_seqs = []
uniq_sseqids.each do |sseqid|
	positions_covered_by_this_sseq = []
	sseqids.size.times do |x|
		if sseqids[x] == sseqid
			([qstarts[x],qends[x]].min).upto([qstarts[x],qends[x]].max) do |x|
				positions_covered_by_this_sseq << x
			end
		end
	end
	n_positions_covered_per_seqs << positions_covered_by_this_sseq.uniq.size
end
best_uniq_sseqids = []
n_positions_covered_by_best_uniq_sseqids = []
uniq_sseqids.size.times do |x|
	if n_positions_covered_per_seqs[x] >= n_positions_covered_per_seqs.sort.reverse[n_scaffolds_for_plot-1]
		best_uniq_sseqids << uniq_sseqids[x]
		n_positions_covered_by_best_uniq_sseqids << n_positions_covered_per_seqs[x]
	end
end
sorted = false
until sorted
	sorted = true
	(best_uniq_sseqids.size-1).times do |x|
		if n_positions_covered_by_best_uniq_sseqids[x] < n_positions_covered_by_best_uniq_sseqids[x+1]
			n_positions_covered_by_best_uniq_sseqids[x],n_positions_covered_by_best_uniq_sseqids[x+1] = n_positions_covered_by_best_uniq_sseqids[x+1],n_positions_covered_by_best_uniq_sseqids[x]
			best_uniq_sseqids[x],best_uniq_sseqids[x+1] = best_uniq_sseqids[x+1],best_uniq_sseqids[x]
			sorted = false
		end
	end
end

# Determine the main direction of each of the best unique sseqids.
main_dirs_of_best_uniq_sseqids = []
best_uniq_sseqids.each do |sseqid|
	dirs_of_uniq_sseqid = []
	sseqids.size.times do |x|
		if sseqids[x] == sseqid
			dirs_of_uniq_sseqid << dirs[x]
		end
	end
	if dirs_of_uniq_sseqid.count("rev") > dirs_of_uniq_sseqid.count("for")
		main_dirs_of_best_uniq_sseqids << "rev"
	else
		main_dirs_of_best_uniq_sseqids << "for"
	end
end

# Initialize the SVG graph.
svg_width = 180
svg_margin = 5
tick_length = 0.5
font_size = 2.4705882353
font_spacer = 0.5
window_width = svg_width-2*svg_margin
scaffold_colors = ['859900', 'b58900', '2aa198', '268bd2', '6c71c4', 'd33682', 'dc322f', 'cb4b16', '002b36', '586e75', '839496', 'a2aca8', 'c6c6bc']
pidents_colors = ['561a44', '6f1642', '7c1441', '8e113f', 'a70e3e', 'b30c3d', 'c50a3c', 'de2d3c', 'ea3e3c', 'fc583c', 'fd8835', 'fd9f32', 'fdc22d']
ld_circle_color = "000000"
scaffold_id_field_width = 25
plot_width = window_width - scaffold_id_field_width
chromosome_window_size = 1000000
ld_window_size = 250000
ld_field_height = (0.5*ld_window_size/chromosome_window_size.to_f) * plot_width
ld_circle_r = 0.15
r2_threshold = 0.8
plot_spacer = 2
scaffold_spacer = 5
scaffold_height = 5
scaffold_field_height = (n_scaffolds_for_plot+1) * scaffold_spacer
svg_height = 2*svg_margin+ld_field_height+plot_spacer+scaffold_field_height+tick_length+font_spacer
window_height = svg_height-2*svg_margin
svg = SVG.new(svg_width,svg_height)

# Find the start and end positions of the window on the chromosome.
name_array = filtered_blast_output_file_name.split("/").last.split("_")
lg = name_array[0]
window_start = name_array[1].to_f
window_end = name_array[2].to_f

# Read the ld file and add lines indicating all SNP positions.
ld_file = File.open(ld_file_name)
ld_lines = ld_file.readlines
snp1_positions = []
ld_lines[1..-1].each do |l|
	line_ary = l.split
	snp1_pos = line_ary[1].to_i
	snp1_positions << snp1_pos
end
snp1_positions.uniq.each do |snp1_pos|
	if snp1_pos > window_start*1000000 and snp1_pos < window_end*1000000
		x1 = svg_margin+scaffold_id_field_width+((snp1_pos-window_start*1000000)/chromosome_window_size.to_f)*(window_width-scaffold_id_field_width)
		y1 = svg_margin+ld_field_height
		line_span = ([0.5*ld_window_size,window_end*1000000-snp1_pos].min/chromosome_window_size.to_f) * plot_width
		x2 = x1 + line_span
		y2 = svg_margin + ld_field_height - line_span
		svg.add_element(Line.new(x1,y1,x2,y2,"#EFEFEF",0.1,nil))
		line_span = ([0.5*ld_window_size,snp1_pos-window_start*1000000].min/chromosome_window_size.to_f) * plot_width
		x2 = x1 - line_span
		y2 = svg_margin + ld_field_height - line_span
		svg.add_element(Line.new(x1,y1,x2,y2,"#EFEFEF",0.1,nil))
	elsif snp1_pos < window_start*1000000 and snp1_pos > window_start*1000000 - 0.5*ld_window_size
		x1 = svg_margin+scaffold_id_field_width
		line_span = ((0.5*ld_window_size-(window_start*1000000-snp1_pos))/chromosome_window_size.to_f) * plot_width
		y1 = svg_margin + line_span
		x2 = x1 + line_span
		y2 = svg_margin
		svg.add_element(Line.new(x1,y1,x2,y2,"#EFEFEF",0.1,nil))		
	elsif snp1_pos > window_end*1000000 and snp1_pos < window_end*1000000 + 0.5*ld_window_size
		x1 = svg_margin+scaffold_id_field_width+plot_width
		line_span = ((0.5*ld_window_size-(snp1_pos-window_end*1000000))/chromosome_window_size.to_f) * plot_width
		y1 = svg_margin + line_span
		x2 = x1 - line_span
		y2 = svg_margin
		svg.add_element(Line.new(x1,y1,x2,y2,"#EFEFEF",0.1,nil))
	end
end

# Add circles to indicate pairs of SNPs that are out of ld.
ld_lines[1..-1].each do |l|
	line_ary = l.split
	snp1_pos = line_ary[1].to_i
	snp2_pos = line_ary[4].to_i
	r2 = line_ary[6].to_f
	if snp2_pos - snp1_pos < ld_window_size
		if [snp1_pos,snp2_pos].min > window_start*1000000 - 0.5*ld_window_size and [snp1_pos,snp2_pos].max < window_end*1000000 + 0.5*ld_window_size
			if (snp1_pos+snp2_pos)/2.0 > window_start*1000000 and (snp1_pos+snp2_pos)/2.0 < window_end*1000000
				x = svg_margin+scaffold_id_field_width + ((((snp1_pos+snp2_pos)/2.0)-window_start*1000000)/chromosome_window_size.to_f) * plot_width
				y = svg_margin + ld_field_height - ((snp2_pos-snp1_pos)/2.0)/chromosome_window_size.to_f * plot_width
				if r2 > r2_threshold
					svg.add_element(Circle.new(x,y,ld_circle_r,"##{ld_circle_color}","none","none",((r2-r2_threshold)/(1.0-r2_threshold))))
				end
			end
		end
	end
end

# Read the features file.
target_region_features_file = File.open(target_region_features_file_name)
target_region_features_lines = target_region_features_file.readlines
feature_from = []
feature_to = []
feature_missing = []
feature_repetitive = []
target_region_features_lines.each do |l|
	line_ary = l.split
	feature_from << line_ary[0].to_i
	feature_to << line_ary[1].to_i
	feature_missing << line_ary[2]
	feature_repetitive << line_ary[3].strip
end

# Add rectangles for the missing data and repetitive regions.
feature_from.size.times do |x|
	if feature_missing[x] == "true" or feature_repetitive[x] == "true"
		rectangle_x = svg_margin + scaffold_id_field_width + (feature_from[x]/chromosome_window_size.to_f)*(window_width-scaffold_id_field_width)
		rectangle_y = svg_margin+ld_field_height+plot_spacer
		rectangle_width = ((feature_to[x]-feature_from[x])/chromosome_window_size.to_f)*(window_width-scaffold_id_field_width)
		rectangle_height = scaffold_field_height
	end
	if feature_missing[x] == "true"
		svg.add_element(Rectangle.new(rectangle_x,rectangle_y,rectangle_width,rectangle_height,"gray","none","none","1.0"))
	elsif feature_repetitive[x] == "true"
		svg.add_element(Rectangle.new(rectangle_x,rectangle_y,rectangle_width,rectangle_height,"#E0E0E0","none","none","1.0"))
	end
end

# Add rectangles for the forward BLAST hits.
count = 0.5
best_uniq_sseqids.size.times do |z|
	sseqid = best_uniq_sseqids[z]
	text_x = svg_margin
	text_y = svg_margin + ld_field_height + plot_spacer + scaffold_spacer*count + 0.5*scaffold_height
	svg.add_element(Text.new(text_x,text_y,font_size,"#{sseqid.chomp("_pilon")}","start","mathematical"))
	sseqids.size.times do |x|
		if sseqids[x] == sseqid
			rectangle_x = svg_margin + scaffold_id_field_width + (qstarts[x]/chromosome_window_size.to_f)*(window_width-scaffold_id_field_width)
			rectangle_width = ((qends[x]-qstarts[x])/chromosome_window_size.to_f)*(window_width-scaffold_id_field_width)
			if dirs[x] == "for" # main_dirs_of_best_uniq_sseqids[z]
				rectangle_height = scaffold_height
				rectangle_y = svg_margin + ld_field_height + plot_spacer + scaffold_spacer*count
				if pidents[x] > 99
					rectangle_color = "##{pidents_colors[-1]}"
				elsif pidents[x] > 98
					rectangle_color = "##{pidents_colors[-2]}"
				elsif pidents[x] > 97
					rectangle_color = "##{pidents_colors[-3]}"
				elsif pidents[x] > 96
					rectangle_color = "##{pidents_colors[-4]}"
				elsif pidents[x] > 95
					rectangle_color = "##{pidents_colors[-5]}"
				elsif pidents[x] > 94
					rectangle_color = "##{pidents_colors[-6]}"
				elsif pidents[x] > 93
					rectangle_color = "##{pidents_colors[-7]}"
				elsif pidents[x] > 92
					rectangle_color = "##{pidents_colors[-8]}"
				elsif pidents[x] > 91
					rectangle_color = "##{pidents_colors[-9]}"
				elsif pidents[x] > 90
					rectangle_color = "##{pidents_colors[-10]}"
				elsif pidents[x] > 89
					rectangle_color = "##{pidents_colors[-11]}"
				elsif pidents[x] > 89
					rectangle_color = "##{pidents_colors[-12]}"
				else
					rectangle_color = "##{pidents_colors[-13]}"
				end
				svg.add_element(Rectangle.new(rectangle_x,rectangle_y,rectangle_width,rectangle_height,rectangle_color,"none","none","1.0"))
			end
		end
	end
	count += 1
	break if count > n_scaffolds_for_plot
end

# Add rectangles for the reverse BLAST hits.
count = 0.5
best_uniq_sseqids.size.times do |z|
	sseqid = best_uniq_sseqids[z]
	text_x = svg_margin
	text_y = svg_margin + ld_field_height + plot_spacer + scaffold_spacer*count + 0.5*scaffold_height
	svg.add_element(Text.new(text_x,text_y,font_size,"#{sseqid.chomp("_pilon")}","start","mathematical"))
	sseqids.size.times do |x|
		if sseqids[x] == sseqid
			rectangle_x = svg_margin + scaffold_id_field_width + (qstarts[x]/chromosome_window_size.to_f)*(window_width-scaffold_id_field_width)
			rectangle_width = ((qends[x]-qstarts[x])/chromosome_window_size.to_f)*(window_width-scaffold_id_field_width)
			if dirs[x] == "rev"
				rectangle_height = 0.5 * scaffold_height
				rectangle_y = svg_margin + ld_field_height + plot_spacer + scaffold_spacer*count + 0.25 * scaffold_height
				if pidents[x] > 99
					rectangle_color = "##{pidents_colors[-1]}"
				elsif pidents[x] > 98
					rectangle_color = "##{pidents_colors[-2]}"
				elsif pidents[x] > 97
					rectangle_color = "##{pidents_colors[-3]}"
				elsif pidents[x] > 96
					rectangle_color = "##{pidents_colors[-4]}"
				elsif pidents[x] > 95
					rectangle_color = "##{pidents_colors[-5]}"
				elsif pidents[x] > 94
					rectangle_color = "##{pidents_colors[-6]}"
				elsif pidents[x] > 93
					rectangle_color = "##{pidents_colors[-7]}"
				elsif pidents[x] > 92
					rectangle_color = "##{pidents_colors[-8]}"
				elsif pidents[x] > 91
					rectangle_color = "##{pidents_colors[-9]}"
				elsif pidents[x] > 90
					rectangle_color = "##{pidents_colors[-10]}"
				elsif pidents[x] > 89
					rectangle_color = "##{pidents_colors[-11]}"
				elsif pidents[x] > 89
					rectangle_color = "##{pidents_colors[-12]}"
				else
					rectangle_color = "##{pidents_colors[-13]}"
				end
				svg.add_element(Rectangle.new(rectangle_x,rectangle_y,rectangle_width,rectangle_height,rectangle_color,"none","none","1.0"))
			end
		end
	end
	count += 1
	break if count > n_scaffolds_for_plot
end

# Read the file with anomalous link information (from the BreakDancer analysis).
anomalous_links_file = File.open(anomalous_links_file_name)
anomalous_links_lines = anomalous_links_file.readlines
anomalous_link_lgs = []
anomalous_link_pos_froms = []
anomalous_link_pos_tos = []
anomalous_links_lines.each do |l|
	line_ary = l.split
	anomalous_link_lgs << line_ary[2]
	anomalous_link_pos_froms << line_ary[3].to_i
	anomalous_link_pos_tos << line_ary[4].to_i
end

# Add lines for anomalous links.
anomalous_link_lgs.size.times do |x|
	plot_this_link = false
	if anomalous_link_lgs[x].downcase.sub("lg0","lg") == lg.downcase.sub("lg0","lg")
		if anomalous_link_pos_froms[x] > window_start*1000000 and anomalous_link_pos_froms[x] < window_end*1000000
			plot_this_link = true
			line_x = svg_margin+scaffold_id_field_width + ((anomalous_link_pos_froms[x]-window_start*1000000)/chromosome_window_size.to_f) * plot_width
		elsif anomalous_link_pos_tos[x] > window_start*1000000 and anomalous_link_pos_tos[x] < window_end*1000000
			plot_this_link = true
			line_x = svg_margin+scaffold_id_field_width + ((anomalous_link_pos_tos[x]-window_start*1000000)/chromosome_window_size.to_f) * plot_width
		end
	end
	if plot_this_link
		line_y1 = svg_margin+ld_field_height+plot_spacer
		line_y2 = svg_margin+ld_field_height+plot_spacer+scaffold_field_height
		svg.add_element(Line.new(line_x,line_y1,line_x,line_y2,nil,nil,nil))
	end
end

# Add lines to vertically separate contigs.
count = 1.5
(n_scaffolds_for_plot-1).times do |x|
	line_y = svg_margin+ld_field_height+plot_spacer+scaffold_spacer*count
	svg.add_element(Line.new(svg_margin+scaffold_id_field_width,line_y,window_width+svg_margin,line_y,nil,0.1,nil))
	count += 1
end

# Add frames to the SVG.
svg.add_element(Rectangle.new(svg_margin+scaffold_id_field_width,svg_margin,window_width-scaffold_id_field_width,ld_field_height,nil,"303030",0.1,nil))
svg.add_element(Rectangle.new(svg_margin+scaffold_id_field_width,svg_margin+ld_field_height+plot_spacer,window_width-scaffold_id_field_width,scaffold_field_height,nil,"black",0.25,nil))
# svg.add_element(Line.new(svg_margin+scaffold_id_field_width,svg_margin,window_width+svg_margin,svg_margin,nil,0.5,nil))
svg.add_element(Line.new(svg_margin+scaffold_id_field_width,svg_margin+ld_field_height,window_width+svg_margin,svg_margin+ld_field_height,nil,0.1,nil))
svg.add_element(Line.new(svg_margin+scaffold_id_field_width,svg_margin+ld_field_height+plot_spacer,window_width+svg_margin,svg_margin+ld_field_height+plot_spacer,nil,0.1,nil))
svg.add_element(Line.new(svg_margin+scaffold_id_field_width,svg_margin+ld_field_height+plot_spacer+scaffold_field_height,window_width+svg_margin,svg_margin+ld_field_height+plot_spacer+scaffold_field_height,nil,0.1,nil))

# Add ticks to the frame.
tick_x = svg_margin + scaffold_id_field_width
svg.add_element(Line.new(tick_x,svg_margin+ld_field_height+plot_spacer+scaffold_field_height,tick_x,svg_margin+ld_field_height+plot_spacer+scaffold_field_height+tick_length,nil,0.1,nil))
svg.add_element(Text.new(tick_x,svg_margin+ld_field_height+plot_spacer+scaffold_field_height+tick_length+font_spacer,font_size,"#{window_start}",nil,"hanging"))
n_ticks = 10
n_ticks.times do |x|
	tick_x = svg_margin + scaffold_id_field_width + (x+1)*((window_width-scaffold_id_field_width)/n_ticks.to_f)
	svg.add_element(Line.new(tick_x,svg_margin+ld_field_height+plot_spacer+scaffold_field_height,tick_x,svg_margin+ld_field_height+plot_spacer+scaffold_field_height+tick_length,nil,0.1,nil))
	svg.add_element(Text.new(tick_x,svg_margin+ld_field_height+plot_spacer+scaffold_field_height+tick_length+font_spacer,font_size,"#{window_start+(window_end-window_start)*((x+1)/n_ticks.to_f)}",nil,"hanging"))
end

# Write the SVG graph to the plot output file.
plot_output_file = File.open(plot_output_file_name,"w")
plot_output_file.write(svg.to_s)
