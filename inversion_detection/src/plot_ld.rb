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
ld_file_name = ARGV[0]
plot_output_file_name = ARGV[1]

# Initialize the SVG graph.
svg_width = 180
svg_margin = 5
tick_length = 0.5
font_size = 2.4705882353
font_spacer = 0.5
window_width = svg_width-2*svg_margin
ld_circle_color = "000000"
plot_width = window_width
chromosome_window_size = 50000000
ld_window_size = 50000000
ld_field_height = (0.5*ld_window_size/chromosome_window_size.to_f) * plot_width
ld_circle_r = 0.15
r2_threshold = 0.5
plot_spacer = 2
svg_height = 2*svg_margin+ld_field_height+plot_spacer+tick_length+font_spacer
window_height = svg_height-2*svg_margin
svg = SVG.new(svg_width,svg_height)

# Find the start and end positions of the window on the chromosome.
window_start = 0
window_end = 50

# Read the ld file and add lines indicating all SNP positions.
ld_file = File.open(ld_file_name)
ld_lines = ld_file.readlines

# Add circles to indicate pairs of SNPs that are out of ld.
ld_lines[1..-1].each do |l|
	line_ary = l.split
	snp1_pos = line_ary[1].to_i
	snp2_pos = line_ary[4].to_i
	r2 = line_ary[6].to_f
	if snp2_pos - snp1_pos < ld_window_size
		if [snp1_pos,snp2_pos].min > window_start*1000000 - 0.5*ld_window_size and [snp1_pos,snp2_pos].max < window_end*1000000 + 0.5*ld_window_size
			if (snp1_pos+snp2_pos)/2.0 > window_start*1000000 and (snp1_pos+snp2_pos)/2.0 < window_end*1000000
				x = svg_margin + ((((snp1_pos+snp2_pos)/2.0)-window_start*1000000)/chromosome_window_size.to_f) * plot_width
				y = svg_margin + ld_field_height - ((snp2_pos-snp1_pos)/2.0)/chromosome_window_size.to_f * plot_width
				if r2 > r2_threshold
					svg.add_element(Circle.new(x,y,ld_circle_r,"##{ld_circle_color}","none","none",((r2-r2_threshold)/(1.0-r2_threshold))))
				end
			end
		end
	end
end

# Add frames to the SVG.
svg.add_element(Rectangle.new(svg_margin,svg_margin,window_width,ld_field_height,nil,"#303030",0.1,nil))
svg.add_element(Line.new(svg_margin,svg_margin+ld_field_height,window_width+svg_margin,svg_margin+ld_field_height,nil,0.1,nil))
svg.add_element(Line.new(svg_margin,svg_margin+ld_field_height+plot_spacer,window_width+svg_margin,svg_margin+ld_field_height+plot_spacer,nil,0.1,nil))
svg.add_element(Line.new(svg_margin,svg_margin+ld_field_height+plot_spacer,window_width+svg_margin,svg_margin+ld_field_height+plot_spacer,nil,0.1,nil))

# Add ticks to the frame.
tick_x = svg_margin
svg.add_element(Line.new(tick_x,svg_margin+ld_field_height+plot_spacer,tick_x,svg_margin+ld_field_height+plot_spacer+tick_length,nil,0.1,nil))
svg.add_element(Text.new(tick_x,svg_margin+ld_field_height+plot_spacer+tick_length+font_spacer,font_size,"#{window_start}",nil,"hanging"))
n_ticks = 50
n_ticks.times do |x|
	tick_x = svg_margin + (x+1)*((window_width)/n_ticks.to_f)
	svg.add_element(Line.new(tick_x,svg_margin+ld_field_height+plot_spacer,tick_x,svg_margin+ld_field_height+plot_spacer+tick_length,nil,0.1,nil))
	svg.add_element(Text.new(tick_x,svg_margin+ld_field_height+plot_spacer+tick_length+font_spacer,font_size,"#{window_start+(window_end-window_start)*((x+1)/n_ticks.to_f)}",nil,"hanging"))
end

# Write the SVG graph to the plot output file.
plot_output_file = File.open(plot_output_file_name,"w")
plot_output_file.write(svg.to_s)
