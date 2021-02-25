# m_matschiner Wed Jun 14 14:15:13 CEST 2017

# Get the command-line argument.
input_file_name = ARGV[0]

# Read the input file.
input_file = File.open(input_file_name)
input_lines = input_file.readlines

# For each chromosome, extract the right lines from the input.
1.upto(23) do |x|

	# Initiate a variable that declares whether or not lines for the right linkage group are currently read.
	in_lg = false

	# Get the name of this linkage group.
	lg = "lg#{x.to_s.rjust(2).sub(" ","0")}"

	# Initiate the output for this chromosome.
	out_string = ""

	# Extract lines from the input.
	input_lines.each do |l|
		if l[0..1] == "lg"
			if l[0..3] == lg
				in_lg = true
			else
				in_lg = false
			end
		else
			if in_lg
				output = true
				output = false if l.strip == ""
				output = false if l.include?("continued")
				output = false if l.include?("elsewhere")
				output = false if l.include?("breakdancer_breakpoint")
				out_string << l if output
			end
		end
	end

	# Define the file name for this output file.
	outfile_name = "#{input_file_name.chomp(".txt")}_#{lg}.txt"

	# Write the output file.
	outfile = File.open(outfile_name,"w")
	outfile.write(out_string)
end
