# m_matschiner Wed Apr 4 13:58:50 CEST 2018

# Get the command-line argument.
rm_input_file_name = ARGV[0]
trf_input_file_name = rm_input_file_name.sub(".rm.",".trf.")
output_file_name = rm_input_file_name.sub(".rm.",".")

# Read the input file with repeatmasker's repeat mask.
rm_input_file = File.open(rm_input_file_name)
rm_lines = rm_input_file.read.lines
rm_ids = []
rm_seqs= []
rm_lines.each do |l|
  if l[0] == ">"
    rm_ids << l.strip[1..-1]
    rm_seqs << ""
  elsif l.strip != ""
    rm_seqs.last << l.strip
  end
end

# Read the input file with trf's repeat mask.
trf_input_file = File.open(trf_input_file_name)
trf_lines = trf_input_file.read.lines
trf_ids = []
trf_seqs = []
trf_lines.each do |l|
  if l[0]== ">"
    trf_ids << l.strip[1..-1]
    trf_seqs << ""
  elsif l.strip != ""
    trf_seqs.last << l.strip
  end
end

# Make sure the sequence ids are identical.
unless rm_ids == trf_ids
  puts "ERROR: The sequence ids differ between the two files!"
  exit 1
end

# Prepare the output sequences.
out_seqs = []
rm_seqs.size.times do |x|
  unless rm_seqs[x].size == trf_seqs[x].size
    puts "ERROR: The sequences with id (#{rm_ids[x]}) differ in length between the two files!"
    exit 1
  end
  out_seq = ""
  rm_seqs[x].size.times do |pos|
    if trf_seqs[x][pos] == "N"
      out_seq << rm_seqs[x][pos].downcase
    else
      out_seq << rm_seqs[x][pos]
    end
  end
  out_seqs << out_seq
end

# Prepare the output string.
out_string = ""
rm_ids.size.times do |x|
  out_string << ">#{rm_ids[x]}\n"
  out_seq = out_seqs[x]
  while out_seq.size > 0
    out_string << "#{out_seq.slice!(0..59)}\n"
  end
end

# Write the output file.
output_file = File.open(output_file_name,"w")
output_file.write(out_string)
