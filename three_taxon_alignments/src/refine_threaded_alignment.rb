# m_matschiner Sat Jun 30 00:08:56 CEST 2018

# Get the command-line arguments.
alignment_file_name = ARGV[0]
refined_alignment_file_name = ARGV[1]

# Read the alignment file.
print "Reading the alignment file #{alignment_file_name}..."
STDOUT.flush
alignment_file = File.open(alignment_file_name)
alignment_lines = alignment_file.readlines
alignment_ids = []
alignment_seqs = []
alignment_lines.each do |l|
  if l.strip[0] == ">"
    alignment_ids << l.strip[1..-1]
    alignment_seqs << ""
  elsif l.strip != ""
    alignment_seqs.last << l.strip
  end
end
alignment_seqs_as_arrays = []
alignment_seqs.each do |s|
  alignment_seq_as_array = []
  s.size.times do |pos|
    alignment_seq_as_array << s[pos]
  end
  alignment_seqs_as_arrays << alignment_seq_as_array
end
puts " done."
STDOUT.flush

# Get the alignment length.
alignment_length = alignment_seqs_as_arrays[0].size

# Set the window size.
window_size = 1000
window_start = 0

# Initialize an array of refined sequences.
refined_alignment_seqs = []
alignment_seqs.size.times { refined_alignment_seqs << "" }

# Extract windows, write them to a temporary file and run mafft with it.
print "Refining the alignment in windows with MAFFT..."
STDOUT.flush
while window_start < alignment_length

  # Extract windows.
  window_end = window_start + window_size - 1
  window_end = alignment_length - 1 if window_end > alignment_length
  window_seqs = []
  alignment_seqs_as_arrays.each do |s|
    window_seq_as_array = s[window_start..window_end]
    window_seq = ""
    window_seq_as_array.each { |i| window_seq << i }
    window_seqs << window_seq
  end

  # Feedback.
  puts "Current window: #{window_start}-#{window_end}"
  STDOUT.flush

  # Prepare an output string for a temporary fasta file.
  tmp_alignment_string = ""
  alignment_ids.size.times do |x|
    tmp_alignment_string << ">#{alignment_ids[x]}\n"
    tmp_alignment_string << "#{window_seqs[x]}\n"
  end

  # Write the temporary fasta file.
  tmp_alignment_file_name = "window.fasta"
  tmp_alignment_file = File.open(tmp_alignment_file_name, "w")
  tmp_alignment_file.write(tmp_alignment_string)
  tmp_alignment_file.close

  # Run mafft to realign the alignment in the temporary file.
  mafft_alignment_file_name = "window_aligned.fasta"
  system( "mafft --quiet --auto #{tmp_alignment_file_name} > #{mafft_alignment_file_name}" )

  # Read the realigned output of mafft.
  mafft_alignment_file = File.open(mafft_alignment_file_name)
  mafft_alignment_lines = mafft_alignment_file.readlines
  mafft_alignment_ids = []
  mafft_alignment_seqs = []
  mafft_alignment_lines.each do |l|
    if l.strip[0] == ">"
      mafft_alignment_ids << l.strip[1..-1]
      mafft_alignment_seqs << ""
    elsif l.strip != ""
      mafft_alignment_seqs.last << l.strip
    end
  end

  # Make sure the order of ids is the same in the mafft output and the original fasta file.
  if mafft_alignment_ids != alignment_ids
    puts "ERROR: The order of IDs is different between the MAFFT output file and the original fasta file!"
    puts "  Original alignment IDs:"
    alignment_ids.each do |i|
      puts "  #{i}"
    end
    puts "  MAFFT alignment IDs:"
    mafft_alignment_ids.each do |i|
      puts "  #{i}"
    end
    exit 1
  end

  # Add the refined sequences from this window to the overall refined alignment.
  mafft_alignment_seqs.size.times do |x|
    refined_alignment_seqs[x] << mafft_alignment_seqs[x]
  end

  # Move to the next window.
  window_start = window_end + 1

end
puts " done."
STDOUT.flush

# Prepare an output string for the refined alignment.
print "Preparing the output..."
STDOUT.flush
refined_alignment_string = ""
mafft_alignment_ids.size.times do |x|
  refined_alignment_string << ">#{mafft_alignment_ids[x]}\n"
  refined_alignment_string << "#{refined_alignment_seqs[x]}\n"
end
puts " done."
STDOUT.flush

# Write the refined alignment to a new file.
refined_alignment_file = File.open(refined_alignment_file_name, "w")
refined_alignment_file.write(refined_alignment_string)
puts "Wrote file #{refined_alignment_file_name}."
