# m_matschiner Wed Jun 27 10:36:48 CEST 2018

# Get the command-line arguments.
merged_alignment_file_name = ARGV[0]
threaded_alignment_file_name = ARGV[1]

# Read the input file.
print "Reading file #{merged_alignment_file_name}..." 
merged_alignment_file = File.open(merged_alignment_file_name)
merged_alignment_lines = merged_alignment_file.readlines
merged_alignment_ids = []
merged_alignment_seqs = []
new_merged_alignment_seqs = []
merged_alignment_lines.each do |l|
    if l[0] == ">"
        merged_alignment_ids << l.strip[1..-1]
    elsif l.strip != ""
        merged_alignment_seq = []
        l.strip!
        l.size.times do |x|
            merged_alignment_seq << l[x]
        end
        merged_alignment_seqs << merged_alignment_seq
        new_merged_alignment_seqs << []
    end
end
puts " done."
STDOUT.flush

# Get the indices of the reference.
ref_indices = []
merged_alignment_ids.size.times do |x|
    if merged_alignment_ids[x] == "reference"
        ref_indices << x
    end
end

# Feedback.
print "Threading alignment..."
STDOUT.flush

# Set up a loop that runs until threading is complete.
threading_complete = false
positions = []
merged_alignment_ids.size.times { |x| positions << 0 }
until threading_complete

    # Read the nucleotides of the reference at this position, make adjustments if any are nil.
    ref_nucleotides_at_this_pos = []
    ref_indices.each { |i| ref_nucleotides_at_this_pos << merged_alignment_seqs[i][positions[i]] }
    ref_nucleotides_at_this_pos_uniq = ref_nucleotides_at_this_pos.uniq
    if ref_nucleotides_at_this_pos_uniq == [nil]

        threading_complete = true

    else

        # Check if any but not all of the reference nucleotides are gaps.
        gap_at_this_pos = false
        gap_at_this_pos = true if ref_nucleotides_at_this_pos.include?("-")
        nucleotide_at_this_pos = false
        nucleotide_at_this_pos = true if ref_nucleotides_at_this_pos.include?("A") or ref_nucleotides_at_this_pos.include?("a")
        nucleotide_at_this_pos = true if ref_nucleotides_at_this_pos.include?("C") or ref_nucleotides_at_this_pos.include?("c")
        nucleotide_at_this_pos = true if ref_nucleotides_at_this_pos.include?("G") or ref_nucleotides_at_this_pos.include?("g")
        nucleotide_at_this_pos = true if ref_nucleotides_at_this_pos.include?("T") or ref_nucleotides_at_this_pos.include?("t")
        nucleotide_at_this_pos = true if ref_nucleotides_at_this_pos.include?("N") or ref_nucleotides_at_this_pos.include?("n")
        if gap_at_this_pos and nucleotide_at_this_pos
            
            # If some but not all reference seqs have gaps, add to the new seqs, but move on only for those that had gaps.
            ref_indices.each do |i|
                if merged_alignment_seqs[i][positions[i]] == "-"
                    new_merged_alignment_seqs[i] << "-"
                    new_merged_alignment_seqs[i+1] << merged_alignment_seqs[i+1][positions[i+1]]
                    positions[i] += 1
                    positions[i+1] += 1
                else
                    new_merged_alignment_seqs[i] << "-"
                    new_merged_alignment_seqs[i+1] << "-"
                end
            end
        else

            # If all or none of the reference seqs have gaps, just add the nucleotides to the new seqs and move on.
            ref_indices.each do |i|
                if merged_alignment_seqs[i][positions[i]] == nil
                    new_merged_alignment_seqs[i] << "-"
                    new_merged_alignment_seqs[i+1] << "-"
                else
                    new_merged_alignment_seqs[i] << merged_alignment_seqs[i][positions[i]]
                    new_merged_alignment_seqs[i+1] << merged_alignment_seqs[i+1][positions[i+1]]
                    positions[i] += 1
                    positions[i+1] += 1
                end
            end
        end

    end

end # until threading_complete

# Feedback.
puts " done."
print "Checking that sequences are unchanged except gap placement..."
STDOUT.flush

# Make sure that without gaps, all new seqs would be identical to their older versions.
merged_alignment_ids.size.times do |x|
    tmp_old = ""
    merged_alignment_seqs[x].size.times do |pos|
        tmp_old <<  merged_alignment_seqs[x][pos] unless merged_alignment_seqs[x][pos] == "-"
    end
    tmp_new = ""
    new_merged_alignment_seqs[x].size.times do |pos|
        tmp_new << new_merged_alignment_seqs[x][pos] unless new_merged_alignment_seqs[x][pos] == "-"
    end
    unless tmp_old == tmp_new
        puts "ERROR: Sequences appear to have changed!"
        exit 1
    end
end

# Feedback.
puts " done."
print "Checking that reference sequences are now all identical..."
STDOUT.flush

# Replace the old seqs with the new ones.
merged_alignment_ids.size.times do |x|
    merged_alignment_seqs[x] = new_merged_alignment_seqs[x]
end

# Check if all reference seqs are now identical.
ref_seqs = []
ref_indices[1..-1].each do |i|
    ref_seqs << merged_alignment_seqs[i]
end
max_ref_seq_length = ref_seqs[0].size
ref_seqs[1..-1].each do |r|
    max_ref_seq_length = r.size if r.size > max_ref_seq_length
end
first_pos_with_difference = max_ref_seq_length
max_ref_seq_length.times do |pos|
    ref_seq_nucleotides_at_this_pos = []
    ref_seqs.each do |r|
         ref_seq_nucleotides_at_this_pos << r[pos]
    end
    if ref_seq_nucleotides_at_this_pos.uniq.size > 1
        first_pos_with_difference = pos
        break
    end
end
if first_pos_with_difference < max_ref_seq_length
    puts " done. Reference sequences differ from position #{first_pos_with_difference+1} (of #{max_ref_seq_length}) on."
else
    puts " done. (#{first_pos_with_difference}, #{max_ref_seq_length})"
end

# Feedback.
print "Preparing the output..."
STDOUT.flush

# Write the threaded alignment to a new file.
threaded_alignment_file = File.open(threaded_alignment_file_name,"w")
merged_alignment_ids.size.times do |x|
    threaded_alignment_file.write(">#{merged_alignment_ids[x]}\n")
    merged_alignment_seq_as_str = ""
    merged_alignment_seqs[x].each do |nuc|
        merged_alignment_seq_as_str << nuc
    end
    threaded_alignment_file.write("#{merged_alignment_seq_as_str}\n")
end

# Feedback.
puts " done."
puts "Wrote file #{threaded_alignment_file_name}."
