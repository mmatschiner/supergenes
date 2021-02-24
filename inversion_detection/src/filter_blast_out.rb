# m_matschiner Mon Jan 16 23:15:49 CET 2017

# Get the command line arguments.
blast_output_file_name = ARGV[0]

# Read the BLAST output file.
blast_output_file = File.open(blast_output_file_name)
blast_output_lines = blast_output_file.readlines
qseqids = []
sseqids = []
pidents = []
lengths = []
mismatchs = []
gapopens = []
qstarts = []
qends = []
sstarts = []
sends = []
evalues = []
bitscores = []
qseqs = []
blast_output_lines.each do |l|
	line_ary = l.split
	qseqids << line_ary[0]
	sseqids << line_ary[1]
	pidents << line_ary[2]
	lengths << line_ary[3].to_i
	mismatchs << line_ary[4]
	gapopens << line_ary[5]
	qstarts << line_ary[6].to_i
	qends << line_ary[7].to_i
	sstarts << line_ary[8].to_i
	sends << line_ary[9].to_i
	evalues << line_ary[10].to_f
	bitscores << line_ary[11]
	qseqs << line_ary[12]
end

# Add the direction of the hits.
dirs = []
qseqids.size.times do |x|
	if qstarts[x] < qends[x] and sstarts[x] < sends[x]
		dirs << "for"
	elsif qstarts[x] < qends[x] and sstarts[x] > sends[x]
		dirs << "rev"
	else
		raise "ERROR: Unexpected constellation of starts and ends."
	end
end

# Quantify the repetitiveness of reads.
repetitivenesses = []
count = 0
qseqs.each do |seq|
	count += 1
	# To speed up the process, look at maximally the first 250 bp.
	upper_end = [250,seq.size-1].min
	n_diagonals = upper_end
	1.upto(upper_end) do |x|
		(x+1).upto(upper_end) do |y|
			n_diagonals += 2 if seq[x] == seq[y] and seq[x-1] == seq[y-1]
		end
	end
	repetitiveness = n_diagonals/((upper_end)**2).to_f
	repetitivenesses << repetitiveness
end

# Remove repetitive and particularly short sequences.
sseqids.size.times do |x|
	remove_hit = false
	if (qends[x]-qstarts[x]).abs < 100
		remove_hit = true
	elsif (qends[x]-qstarts[x]).abs < 500
		remove_hit = true if repetitivenesses[x] > 0.1
	elsif (qends[x]-qstarts[x]).abs < 2000
		remove_hit = true if repetitivenesses[x] > 0.2
	end
	if remove_hit
		qseqids[x] = nil
		sseqids[x] = nil
        pidents[x] = nil
        lengths[x] = nil
        mismatchs[x] = nil
        gapopens[x] = nil
        qstarts[x] = nil
        qends[x] = nil
        sstarts[x] = nil
        sends[x] = nil
        evalues[x] = nil
        bitscores[x] = nil
        qseqs[x] = nil
        dirs[x] = nil
        repetitivenesses[x] = nil
     end
end
qseqids.compact!
sseqids.compact!
pidents.compact!
lengths.compact!
mismatchs.compact!
gapopens.compact!
qstarts.compact!
qends.compact!
sstarts.compact!
sends.compact!
evalues.compact!
bitscores.compact!
qseqs.compact!
dirs.compact!
repetitivenesses.compact!

# Prepare a filtered version of the BLAST output.
filtered_string = ""
qseqids.size.times do |x|
	filtered_string << "#{qseqids[x]}\t"
	filtered_string << "#{sseqids[x]}\t"
	filtered_string << "#{pidents[x]}\t"
	filtered_string << "#{lengths[x]}\t"
	filtered_string << "#{mismatchs[x]}\t"
	filtered_string << "#{gapopens[x]}\t"
	filtered_string << "#{qstarts[x]}\t"
	filtered_string << "#{qends[x]}\t"
	filtered_string << "#{sstarts[x]}\t"
	filtered_string << "#{sends[x]}\t"
	filtered_string << "#{evalues[x]}\t"
	filtered_string << "#{bitscores[x]}\t"
	filtered_string << "#{dirs[x]}\t"
	filtered_string << "#{repetitivenesses[x]}\n"
end

# Write a filtered list of the BLAST output.
filtered_blast_output_file_name = "#{blast_output_file_name.chomp(".out")}_filtered.out"
filtered_blast_output_file = File.open(filtered_blast_output_file_name,"w")
filtered_blast_output_file.write(filtered_string)
