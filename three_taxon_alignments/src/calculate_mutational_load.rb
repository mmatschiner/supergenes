# m_matschiner Mon Dec 7 16:08:59 CET 2020

# Define a class for exon information.
class Exon
	attr_reader :gene_id, :strand, :from, :to
	def initialize(from, to, strand, gene_id)
		@from = from
		@to = to
		@strand = strand
		@gene_id = gene_id
	end
end

# Define a class for gene information.
class Gene
	attr_reader :id, :center, :aa_seqs
	def initialize(id, center, aa_seqs)
		@id = id
		@center = center
		@aa_seqs = aa_seqs
	end
end

# Get the command-line arguments.
fasta_file_name = ARGV[0]
lg_id = ARGV[1].upcase
gff_file_name = ARGV[2]
window_size = ARGV[3].to_i
out_file_name = ARGV[4]

# Read the alignment file.
fasta_file = File.open(fasta_file_name)
fasta_lines = fasta_file.readlines
fasta_ids = []
fasta_seqs = []
3.times do |x|
	fasta_ids << fasta_lines[x*2][1..-1].strip.upcase
	fasta_seqs << fasta_lines[x*2+1].strip.upcase
end
screen_out = "Read file #{fasta_file_name}."
puts screen_out
STDOUT.flush

# Write screen output.
screen_out = "Sequence IDs are"
fasta_ids.each do |i|
	screen_out << " #{i},"
end
screen_out.chomp!(",")
screen_out << ".\n"
screen_out << "The first sequence with ID #{fasta_ids[0]} is assumed to be the reference."
puts screen_out
STDOUT.flush

# Read the annotation file.
exons = []
gene_ids = []
chr_length = nil
gff_file = File.open(gff_file_name)
gff_file.each do |l|
	line_ary = l.split
	if line_ary[0] == lg_id
		if line_ary[2] == "CDS"
			from = line_ary[3].to_i
			to = line_ary[4].to_i
			strand = line_ary[6]
			gene_id = line_ary[8].split(":")[0].split("=")[1]
			exons << Exon.new(from, to, strand, gene_id)
			gene_ids << gene_id
		else
			if line_ary[2] == "contig"
				chr_length = line_ary[4].to_i
			end
		end
	end
end
# Make sure a chromosome length is found.
if chr_length == nil
	puts "ERROR: Chromosome length could not be determined."
	exit 1
end
# Make sure the chromosome length matches the alignment length.
unless fasta_seqs[0].size == chr_length
	puts "ERROR: Chromosome length differs from alignment length!"
	puts "    Perhaps the annotation does not fit to this sequence alignment."
	exit 1
end

# Reduce gene ids to unique ones.
uniq_gene_ids = gene_ids.uniq

# Compose the nucleotide sequences.
genes = []
uniq_gene_ids.each do |g|
	gene_seqs = []
	exons_this_gene = []
	exons.each do |e|
		if e.gene_id == g
			exons_this_gene << e
		end
	end
	# Check that all exons of this gene have the same strand.
	all_exons_same_strand = true
	gene_strand = exons_this_gene[0].strand
	exons_this_gene.each do |e|
		all_exons_same_strand = false unless e.strand == gene_strand
	end
	unless all_exons_same_strand
		puts "ERROR: Different strands for gene #{g}"
		exit 1
	end
	fasta_seqs.size.times do |x|
		gene_seqs[x] = ""
		exons_this_gene.each do |e|
			gene_seqs[x] << fasta_seqs[x][(e.from-1)..(e.to-1)]
		end
	end
	if gene_strand == "-"
		# Make the reverse complement.
		gene_seqs.size.times do |x|
			rev_comp_seq = ""
			(gene_seqs[x].size-1).downto(0) do |pos|
				if gene_seqs[x][pos] == "A"
					rev_comp_seq << "T"
				elsif gene_seqs[x][pos] == "C"
					rev_comp_seq << "G"
				elsif gene_seqs[x][pos] == "G"
					rev_comp_seq << "C"
				elsif gene_seqs[x][pos] == "T"
					rev_comp_seq << "A"
				elsif gene_seqs[x][pos] == "N"
					rev_comp_seq << "N"
				elsif gene_seqs[x][pos] == "-"
					rev_comp_seq << "-"
				else
					puts "ERROR: Unknown nucleotide #{gene_seqs[x][pos]} in sequence #{x} for gene #{g}!"
					exit 1
				end
			end
			gene_seqs[x] = rev_comp_seq
		end
	elsif gene_strand != "+"
		puts "ERROR: Unknown strand #{gene_strand} for gene #{g}!"
		exit 1
	end
	# Make sure that all gene sequences have a length that is a multiple of three.
	unless (gene_seqs[0].size/3)*3 == gene_seqs[0].size
		puts "ERROR: Sequence length of gene #{g} (#{gene_seqs[0].size}) is not a multiple of 3!"
		exit 1
	end
	# Only use genes with sequences that start with a start codon.
	if gene_seqs[0][0..2] == "ATG"
		# Only use genes if all sequences have at least some 12 nonmissing sites.
		n_non_missing_sites = []
		gene_seqs.each do |s|
			n_non_missing_sites << s.count("A") + s.count("C") + s.count("G") + s.count("T")
		end
		if n_non_missing_sites.min > 12
			# Get the first and last position of exons of this gene, and calculate the center.
			gene_first = 100000000
			gene_last = 0
			exons_this_gene.each do |e|
				gene_first = e.from if e.from < gene_first
				gene_last = e.to if e.to > gene_last
			end
			if gene_first == 100000000 or gene_last == 0
				puts "ERROR: Gene location not found!"
				exit 1
			end
			gene_center = (gene_first + gene_last)/2.0
			# Translate sequences into amino-acid sequences.
			aa_seqs = []
			gene_seqs.each do |s|
				aa_seq = ""
				0.upto((s.size/3)-1) do |aa_pos|
					codon = s[(aa_pos*3)..((aa_pos*3)+2)]
					if codon == "TTT" or codon == "TTC"
						aa_seq << "P"
					elsif codon == "TTA" or codon == "TTG" or codon == "CTT" or codon == "CTC" or codon == "CTA" or codon == "CTG"
						aa_seq << "L"
					elsif codon == "ATT" or codon == "ATC" or codon == "ATA"
						aa_seq << "I"
					elsif codon == "ATG"
						aa_seq << "M"
					elsif codon == "GTT" or codon == "GTC" or codon == "GTA" or codon == "GTG"
						aa_seq << "V"
					elsif codon == "TCT" or codon == "TCC" or codon == "TCA" or codon == "TCG"
						aa_seq << "S"
					elsif codon == "CCT" or codon == "CCC" or codon == "CCA" or codon == "CCG"
						aa_seq << "P"
					elsif codon == "ACT" or codon == "ACC" or codon == "ACA" or codon == "ACG"
						aa_seq << "T"
					elsif codon == "GCT" or codon == "GCC" or codon == "GCA" or codon == "GCG"
						aa_seq << "A"
					elsif codon == "TAT" or codon == "TAC"
						aa_seq << "Y"
					elsif codon == "TAA" or codon == "TAG"
						aa_seq << "*"
					elsif codon == "CAT" or codon == "CAC"
						aa_seq << "H"
					elsif codon == "CAA" or codon == "CAG"
						aa_seq << "Q"
					elsif codon == "AAT" or codon == "AAC"
						aa_seq << "N"
					elsif codon == "AAA" or codon == "AAG"
						aa_seq << "K"
					elsif codon == "GAT" or codon == "GAC"
						aa_seq << "D"
					elsif codon == "GAA" or codon == "GAG"
						aa_seq << "E"
					elsif codon == "TGT" or codon == "TGC"
						aa_seq << "C"
					elsif codon == "TGA"
						aa_seq << "*"
					elsif codon == "TGG"
						aa_seq << "W"
					elsif codon == "CGT" or codon == "CGC" or codon == "CGA" or codon == "CGG"
						aa_seq << "R"
					elsif codon == "AGT" or codon == "AGC"
						aa_seq << "S"
					elsif codon == "AGA" or codon == "AGG"
						aa_seq << "R"
					elsif codon == "GGT" or codon == "GGC" or codon == "GGA" or codon == "GGG"
						aa_seq << "G"
					elsif codon.include?("N") or codon.include?("-")
						aa_seq << "X"
					else
						puts "Unexpected codon #{codon}!"
						exit 1
					end
				end
				aa_seqs << aa_seq
			end
			genes << Gene.new(g, gene_center, aa_seqs)
		end
	end
end

# Check if any gene sequences have stop codons.
n_genes_with_stop_codons_in_ref = 0
genes.each do |g|
	n_genes_with_stop_codons_in_ref += 1 if g.aa_seqs[0][0..-2].include?("*")
end
if n_genes_with_stop_codons_in_ref > 0
	puts "ERROR: #{n_genes_with_stop_codons_in_ref} genes contain stop codons in the reference!"
	exit 1
end

# Analyze genes in sliding windows.
file_out = "chr\tstart\tend"
fasta_ids.size.times do |x|
	file_out << "\tn_aas_#{fasta_ids[x]}"
end
fasta_ids.size.times do |x|
	file_out << "\tn_stops_#{fasta_ids[x]}"
end
fasta_ids.size.times do |x|
	fasta_ids.size.times do |y|
		if y > x
			file_out << "\tn_aas_unchanged_#{fasta_ids[x]}_#{fasta_ids[y]}\tn_aas_changed_#{fasta_ids[x]}_#{fasta_ids[y]}"
		end
	end
end
file_out << "\n"
window_from = 1
window_to = window_from + window_size - 1
until window_to > chr_length
	n_aas_this_window = []
	n_stop_codons_this_window = []
	n_aas_changed_this_window = []
	n_aas_unchanged_this_window = []
	fasta_ids.size.times do |x|
		n_aas_this_window[x] = 0
		n_stop_codons_this_window[x] = 0
		n_aas_changed_this_window[x] = []
		n_aas_unchanged_this_window[x] = []
		fasta_ids.size.times do |y|
			if y > x
				n_aas_changed_this_window[x][y] = 0
				n_aas_unchanged_this_window[x][y] = 0
			end
		end
	end
	genes.each do |g|
		if g.center >= window_from and g.center <= window_to
			fasta_ids.size.times do |x|
				(g.aa_seqs[0].size-1).times do |pos| # the very last position is skipped, it may be a stop codon.
					n_aas_this_window[x] += 1 if g.aa_seqs[x][pos] != "X"
					n_stop_codons_this_window[x] += 1 if g.aa_seqs[x][pos] == "*"
					fasta_ids.size.times do |y|
						if y > x
							if g.aa_seqs[x][pos] != "X" and g.aa_seqs[y][pos] != "X"
								if g.aa_seqs[x][pos] == g.aa_seqs[y][pos]
									n_aas_unchanged_this_window[x][y] += 1
								else
									n_aas_changed_this_window[x][y] += 1
								end
							end
						end
					end
				end
			end
		end
	end
	file_out << "#{lg_id.upcase}\t#{window_from}\t#{window_to}"
	fasta_ids.size.times do |x|
		file_out << "\t#{n_aas_this_window[x]}"
	end
	fasta_ids.size.times do |x|
		file_out << "\t#{n_stop_codons_this_window[x]}"
	end
	fasta_ids.size.times do |x|
		fasta_ids.size.times do |y|
			if y > x
				file_out << "\t#{n_aas_unchanged_this_window[x][y]}\t#{n_aas_changed_this_window[x][y]}"
			end
		end
	end
	file_out << "\n"
	window_from = window_from + window_size
	window_to = window_from + window_size - 1
end

# Write the output file.
out_file = File.open(out_file_name,"w")
out_file.write(file_out)

# Write screen output.
puts "Wrote file #{out_file_name}."
