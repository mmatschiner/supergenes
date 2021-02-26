# m_matschiner Mon Feb 22 11:44:18 CET 2021

# See https://en.wikipedia.org/wiki/Binomial_coefficient
def n_choose_k(n, k)
	if k > n
		puts "ERROR: k should never be larger than n!"
		exit 1
	elsif k == 0
		return 1
	elsif k > n/2
		return n_choose_k(n, n-k)
	end
	r = 1
	1.upto(k) do |j|
		r = r * (n + 1 - j) / j
	end
	r
end

# Get the command-line argument.
vcf_file_name = ARGV[0]

# Read the vcf file line by line.
sample_ids = []
biases_per_sample = [] 
File.open(vcf_file_name).each do |l|
	if l[0..1] == "##"
		next
	elsif l[0] == "#"
		sample_ids = l.split[9..-1]
		sample_ids.size.times do
			biases_per_sample << []
		end
	else
		if sample_ids == []
			puts l
			puts "ERROR: Sample IDs not found!"
			exit 1
		end
		line_ary = l.split
		sample_ids.size.times do |x|
			gt = line_ary[x+9].split(":")[0]
			if gt == "0/1"
				ad = line_ary[x+9].split(":")[1]
				n_ref_alleles = ad.split(",")[0].to_i
				n_alt_alleles = ad.split(",")[1].to_i
				if n_ref_alleles > n_alt_alleles
					biases_per_sample[x] << "ref"
				elsif n_alt_alleles > n_ref_alleles
					biases_per_sample[x] << "alt"
				else
					biases_per_sample[x] << "equ"
				end
			end
		end
	end
end

# Report, per sample, the numbers of heterozygous sites with larger numbers of reference alleles, larger numbers of alternate alleles, and equal numbers of both alleles.
outstring = "sample_id	ref_bias	no_bias	alt_bias	p\n"
sample_ids.size.times do |x|

	# Get biases for this sample.
	ref_bias = biases_per_sample[x].count("ref")
	alt_bias = biases_per_sample[x].count("alt")
	no_bias = biases_per_sample[x].count("equ")

	# Estimate the p value by bootstrapping, this might be faster than the binomial test.
	bs_more_extreme = 0
	1000.times do
		ref_bias_bs = 0
		alt_bias_bs = 0
		(ref_bias + alt_bias).times do
			if rand > 0.5
				ref_bias_bs += 1
			else
				alt_bias_bs += 1
			end
		end
		bs_more_extreme += 1 if ref_bias_bs - alt_bias_bs > ref_bias - alt_bias
	end
	pr = bs_more_extreme / 1000.0


	# # Test if the imbalance between discordant trees is greater than expected by chance.
	# pr = 0
	# n = ref_bias + alt_bias
	# k = ([ref_bias,alt_bias].min+1)
	# puts "n: #{n}, k: #{k}"
	# k.times do |i|
	# 	pr += Math.exp(Math.log(n_choose_k(n,i)) + Math.log(0.5**i) + Math.log((1-0.5)**(n-i)))
	# end
	# pr = pr * 2 # Make the test two-tailed instead of one-tailed.
	# pr = 1.0 if pr > 1.0 # This is possible when both alternative topologies are equally frequent.

	outstring << "#{sample_ids[x]}	#{ref_bias}	#{no_bias}	#{alt_bias}	#{('%.3f' % pr).rjust(12)}\n"
end
puts outstring
