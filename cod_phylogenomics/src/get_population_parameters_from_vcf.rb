# m_matschiner Wed Feb 20 11:34:09 CET 2019

# Get the command-line arguments.
vcf_file_name = ARGV[0]
callable_genome_size = ARGV[1].to_i
table_file_name = ARGV[2]

# Read the input vcf file.
vcf_file = File.open(vcf_file_name)
vcf_lines = vcf_file.readlines

# Read the vcf file.
print "Reading the vcf file #{vcf_file_name}..."
genotypes_per_site = []
vcf_lines_without_header = []
vcf_lines.each do |l|
	unless l[0] == "#"
		vcf_lines_without_header << l
		line_ary = l.split
		genotypes_at_this_site = []
		line_ary[9..-1].each do |gt_with_info|
			genotypes_at_this_site << gt_with_info.split(":")[0].sub("|","/")
		end
		# If any genotypes are half-called, replace them with missing data.
		genotypes_at_this_site.size.times do |x|
			if genotypes_at_this_site[x].count(".") == 1
				genotypes_at_this_site[x] = "./."
			end
		end
		# If there are more 1s than 0s, swap ref and alt1 alleles.
		n_zeros_at_this_site = 0
		n_ones_at_this_site = 0
		genotypes_at_this_site.each do |gt|
			n_zeros_at_this_site += gt.count("0")
			n_ones_at_this_site += gt.count("1")
		end
		if n_ones_at_this_site > n_zeros_at_this_site
			replaced_genotypes_at_this_site = []
			genotypes_at_this_site.each do |gt|
				replaced_genotype_at_this_site = gt.gsub("0","X").gsub("1","0").gsub("X","1")
				if replaced_genotype_at_this_site == "./0" or replaced_genotype_at_this_site == "0/."
					puts "ERROR!"
					puts l
					exit
				end
				replaced_genotype_at_this_site = "0/1" if replaced_genotype_at_this_site == "1/0"
				replaced_genotypes_at_this_site << replaced_genotype_at_this_site
			end
			genotypes_per_site << replaced_genotypes_at_this_site
		else
			genotypes_per_site << genotypes_at_this_site
		end
	end
end
n_samples = genotypes_per_site[0].size
puts " done."

# Analyze the genotypes.
print "Analyzing genotypes..."
site_variable_flags = []
site_biallelic_flags = []
n_gts_heterozygous = 0
pi_per_site = []
genotypes_per_site.each do |gts|

	# Get the non-missing genotypes at this site and determine if it is variable.
	site_variable = false
	uniq_non_missing_gts_at_this_site = []
	gts.uniq.each do |gt|
		uniq_non_missing_gts_at_this_site << gt unless gt == "./."
	end
	if uniq_non_missing_gts_at_this_site.size > 1
		site_variable = true
	elsif uniq_non_missing_gts_at_this_site.size == 1
		alleles_at_this_site = uniq_non_missing_gts_at_this_site[0].split("/")
		site_variable = true if alleles_at_this_site[0] != alleles_at_this_site[1]
	end
	site_variable_flags << site_variable

	# Determine if the site is biallelic.
	site_biallelic = false
	if site_variable
		site_biallelic = true
		alleles_at_this_site = []
		uniq_non_missing_gts_at_this_site.each do |gt|
			gt_ary = gt.split("/")
			alleles_at_this_site << gt_ary[0]
			alleles_at_this_site << gt_ary[1]
		end
		site_biallelic = false if alleles_at_this_site.uniq.size > 2
	end
	site_biallelic_flags << site_biallelic

	# Determine the numbers of heterozygous genotypes at this site.
	gts.each do |gt|
		unless gt.include?(".")
			gt_ary = gt.split("/")
			n_gts_heterozygous += 1 if gt_ary[0] != gt_ary[1]
		end
	end

	# Determine nucleotide diversity for biallelic sites.
	if site_biallelic
		n_zeros = 0
		n_ones = 0
		gts.each do |gt|
			n_zeros += gt.count("0")
			n_ones += gt.count("1")
		end
		# Calculate pi for this site according to Ruegg et al. (2014).
		pi_nominator = n_zeros * n_ones
		# Using the multiplicative formulat of https://en.wikipedia.org/wiki/Binomial_coefficient#Computing_the_value_of_binomial_coefficients,
		# with n = (n_zeros + n_ones) and k = 2 as described in Ruegg et al. (2014):
		# pi_denominator = (((n_zeros + n_ones) + 1 - 1) / 1 ) * (((n_zeros + n_ones) + 1 - 2) / 2.0 )
		pi_denominator = (n_zeros + n_ones) * (n_zeros + n_ones - 1) / 2.0
		pi_at_this_site = pi_nominator / pi_denominator.to_f
		pi_per_site << pi_at_this_site
	elsif site_variable # not biallelic but variable.
		pi_per_site << nil
	elsif uniq_non_missing_gts_at_this_site.size > 0 # not variable but not completely missing.
		pi_per_site << 0
	else
		pi_per_site << nil # completely missing.
	end

end
puts " done."

# Calculate summary statistics.
heterozygosity = n_gts_heterozygous / (n_samples * callable_genome_size).to_f
pi_sum = 0
pi_per_site.size.times do |x|
	unless pi_per_site[x] == nil
		pi_sum += pi_per_site[x]
	end
end
pi = pi_sum / callable_genome_size.to_f

# Make a pool of all non-missing genotypes.
print "Collecting all non-missing genotypes..."
n_missing_genotypes = 0
non_missing_genotypes = []
genotypes_per_site.each do |gts|
	gts.each do |gt|
		if gt.include?(".")
			n_missing_genotypes += 1
		else
			non_missing_genotypes << gt
		end
	end
end
puts " done."

# Make a new set of genotypes in which missing data is replaced randomly with genotypes sampled from the total pool.
print "Replacing missing data with genotypes from the total pool..."
rp_genotypes_per_site = []
genotypes_per_site.each do |gts_at_a_site|
	rp_gts_at_this_site = []
	gts_at_a_site.each do |gt|
		if gt.include?(".")
			rp_gts_at_this_site << non_missing_genotypes.sample
		else
			rp_gts_at_this_site << gt
		end
	end
	rp_genotypes_per_site << rp_gts_at_this_site
end
puts " done."

# Analyze the genotypes.
print "Analyzing replaced genotypes..."
rp_site_variable_flags = []
rp_site_biallelic_flags = []
n_rp_gts_heterozygous = 0
rp_pi_per_site = []
rp_genotypes_per_site.each do |rp_gts|

	# Get the non-missing genotypes at this site and determine if it is variable.
	rp_site_variable = false
	uniq_rp_gts_at_this_site = rp_gts.uniq
	if uniq_rp_gts_at_this_site.size > 1
		rp_site_variable = true
	elsif uniq_rp_gts_at_this_site.size == 1
		rp_alleles_at_this_site = uniq_rp_gts_at_this_site[0].split("/")
		rp_site_variable = true if rp_alleles_at_this_site[0] != rp_alleles_at_this_site[1]
	end
	rp_site_variable_flags << rp_site_variable

	# Determine if the site is biallelic.
	rp_site_biallelic = false
	if rp_site_variable
		rp_site_biallelic = true
		rp_alleles_at_this_site = []
		uniq_rp_gts_at_this_site.each do |rp_gt|
			rp_gt_ary = rp_gt.split("/")
			rp_alleles_at_this_site << rp_gt_ary[0]
			rp_alleles_at_this_site << rp_gt_ary[1]
		end
		rp_site_biallelic = false if rp_alleles_at_this_site.uniq.size > 2
	end
	rp_site_biallelic_flags << rp_site_biallelic

	# Determine the numbers of heterozygous genotypes at this site.
	rp_gts.each do |rp_gt|
		rp_gt_ary = rp_gt.split("/")
		n_rp_gts_heterozygous += 1 if rp_gt_ary[0] != rp_gt_ary[1]
	end

	# Determine nucleotide diversity for biallelic sites.
	if rp_site_biallelic
		n_zeros = 0
		n_ones = 0
		rp_gts.each do |rp_gt|
			n_zeros += rp_gt.count("0")
			n_ones += rp_gt.count("1")
		end
		# Calculate pi for this site according to Ruegg et al. (2014).
		pi_nominator = n_zeros * n_ones
		# Using the multiplicative formulat of https://en.wikipedia.org/wiki/Binomial_coefficient#Computing_the_value_of_binomial_coefficients,
		# with n = (n_zeros + n_ones) and k = 2 as described in Ruegg et al. (2014):
		# pi_denominator = (((n_zeros + n_ones) + 1 - 1) / 1 ) * (((n_zeros + n_ones) + 1 - 2) / 2.0 )
		pi_denominator = (n_zeros + n_ones) * (n_zeros + n_ones - 1) / 2.0
		pi_at_this_site = pi_nominator / pi_denominator.to_f
		rp_pi_per_site << pi_at_this_site
	elsif rp_site_variable # not biallelic but variable.
		rp_pi_per_site << nil
	else # not variable.
		rp_pi_per_site << 0
	end

end
puts " done."

# Calculate summary statistics.
rp_heterozygosity = n_rp_gts_heterozygous / (n_samples * callable_genome_size).to_f
rp_pi_sum = 0
rp_pi_per_site.size.times do |x|
	unless rp_pi_per_site[x] == nil
		rp_pi_sum += rp_pi_per_site[x]
	end
end
rp_pi = rp_pi_sum / callable_genome_size.to_f

# Calculate Watterson's estimator of theta assuming that all missing data is invariant.
n_variants = site_variable_flags.count(true)
theta_denominator = 0
1.upto(2*n_samples - 1) do |i|
	theta_denominator += 1/i.to_f
end
theta_nominator = n_variants
theta = theta_nominator / theta_denominator
theta_per_site = theta / callable_genome_size.to_f

# Calculate Watterson's estimator of theta after replacing missing data.
rp_n_variants = rp_site_variable_flags.count(true)
rp_theta_nominator = rp_n_variants
rp_theta = rp_theta_nominator / theta_denominator
rp_theta_per_site = rp_theta / callable_genome_size.to_f

# Prepare the output string.
outstring = ""
outstring << "Number of sites: #{genotypes_per_site.size}\n"
outstring << "Number of samples: #{n_samples}\n"
outstring << "Callable genome size: #{callable_genome_size}\n"
outstring << "Completeness: #{1-(n_missing_genotypes/(genotypes_per_site.size * n_samples).to_f)}\n"
outstring << "Number of variable sites: #{site_variable_flags.count(true)}\n"
outstring << "Number of biallelic sites: #{site_biallelic_flags.count(true)}\n"
outstring << "Heterozygosity: #{heterozygosity}\n"
outstring << "Nucleotide diversity: #{pi}\n"
outstring << "Watterson's estimator of Theta: #{theta_per_site}\n"
outstring << "Number of variable sites after replacing missing data: #{rp_site_variable_flags.count(true)}\n"
outstring << "Number of biallelic sites after replacing missing data: #{rp_site_biallelic_flags.count(true)}\n"
outstring << "Heterozygosity after replacing missing data: #{rp_heterozygosity}\n"
outstring << "Nucleotide diversity after replacing missing data: #{rp_pi}\n"
outstring << "Watterson's estimator of Theta after replacing missing data: #{rp_theta_per_site}\n"

# Write the output string.
table_file = File.open(table_file_name, "w")
table_file.write(outstring)
puts "Wrote file #{table_file_name}"
