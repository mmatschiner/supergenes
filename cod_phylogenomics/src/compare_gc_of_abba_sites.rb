# m_matschiner Fri Aug 28 17:27:24 CEST 2020

module Enumerable
	def sum
		self.inject(0){|accum, i| accum + i }
	end
	def mean
		if self.length == 0
			nil
		else
			self.sum/self.length.to_f
		end
	end
	def sample_variance
		if self.length == 0
			nil
		else
			m = self.mean
			sum = self.inject(0){|accum, i| accum +(i-m)**2 }
			sum/(self.length - 1).to_f
		end
	end
	def standard_deviation
		if self.length == 0
			nil
		else
			return Math.sqrt(self.sample_variance)
		end
	end
end

# Get the command-line arguments.
vcf_file_name = ARGV[0]
p1_string = ARGV[1]
p2_string = ARGV[2]
p3_string = ARGV[3]
gc_bbaa_file_name = ARGV[4]
gc_abba_file_name = ARGV[5]
gc_baba_file_name = ARGV[6]

# Parse the taxon strings.
p1s = p1_string.split(",")
p2s = p2_string.split(",")
p3s = p3_string.split(",")

puts "p1s: #{p1s}"
puts "p2s: #{p2s}"
puts "p3s: #{p3s}"

# Read the vcf file.
vcf_file = File.open(vcf_file_name)
vcf_lines = vcf_file.readlines
header_ary = []
gcs = []
# gc_bbaa_weights = []
# gc_abba_weights = []
# gc_baba_weights = []
gc_bbaa = []
gc_abba = []
gc_baba = []
c_bbaas = []
c_abbas = []
c_babas = []
c_aaaas = []
vcf_lines.each do |l|
	if l[0..1] == "##"
		next
	elsif l[0] == "#"
		header_ary = l.split
		p1s.each do |i|
			unless header_ary.include?(i)
				puts "ERROR: Sample #{i} not found in header!"
				exit 1
			end
		end
		p2s.each do |i|
			unless header_ary.include?(i)
				puts "ERROR: Sample #{i} not found in header!"
				exit 1
			end
		end
		p3s.each do |i|
			unless header_ary.include?(i)
				puts "ERROR: Sample #{i} not found in header!"
				exit 1
			end
		end
	else
		# Get counts of bbaa, abba, and baba patterns.
		p1_alleles = []
		p2_alleles = []
		p3_alleles = []
		line_ary = l.split
		p1s.each {|i| p1_alleles << line_ary[header_ary.index(i)].split(":")[0].sub("|","/").split("/")}
		p2s.each {|i| p2_alleles << line_ary[header_ary.index(i)].split(":")[0].sub("|","/").split("/")}
		p3s.each {|i| p3_alleles << line_ary[header_ary.index(i)].split(":")[0].sub("|","/").split("/")}
		p1_alleles.flatten!
		p2_alleles.flatten!
		p3_alleles.flatten!
		p1_alleles.delete(".")
		p2_alleles.delete(".")
		p3_alleles.delete(".")
		n_bbaa = 0
		n_abba = 0
		n_baba = 0
		n_aaaa = 0
		p1_alleles.each do |p1|
			p2_alleles.each do |p2|
				p3_alleles.each do |p3|
					if p1 == p2 and p2 != p3
						n_bbaa += 1
					elsif p2 == p3 and p1 != p2
						n_abba += 1
					elsif p1 == p3 and p1 != p2
						n_baba += 1
					elsif p1 == p2 and p2 == p3
						n_aaaa += 1
					else
						puts "ERROR: Unexpected case (#{p1}, #{p2}, #{p3})!"
					end
				end
			end
		end
		next if n_bbaa + n_abba + n_baba + n_aaaa == 0
		c_bbaa = n_bbaa / (n_bbaa + n_abba + n_baba + n_aaaa).to_f
		c_abba = n_abba / (n_bbaa + n_abba + n_baba + n_aaaa).to_f
		c_baba = n_baba / (n_bbaa + n_abba + n_baba + n_aaaa).to_f
		c_aaaa = n_aaaa / (n_bbaa + n_abba + n_baba + n_aaaa).to_f
		c_bbaas << c_bbaa
		c_abbas << c_abba
		c_babas << c_baba
		c_aaaas << c_aaaa

		# Get the gc_content.
		ref_alt_alleles = [line_ary[3],line_ary[4]]
		# gcs << ref_alt_alleles.count("C") + ref_alt_alleles.count("G") / 2.0
		gc = (ref_alt_alleles.count("C") + ref_alt_alleles.count("G")) / 2.0
		# gc_bbaa_weights << c_bbaa
		# gc_abba_weights << c_abba
		# gc_baba_weights << c_baba
		if c_bbaa > c_abba and c_bbaa > c_baba
			gc_bbaa << gc
		elsif c_abba > c_bbaa and c_abba > c_baba
			gc_abba << gc
		elsif c_baba > c_bbaa and c_baba > c_abba
			gc_baba << gc
		end
	end
end

puts "c_bbaa: #{c_bbaas.sum}"
puts "c_abba: #{c_abbas.sum}"
puts "c_baba: #{c_babas.sum}"
puts "c_aaaa: #{c_aaaas.sum}"

puts "gc_bbaa.mean: #{gc_bbaa.mean}"
puts "gc_bbaa.standard_deviation: #{gc_bbaa.standard_deviation}"
puts "gc_abba.mean: #{gc_abba.mean}"
puts "gc_abba.standard_deviation: #{gc_abba.standard_deviation}"
puts "gc_baba.mean: #{gc_baba.mean}"
puts "gc_baba.standard_deviation: #{gc_baba.standard_deviation}"

gc_bbaa_string = ""
gc_bbaa.each {|gc| gc_bbaa_string << "#{gc}\n"}
gc_bbaa_file = File.open(gc_bbaa_file_name, "w")
gc_bbaa_file.write(gc_bbaa_string)

gc_abba_string = ""
gc_abba.each {|gc| gc_abba_string << "#{gc}\n"}
gc_abba_file = File.open(gc_abba_file_name, "w")
gc_abba_file.write(gc_abba_string)

gc_baba_string = ""
gc_baba.each {|gc| gc_baba_string << "#{gc}\n"}
gc_baba_file = File.open(gc_baba_file_name, "w")
gc_baba_file.write(gc_baba_string)

# gc_bbaa_summed = 0
# gc_bbaa_total_weight = 0
# gc_bbaa_weights.size.times do |x|
# 	gc_bbaa_summed += gcs[x] * gc_bbaa_weights[x]
# 	gc_bbaa_total_weight += gc_bbaa_weights[x]
# end
# puts gc_bbaa_summed/gc_bbaa_total_weight
# gc_abba_summed = 0
# gc_abba_total_weight = 0
# gc_abba_weights.size.times do |x|
# 	gc_abba_summed += gcs[x] * gc_abba_weights[x]
# 	gc_abba_total_weight += gc_abba_weights[x]
# end
# puts gc_abba_summed/gc_abba_total_weight
# gc_baba_summed = 0
# gc_baba_total_weight = 0
# gc_baba_weights.size.times do |x|
# 	gc_baba_summed += gcs[x] * gc_baba_weights[x]
# 	gc_baba_total_weight += gc_baba_weights[x]
# end
# puts gc_baba_summed/gc_baba_total_weight
