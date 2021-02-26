# Michael Matschiner, 2015-06-04
# 
# This script calculates the mean coverage from
# files produced with get_coverage_distribution.sh
# One or more filenames of coverage distribution files
# should be given as command line arguments.
# Start e.g with
# ruby get_coverage_stats.rb *.covdist.txt

ids = []
id_max_length = 0
medians = []
means = []
ARGV.each do |a|
	id = "#{a.sub(".merged.sorted.dedup.realn.covdist.txt","").split("/").last}"
	ids << id
	id_max_length = id.size if id.size > id_max_length
	sum = 0
	cumulative_freq = 0
	lines = File.open(a).readlines
	lines.each do |l|
		cov = l.split[1].to_i
		freq = l.split[2].to_i
		sum += cov * freq
		cumulative_freq += freq
	end
	total_freq = cumulative_freq
	means << (sum/total_freq.to_f).round(3)
	cumulative_freq = 0
	lines.each do |l|
		cov = l.split[1].to_i
		freq = l.split[2].to_i
		cumulative_freq += freq
		if cumulative_freq > total_freq/2.0
			medians << cov
			break
		end
	end
end

outstring = "#{"id".ljust(id_max_length+3)}#{"mean".rjust(10)}#{"median".rjust(10)}\n"
ids.size.times do |x|
	outstring << "#{ids[x].ljust(id_max_length+3)}#{("%5.2f" % means[x]).rjust(10)}#{medians[x].to_s.rjust(10)}\n"
end

puts outstring