# m_matschiner Fri Jul 10 13:47:52 CEST 2020

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
p_bottleneck_file_name = ARGV[0]
pi_file_name = ARGV[1]

# Read both files.
p_bottleneck_file = File.open(p_bottleneck_file_name)
p_bottleneck_lines = p_bottleneck_file.readlines[1..-1]
pi_file = File.open(pi_file_name)
pi_lines = pi_file.readlines[1..-1]

# Get the unique population sizes and bottleneck times.
bottleneck_times = []
pop_sizes = []
p_bottleneck_lines.each do |l|
	line_ary = l.split
	bottleneck_times << line_ary[0].to_i
	pop_sizes << line_ary[1].to_i
end
uniq_bottleneck_times = bottleneck_times.uniq
uniq_pop_sizes = pop_sizes.uniq

# Write summaries for each combination of population size and bottleneck time.
puts "Ne\tt_bottleneck\tp_bottleneck\tpi1\tpi2"
uniq_pop_sizes.each do |uniq_pop_size|
	uniq_bottleneck_times.each do |uniq_bottleneck_time|
		ps_after_bottleneck = []
		p_bottleneck_lines.each do |l|
			line_ary = l.split
			if line_ary[0].to_i == uniq_bottleneck_time and line_ary[1].to_i == uniq_pop_size
				ps_after_bottleneck << 1-line_ary[4].to_f
			end
		end
		pi1s = []
		pi2s = []
		pi_lines.each do |l|
			line_ary = l.split
			if line_ary[0].to_i == uniq_bottleneck_time and line_ary[1].to_i == uniq_pop_size
				pi1s << line_ary[4].to_f
				pi2s << line_ary[5].to_f
			end
		end
		puts "#{uniq_pop_size}\t#{uniq_bottleneck_time}\t#{'%.2f' % ps_after_bottleneck.mean} (#{'%.2f' % ps_after_bottleneck.standard_deviation})\t#{'%.3f' % (1000*pi2s.mean)} (#{'%.3f' % (1000*pi2s.standard_deviation)})\t#{'%.3f' % (1000*pi1s.mean)} (#{'%.3f' % (1000*pi1s.standard_deviation)})"
	end
end
