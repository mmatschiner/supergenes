# m_matschiner Sat Oct 9 10:20:02 CEST 2021

pop_size = ARGV[0].to_i
n_simulations = ARGV[1].to_i
output_filename = ARGV[2]

n_generations_to_fixation = []

n_simulations.times do
	alleles = []
	pop_size.times {alleles << 0}
	alleles[0] = 1
	fixed = false
	n_generations = 0
	until fixed do
		n_generations += 1
		new_alleles = []
		pop_size.times {new_alleles << alleles.sample}
		alleles = new_alleles
		fixed = true if alleles.count(0) == 0
		fixed = true if alleles.count(1) == 0
	end
	n_generations_to_fixation << n_generations
end
puts n_generations_to_fixation