# m_matschiner Mon Oct 15 19:49:12 CEST 2018

# Load required libraries.
require 'optparse'

# Feedback.
puts ""
puts "aim_prep.rb"
puts ""
puts "----------------------------------------------------------------------------------------"
puts ""

# Define default options.
options = {}
options[:dir] = nil
options[:outgroup] = nil
options[:age] = nil
options[:length] = 1000000000
options[:xml] = "aim.xml"

# Get the command line options.
ARGV << '-h' if ARGV.empty?
opt_parser = OptionParser.new do |opt|
	opt.banner = "Usage: ruby #{$0} [OPTIONS]"
	opt.separator  ""
	opt.separator  "Example:"
	opt.separator  "ruby #{$0} -d alignments -o spc1 -a 10.0 -l #{options[:length]} -x #{options[:xml]}"
	opt.separator  ""
	opt.separator  "Options:"
	opt.on("-d","--dir NAME","Directory with alignment files (default: none).") {|d| options[:dir] = d}
	opt.on("-o","--outgroup NAME","Name of outgroup species (default: none).") {|o| options[:outgroup] = o}
	opt.on("-a","--age STRING","Root age of the species tree (default: none).") {|a| options[:age] = a}
	opt.on("-l","--length INT","Number of MCMC generations (default: #{options[:length]}).") {|l| options[:length] = l}
	opt.on("-x","--xml FILENAME","Output file in XML format (default: #{options[:xml]}).") {|x| options[:xml] = x}
	opt.on("-h","--help","Print this help text.") {
		puts opt_parser
		exit(0)
	}
	opt.separator  ""
end
opt_parser.parse!

# Get the command-line arguments.
align_dir = options[:dir]
xml_file_name = options[:xml]
xml_id = xml_file_name.split("/").last.chomp(".xml")
outgroup_id = options[:outgroup]
root_age = options[:age]
chain_length = options[:length].to_i
log_every = chain_length/1000

# Get alignment file names and format.
align_file_names = []
align_format = ""
align_ids = []
align_lengths = []
Dir.entries(align_dir).each do |e|
	if e.match(/.nex/)
		align_file_names << e
		align_ids << e.split("/").last.chomp(".nex")
		if align_format == ""
			align_format = "nex"
		elsif align_format != "nex"
			puts "ERROR: Found files with mixed formats: #{align_format}, nex!"
			exit(1)
		end
	elsif e.match(/.fasta/)
		align_file_names << e
		align_ids << e.split("/").last.chomp(".fasta")
		if align_format == ""
			align_format = "fas"
		elsif align_format != "fas"
			puts "ERROR: Found files with mixed formats: #{align_format}, fas!"
			exit(1)
		end
	elsif e.match(/.phy/)
		align_file_names << e
		align_ids << e.split("/").last.chomp(".phy")
		if align_format == ""
			align_format = "phy"
		elsif align_format != "phy"
			puts "ERROR: Found files with mixed formats: #{align_format}, phy!"
			exit(1)
		end
	end
end

# Initiate the xml string.
first_ids = []
xml_string = ""
xml_string << "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"?><beast beautitemplate='AIM' beautistatus='noAutoSetClockRate' namespace=\"beast.core:beast.evolution.alignment:beast.evolution.tree.coalescent:beast.core.util:beast.evolution.nuc:beast.evolution.operators:beast.evolution.sitemodel:beast.evolution.substitutionmodel:beast.evolution.likelihood\" required=\"BEAST v2.5.0:starbeast2 v0.15.0\" version=\"2.5\">\n\n"
align_file_names.size.times do |x|
	# Get the ids and sequences from the alignment.
	align_file = File.open("#{align_dir}/#{align_file_names[x]}")
	align_lines = align_file.readlines
	ids = []
	seqs = []
	if align_format == "nex"
		in_matrix = false
		align_lines.each do |l|
			if l.match(/matrix/)
				in_matrix = true
			elsif l.match(/end;/) or l.match(/END;/) or l.match(/End;/)
				in_matrix = false
			elsif l.strip != ""
				ids << l.strip.split[0]
				seqs << l.strip.split[1]
			end
		end
	elsif align_format == "fas"
		align_lines.each do |l|
			if l[0] == ">"
				ids << l[1..-1].strip
				seqs << ""
			elsif l.strip != ""
				seqs.last << l.strip
			end
		end
	elsif align_format == "phy"
		align_lines[1..-1].each do |l|
			ids << l.strip.split[0]
			seqs << l.strip.split[1]
		end
	end
	first_ids = ids if first_ids == []
	unless first_ids == ids
		puts "ERROR: Taxon IDs seem to differ between alignment files!"
		exit(1)
	end
	xml_string << "\t<data id=\"#{align_ids[x]}\" name=\"alignment\">\n"
	ids.size.times do |y|
		xml_string << "\t\t<sequence id=\"#{align_ids[x]}_#{ids[y]}\" taxon=\"#{ids[y]}\" totalcount=\"4\" value=\"#{seqs[y]}\"/>\n"
	end
	xml_string << "\t</data>\n"
	align_lengths << seqs[0].size
end
ids = first_ids

# Ensure that the outgroup is among the ids.
unless ids.include?(outgroup_id)
	puts "ERROR: The outgroup ID #{outgroup_id} is not included in alignment IDs!"
	exit(1)
end

# Add maps to the xml.
xml_string << "\n"
xml_string << "\t<map name=\"Uniform\" >beast.math.distributions.Uniform</map>\n"
xml_string << "\t<map name=\"Exponential\" >beast.math.distributions.Exponential</map>\n"
xml_string << "\t<map name=\"LogNormal\" >beast.math.distributions.LogNormalDistributionModel</map>\n"
xml_string << "\t<map name=\"Normal\" >beast.math.distributions.Normal</map>\n"
xml_string << "\t<map name=\"Beta\" >beast.math.distributions.Beta</map>\n"
xml_string << "\t<map name=\"Gamma\" >beast.math.distributions.Gamma</map>\n"
xml_string << "\t<map name=\"LaplaceDistribution\" >beast.math.distributions.LaplaceDistribution</map>\n"
xml_string << "\t<map name=\"prior\" >beast.math.distributions.Prior</map>\n"
xml_string << "\t<map name=\"InverseGamma\" >beast.math.distributions.InverseGamma</map>\n"
xml_string << "\t<map name=\"OneOnX\" >beast.math.distributions.OneOnX</map>\n"
xml_string << "\n"

# Initiate the run element.
xml_string << "\t<run id=\"mcmc\" spec=\"beast.core.MCMC\" chainLength=\"#{chain_length}\">\n"

# Write the state element.
xml_string << "\t\t<state id=\"state\" storeEvery=\"#{log_every}\">\n"
xml_string << "\t\t\t<stateNode id=\"Tree.t:Species\" spec=\"starbeast2.SpeciesTree\">\n"
xml_string << "\t\t\t\t<taxonset id=\"taxonsuperset\" spec=\"starbeast2.StarBeastTaxonSet\">\n"
ids.size.times do |x|
	xml_string << "\t\t\t\t\t<taxon id=\"#{ids[x]}_spc\" spec=\"TaxonSet\">\n"
	xml_string << "\t\t\t\t\t\t<taxon id=\"#{ids[x]}\" spec=\"Taxon\"/>\n"
	xml_string << "\t\t\t\t\t</taxon>\n"
end
xml_string << "\t\t\t\t</taxonset>\n"
xml_string << "\t\t\t</stateNode>\n"
xml_string << "\t\t\t<parameter id=\"popSizes.Species\" name=\"stateNode\" lower=\"0.001\" upper=\"100\">1</parameter>\n"
xml_string << "\t\t\t<parameter id=\"migRates.Species\" dimension=\"1\" upper=\"0.5\" name=\"stateNode\">0.1</parameter>\n"
xml_string << "\t\t\t<stateNode id=\"migIndicators.Species\" spec=\"parameter.BooleanParameter\" dimension=\"0\">false</stateNode>\n"
xml_string << "\t\t\t<parameter id=\"migMean.Species\" lower=\"0.1\" upper=\"100.0\" name=\"stateNode\">1</parameter>\n"
xml_string << "\t\t\t<parameter id=\"popMean.Species\"  name=\"stateNode\">0.005</parameter>\n"
xml_string << "\t\t\t<parameter id=\"speciationRate.t:Species\" lower=\"0.0\" name=\"stateNode\">1.0</parameter>\n"
xml_string << "\t\t\t<parameter id=\"strictClockRate\" lower=\"0.0\" name=\"stateNode\">1</parameter>\n"
align_ids.size.times do |x|
	xml_string << "\t\t\t<tree id=\"Tree.t:#{align_ids[x]}\" name=\"stateNode\">\n"
	xml_string << "\t\t\t\t<taxonset id=\"TaxonSet.#{align_ids[x]}\" spec=\"TaxonSet\">\n"
	xml_string << "\t\t\t\t\t<alignment idref=\"#{align_ids[x]}\"/>\n"
	xml_string << "\t\t\t\t</taxonset>\n"
	xml_string << "\t\t\t</tree>\n"
end
align_ids.size.times do |x|
	xml_string << "\t\t\t<parameter id=\"kappa.s:#{align_ids[x]}\" lower=\"0.0\" name=\"stateNode\">2.0</parameter>\n"
end
align_ids.size.times do |x|
	xml_string << "\t\t\t<parameter id=\"mutationRate.s:#{align_ids[x]}\" lower=\"0.0\" name=\"stateNode\">1.0</parameter>\n"
end
xml_string << "\t\t</state>\n"
xml_string << "\n"

# Write the init element.
xml_string << "\t\t<init id=\"SBI\" spec=\"starbeast2.StarBeastInitializer\" birthRate=\"@speciationRate.t:Species\" estimate=\"false\" speciesTree=\"@Tree.t:Species\">\n"
align_ids.size.times do |x|
	xml_string << "\t\t\t<geneTree idref=\"Tree.t:#{align_ids[x]}\"/>\n"
end
xml_string << "\t\t\t<populationModel id=\"popModelAIM.Species\" spec=\"starbeast2.ConstantWithGeneFlow\" Ne=\"@popSizes.Species\" NeMean=\"@popMean.Species\" indicator=\"@migIndicators.Species\" m=\"@migRates.Species\">\n"
xml_string << "\t\t\t\t<migrationModel id=\"migModel.Species\" spec=\"starbeast2.Overlap\" minimalBranchLength=\"0.00001\" exclude=\"#{outgroup_id}_spc\" effectiveMigrants=\"@migMean.Species\" speciesTree=\"@Tree.t:Species\"/>\n"
xml_string << "\t\t\t</populationModel>\n"
xml_string << "\t\t</init>\n"
xml_string << "\n"

# Write the posterior element.
xml_string << "\t\t<distribution id=\"posterior\" spec=\"util.CompoundDistribution\">\n"
xml_string << "\t\t\t<distribution id=\"speciescoalescent\" spec=\"util.CompoundDistribution\">\n"
align_ids.size.times do |x|
	xml_string << "\t\t\t\t<distribution id=\"geneTree.t:#{align_ids[x]}\" spec=\"starbeast2.GeneTreeWithMigration\" ploidy=\"2.0\" populationModel=\"@popModelAIM.Species\" tree=\"@Tree.t:#{align_ids[x]}\"/>\n"
end
xml_string << "\t\t\t</distribution>\n"
xml_string << "\n"
xml_string << "\t\t\t<distribution id=\"prior\" spec=\"util.CompoundDistribution\">\n"
xml_string << "\t\t\t\t<distribution id=\"all.prior\" monophyletic=\"true\"  spec=\"beast.math.distributions.MRCAPrior\" tree=\"@Tree.t:Species\">\n"
xml_string << "\t\t\t\t\t<taxonset id=\"all\" spec=\"TaxonSet\">\n"
ids.size.times do |x|
	xml_string << "\t\t\t\t\t\t<taxon idref=\"#{ids[x]}_spc\" spec=\"Taxon\"/>\n"
end
xml_string << "\t\t\t\t\t</taxonset>\n"
if root_age.match(/lognormal\(([0-9\.]+),([0-9\.]+)\)/)
	xml_string << "\t\t\t\t\t<LogNormal id=\"root_age\" meanInRealSpace=\"true\" name=\"distr\">\n"
	xml_string << "\t\t\t\t\t\t<parameter estimate=\"false\" name=\"M\">#{$1}</parameter>\n"
	xml_string << "\t\t\t\t\t\t<parameter estimate=\"false\" name=\"S\">#{$2}</parameter>\n"
	xml_string << "\t\t\t\t\t</LogNormal>\n"
elsif root_age.match(/[0-9\.]+/)
	xml_string << "\t\t\t\t\t<Uniform id=\"root_age\" name=\"distr\" lower=\"#{root_age.to_f-0.01}\" upper=\"#{root_age.to_f+0.01}\"/>\n"
else
	puts "ERROR: Root age could not be interpreted!"
	exit(1)
end
xml_string << "\t\t\t\t</distribution>\n"
xml_string << "\n"
xml_string << "\t\t\t\t<distribution id=\"ingroup.prior\" monophyletic=\"true\" spec=\"beast.math.distributions.MRCAPrior\" tree=\"@Tree.t:Species\">\n"
xml_string << "\t\t\t\t\t<taxonset id=\"ingroup\" spec=\"TaxonSet\">\n"
ids.size.times do |x|
	xml_string << "\t\t\t\t\t\t<taxon idref=\"#{ids[x]}_spc\" spec=\"Taxon\"/>\n" unless ids[x] == outgroup_id
end
xml_string << "\t\t\t\t\t</taxonset>\n"
xml_string << "\t\t\t\t</distribution>\n"
xml_string << "\n"
xml_string << "\t\t\t\t<distribution id=\"YuleModel.t:Species\" spec=\"beast.evolution.speciation.YuleModel\" birthDiffRate=\"@speciationRate.t:Species\" tree=\"@Tree.t:Species\"/>\n"
align_ids.size.times do |x|
	xml_string << "\t\t\t\t<prior id=\"KappaPrior.s:#{align_ids[x]}\" name=\"distribution\" x=\"@kappa.s:#{align_ids[x]}\">\n"
	xml_string << "\t\t\t\t\t<LogNormal id=\"LogNormalDistributionModel.kappa.#{align_ids[x]}\" name=\"distr\">\n"
	xml_string << "\t\t\t\t\t\t<parameter id=\"RealParameter.M.#{align_ids[x]}\" estimate=\"false\" name=\"M\">1.0</parameter>\n"
	xml_string << "\t\t\t\t\t\t<parameter id=\"RealParameter.S.#{align_ids[x]}\" estimate=\"false\" name=\"S\">1.25</parameter>\n"
	xml_string << "\t\t\t\t\t</LogNormal>\n"
	xml_string << "\t\t\t\t</prior>\n"
end
xml_string << "\t\t\t\t<prior id=\"migIndicatorSumPrior.Species\" name=\"distribution\">\n"
xml_string << "\t\t\t\t\t<x id=\"migIndicatorSum.Species\" spec=\"util.Sum\">\n"
xml_string << "\t\t\t\t\t\t<arg idref=\"migIndicators.Species\"/>\n"
xml_string << "\t\t\t\t\t</x>\n"
xml_string << "\t\t\t\t\t<distr id=\"Poisson.0\" spec=\"beast.math.distributions.Poisson\">\n"
xml_string << "\t\t\t\t\t\t<parameter id=\"RealParameter.105\" name=\"lambda\">0.695</parameter>\n" # XXX TODO Check the origin of this number.
xml_string << "\t\t\t\t\t</distr>\n"
xml_string << "\t\t\t\t</prior>\n"
xml_string << "\t\t\t\t<prior id=\"migIndicatorPrior.Species\" name=\"distribution\" x=\"@migIndicators.Species\">\n"
xml_string << "\t\t\t\t\t<Uniform id=\"Uniform.52\" name=\"distr\" upper=\"Infinity\"/>\n"
xml_string << "\t\t\t\t</prior>\n"
xml_string << "\n"
xml_string << "\t\t\t\t<prior id=\"migMeanPrior.Species\" name=\"distribution\" x=\"@migMean.Species\">\n"
xml_string << "\t\t\t\t\t<OneOnX id=\"inverseUniform.1\" name=\"distr\"/>\n"
xml_string << "\t\t\t\t</prior>\n"
xml_string << "\t\t\t\t<prior id=\"popMeanPrior.Species\" name=\"distribution\" x=\"@popMean.Species\">\n"
xml_string << "\t\t\t\t\t<OneOnX id=\"OneOnX.0\" name=\"distr\"/>\n"
xml_string << "\t\t\t\t</prior>\n"
xml_string << "\n"
xml_string << "\t\t\t\t<prior id=\"migRatesPrior.Species\" name=\"distribution\" x=\"@migRates.Species\">\n"
xml_string << "\t\t\t\t\t<Exponential id=\"Exponential.0\" name=\"distr\" mean=\"0.05\"/>\n"
xml_string << "\t\t\t\t</prior>\n"
xml_string << "\t\t\t\t<prior id=\"popSizesPrior.Species\" name=\"distribution\" x=\"@popSizes.Species\">\n"
xml_string << "\t\t\t\t\t<LogNormal id=\"LogNormal.0000\" name=\"distr\" M=\"1\" S=\"0.5\" meanInRealSpace=\"true\"/>\n"
xml_string << "\t\t\t\t</prior>\n"
xml_string << "\n"
xml_string << "\t\t\t\t<prior id=\"speciationRatePrior.t:Species\" name=\"distribution\" x=\"@speciationRate.t:Species\">\n"
xml_string << "\t\t\t\t\t<Uniform id=\"Uniform.0\" name=\"distr\" upper=\"10000.0\"/>\n"
xml_string << "\t\t\t\t</prior>\n"
xml_string << "\t\t\t</distribution>\n"
xml_string << "\n"
xml_string << "\t\t\t<distribution id=\"likelihood\" spec=\"util.CompoundDistribution\">\n"
align_ids.size.times do |x|
	xml_string << "\t\t\t\t<distribution id=\"treeLikelihood.#{align_ids[x]}\" spec=\"TreeLikelihood\" data=\"@#{align_ids[x]}\" tree=\"@Tree.t:#{align_ids[x]}\" useAmbiguities=\"true\">\n"
	xml_string << "\t\t\t\t\t<siteModel id=\"SiteModel.s:#{align_ids[x]}\" spec=\"SiteModel\" mutationRate=\"@mutationRate.s:#{align_ids[x]}\">\n"
	xml_string << "\t\t\t\t\t\t<parameter id=\"gammaShape.s:#{align_ids[x]}\" estimate=\"false\" name=\"shape\">1.0</parameter>\n"
	xml_string << "\t\t\t\t\t\t<parameter id=\"proportionInvariant.s:#{align_ids[x]}\" estimate=\"false\" lower=\"0.0\" name=\"proportionInvariant\" upper=\"1.0\">0.0</parameter>\n"
	xml_string << "\t\t\t\t\t\t<substModel id=\"hky.s:#{align_ids[x]}\" spec=\"HKY\" kappa=\"@kappa.s:#{align_ids[x]}\">\n"
	xml_string << "\t\t\t\t\t\t\t<frequencies id=\"empiricalFreqs.s:#{align_ids[x]}\" spec=\"Frequencies\" data=\"@#{align_ids[x]}\"/>\n"
	xml_string << "\t\t\t\t\t\t</substModel>\n"
	xml_string << "\t\t\t\t\t</siteModel>\n"
	xml_string << "\t\t\t\t\t<branchRateModel id=\"StrictClock.c:#{align_ids[x]}\" spec=\"beast.evolution.branchratemodel.StrictClockModel\">\n"
	xml_string << "\t\t\t\t\t\t<parameter idref=\"strictClockRate\" name=\"clock.rate\"/>\n"
	xml_string << "\t\t\t\t\t</branchRateModel>\n"
	xml_string << "\t\t\t\t</distribution>\n"
end
xml_string << "\t\t\t</distribution>\n"
xml_string << "\t\t</distribution>\n"
xml_string << "\n"

# Write the operators.
xml_string << "\t\t<operator id=\"ClockScale.Species\" spec=\"ScaleOperator\" parameter=\"@strictClockRate\" scaleFactor=\"0.8\" weight=\"100.0\"/>\n"
xml_string << "\t\t<operator id=\"popSizesSwap.Species\" spec=\"starbeast2.RealCycle\" k=\"2\" optimise=\"false\" parameter=\"@popSizes.Species\" weight=\"10.0\"/>\n"
xml_string << "\t\t<operator id=\"popSizesScale.Species\" scaleAll=\"true\" scaleAllIndependently=\"true\" spec=\"ScaleOperator\" parameter=\"@popSizes.Species\" scaleFactor=\"0.5\" weight=\"10.0\"/>\n"
xml_string << "\t\t<operator id=\"popMeanScale.Species\" spec=\"ScaleOperator\" parameter=\"@popMean.Species\" scaleFactor=\"0.75\" weight=\"10.0\"/>\n"
xml_string << "\t\t<operator id=\"PopScaleUpDown:Species\" spec=\"UpDownOperator\" scaleFactor=\"0.75\" weight=\"10.0\">\n"
xml_string << "\t\t\t<up idref=\"popSizes.Species\"/>\n"
xml_string << "\t\t\t<down idref=\"popMean.Species\"/>\n"
xml_string << "\t\t</operator>\n"
xml_string << "\t\t<operator id=\"migRatesScale.Species\" scaleAll=\"true\" scaleAllIndependently=\"true\" spec=\"ScaleOperator\" parameter=\"@migRates.Species\" scaleFactor=\"0.8\" weight=\"30.0\"/>\n"
xml_string << "\t\t<operator id=\"migIndicatorFlipper.c:$(n)\" spec=\"BitFlipOperator\" parameter=\"@migIndicators.Species\" weight=\"10.0\"/>\n"
xml_string << "\t\t<operator id=\"coordinatedUniform.t:Species\" spec=\"starbeast2.CoordinatedUniform\" speciesTree=\"@Tree.t:Species\" weight=\"30.0\">\n"
align_ids.size.times do |x|
	xml_string << "\t\t\t<geneTree idref=\"Tree.t:#{align_ids[x]}\"/>\n"
end
xml_string << "\t\t</operator>\n"
xml_string << "\t\t<operator id=\"coordinatedExponential.t:Species\" spec=\"starbeast2.CoordinatedExponential\" speciesTree=\"@Tree.t:Species\" weight=\"30.0\">\n"
align_ids.size.times do |x|
	xml_string << "\t\t\t<geneTree idref=\"Tree.t:#{align_ids[x]}\"/>\n"
end
xml_string << "\t\t</operator>\n"
xml_string << "\t\t<operator id=\"SubtreeSlideAndSwap.t:Species\" spec=\"starbeast2.SubtreeSlideAndSwap\" size=\"0.5\" Ne=\"@popSizes.Species\" tree=\"@Tree.t:Species\" weight=\"30.0\"/>\n"
xml_string << "\t\t<operator id=\"WilsonBalding.t:Species\" spec=\"WilsonBalding\" tree=\"@Tree.t:Species\" weight=\"30.0\"/>\n"
xml_string << "\t\t<operator id=\"Wide.t:Species\" spec=\"Exchange\" isNarrow=\"false\" tree=\"@Tree.t:Species\" weight=\"30.0\"/>\n"
xml_string << "\t\t<operator id=\"Narrow.t:Species\" spec=\"Exchange\" tree=\"@Tree.t:Species\" weight=\"30.0\"/>\n"
xml_string << "\t\t<operator id=\"UniformOperator.t:Species\" spec=\"Uniform\" tree=\"@Tree.t:Species\" weight=\"30.0\"/>\n"
xml_string << "\t\t<operator id=\"TreeRootScaler.t:Species\" spec=\"ScaleOperator\" rootOnly=\"true\" scaleFactor=\"0.7\" tree=\"@Tree.t:Species\" weight=\"30.0\"/>\n"
xml_string << "\t\t<operator id=\"TreeScaler.t:Species\" spec=\"ScaleOperator\" scaleFactor=\"0.95\" tree=\"@Tree.t:Species\" weight=\"30.0\"/>\n"
xml_string << "\t\t<operator id=\"speciationRateScale.t:Species\" spec=\"ScaleOperator\" parameter=\"@speciationRate.t:Species\" scaleFactor=\"0.5\" weight=\"30.0\"/>\n"
xml_string << "\t\t<operator id=\"updownAll:Species\" spec=\"UpDownOperator\" scaleFactor=\"0.75\" weight=\"30.0\">\n"
xml_string << "\t\t\t<up idref=\"speciationRate.t:Species\"/>\n"
xml_string << "\t\t\t<down idref=\"popMean.Species\"/>\n"
xml_string << "\t\t\t<up idref=\"strictClockRate\"/>\n"
align_ids.size.times do |x|
	xml_string << "\t\t\t<down idref=\"Tree.t:#{align_ids[x]}\"/>\n"
end
xml_string << "\t\t</operator>\n"
align_ids.size.times do |x|
	xml_string << "\t\t<operator id=\"TreeScaler.t:#{align_ids[x]}\" spec=\"ScaleOperator\" scaleFactor=\"0.95\" tree=\"@Tree.t:#{align_ids[x]}\" weight=\"3.0\"/>\n"
	xml_string << "\t\t<operator id=\"TreeRootScaler.t:#{align_ids[x]}\" spec=\"ScaleOperator\" rootOnly=\"true\" scaleFactor=\"0.7\" tree=\"@Tree.t:#{align_ids[x]}\" weight=\"3.0\"/>\n"
	xml_string << "\t\t<operator id=\"UniformOperator.t:#{align_ids[x]}\" spec=\"Uniform\" tree=\"@Tree.t:#{align_ids[x]}\" weight=\"15.0\"/>\n"
	xml_string << "\t\t<operator id=\"SubtreeSlide.t:#{align_ids[x]}\" spec=\"SubtreeSlide\" size=\"0.002\" tree=\"@Tree.t:#{align_ids[x]}\" weight=\"15.0\"/>\n"
	xml_string << "\t\t<operator id=\"Narrow.t:#{align_ids[x]}\" spec=\"Exchange\" tree=\"@Tree.t:#{align_ids[x]}\" weight=\"15.0\"/>\n"
	xml_string << "\t\t<operator id=\"Wide.t:#{align_ids[x]}\" spec=\"Exchange\" isNarrow=\"false\" tree=\"@Tree.t:#{align_ids[x]}\" weight=\"15.0\"/>\n"
	xml_string << "\t\t<operator id=\"WilsonBalding.t:#{align_ids[x]}\" spec=\"WilsonBalding\" tree=\"@Tree.t:#{align_ids[x]}\" weight=\"15.0\"/>\n"
end
xml_string << "\t\t<operator id=\"FixMeanMutationRatesOperator\" spec=\"DeltaExchangeOperator\" delta=\"0.75\" weight=\"100.0\">\n"
align_ids.size.times do |x|
	xml_string << "\t\t\t<parameter idref=\"mutationRate.s:#{align_ids[x]}\"/>\n"
end
xml_string << "\t\t\t<weightvector id=\"weightparameter\" spec=\"parameter.IntegerParameter\" dimension=\"#{align_ids.size}\" estimate=\"false\" lower=\"0\" upper=\"0\">"
align_ids.size.times do |x|
	xml_string << "\t#{align_lengths[x]}"
end
xml_string << "</weightvector>\n"
xml_string << "\t\t</operator>\n"
align_ids.size.times do |x|
	xml_string << "\t\t<operator id=\"KappaScaler.s:#{align_ids[x]}\" spec=\"ScaleOperator\" parameter=\"@kappa.s:#{align_ids[x]}\" scaleFactor=\"0.75\" weight=\"1.0\"/>\n"
end

# Write the loggers.
xml_string << "\t\t<logger id=\"tracelog1\" fileName=\"#{xml_id}.log\" logEvery=\"#{log_every}\" sort=\"smart\">\n"
xml_string << "\t\t\t<log idref=\"posterior\"/>\n"
xml_string << "\t\t\t<log idref=\"likelihood\"/>\n"
xml_string << "\t\t\t<log idref=\"prior\"/>\n"
xml_string << "\t\t\t<log idref=\"speciescoalescent\"/>\n"
xml_string << "\t\t\t<log idref=\"speciationRate.t:Species\"/>\n"
xml_string << "\t\t\t<log idref=\"strictClockRate\"/>\n"
xml_string << "\t\t\t<log idref=\"YuleModel.t:Species\"/>\n"
xml_string << "\t\t\t<log idref=\"popMean.Species\"/>\n"
xml_string << "\t\t\t<log idref=\"migMean.Species\"/>\n"
xml_string << "\t\t\t<log idref=\"migRates.Species\"/>\n"
xml_string << "\t\t\t<log idref=\"migIndicatorSum.Species\"/>\n"
xml_string << "\t\t\t<log idref=\"migIndicators.Species\"/>\n"
xml_string << "\t\t\t<log idref=\"popSizes.Species\"/>\n"
xml_string << "\t\t\t<log id=\"TreeHeight.Species\" spec=\"beast.evolution.tree.TreeHeightLogger\" tree=\"@Tree.t:Species\"/>\n"
xml_string << "\t\t\t<log id=\"TreeLength.Species\" spec=\"starbeast2.TreeLengthLogger\" tree=\"@Tree.t:Species\"/>\n"
xml_string << "\t\t</logger>\n"
xml_string << "\n"
xml_string << "\t\t<logger id=\"tracelog2\" fileName=\"#{xml_id}_lik.log\" logEvery=\"#{log_every}\" sort=\"smart\">\n"
xml_string << "\t\t\t<log idref=\"posterior\"/>\n"
xml_string << "\t\t\t<log idref=\"likelihood\"/>\n"
xml_string << "\t\t\t<log idref=\"prior\"/>\n"
align_ids.size.times do |x|
	xml_string << "\t\t\t<log idref=\"treeLikelihood.#{align_ids[x]}\"/>\n"
end
xml_string << "\t\t</logger>\n"
xml_string << "\n"
xml_string << "\t\t<logger id=\"tracelog3\" fileName=\"#{xml_id}_heights.log\" logEvery=\"#{log_every}\" sort=\"smart\">\n"
xml_string << "\t\t\t<log idref=\"posterior\"/>\n"
xml_string << "\t\t\t<log idref=\"likelihood\"/>\n"
xml_string << "\t\t\t<log idref=\"prior\"/>\n"
align_ids.size.times do |x|
	xml_string << "\t\t\t<log id=\"TreeHeight.t:#{align_ids[x]}\" spec=\"beast.evolution.tree.TreeHeightLogger\" tree=\"@Tree.t:#{align_ids[x]}\"/>\n"
end
xml_string << "\t\t</logger>\n"
xml_string << "\n"
xml_string << "\t\t<logger id=\"tracelog4\" fileName=\"#{xml_id}_kappa.log\" logEvery=\"#{log_every}\" sort=\"smart\">\n"
xml_string << "\t\t\t<log idref=\"posterior\"/>\n"
xml_string << "\t\t\t<log idref=\"likelihood\"/>\n"
xml_string << "\t\t\t<log idref=\"prior\"/>\n"
align_ids.size.times do |x|
	xml_string << "\t\t\t<log idref=\"kappa.s:#{align_ids[x]}\"/>\n"
end
xml_string << "\t\t</logger>\n"
xml_string << "\n"
xml_string << "\t\t<logger id=\"tracelog5\" fileName=\"#{xml_id}_mutation.log\" logEvery=\"#{log_every}\" sort=\"smart\">\n"
xml_string << "\t\t\t<log idref=\"posterior\"/>\n"
xml_string << "\t\t\t<log idref=\"likelihood\"/>\n"
xml_string << "\t\t\t<log idref=\"prior\"/>\n"
align_ids.size.times do |x|
	xml_string << "\t\t\t<log idref=\"mutationRate.s:#{align_ids[x]}\"/>\n"
end
xml_string << "\t\t</logger>\n"
xml_string << "\n"
xml_string << "\t\t<logger id=\"speciesTreeLogger\" fileName=\"#{xml_id}_species.trees\" logEvery=\"#{log_every}\" mode=\"tree\">\n"
xml_string << "\t\t\t<log id=\"SpeciesTreeLoggerX\" spec=\"starbeast2.SpeciesTreeLoggerWithGeneFlow\" populationModel=\"@popModelAIM.Species\"/>\n"
xml_string << "\t\t</logger>\n"
xml_string << "\n"
xml_string << "\t\t<logger id=\"screenlog\" logEvery=\"#{(0.1*log_every).round}\">\n"
xml_string << "\t\t\t<log idref=\"posterior\"/>\n"
xml_string << "\t\t\t<log id=\"ESS.0\" spec=\"util.ESS\" arg=\"@posterior\"/>\n"
xml_string << "\t\t\t<log idref=\"likelihood\"/>\n"
xml_string << "\t\t\t<log idref=\"prior\"/>\n"
xml_string << "\t\t</logger>\n"
xml_string << "\n"
align_ids.size.times do |x|
	xml_string << "\t\t<logger id=\"treelog.t:#{align_ids[x]}\" fileName=\"#{xml_id}_#{align_ids[x]}.trees\" logEvery=\"#{log_every}\" mode=\"tree\">\n"
	xml_string << "\t\t\t<log id=\"TreeWithMetaDataLogger.t:#{align_ids[x]}\" spec=\"beast.evolution.tree.TreeWithMetaDataLogger\" tree=\"@Tree.t:#{align_ids[x]}\"/>\n"
	xml_string << "\t\t</logger>\n"
end
xml_string << "\n"

# Complete the run element.
xml_string << "\t</run>\n"

# Complete the xml string.
xml_string << "</beast>\n"

# Write the xml file for aim.
xml_file = File.open(xml_file_name, "w")
xml_file.write(xml_string)
