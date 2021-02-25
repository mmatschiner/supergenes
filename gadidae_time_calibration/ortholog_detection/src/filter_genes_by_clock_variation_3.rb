# m_matschiner Tue Jul 24 13:02:25 CEST 2018

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
    def hpd_lower(proportion)
        raise "The interval should be between 0 and 1!" if proportion >= 1 or proportion <= 0
        sorted_array = self.sort
        hpd_index = 0
        min_range = sorted_array[-1]
        diff = (proportion*self.size).round
        (self.size-diff).times do |i|
            min_value = sorted_array[i]
            max_value = sorted_array[i+diff-1]
            range = max_value - min_value
            if range < min_range
                min_range = range
                hpd_index = i
            end
        end
        sorted_array[hpd_index]
    end
    def hpd_upper(proportion)
        raise "The interval should be between 0 and 1!" if proportion >= 1 or proportion <= 0
        sorted_array = self.sort
        hpd_index = 0
        min_range = sorted_array[-1]
        diff = (proportion*self.size).round
        (self.size-diff).times do |i|
            min_value = sorted_array[i]
            max_value = sorted_array[i+diff-1]
            range = max_value - min_value
            if range < min_range
                min_range = range
                hpd_index = i
            end
        end
        sorted_array[hpd_index+diff-1]
    end
end

# Get the command line arguments.
alignment_directory_in = ARGV[0].chomp("/")
alignment_directory_out = ARGV[1].chomp("/")
filtered_nuclear_exons_info_file_name = ARGV[2]
selected_nuclear_exons_info_file_name = ARGV[3]
coefficients_of_variation_max_mean = ARGV[4].to_f
coefficients_of_variation_max_upper = ARGV[5].to_f

# Create the output directory if it does not exist yet.
unless Dir.exists?(alignment_directory_out)
    Dir.mkdir(alignment_directory_out)
end

# Get the names of the directories with exon alignments for each gene.
alignment_directory_entries = Dir.entries(alignment_directory_in)
gene_alignment_directories = []
alignment_directory_entries.each {|e| gene_alignment_directories << e if e[0..2] == "ENS"}

# Read the nuclear_queries_exons.txt file from the info directory.
filtered_nuclear_exons_info_file = File.open(filtered_nuclear_exons_info_file_name)
filtered_nuclear_exons_info_lines = filtered_nuclear_exons_info_file.readlines
selected_nuclear_exons_info_lines = []

# For each of the alignment directories, do the following.
is_first_dir = true
nexus_ids_in_first_alignment = []
gene_alignment_directories.each do |dir|

    # Read the BEAST log file.
    log_file_name = "#{alignment_directory_in}/#{dir}/#{dir}.log"
    print "Analysing log file #{log_file_name.split("/").last}..."
    log_file = File.open(log_file_name)
    log_lines = log_file.readlines
    coefficients_of_variation = []
    mutation_rates = []
    out_of_comments = false
    coefficients_of_variation_index = nil
    mutation_rate_index = nil
    log_lines.each do |l|
        line_ary = l.split
        if line_ary[0].downcase == "sample"
            out_of_comments = true
            line_ary.size.times do |x|
                coefficients_of_variation_index = x if line_ary[x].match(/ucldStdev/)
                mutation_rate_index = x if line_ary[x].match(/mutationRate/)
            end
        elsif out_of_comments
            unless l[0] == "#"
                if coefficients_of_variation_index == 0
                    puts "ERROR: Column for coefficient of variation could not be found!"
                    exit 1
                end
                coefficients_of_variation << line_ary[coefficients_of_variation_index].to_f
                if mutation_rate_index == 0
                    puts "ERROR: Column for mutation rate could not be found!"
                    exit 1
                end
                mutation_rates << line_ary[mutation_rate_index].to_f
            end
        end
    end
    if coefficients_of_variation.size == 0
        puts "ERROR: Coefficient of variation could not be read!"
        exit 1
    end
    coefficients_of_variation = coefficients_of_variation[1..-1]
    number_of_burnin_samples = (coefficients_of_variation.size)/5
    coefficients_of_variation_wo_burnin = coefficients_of_variation[number_of_burnin_samples+1..-1]
    coefficients_of_variation_wo_burnin_mean = coefficients_of_variation_wo_burnin.mean
    coefficients_of_variation_wo_burnin_hpd_upper = coefficients_of_variation_wo_burnin.hpd_upper(0.95)
    coefficients_of_variation_wo_burnin_hpd_lower = coefficients_of_variation_wo_burnin.hpd_lower(0.95)
    mutation_rates_wo_burnin = mutation_rates[number_of_burnin_samples+1..-1]
    mutation_rates_wo_burnin_mean = mutation_rates_wo_burnin.mean

    if coefficients_of_variation_wo_burnin_mean < coefficients_of_variation_max_mean and coefficients_of_variation_wo_burnin_hpd_upper < coefficients_of_variation_max_upper

        # Get the ids and seqs from the Nexus format alignment.
        nexus_file = File.open("#{alignment_directory_in}/#{dir}/#{dir}.nex")
        nexus_lines = nexus_file.readlines
        nexus_ids_in_this_alignment = []
        nexus_seqs_in_this_alignment = []
        in_matrix = false
        nexus_lines.each do |l|
            if l.strip.downcase == "matrix"
                in_matrix = true
            elsif l.strip == ";"
                in_matrix = false
            elsif l.strip != "" and in_matrix
                nexus_ids_in_this_alignment << l.split[0].strip
                nexus_seqs_in_this_alignment << l.split[1].strip
            end
        end

        # Store the nexus ids of the first alignment.
        if is_first_dir
            nexus_ids_in_first_alignment = nexus_ids_in_this_alignment
        else
            if nexus_ids_in_first_alignment != nexus_ids_in_this_alignment
                puts "ERROR: Nexus ids differ between alignments!"
                puts "nexus_ids_in_first_alignment: #{nexus_ids_in_first_alignment}"
                puts
                puts "nexus_ids_in_this_alignment: #{nexus_ids_in_this_alignment}"
                puts
                exit 1
            end
        end

        # Prepare the string for a new nexus file.
        new_nexus_string = "#nexus\n"
        new_nexus_string << "\n"
        new_nexus_string << "begin data;\n"
        new_nexus_string << "dimensions  ntax=#{nexus_ids_in_this_alignment.size} nchar=#{nexus_seqs_in_this_alignment[0].size};\n"
        new_nexus_string << "format datatype=DNA gap=- missing=?;\n"
        new_nexus_string << "matrix\n"
        nexus_seqs_in_this_alignment.size.times do |x|
            new_nexus_string << "#{nexus_ids_in_this_alignment[x].ljust(12)}#{nexus_seqs_in_this_alignment[x]}\n"
        end
        new_nexus_string << ";\n"
        new_nexus_string << "end;\n"

        # Write the new nexus file.
        new_nexus_file_name = "#{dir}.nex"
        new_nexus_file = File.open("#{alignment_directory_out}/#{new_nexus_file_name}","w")
        new_nexus_file.write(new_nexus_string)

        # Add to the selected nuclear exons info lines.
        filtered_nuclear_exons_info_lines.each do |l|
            selected_nuclear_exons_info_lines << l if l.include?(dir)
        end
        is_first_dir = false

    end

    # Feedback.
    puts " done. Coefficient of variation: #{coefficients_of_variation_wo_burnin_mean.round(3)} (#{coefficients_of_variation_wo_burnin_hpd_lower.round(3)}-#{coefficients_of_variation_wo_burnin_hpd_upper.round(3)}); mutation rate: #{mutation_rates_wo_burnin_mean.round(5)}."

end

# Prepare the string for the selected nuclear exons file.
selected_nuclear_exons_info_string = "exon_id\tgene_id\ttranscript_id\ttranslation\tlength\tthreshold_bitscore\tcorrect_bitscores\tincorrect_bitscores\tmin(correct_bitscores)-max(incorrect_bitscores)\n"
selected_nuclear_exons_info_lines.each {|l| selected_nuclear_exons_info_string << l}

# Write the selected nuclear exons file to the info directory.
selected_nuclear_exons_info_file = File.open(selected_nuclear_exons_info_file_name,"w")
selected_nuclear_exons_info_file.write(selected_nuclear_exons_info_string)
selected_nuclear_exons_info_file.close
