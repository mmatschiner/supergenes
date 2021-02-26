# michaelm Wed Mar 8 00:08:44 CET 2017

# Make the log directory.
mkdir -p ../log/misc

# Set the data directory.
data_dir="../res/mapping/"

# Set the reference.
ref=../data/assemblies/gadMor2.fasta

# Read the list of bam files in the data directory.
bams_with_relative_path=`ls ${data_dir}*.bam`

# Produce a list of unique sample identifiers.
truncated_bams_with_relative_path=()
for bam_with_relative_path in ${bams_with_relative_path[@]}
do
    split_ary=(${bam_with_relative_path//_/ })
    split_ary_size=(${#split_ary[@]})
    trim_part=${split_ary[$split_ary_size-1]}
    truncated_bam_with_relative_path=${bam_with_relative_path%_$trim_part}
    truncated_bams_with_relative_path+=($truncated_bam_with_relative_path)
done
unique_truncated_bams_with_relative_path=$(echo "${truncated_bams_with_relative_path[@]}" | tr ' ' '\n' | sort -u | tr '\n' ' ')

# For each item in that list, start merging with all bam files that match that item.
for unique_truncated_bam_with_relative_path in ${unique_truncated_bams_with_relative_path[@]}
do
    if ls ${unique_truncated_bam_with_relative_path}_*.bam 1> /dev/null 2>&1
    then
        if [ ! -f ${unique_truncated_bam_with_relative_path}.bam ]
        then
            infiles=(`ls ${unique_truncated_bam_with_relative_path}_*.bam`)
            number_of_infiles=(${#infiles[@]})
            # If the individual has sequence data for only a single library, just give it a new name.
            if [ ${number_of_infiles} -eq 1 ]
            then
                mv -v ${infiles[0]} ${unique_truncated_bam_with_relative_path}.bam
                mv -v ${infiles[0]}.bai ${unique_truncated_bam_with_relative_path}.bam.bai

                # If not, merge all sequence files for this individual.
                else
                unique_truncated_bam=`basename ${unique_truncated_bam_with_relative_path}`
                out=../log/misc/merge.${unique_truncated_bam}.out
                log=../log/misc/merge.${unique_truncated_bam}.log
                rm -f ${out}
                if [ ${unique_truncated_bam} == "melaeg" ]
                then
                    sbatch -o $out --partition=hugemem --time=100:00:00 --mem-per-cpu=25G merge.slurm ${ref} ${log} ${unique_truncated_bam_with_relative_path}*.bam
                else
                    sbatch -o $out merge.slurm ${ref} ${log} ${unique_truncated_bam_with_relative_path}*.bam
                fi
            fi
        fi
    fi
done
