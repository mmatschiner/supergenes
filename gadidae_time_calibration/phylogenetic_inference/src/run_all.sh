# m_matschiner Tue Aug 7 12:24:03 CEST 2018

# Copy the alignments.
bash get_alignments.sh

# Run raxml for all genes.
bash run_raxml_per_gene.sh

# Get the bootstrap node support values.
bash get_node_support_of_raxml_trees.sh

# Make a directory only with alignments that had greater than 65 node support.
bash select_alignments.sh

# Prepare beast analyses.
bash prepare_beast_analyses.sh

# Run beast analyses.
bash run_beast_analyses.sh

# Summarize beast analyses.
bash summarize_beast_analyses.sh