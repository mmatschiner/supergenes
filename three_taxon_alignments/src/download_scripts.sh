# m_matschiner Thu Feb 25 16:19:24 CET 2021

# Make bin directory.
mkdir -p ../bin/kent

# Download scripts by jim kent.
cd ../bin/kent
rsync -aP rsync://hgdownload.soe.ucsc.edu/genome/admin/exe/linux.x86_64/ ./
cd -
