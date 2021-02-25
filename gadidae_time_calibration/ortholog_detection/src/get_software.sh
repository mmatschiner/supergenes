# m_matschiner Thu Feb 25 19:39:47 CET 2021

# Download bmge.
cd ../bin
wget ftp://ftp.pasteur.fr:21/pub/gensoft/projects/BMGE/BMGE-1.0.tar.gz
tar -xzf BMGE-1.0.tar.gz
mv BMGE-1.0/BMGE.jar .
rm -r BMGE-1.0

# Uncompress concaterpillar and raxmlhpc.
tar -xzf tools.tgz
rm tools.tgz
mv tools/* .
rm -r tools