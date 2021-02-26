# m_matschiner Tue Oct 9 16:32:32 CEST 2018

# Make a nobackup directory.
mkdir -p ../res/windows/5000bp/nobackup

# Move info files without alignment files into the nobackup directory.
for info in ../res/windows/5000bp/LG*.info.txt
do
    align=${info%.info.txt}.phy
    if [ ! -f ${align} ]
    then
        mv ${info} ../res/windows/5000bp/nobackup
    fi
done