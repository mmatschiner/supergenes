# m_matschiner Thu Feb 25 17:30:40 CET 2021

# Make data directories.
mkdir -p ../data/ped

# Download snp datasets for lgs 1, 2, 7, and 12 in ped format from the zenodo repository.
cd ../data/ped
wget https://zenodo.org/record/4560275/files/gadus_morhua_supergene_lg01_plink.ped
wget https://zenodo.org/record/4560275/files/gadus_morhua_supergene_lg01_plink.map
wget https://zenodo.org/record/4560275/files/gadus_morhua_supergene_lg02_plink.ped
wget https://zenodo.org/record/4560275/files/gadus_morhua_supergene_lg02_plink.map
wget https://zenodo.org/record/4560275/files/gadus_morhua_supergene_lg07_plink.ped
wget https://zenodo.org/record/4560275/files/gadus_morhua_supergene_lg07_plink.map
wget https://zenodo.org/record/4560275/files/gadus_morhua_supergene_lg12_plink.ped
wget https://zenodo.org/record/4560275/files/gadus_morhua_supergene_lg12_plink.map
cd -