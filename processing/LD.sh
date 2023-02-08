conda activate geno

#1
plink --bfile ../../BigData/WiDivGeno/RNAseqagpv4/plink/WiDiv752 --chr 3 --from-bp 175652936 --to-bp 178038945 --out data/work/loci/psbs --make-bed

plink --bfile data/work/loci/psbs --r2 --ld-snp rs3_176152936 --ld-window-kb 4000 --ld-window 99999 --ld-window-r2 0 --out data/work/loci/LD_rs3_176152936
