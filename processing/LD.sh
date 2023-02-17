conda activate geno

#1
plink --bfile ../../BigData/Maize/v4/Genotypes/WiDiv_RNAseq/plink/WiDiv752 --chr 3 --from-bp 175652936 --to-bp 178038945 --out ../Data/Maize_NPQ_Natural_Variation/data/LD/rs3_176152936 --make-bed

plink --bfile ../Data/Maize_NPQ_Natural_Variation/data/LD/rs3_176152936 --r2 --ld-snp rs3_176152936 --ld-window-kb 4000 --ld-window 99999 --ld-window-r2 0 \
--out ../Data/Maize_NPQ_Natural_Variation/data/LD/LD_rs3_176152936

plink --bfile ../../BigData/Maize/v4/Genotypes/WiDiv_RNAseq/plink/WiDiv752 --chr 2 --from-bp 181882835 --to-bp 182882835 \
--out ../Data/Maize_NPQ_Natural_Variation/data/LD/rs2_182382835 --make-bed

plink --bfile ../Data/Maize_NPQ_Natural_Variation/data/LD/rs2_182382835 --r2 --ld-snp rs2_182382835 --ld-window-kb 4000 --ld-window 99999 --ld-window-r2 0 \
--out ../Data/Maize_NPQ_Natural_Variation/data/LD/LD_rs2_182382835
