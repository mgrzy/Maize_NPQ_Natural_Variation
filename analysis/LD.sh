conda activate geno

### PSI-PSBS
plink --bfile ../../BigData/Maize/v4/Genotypes/WiDiv_RNAseq/plink/WiDiv752 --chr 3 --from-bp 175652936 \
--to-bp 178038945 --out ../Data/Maize_NPQ_Natural_Variation/data/LD/rs3_176152936 --make-bed

plink --bfile ../Data/Maize_NPQ_Natural_Variation/data/LD/rs3_176152936 --r2 --ld-snp rs3_176152936 \
--ld-window-kb 4000 --ld-window 99999 --ld-window-r2 0 \
--out ../Data/Maize_NPQ_Natural_Variation/data/LD/LD_rs3_176152936

###  IRM1 
plink --bfile ../../BigData/Maize/v4/Genotypes/WiDiv_RNAseq/plink/WiDiv752 --chr 2 --from-bp 182332835 \
--to-bp 182432835 --out ../Data/Maize_NPQ_Natural_Variation/data/LD/rs2_182382835 --make-bed

plink --bfile ../Data/Maize_NPQ_Natural_Variation/data/LD/rs2_182382835 --r2 --ld-snp rs2_182382835 \
--ld-window-kb 4000 --ld-window 99999 --ld-window-r2 0 \
--out ../Data/Maize_NPQ_Natural_Variation/data/LD/LD_rs2_182382835

###TRX-Y1
plink --bfile ../../BigData/Maize/v4/Genotypes/WiDiv_RNAseq/plink/WiDiv752 --chr 3 --from-bp 147481265 \
--to-bp 147581265 --out ../Data/Maize_NPQ_Natural_Variation/data/LD/rs3_147531265 --make-bed

plink --bfile ../Data/Maize_NPQ_Natural_Variation/data/LD/rs3_147531265 --r2 --ld-snp rs3_147531265 \
--ld-window-kb 4000 --ld-window 99999 --ld-window-r2 0 \
--out ../Data/Maize_NPQ_Natural_Variation/data/LD/LD_rs3_147531265

### ACHT3
plink --bfile ../../BigData/Maize/v4/Genotypes/WiDiv_RNAseq/plink/WiDiv752 --chr 7 --from-bp 179231748 \
--to-bp 179331748 --out ../Data/Maize_NPQ_Natural_Variation/data/LD/rs7_179281748 --make-bed

plink --bfile ../Data/Maize_NPQ_Natural_Variation/data/LD/rs7_179281748 --r2 --ld-snp rs7_179281748 \
--ld-window-kb 4000 --ld-window 99999 --ld-window-r2 0 \
--out ../Data/Maize_NPQ_Natural_Variation/data/LD/LD_rs7_179281748

### PMI1
plink --bfile ../../BigData/Maize/v4/Genotypes/WiDiv_RNAseq/plink/WiDiv752 --chr 1 --from-bp 86014177 \
--to-bp 86114177 --out ../Data/Maize_NPQ_Natural_Variation/data/LD/rs1_86064177 --make-bed

plink --bfile ../Data/Maize_NPQ_Natural_Variation/data/LD/rs1_86064177 --r2 --ld-snp rs1_86064177 \
--ld-window-kb 4000 --ld-window 99999 --ld-window-r2 0 \
--out ../Data/Maize_NPQ_Natural_Variation/data/LD/LD_rs1_86064177

### OEP37
plink --bfile ../../BigData/Maize/v4/Genotypes/WiDiv_RNAseq/plink/WiDiv752 --chr 5 --from-bp 187809720 \
--to-bp 187909720 --out ../Data/Maize_NPQ_Natural_Variation/data/LD/rs5_187859720 --make-bed

plink --bfile ../Data/Maize_NPQ_Natural_Variation/data/LD/rs5_187859720 --r2 --ld-snp rs5_187859720 \
--ld-window-kb 4000 --ld-window 99999 --ld-window-r2 0 \
--out ../Data/Maize_NPQ_Natural_Variation/data/LD/LD_rs5_187859720

