conda activate geno

### PsbS (Zm00001eb146510)
bcftools view ../../BigData/Maize/v5/Genotypes/All/VCF/imputed/chr_3_imputed.vcf.gz \
-S ../../BigData/Maize/v5/Genotypes/WiDiv/rMVP/Bugeater.geno.ind -i "MAF>0.05" -r chr3:179564407-179566950 \
-Oz -o ../Data/Maize_NPQ_Natural_Variation/data/genomic/VCF/PsbS.vcf.gz

bcftools index ../Data/Maize_NPQ_Natural_Variation/data/genomic/VCF/PsbS.vcf.gz -t

### PSI3 (Zm00001eb146270)
bcftools view ../../BigData/Maize/v5/Genotypes/All/VCF/imputed/chr_3_imputed.vcf.gz \
-S ../../BigData/Maize/v5/Genotypes/WiDiv/rMVP/Bugeater.geno.ind -i "MAF>0.05" -r chr3:178404716-178415602 \
-Oz -o ../Data/Maize_NPQ_Natural_Variation/data/genomic/VCF/PSI3.vcf.gz

bcftools index ../Data/Maize_NPQ_Natural_Variation/data/genomic/VCF/PSI3.vcf.gz -t

### IRM1 (Zm00001eb098640)
bcftools view ../../BigData/Maize/v5/Genotypes/All/VCF/imputed/chr_2_imputed.vcf.gz \
-S ../../BigData/Maize/v5/Genotypes/WiDiv/rMVP/Bugeater.geno.ind -i "MAF>0.05" -r chr2:181981412-181982265 \
-Oz -o ../Data/Maize_NPQ_Natural_Variation/data/genomic/VCF/IRM1.vcf.gz

bcftools index ../Data/Maize_NPQ_Natural_Variation/data/genomic/VCF/IRM1.vcf.gz -t

### TRX1 (Zm00001eb330920)
bcftools view ../../BigData/Maize/v5/Genotypes/All/VCF/imputed/chr_3_imputed.vcf.gz \
-S ../../BigData/Maize/v5/Genotypes/WiDiv/rMVP/Bugeater.geno.ind -i "MAF>0.05" -r chr3:148313086-148315944 \
-Oz -o ../Data/Maize_NPQ_Natural_Variation/data/genomic/VCF/TRX1.vcf.gz

bcftools index ../Data/Maize_NPQ_Natural_Variation/data/genomic/VCF/TRX1.vcf.gz -t

### ACHT3 (Zm00001d022518)
bcftools view ../../BigData/Maize/v5/Genotypes/All/VCF/imputed/chr_7_imputed.vcf.gz \
-S ../../BigData/Maize/v5/Genotypes/WiDiv/rMVP/Bugeater.geno.ind -i "MAF>0.05" -r chr7:182741326-182743238 \
-Oz -o ../Data/Maize_NPQ_Natural_Variation/data/genomic/VCF/ACHT3.vcf.gz

bcftools index ../Data/Maize_NPQ_Natural_Variation/data/genomic/VCF/ACHT3.vcf.gz -t

### PMI (Zm00001eb022210)
bcftools view ../../BigData/Maize/v5/Genotypes/All/VCF/imputed/chr_1_imputed.vcf.gz \
-S ../../BigData/Maize/v5/Genotypes/WiDiv/rMVP/Bugeater.geno.ind -i "MAF>0.05" -r chr1:85240151-85243577 \
-Oz -o ../Data/Maize_NPQ_Natural_Variation/data/genomic/VCF/PMI1.vcf.gz

bcftools index ../Data/Maize_NPQ_Natural_Variation/data/genomic/VCF/PMI1.vcf.gz -t

### OEP37 (Zm00001eb246750)
bcftools view ../../BigData/Maize/v5/Genotypes/All/VCF/imputed/chr_5_imputed.vcf.gz \
-S ../../BigData/Maize/v5/Genotypes/WiDiv/rMVP/Bugeater.geno.ind -i "MAF>0.05" -r chr5:186951702-186956396 \
-Oz -o ../Data/Maize_NPQ_Natural_Variation/data/genomic/VCF/OEP37.vcf.gz

bcftools index ../Data/Maize_NPQ_Natural_Variation/data/genomic/VCF/OEP37.vcf.gz -t

