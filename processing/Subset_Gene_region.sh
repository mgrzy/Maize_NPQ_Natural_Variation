conda activate geno

bcftools view ../../BigData/Maize/v5/Genotypes/All/VCF/imputed/chr_3_imputed.vcf.gz \
-S ../../BigData/Maize/v5/Genotypes/WiDiv/rMVP/Bugeater.geno.ind -i "MAF>0.05" -r chr3:179564407-179566950 \
-Oz -o ../Data/Maize_NPQ_Natural_Variation/data/genomic/VCF/PsbS.vcf.gz

bcftools index ../Data/Maize_NPQ_Natural_Variation/data/genomic/VCF/PsbS.vcf.gz -t

bcftools view ../../BigData/Maize/v5/Genotypes/All/VCF/imputed/chr_3_imputed.vcf.gz \
-S ../../BigData/Maize/v5/Genotypes/WiDiv/rMVP/Bugeater.geno.ind -i "MAF>0.05" -r chr3:178404716-178415602 \
-Oz -o ../Data/Maize_NPQ_Natural_Variation/data/genomic/VCF/PSI3.vcf.gz

bcftools index ../Data/Maize_NPQ_Natural_Variation/data/genomic/VCF/PSI3.vcf.gz -t
