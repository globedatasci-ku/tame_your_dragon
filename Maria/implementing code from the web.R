# Tip 2 - implementing code for PCA plot
setwd("/xxx/pca")
library("SNPRelate")
vcf.fn<-"~/xxx/tmp.vcf"
snpgdsVCF2GDS(vcf.fn, "ccm.gds",  method="biallelic.only")
genofile <- openfn.gds("ccm.gds")
ccm_pca<-snpgdsPCA(genofile)
plot(ccm_pca$eigenvect[,1],ccm_pca$eigenvect[,2] ,col=as.numeric(substr(ccm_pca$sample, 1,3) == 'CCM')+3, pch=2)