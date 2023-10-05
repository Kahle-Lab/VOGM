library(devtools)
library(denovolyzeR)
library(reshape)
library(dplyr)
options(stringsAsFactors = F)

setwd("/Users/shujuanzhao/Dropbox/shujuan/Statistical Analysis")
raw_case = read.table("114-DeNovo-Input.txt",sep="\t",header=TRUE,stringsAsFactors = FALSE)
raw_case$Cases10GeneName <- toupper(raw_case$Cases10GeneName)


cases = raw_case
case_num = 114

#ctrls = raw_control
#ctrls_num = 1798

HBE_genes <- read.table(file='HBE Genes.txt',sep="\t",header=TRUE,stringsAsFactors = FALSE)
HBE_genes  <- toupper(HBE_genes[[2]])
index <- which(HBE_genes=="N/A")
HBE_genes <- HBE_genes[-index]
HBE_genes <- unique(HBE_genes)

Int_genes <- read.table(file="LoF-Intolerant-Genes-gnomAD2.1.1.txt",header=TRUE,stringsAsFactors = FALSE,sep="\t")
Int_genes <- toupper(Int_genes[[2]])
index <- which(Int_genes=="N/A")
if(length(index) != 0){
	Int_genes <- Int_genes[-index]
}

# Intolerant HBE genes
HBE_Int_genes <- intersect(HBE_genes, Int_genes)


## Read in and reformat probability tables for both cases and controls
## code for reformatting the table ##
cases10 = read.delim("1210_hs37d5_coding_idt_med_v2_spikein_padded_Mar2018_19347-unique-none-0-gene_mpc-modified.txt",sep="\t",stringsAsFactors = F)
#controls10 = read.delim("Padded_control_adj_19347-unique-gene_mpc-modified.txt",sep="\t",stringsAsFactors = F)

reformat_pDNM = function(x, mis_filter="Mis_mpc2_or_MetaSVM"){  # MetaSVM only
  names(x)[names(x)==mis_filter] <- "misD"
  x$prot <- x[,names(x) %in% c("lof","mis")] %>% apply(MARGIN=1, sum)
  x$protD <- x[,names(x) %in% c("lof","misD")] %>% apply(MARGIN=1,sum)
  x <- melt(x)                                          #use anything that is 'chr' as id variables
  names(x)[names(x)=="variable"] <- "class"             #change whichever column with the name 'variable' to 'class'
  return(x)
}
unique(reformat_pDNM(cases10)[,"class"])
pDNM_cases10 <- reformat_pDNM(cases10)

## Gene level significance using the column denovo_call_metaSVM_MPC
# Cases
casesByGene = denovolyzeByGene(cases$Cases10GeneName,cases$denovo_call_metaSVM_MPC,case_num, geneId="gene",probTable=pDNM_cases10, includeGenes="all", signifP=3, roundExpected =15, includeClasses=c("protD", "lof", "misD","prot"))
casesByGene = cbind.data.frame(rownames(casesByGene),casesByGene)
colnames(casesByGene)[1] = 'GeneOrder'
write.table(file="Gene_Level_Significance_MetaSVM_MPC-DMis_114Cases.txt",casesByGene,col.names=T,row.names=F,sep="\t",append=F,quote=F)
##### Generate enrichment analysis for each functional class
# All genes in 114 cases using the column denovo_call_metaSVM_MPC
denovolyzeByClass(genes=cases$Cases10GeneName,classes= cases$denovo_call_metaSVM_MPC,geneId="gene", nsamples=case_num, probTable=pDNM_cases10, includeClasses=c("all", "syn", "mis", "protD", "misD","lof", "prot"))
write.table(file="All_gene_enrichment_114Cases.txt",denovolyzeByClass(genes=cases$Cases10GeneName,classes= cases$denovo_call_metaSVM_MPC,geneId="gene", nsamples=case_num, probTable=pDNM_cases10, includeClasses=c("all", "syn", "mis", "protD", "misD","lof", "prot")),col.names=T,row.names=F,sep="\t",append=F,quote=F)
# HBE genes in 114 cases using the column denovo_call_metaSVM_MPC
denovolyzeByClass(genes=cases$Cases10GeneName,classes= cases$denovo_call_metaSVM_MPC,geneId="gene", nsamples=case_num, probTable=pDNM_cases10, includeGenes=HBE_genes, includeClasses=c("all", "syn", "mis", "protD", "misD","lof", "prot"))
write.table(file="HBE_gene_enrichment_114Cases.txt",denovolyzeByClass(genes=cases$Cases10GeneName,classes= cases$denovo_call_metaSVM_MPC,geneId="gene", nsamples=case_num, probTable=pDNM_cases10, includeGenes=HBE_genes, includeClasses=c("all", "syn", "mis", "protD", "misD","lof", "prot")),col.names=T,row.names=F,sep="\t",append=F,quote=F)
# LoF_Intolerant genes in 114 cases using the column denovo_call_metaSVM_MPC
denovolyzeByClass(genes=cases$Cases10GeneName,classes= cases$denovo_call_metaSVM_MPC,geneId="gene", nsamples=case_num, probTable=pDNM_cases10, includeGenes=Int_genes, includeClasses=c("all", "syn", "mis", "protD", "misD","lof", "prot"))
write.table(file="LoF_Intolerant_gene_enrichment_114Cases.txt",denovolyzeByClass(genes=cases$Cases10GeneName,classes= cases$denovo_call_metaSVM_MPC,geneId="gene", nsamples=case_num, probTable=pDNM_cases10, includeGenes=Int_genes, includeClasses=c("all", "syn", "mis", "protD", "misD","lof", "prot")),col.names=T,row.names=F,sep="\t",append=F,quote=F)
# HBE_LoF_Intolerant genes in 114 cases using the column denovo_call_metaSVM_MPC
denovolyzeByClass(genes=cases$Cases10GeneName,classes= cases$denovo_call_metaSVM_MPC,geneId="gene", nsamples=case_num, probTable=pDNM_cases10, includeGenes=HBE_Int_genes, includeClasses=c("all", "syn", "mis", "protD", "misD","lof", "prot"))
write.table(file="HBE_LoF_Intolerant_gene_enrichment_114Cases.txt",denovolyzeByClass(genes=cases$Cases10GeneName,classes= cases$denovo_call_metaSVM_MPC,geneId="gene", nsamples=case_num, probTable=pDNM_cases10, includeGenes=HBE_Int_genes, includeClasses=c("all", "syn", "mis", "protD", "misD","lof", "prot")),col.names=T,row.names=F,sep="\t",append=F,quote=F)

##### Generate enrichment analysis for each functional class

# All genes in 1798 controls using the column denovo_call_metaSVM_CADD
denovolyzeByClass(genes=ctrls$Ctrls10GeneName,classes= ctrls$denovo_call,geneId="gene", nsamples=ctrls_num, probTable=pDNM_controls10, includeClasses=c("all", "syn", "mis", "protD", "misD","lof", "prot"))

# HBE genes in 1798 controls using the column denovo_call_metaSVM_CADD
denovolyzeByClass(genes=ctrls$Ctrls10GeneName,classes= ctrls$denovo_call,geneId="gene", nsamples=ctrls_num, probTable=pDNM_controls10, includeGenes=HBE_genes, includeClasses=c("all", "syn", "mis", "protD", "misD","lof", "prot"))

# LoF_Intolerant genes in 1798 controls using the column denovo_call_metaSVM_CADD
denovolyzeByClass(genes=ctrls$Ctrls10GeneName,classes= ctrls$denovo_call,geneId="gene", nsamples=ctrls_num, probTable=pDNM_controls10, includeGenes=Int_genes, includeClasses=c("all", "syn", "mis", "protD", "misD","lof", "prot"))

# HBE_LoF_Intolerant genes in 1798 controls using the column denovo_call_metaSVM_CADD
denovolyzeByClass(genes=ctrls$Ctrls10GeneName,classes= ctrls$denovo_call,geneId="gene", nsamples=ctrls_num, probTable=pDNM_controls10, includeGenes=HBE_Int_genes, includeClasses=c("all", "syn", "mis", "protD", "misD","lof", "prot"))



###############################################################################
# Generate the qq plot for all genes using protein-altering de novo mutations
###############################################################################

data  <- read.table(file="Gene_Level_Significance_MetaSVM_MPC-DMis_114Cases.txt",header=TRUE,sep="\t")

data$min_p <- sapply(1:nrow(data),function(i) min(data[i,c(8,11,14)]) )



pvalue <- append(data$min_p,rep(1, 19347 - nrow(data)))

observed <- sort(as.numeric(pvalue), decreasing = FALSE)
lobs <- -log10(observed)

expected <- c(1:length(observed))
lexp <- -(log10(expected / (length(expected)+1)))

pdf("qqplot_All-Genes_min-p_Poisson-P-values_114-Cases_MetaSVM_MPC-DMis.pdf", width=6, height=6)
plot(c(0,16), c(0,16), col="red", lwd=3, type="l", xlab="Expected (-logP)", ylab="Observed (-logP)", xlim=c(0,16), ylim=c(0,16), las=1, xaxs="i", yaxs="i", bty="l")
points(lexp, lobs, pch=23, cex=.4, bg="black") 
dev.off()
