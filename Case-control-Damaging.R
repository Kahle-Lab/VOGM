setwd("/Users/shujuanzhao/Dropbox/VOGM/case_control/v1/")
library(data.table)
options(stringsAsFactors = F)

# Read in the synonymous mutations in cases
data <- as.data.frame(fread("114VOGM_LoF_MetaSVM-D_MPC-2_BravoAF_0.00005.txt"))

# Read in the count table for synonymous mutations in gnomAD
gnomAD.damaging <- as.data.frame(fread("gnomAD_WGS_WES_0.00005_bravo_0.005_Incohort_MetaSVMD_MPC2_LoF_Eth_All_MQ40_method_sum_Data_Non_Topmed_gene_count.txt"))
#gnomAD_WGS_WES_0.00005_bravo_0.005_Incohort_DmisMPCLoF_Eth_All_MQ40_method_sum_Data_Non_Topmed_sum.txt, if this is the case, no need to remove header, see below
gnomAD.damaging$Gene <- toupper(gnomAD.damaging$Gene)

# Remove headers (due to Weilai's script only combining files from chr1-22)
index <- which(gnomAD.damaging$Gene == "GENE")
gnomAD.damaging <- gnomAD.damaging[-index,]

cases <- data
cases$gnomAD_gene <- toupper(cases$gnomAD_gene)
cases$Case10Gene <- toupper(cases$Case10Gene)

####Remove MQ<40
index <- which(as.numeric(cases$MQ) >= 40)
cases <- cases[index,]

# Only include LoF + MetaSVM-Damaging missense
class <- c("stopgain","stoploss","frameshift_insertion","frameshift_deletion",".")
index <- which( (cases$ExonicFunc.refGene %in% class) | ((cases$ExonicFunc.refGene == "nonsynonymous_SNV") & (cases$MetaSVM_pred == "D")) | ((cases$ExonicFunc.refGene == "nonsynonymous_SNV") & (as.numeric(cases$MPC_score >= 2))) )
cases <- cases[index,]

# Calculate burden
table.cases = table(cases$gnomAD_gene)
out.table.cases =  data.frame(names(table.cases),as.numeric(table.cases),stringsAsFactors=F)
##remove NAs, already removed
#out.table.cases <- out.table.cases[- which(out.table.cases[[1]] == "N/A"),]
cases10gene <- sapply( 1:nrow(out.table.cases), function(i) cases$Case10Gene[which(cases$gnomAD_gene == out.table.cases[[1]][i])[1]] )
out.table.cases$Case10Gene <- cases10gene


# Read in the file with all genes on the mutability table
mut_genes <- as.data.frame(fread("/Users/shujuanzhao/Dropbox/VOGM/case_control/v1/Genes-in-MutabilityTable.txt"))[[1]]
mut_genes <- toupper(mut_genes)

# Figure out which gene does not belong to "out.table.cases" and fill out the remaining genes on the mutability table
index <- which(!mut_genes %in% out.table.cases$Case10Gene)
add.data.frame <- data.frame( mut_genes[index], rep(0,length(index)), mut_genes[index] )
out.table.cases <- rbind(out.table.cases, setNames(add.data.frame,names(out.table.cases)) )
colnames(out.table.cases)[1:2] <- c("gnomAD_gene", "count")

first.line = c("Gene","N_alt_cases","N_ref_cases","N_alt_ctrls","N_ref_ctrls","OR","L95","U95","P")

output.file1 = "Genome_wide_Burden_114Cases_versus_gnomAD_OnlyPASS_DamagingVariants_MAF0.00005.txt"
write.table(t(first.line),file=output.file1,col.names=F,row.names=F,sep="\t",quote=F)

####VOGM cases
N.cases <- 114

# loop through each gene
for(i in 1:nrow(out.table.cases)){
  gene <- out.table.cases[i,1]
  n.alt.cases = as.numeric(out.table.cases$count[i])
  n.ref.cases = (2*N.cases) - n.alt.cases

  # in control
  index <-  which(gnomAD.damaging$Gene %in% gene)
  if(length(index) != 0){
    n.alt.ctrls = as.numeric(gnomAD.damaging$Het_AC[index])
    n.ref.ctrls = as.numeric(gnomAD.damaging$Het_AN[index]) - as.numeric(gnomAD.damaging$Het_AC[index])

  } else {
    n.alt.ctrls = 0
    n.ref.ctrls = 2*135743 # non-TopMed (v2)
  }

  model = fisher.test( matrix(c( n.alt.cases, n.alt.ctrls, n.ref.cases, n.ref.ctrls),2) , alternative = "greater")
  output = c( gene, n.alt.cases, n.ref.cases, n.alt.ctrls, n.ref.ctrls, model$estimate,
               model$conf.int[1], model$conf.int[2], model$p.value )
  write.table(t(output), file=output.file1, col.names=F, row.names=F, sep="\t", append=T, quote=F)
}
# Read in the table
data <- as.data.frame(fread(output.file1))


# Q-Q Plot
observed <- sort(as.numeric(data$P), decreasing = FALSE)
lobs <- -log10(observed)

expected <- c(1:length(observed))
lexp <- -(log10(expected / (length(expected)+1)))

r=cor(lobs,lexp)

pdf("qqplot_Damaging_114Cases_versus_gnomAD_OnlyPASS_MAF0.00005.pdf", width=6, height=6)
plot(c(0,10), c(0,10), col="red", lwd=3, type="l", xlab="Expected (-logP)", ylab="Observed (-logP)", cex=0.5, cex.lab=1.4, cex.axis=1.4, cex.main=1.4,xlim=c(0,10), ylim=c(0,10), las=1, xaxs="i", yaxs="i", bty="l")
points(lexp, lobs, pch=23, cex=.4, bg="black")
dev.off()
