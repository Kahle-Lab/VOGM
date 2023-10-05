library(data.table)
options(stringsAsFactors = F)

# Use mutability for binomial test
setwd("/Users/shujuanzhao/Dropbox/VOGM/Transmitted_mutations/Statistical_Analysis/")
library(dplyr)
library(reshape)
library(devtools)
input = "case"

reformat_pDNM = function(x, mis_filter="Mis_mpc2_or_MetaSVM") # MetaSVM + MPC-2
{
  names(x)[names(x)==mis_filter] <- "misD"
  x$prot <- x[,names(x) %in% c("lof","mis")] %>% apply(MARGIN=1, sum)
  x$protD = x[,names(x) %in% c("lof","misD")] %>% apply(MARGIN=1,sum)
  x <- melt(x)                                          #use anything that is 'chr' as id variables
  names(x)[names(x)=="variable"] <- "class"             #change whichever column with the name 'variable' to 'class'
  return(x)
}

if(input == "case"){
  cases10 = read.delim("0422_hs37d5_coding_idt_med_v2_spikein_padded_Mar2018_adj_19347-unique-gene_mpc-modified.txt",sep="\t",stringsAsFactors = F)
  unique(reformat_pDNM(cases10)[,"class"])
  pDNM<-reformat_pDNM(cases10)
  geneId="gene"
  includeGenes="all"
  probTable=pDNM
  names(probTable)[names(probTable)==geneId] <- "gene"
  probTable$gene <- toupper(as.character(probTable$gene))
  includeGenes <- toupper(as.character(includeGenes))
} else {
  controls10 = read.delim("Padded_control_adj_19347-unique-gene_mpc-modified.txt",sep="\t",stringsAsFactors = F)
  unique(reformat_pDNM(controls10)[,"class"])
  geneId="gene"
  includeGenes="all"
  probTable=pDNM
  names(probTable)[names(probTable)==geneId] <- "gene"
  probTable$gene <- toupper(as.character(probTable$gene))
  includeGenes <- toupper(as.character(includeGenes))
}

# --------------------
# If a list of genes for inclusion is specified, restrict analysis to these genes
probtable_damaging <- probTable[probTable$class=="protD",c("gene","value")]

# Read in the pli scores
PLI = read.table(file="Mutability_pLI_Table.txt",header=T,sep="\t",stringsAsFactors=F)

probtable_damaging$pLI <- unlist(sapply(1:nrow(probtable_damaging), function(k) {
	index <- which(PLI$gene.symbol %in% probtable_damaging$gene[k])
	if(length(index) == 0){
		return (NA)
	} else {
		return (PLI$pLI[index[1]])
	}
}))

### Sum of the probability
#probtable_lof = probtable_lof[-which(is.na(probtable_lof$pLI)),]
#N3 = sum(probtable_lof$value)
#Q.PLI <- quantile(probtable_lof$pLI)[2:4]	# c(25%, 50%, 75%)
#Q.PLI <- quantile(probtable_lof$pLI,na.rm=TRUE)[2:4] # c(25%, 50%, 75%)
Q.PLI <- quantile(probtable_damaging$pLI,na.rm=TRUE)[2:4] # c(25%, 50%, 75%)

### LOF and Dmis
cases <- as.data.frame(fread("114VOGM_LoF_MetaSVM-D_MPC-2_BravoAF_0.00005.txt"))
cases$Case10Gene <- toupper(as.character(cases$Case10Gene))
#cases <- tbl_df(cases)
cases <- tibble::as_tibble(cases)

#####for BarvoAF 0.00002, run the folowing:
ii <- which(colnames(cases) == "bravo")
##convert "." to "0"
na_code <- c(".")
cases[ii] <- apply(cases[ii], 1, function(x){ ifelse(x %in% na_code, 0, x) } )
##filter for bravo AF <= 0.00002
index <- which(as.numeric(cases[[ii]]) <= 0.00002)
cases <- cases[index,]
##convert "0" to "."
cases[ii] <- apply(cases[ii], 1, function(x){ ifelse(x %in% 0, na_code, x) } )

### Select True variants
index <- which(cases$Visualization %in% c("Yes",""))
cases <- cases[index,]

# Only include LoF + MetaSVM-Damaging missense
class <- c("stopgain","stoploss","frameshift_insertion","frameshift_deletion",".")
index <- which( (cases$ExonicFunc.refGene %in% class) | ((cases$ExonicFunc.refGene == "nonsynonymous_SNV") & (cases$MetaSVM_pred == "D")) | ((cases$ExonicFunc.refGene == "nonsynonymous_SNV") & (as.numeric(cases$MPC_score >= 2))) )
cases <- cases[index,]

# Remove duplicated mutations in the same gene within an individual
cases <- mutate(cases, ID.Gene = paste(SAMPLE,Case10Gene,sep=":"))
cases <- distinct(cases, ID.Gene,.keep_all = TRUE)


Q1.genes = probtable_damaging[which(probtable_damaging$pLI >= 0 & probtable_damaging$pLI < Q.PLI[1]),]
Q2.genes = probtable_damaging[which(probtable_damaging$pLI >= Q.PLI[1] & probtable_damaging$pLI < Q.PLI[2]),]
Q3.genes = probtable_damaging[which(probtable_damaging$pLI >= Q.PLI[2] & probtable_damaging$pLI < Q.PLI[3]),]
Q4.genes = probtable_damaging[which(probtable_damaging$pLI >= Q.PLI[3] & probtable_damaging$pLI <= 1),]
Q5.genes = probtable_damaging[which(is.na(probtable_damaging$pLI)),]

cases <- filter(cases, Case10Gene %in% probtable_damaging$gene)
table.cases <- table(cases$Case10Gene)
count.table <- data.frame(names(table.cases),as.numeric(table.cases),stringsAsFactors=F)
colnames(count.table) <- c("gene","count")
#count.table<-count_table(cases)
total.no.mutations <- sum(count.table[,2])
total.no.mutations.Q1 <- sum(count.table[count.table$gene %in% Q1.genes$gene,2])
total.no.mutations.Q2 <- sum(count.table[count.table$gene %in% Q2.genes$gene,2])
total.no.mutations.Q3 <- sum(count.table[count.table$gene %in% Q3.genes$gene,2])
total.no.mutations.Q4 <- sum(count.table[count.table$gene %in% Q4.genes$gene,2])
total.no.mutations.Q5 <- sum(count.table[count.table$gene %in% Q5.genes$gene,2])

# Sum of Q1 genes
N.Q1 = sum(Q1.genes$value)
# Sum of Q2 genes
N.Q2 = sum(Q2.genes$value)
# Sum of Q3 genes
N.Q3 = sum(Q3.genes$value)
# Sum of Q4 genes
N.Q4 = sum(Q4.genes$value)
# Sunm of Q5 genes
N.Q5 = sum(Q5.genes$value)

### Generate Binomial Output
binom.table <- NULL
iterations<-dim(probtable_damaging)[[1]]
pb = txtProgressBar(0, iterations, style=3)
for(a in 1:dim(probtable_damaging)[[1]]) {
  gene.probe <- probtable_damaging$gene[a]
  prob.probe <- probtable_damaging$value[a]
  pli.probe <- probtable_damaging$pLI[a]
  if(length(which(count.table$gene==gene.probe))==0){ # If no mutation
    count <- 0
    if(gene.probe %in% Q1.genes$gene){
      gene.phyper <- binom.test(count, total.no.mutations.Q1, p = prob.probe/N.Q1, alternative = "greater")$p.value
    }
    if(gene.probe %in% Q2.genes$gene){
      gene.phyper <- binom.test(count, total.no.mutations.Q2, p = prob.probe/N.Q2, alternative = "greater")$p.value
    }
    if(gene.probe %in% Q3.genes$gene){
      gene.phyper <- binom.test(count, total.no.mutations.Q3, p = prob.probe/N.Q3, alternative = "greater")$p.value
    }
    if(gene.probe %in% Q4.genes$gene){
      gene.phyper <- binom.test(count, total.no.mutations.Q4, p = prob.probe/N.Q4, alternative = "greater")$p.value
    }
    if(gene.probe %in% Q5.genes$gene){
      gene.phyper <- binom.test(count, total.no.mutations.Q5, p = prob.probe/N.Q5, alternative = "greater")$p.value
    }
  } else{ # If there is at least one mutation
    ii <- which(count.table$gene==gene.probe)
    count <- as.numeric(as.character(count.table$count[ii]))
    if(gene.probe %in% Q1.genes$gene){
      gene.phyper <- binom.test(count, total.no.mutations.Q1, p = prob.probe/N.Q1, alternative = "greater")$p.value
    }
    if(gene.probe %in% Q2.genes$gene){
      gene.phyper <- binom.test(count, total.no.mutations.Q2, p = prob.probe/N.Q2, alternative = "greater")$p.value
    }
    if(gene.probe %in% Q3.genes$gene){
      gene.phyper <- binom.test(count, total.no.mutations.Q3, p = prob.probe/N.Q3, alternative = "greater")$p.value
    }
    if(gene.probe %in% Q4.genes$gene){
      gene.phyper <- binom.test(count, total.no.mutations.Q4, p = prob.probe/N.Q4, alternative = "greater")$p.value
    }
    if(gene.probe %in% Q5.genes$gene){
      gene.phyper <- binom.test(count, total.no.mutations.Q5, p = prob.probe/N.Q5, alternative = "greater")$p.value
    }
  }
  if(gene.probe %in% Q1.genes$gene){ # Calculat the expected transmitted variants
        expected <- (total.no.mutations.Q1*(prob.probe/N.Q1))
  }
  if(gene.probe %in% Q2.genes$gene){
        expected <- (total.no.mutations.Q2*(prob.probe/N.Q2))
  }
  if(gene.probe %in% Q3.genes$gene){
        expected <- (total.no.mutations.Q3*(prob.probe/N.Q3))
  }
  if(gene.probe %in% Q4.genes$gene){
        expected <- (total.no.mutations.Q4*(prob.probe/N.Q4))
  }
	if(gene.probe %in% Q5.genes$gene){
        expected <- (total.no.mutations.Q5*(prob.probe/N.Q5))
  }
  fold.change <- count/expected
  temp <- t(as.matrix(c(gene.probe, prob.probe, pli.probe, count, expected, fold.change, gene.phyper)))
  binom.table <- rbind(binom.table,temp)
  setTxtProgressBar(pb, a)
}
colnames(binom.table) <- c("gene","mutability","pLI","observed","expected","fold_change","binom")

####output
write.table(binom.table,"binomial_transmitted_dominant_LoF_MetaSVM-D_MPC2_All_Genes_114VGAM_MAF_5E-5_SepBypLI_5quartiles_withDNMsRGs.txt",row.names = F, quote = F, sep="\t")
# only selected intolerant genes
binom.table2 <- binom.table
#binom.table2 <- binom.table[which(as.numeric(binom.table[,3]) >= 0.9 & !is.na(binom.table[,3])), ]
total.number.of.genes <- dim(binom.table2)
# Q-Q Plot
observed <- sort(as.numeric(binom.table2[,7]), decreasing = FALSE)
lobs <- -log10(observed)
expected <- c(1:length(observed))
lexp <- -(log10(expected / (length(expected)+1)))
pdf("qqplot_transmitted_dominant_LoF_MetaSVM-D_MPC2_All_Genes_114VGAM_MAF_5E-5_SepBypLI_5quartiles_withDNMsRGs.pdf", width=6, height=6)
plot(c(0,10), c(0,10), col="red", lwd=3, type="l", xlab="Expected (-logP)", ylab="Observed (-logP)", cex=0.5, cex.lab=1.4, cex.axis=1.4, cex.main=1.4,xlim=c(0,10), ylim=c(0,10), las=1, xaxs="i", yaxs="i", bty="l")
points(lexp, lobs, pch=23, cex=.4, bg="black")
dev.off()
