# Load libraries
library("here")
library("ShortRead")
library("dada2")
library("ggplot2")
library("msa")
library("phangorn")
package.version("phangorn")
library("phyloseq")

setwd("[your working directory]")
current_wd <- getwd()

#be sure to create these paths prior to running this script
pathF <- file.path(current_wd,"reads")
filtPath <- file.path(current_wd, "filtered_samples")
fastqFs <- sort(list.files(pathF,pattern="fastq"))

p.qual.f <- plotQualityProfile(pathF[1]) + ggtitle("fwd")
p.qual.f

filterAndTrim(fwd=file.path(pathF,fastqFs),
              filt=file.path(filtPath,fastqFs),
              trimLeft=10,
              truncLen=225,
              maxEE=2,
              truncQ=11,
              maxN=0,
              rm.phix=TRUE,
              compress=TRUE,
              verbose=TRUE,
              multithread=TRUE)

filtFs <- list.files(filtPath,patter="fastq",full.names=TRUE)
sample.names <- sapply(strsplit(basename(filtFs),"\\\\."),`[`,1)
names(filtFs) <- sample.names

set.seed(212)
errF <- learnErrors(filtFs[2],nreads=2e6)

p.err.F <- plotErrors(errF,nominalQ=TRUE)
p.err.F

singles <- vector("list",length(sample.names))
names(singles) <- sample.names
for(sam in sample.names) {
  cat("Processing:",sam,'\\n')
  derepF <- derepFastq(filtFs[[sam]])
  singles[[sam]] <- dada(derepF,err=errF,multithread=TRUE)
}
rm(derepF);

seqtab <- makeSequenceTable(singles)
seqtab.nochim <- removeBimeraDenovo(seqtab,multithread=TRUE)
saveRDS(seqtab.nochim,"seqtab.nochim.rds")

taxa.rdp <- assignTaxonomy(seqtab.nochim,"rdp_files/rdp_train_set_18.fa.gz",multithread=TRUE)
unname(head(taxa.rdp))
colnames(taxa.rdp) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")

taxa.rdp.plus <- addSpecies(taxa.rdp,"rdp_species_assignment_18.fa.gz")
unname(head(taxa.rdp.plus))
colnames(taxa.rdp.plus) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

save.image(here("preprocessing.RData"))

seqs <- getSequences(seqtab.nochim)
names(seqs) <- seqs

mult <- msa(seqs,method="ClustalW",type="dna",order="input")

phang.align <- as.phyDat(mult,type="DNA",names=getSequences(seqtab.nochim))

dm <- dist.ml(phang.align)
treeNJ <- NJ(dm)

fit <- pml(treeNJ,data=phang.align)

fitGTR <- update(fit,k=4,inv=0.2)
fitGTR <- optim.pml(fitGTR,model="GTR",optInv=TRUE,optGamma=TRUE,
                    rearrangement="stochastic",control=pml.control(trace=0))

detach("package:phangorn",unload=TRUE)
save.image(here("preprocessing.RData"))

ps0.rdp <- phyloseq(otu_table(seqtab.nochim,taxa_are_rows=FALSE),tax_table(taxa.rdp.plus),phy_tree(fitGTR$tree))
ps0.rdp

saveRDS(ps0.rdp, file = file.path(current_wd,"data", "rdp_single.RDS"))

save.image(here("preprocessing.RData"))
