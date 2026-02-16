####
#### R script for Ohigashi et al (2025)
#### alignment for Maizuru fish sequences by gene types, and supermatrix
#### 2025.05.21 written by Ohigashi
#### R 4.5.0
####


### load packages
library(seqinr); packageVersion("seqinr")
library(msa); packageVersion("msa")
library(phangorn); packageVersion("phangorn")
library(phylotools); packageVersion("phylotools")
library(ips); packageVersion("ips")


### load data
seq12S_MiFish.U <- readDNAStringSet("05_phylogeny_out/maizuru_selected_seqs_12S_MiFish.U.fasta")
seq12S_Ac12S <- readDNAStringSet("05_phylogeny_out/maizuru_selected_seqs_12S_Ac12S.fasta")
seq12S_AcMDB <- readDNAStringSet("05_phylogeny_out/maizuru_selected_seqs_12S_AcMDB.fasta")
seq16S_Fish16S.FD <- readDNAStringSet("05_phylogeny_out/maizuru_selected_seqs_16S_Fish16S.FD.fasta")
seq16S_Vert.16S <- readDNAStringSet("05_phylogeny_out/maizuru_selected_seqs_16S_Vert.16S.fasta")
seqcytb_Fish2deg <- readDNAStringSet("05_phylogeny_out/maizuru_selected_seqs_cytb_Fish2deg.fasta")
seqcytb_L14735c2 <- readDNAStringSet("05_phylogeny_out/maizuru_selected_seqs_cytb_L14735c2.fasta")
seqENC1 <- readDNAStringSet("05_phylogeny_out/maizuru_selected_seqs_ENC1.fasta")
seqmyh6 <- readDNAStringSet("05_phylogeny_out/maizuru_selected_seqs_myh6.fasta")
seqptr <- readDNAStringSet("05_phylogeny_out/maizuru_selected_seqs_ptr.fasta")
seqtbr1 <- readDNAStringSet("05_phylogeny_out/maizuru_selected_seqs_tbr1.fasta")

# backbone tree
bbtree <- read.tree("05_phylogeny_out/family_tree_wo_akajei.nwk")


### Alignment
# multiple sequence alignment
mult12S_MiFish.U <- msa(seq12S_MiFish.U, method = "ClustalW", type = "dna", order = "aligned", verbose = TRUE)
mult12S_Ac12S <- msa(seq12S_Ac12S, method = "ClustalW", type = "dna", order = "aligned", verbose = TRUE)
mult12S_AcMDB <- msa(seq12S_AcMDB, method = "ClustalW", type = "dna", order = "aligned", verbose = TRUE)
mult16S_Fish16S.FD <- msa(seq16S_Fish16S.FD, method = "ClustalW", type = "dna", order = "aligned", verbose = TRUE)
mult16S_Vert.16S <- msa(seq16S_Vert.16S, method = "ClustalW", type = "dna", order = "aligned", verbose = TRUE)
multcytb_Fish2deg <- msa(seqcytb_Fish2deg, method = "ClustalW", type = "dna", order = "aligned", verbose = TRUE)
multcytb_L14735c2 <- msa(seqcytb_L14735c2, method = "ClustalW", type = "dna", order = "aligned", verbose = TRUE)
multENC1 <- msa(seqENC1, method = "ClustalW", type = "dna", order = "aligned", verbose = TRUE)
multmyh6 <- msa(seqmyh6, method = "ClustalW", type = "dna", order = "aligned", verbose = TRUE)
multptr <- msa(seqptr, method = "ClustalW", type = "dna", order = "aligned", verbose = TRUE)
multtbr1 <- msa(seqtbr1, method = "ClustalW", type = "dna", order = "aligned", verbose = TRUE)


### save alignment as fasta
writeXStringSet(unmasked(mult12S_MiFish.U), file="05_phylogeny_out/msa12S_MiFish.U.fasta")
writeXStringSet(unmasked(mult12S_Ac12S), file="05_phylogeny_out/msa12S_Ac12S.fasta")
writeXStringSet(unmasked(mult12S_AcMDB), file="05_phylogeny_out/msa12S_AcMDB.fasta")
writeXStringSet(unmasked(mult16S_Fish16S.FD), file="05_phylogeny_out/msa16S_Fish16S.FD.fasta")
writeXStringSet(unmasked(mult16S_Vert.16S), file="05_phylogeny_out/msa16S_Vert.16S.fasta")
writeXStringSet(unmasked(multcytb_Fish2deg), file="05_phylogeny_out/msacytb_Fish2deg.fasta")
writeXStringSet(unmasked(multcytb_L14735c2), file="05_phylogeny_out/msacytb_L14735c2.fasta")
writeXStringSet(unmasked(multENC1), file="05_phylogeny_out/msaENC1.fasta")
writeXStringSet(unmasked(multmyh6), file="05_phylogeny_out/msamyh6.fasta")
writeXStringSet(unmasked(multptr), file="05_phylogeny_out/msaptr.fasta")
writeXStringSet(unmasked(multtbr1), file="05_phylogeny_out/msatbr1.fasta")


###### run supermat to make a supermatrix ######
supermat(infiles = c("05_phylogeny_out/msa12S_MiFish.U.fasta",
                     "05_phylogeny_out/msa12S_Ac12S.fasta",
                     "05_phylogeny_out/msa12S_AcMDB.fasta",
                     "05_phylogeny_out/msa16S_Vert.16S.fasta",
                     "05_phylogeny_out/msa16S_Fish16S.FD.fasta",
                     "05_phylogeny_out/msacytb_Fish2deg.fasta",
                     "05_phylogeny_out/msacytb_L14735c2.fasta",
                     "05_phylogeny_out/msaENC1.fasta",
                     "05_phylogeny_out/msamyh6.fasta",
                     "05_phylogeny_out/msaptr.fasta",
                     "05_phylogeny_out/msatbr1.fasta"
                     ),
         outfile = "05_phylogeny_out/supermat.out_w.B.phy",
         partition.file = "05_phylogeny_out/gene_partition_w.B.txt"
         )

# load the supermatrix
smat <- read.dna("05_phylogeny_out/supermat.out_w.B.phy")


###### run RAxML #####
## format the backbone tree
# currently, it includes Family name so exclude them
bbtree$tip.label <- sub("^.*?_", "", bbtree$tip.label)

# set the execute file (installed for the Ushiolab workstation)
exec <- "/home/taka/standard-RAxML/raxmlHPC-PTHREADS-AVX"

# create the tree
setwd("05_phylogeny_out/")
tr <- raxml(smat, m = "GTRGAMMA", f = "a", N = 2, p = 123,
            outgroup = "Hemitrygon_akajei", x = 1000,
            backbone = bbtree,
            exec = exec, threads = 8)
setwd("../")


### save tree
saveRDS(tr, "05_phylogeny_out/RAxML_res_w.B.rds")
write.tree(tr$bestTree, "05_phylogeny_out/RAxML_tree_w.B.nwk")

pdf("05_phylogeny_out/RAxML_tree_w.B.pdf", height = 20, width = 8) 
plot(tr$bestTree)
add.scale.bar(x = 0, y = -1, length = 0.2) 
dev.off()

pdf("05_phylogeny_out/RAxML_tree_w.B_c.pdf", height = 15, width = 15) # circle version
plot(tr$bestTree, type = "fan")
dev.off() 


### save session info
writeLines(capture.output(sessionInfo()),
           # please change 0X or XX below to the script number you used.
           sprintf("00_SessionInfo/05_3_SessionInfo_%s.txt", substr(Sys.time(), 1, 10)))




