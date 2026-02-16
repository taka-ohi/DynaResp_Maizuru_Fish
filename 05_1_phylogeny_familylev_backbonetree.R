####
#### R script for Ohigashi et al (2025)
#### Create a backbone phylogenetic tree for Maizuru fish (just using Family-level topology)
#### 2025.05.21 written by Ohigashi
#### R 4.5.0
#### 


### load packages
library(dplyr); packageVersion("dplyr")
library(ape); packageVersion("ape")


### load data
# taxonomy data of Maizuru fish
tax_sheet <- read.csv("01_dataformatting_out/fishinfo_w_ecology.csv", header = T)

# phylogenetic tree of bony fish
# obtained from Betancur-R et al. (2017) https://github.com/projectdigest/betancur_r-fish-tree
full.tree <- read.tree("Data/12862_2017_958_MOESM2_ESM.tre")


### prepare for making a tree
# extract family from the full tree
family_labels <- sapply(strsplit(full.tree$tip.label, "_"), `[`, 1)

# choose only one label from the same family
unique_family <- unique(family_labels)

# extract family in Maizuru
maizuru_family <- unique_family[unique_family %in% tax_sheet$Family]

# delete tips unrelated to Maizuru, and create a Family-level tree
family_tree <- drop.tip(full.tree, setdiff(full.tree$tip.label, full.tree$tip.label[match(maizuru_family, family_labels)]))

# change the tip names of the Family-level tree to only Family name (no use Genus Species etc.)
family_tree$tip.label <- sapply(strsplit(family_tree$tip.label, "_"), `[`, 1)


### create a tree of Maizuru fish in Family level
# get a list of species by combining Family and Scientific name (Genus+Species)
species_to_add <- tax_sheet %>%
  mutate(label = paste(Family, FishBase_name, sep = "_")) %>%
  pull(label)
# remove akajei (Fish001) because it's not bony
species_to_add <- species_to_add[-1]


# get lists of "Family_Genus_Species" names for each Family
extended_tips <- list()
for (family in unique(tax_sheet$Family)) {
  # get species in the family
  species_in_family <- species_to_add[grep(family, species_to_add)]
  
  # check if there family name in 
  matching_tips <- grep(family, family_tree$tip.label, value = TRUE)
  
  # replace tips with species name
  if (length(matching_tips) > 0) {
    extended_tips[[family]] <- species_in_family
  }
}


# extend the Family-level tree with multiple branches of Species in each Family
new_tree <- family_tree
for (family in names(extended_tips)) {
  # create a sub tree for the family
  sub_tree <- stree(length(extended_tips[[family]]), tip.label = extended_tips[[family]])
  # set branch length of the sub tree (pseudo)
  sub_tree$edge.length <- rep(0.1, nrow(sub_tree$edge))
  
  # bind the sub tree with the family tree
  new_tree <- bind.tree(
    x = new_tree,
    y = sub_tree,
    where = which(new_tree$tip.label == family)
  )
  new_tree <- drop.tip(new_tree, family)
}


### save data
dir.create("05_phylogeny_out")

# save the new tree (w/o akajei)
write.tree(new_tree, file = "05_phylogeny_out/family_tree_wo_akajei.nwk")

# save pdf
pdf("05_phylogeny_out/family_tree_wo_akajei.pdf", height = 20, width = 8)
plot(new_tree)
dev.off()
pdf("05_phylogeny_out/family_tree_wo_akajei_c.pdf", height = 15, width = 15) # circle version
plot(new_tree, type = "fan")
dev.off() 


### save session info
writeLines(capture.output(sessionInfo()),
           # please change 0X or XX below to the script number you used.
           sprintf("00_SessionInfo/05_1_SessionInfo_%s.txt", substr(Sys.time(), 1, 10)))

