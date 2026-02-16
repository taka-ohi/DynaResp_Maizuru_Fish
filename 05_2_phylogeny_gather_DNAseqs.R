####
#### R script for Ohigashi et al (2025)
#### Create a fasta file for Maizuru fish
#### 2025.05.21 written by Ohigashi
#### R 4.5.0
#### NOTE: MitoFish complete + partial sequences will be used in this script (https://mitofish.aori.u-tokyo.ac.jp/download/ )


### load packages
library(seqinr); packageVersion("seqinr")
library(Biostrings); packageVersion("Biostrings")


### load data
# taxonomy data of Maizuru fish
tax_sheet <- read.csv("01_dataformatting_out/fishinfo_w_ecology.csv", header = T)

# mitochondrial sequence data from MitoFish
mitoall <- read.fasta("Data/mito-all", seqtype = "DNA") # download mito-all first!


### extract Maizuru fish sequences
## 1. using MitoFish
selected_seqs <- list()

# set rng seed
set.seed(123)

# make a list for primers
# ref https://besjournals.onlinelibrary.wiley.com/doi/10.1111/2041-210X.13485 
primers <- list(
  "12S_MiFish.U" = c("GTCGGTAAAACTCGTGCCAGC", # MiFish-U-F
                     "CATAGTGGGGTATCTAATCCCAGTTTG"), # MiFish-U-R
  "12S_Ac12S" = c("ACTGGGATTAGATACCCCACTATG", # Ac12S
                  "GAGAGTGACGGGCGGTGT"),
  "12S_AcMDB" = c("GCCTATATACCGCCGTCG", # AcMDB07
                  "GTACACTTACCATGTTACGACTT"),
  "16S_Fish16S" = c("GGTCGCCCCAACCRAAG", # Fish16S
                    "CGAGAAGACCCTWTGGAGCTTIAG"),
  "16S_Vert.16S" = c("AGACGAGAAGACCCYdTGGAGCTT", # Vert-16S
                     "GATCCAACATCGAGGTCGTAA"),
  "16S_Fish16S.FD" = c("GACCCTATGGAGCTTTAGAC", # Fish16SF/D-2R
                       "CGCTGTTATCCCTADRGTAACT"),
  "cytb_L14912" = c("TTCCTAGCCATACAYTAYAC",# L14912/H15149
                    "GGTGGCKCCTCAGAAGGACATTTGKCCYCA"),
  "cytb_L14735c2" = c("AAAAACCACCGTTGTTATTCAACTA", # L14735/H15149c2
                      "GCDCCTCARAATGAYATTTGTCCTCA"),
  "cytb_Fish2deg" = c("ACAACTTCACCCCTGCRAAY", # Fish2degCBL/CBR
                      "GATGGCGTAGGCAAATAGGA"),
  "COI_PS1" = c("ACCTGCCTGCCGTATTTGGYGCYTGRGCCGGRATAGT", # PS1
                "ACGCCACCGAGCCARAARCTYATRTTRTTYATTCG")
)

# loop by tax_sheet$FishBase_name
for (fish_name in tax_sheet$FishBase_name) {
  # get sequence names from mito_all, which includes fish_name
  matching_seqs <- mitoall[grepl(fish_name, names(mitoall))]
  
  # skip if there is no sequences matched
  if (length(matching_seqs) == 0) {
    cat("No sequences found for:", fish_name, "\n")
    next
  }
  
  # amplify by gene
  for (gene in names(primers)) {
    primer <- primers[[gene]]  # get the primer
    amplified_seqs <- list()  # a list for amplified sequences for the primer for the species
    
    for (seq in matching_seqs) {
      # write a temporary file
      temp_fasta <- tempfile(fileext = ".fasta")
      write.fasta(seq, fish_name, temp_fasta)
      
      amp_result <- system2("seqkit", args = c(
        "amplicon",
        "-F", primer[1],
        "-R", primer[2],
        "-m 3",
        "-w 0",
        temp_fasta
      ), stdout = TRUE)
      
      # add the amplified sequences to the list
      if (length(amp_result) > 1) {
        amplified_seqs[[length(amplified_seqs) + 1]] <- amp_result[-1]
      }
    }
    
    # randomly select a sequence
    if (length(amplified_seqs) > 0) {
      selected_seq <- sample(unlist(amplified_seqs), 1)
      
      genelist_name <- paste("selected_seqs", gene, sep = "_")
      selected_seqs[[genelist_name]][[fish_name]] <- selected_seq
    }
  }
  
  # remove for memory use
  rm(matching_seqs)
  
  # print the progress
  cat("Finished processing:", fish_name, "\n")
}


# list detected species
detected_sp <- unique(unlist(lapply(selected_seqs, names)))

# get a list of unlisted fish in the mito_all fasta
unlisted_in_mitoall <- setdiff(tax_sheet$FishBase_name, detected_sp)

# NOTE: there are synonyms for 9 of the 10 unlisted species so I will get their sequences next
# only Petroscirtes_springeri does not have any synonym


## 2. get sequences of the unlisted species that have synonyms in the mito-all
# create a data frame for synonym list
synonym.df <- data.frame(
  FishBase_name = c("Ophichthus_altipennis", "Ostorhinchus_semilineatus", "Roa_modesta",
                    "Goniistius_zonatus", "Goniistius_quadricornis", "Paradiplogrammus_enneactis", 
                    "Calliurichthys_japonicus", "Callionymus_curvicornis", "Takifugu_niphobles"), # from the unlisted_in_mitoall
  synonym = c("Ophichthus_zophistius", "Apogon_semilineatus", "Chaetodon_modestus",
              "Cheilodactylus_zonatus", "Cheilodactylus_quadricornis", "Callionymus_enneactis",
              "Callionymus_japonicus", "Repomucenus_curvicornis", "Takifugu_alboplumbeus"),
  source_URL = c("https://ja.wikipedia.org/wiki/%E3%83%9B%E3%82%BF%E3%83%86%E3%82%A6%E3%83%9F%E3%83%98%E3%83%93",
                 "https://www.fishbase.se/Nomenclature/SynonymsList.php?ID=9241&SynCode=170340&GenusName=Ostorhinchus&SpeciesName=semilineatus",
                 "https://fishbase.se/Nomenclature/SynonymsList.php?ID=11668&SynCode=162911&GenusName=Roa&SpeciesName=modesta",
                 "https://en.wikipedia.org/wiki/Spottedtail_morwong",
                 "https://www.zukan-bouz.com/syu/%E3%83%A6%E3%82%A6%E3%83%80%E3%83%81%E3%82%BF%E3%82%AB%E3%83%8E%E3%83%8F",
                 "https://fishbase.se/Nomenclature/SynonymsList.php?ID=17468&SynCode=51782&GenusName=Callionymus&SpeciesName=enneactis",
                 "https://en.wikipedia.org/wiki/Callionymus_japonicus",
                 "https://fishbase.se/Nomenclature/SynonymsList.php?ID=55180&SynCode=134450&GenusName=Callionymus&SpeciesName=curvicornis",
                 "https://ja.wikipedia.org/wiki/%E3%82%AF%E3%82%B5%E3%83%95%E3%82%B0")
)

# loop to get sequences
for (fish_name in unlisted_in_mitoall) {
  # get synonym name
  synonym <- synonym.df$synonym[synonym.df$FishBase_name == fish_name]
  
  if (length(synonym) == 0) {
    cat("No synonyms found for:", fish_name, "\n")
    next
  }
  
  # get sequence names from mito_all, which includes synonym
  matching_seqs <- mitoall[grepl(synonym, names(mitoall))]
  
  # skip if there is no sequences matched
  if (length(matching_seqs) == 0) {
    cat("No sequences found for:", fish_name, "\n")
    next
  }
  
  # amplify by gene
  for (gene in names(primers)) {
    primer <- primers[[gene]]  # get the primer
    amplified_seqs <- list()  # a list for amplified sequences for the primer for the species
    
    for (seq in matching_seqs) {
      # write a temporary file
      temp_fasta <- tempfile(fileext = ".fasta")
      write.fasta(seq, fish_name, temp_fasta)
      
      amp_result <- system2("seqkit", args = c(
        "amplicon",
        "-F", primer[1],
        "-R", primer[2],
        "-m 3",
        "-w 0",
        temp_fasta
      ), stdout = TRUE)
      
      # add the amplified sequences to the list
      if (length(amp_result) > 1) {
        amplified_seqs[[length(amplified_seqs) + 1]] <- amp_result[-1]
      }
    }
    
    # randomly select a sequence
    if (length(amplified_seqs) > 0) {
      selected_seq <- sample(unlist(amplified_seqs), 1)
      
      genelist_name <- paste("selected_seqs", gene, sep = "_")
      selected_seqs[[genelist_name]][[fish_name]] <- selected_seq
    }
  }
  
  # remove for memory use
  rm(matching_seqs)
  
  # print the progress
  cat("Finished processing:", fish_name, "\n")
}

# list detected species again
detected_sp2 <- unique(unlist(lapply(selected_seqs, names)))
length(detected_sp2) # 111 (no sequence for Petroscirtes_springeri and Paradiplogrammus_enneactis)


## 3. extract from a complete genome in the NCBI
# NOTE: download Callionymus enneactis complete genome https://www.ncbi.nlm.nih.gov/nuccore/NC_024194.1
# and put it under "Data/C.enneactis.fasta"

amp_result_manual <- system2("seqkit", args = c(
  "amplicon",
  "-F", "GTCGGTAAAACTCGTGCCAGC", # MiFish-U
  "-R", "CATAGTGGGGTATCTAATCCCAGTTTG",
  "-m 3",
  "-w 0",
  "Data/C.enneactis.fasta" # Callionymus enneactis https://www.ncbi.nlm.nih.gov/nuccore/NC_024194.1
), stdout = TRUE)
# input this to selected_seqs
selected_seqs[["selected_seqs_12S_MiFish.U"]][["Paradiplogrammus_enneactis"]] <- tolower(amp_result_manual[-1])

detected_sp3 <- unique(unlist(lapply(selected_seqs, names)))
length(detected_sp3) # 112 (no sequence for Petroscirtes_springeri)

##### save the sequences (mitogenome based) #####
# save the sequences by r object
saveRDS(selected_seqs, "05_phylogeny_out/seqs_all.mito.regions.obj")

# save the selected sequences as fasta
for (genelist_name in names(selected_seqs)){
  # get sequence data for the list
  seq_data <- selected_seqs[[genelist_name]]
  seq_names <- names(seq_data)
  
  # output
  write.fasta(
    sequences = seq_data,
    names = seq_names,
    file.out = sprintf("05_phylogeny_out/maizuru_%s.fasta", genelist_name)
  )
}


## 4. deal with Petroscirtes springeri
# accession num. for Blenniidae protein coding genome from Hundt et al. (2014) Mol. Phylogenet. Evol.
blenniidae_ac <- read.csv("Data/Blenniidae_accession.csv", header = T)
# get gene names
genes <- colnames(blenniidae_ac)[-1]  

for (gene in genes) {
  # get accession numbers
  acc_numbers <- blenniidae_ac[[gene]]
  species_names <- blenniidae_ac$Species
  
  # get sequence using rentrez
  sequences <- lapply(acc_numbers, function(acc) {
    if (!is.na(acc)) {
      tryCatch({
        # rentrez
        entrez_fetch(db = "nuccore", id = acc, rettype = "fasta")
      }, error = function(e) {
        warning(sprintf("Failed to fetch sequence for %s", acc))
        return(NA)
      })
    } else {
      return(NA)
    }
  })
  
  fasta_content <- mapply(function(seq, species) {
    if (!is.na(seq)) {
      # replace the header with species name only
      sub("^(>.*?\\n)", paste0(">", species, "\n"), seq)
    } else {
      return("")
    }
  }, sequences, species_names, SIMPLIFY = TRUE)
  
  # output as fasta
  file_name <- sprintf("05_phylogeny_out/maizuru_selected_seqs_%s.fasta", gene)
  writeLines(fasta_content, con = file_name)
  
  cat(sprintf("FASTA file created for %s: %s\n", gene, file_name))
}


### remove unused large object
rm(mitoall)


### save session info
writeLines(capture.output(sessionInfo()),
           # please change 0X or XX below to the script number you used.
           sprintf("00_SessionInfo/05_2_SessionInfo_%s.txt", substr(Sys.time(), 1, 10)))






