# Load required libraries
library(Biostrings)
library(ape)
library(phangorn)

# Load aligned genome sequences
aligned_sequences <- readDNAStringSet(
  "results/aligned_genomes.fasta"
)

# Convert aligned sequences into a character matrix
sequence_matrix <- do.call(
  rbind,
  strsplit(as.character(aligned_sequences), split = "")
)

rownames(sequence_matrix) <- names(aligned_sequences)

# Convert to DNA format
dna_data <- ape::as.DNAbin(sequence_matrix)

# Calculate pairwise genetic distances
distance_matrix <- ape::dist.dna(
  dna_data,
  model = "raw",
  pairwise.deletion = TRUE
)

# Check distance matrix
cat("Missing distances:", sum(is.na(distance_matrix)), "\n")

cat(
  "All distances finite:",
  all(is.finite(as.matrix(distance_matrix))),
  "\n"
)

print(distance_matrix)

# Construct Neighbor-Joining tree
tree <- ape::nj(distance_matrix)

# Load metadata
metadata <- read.csv(
  "data/raw/genome-metadata.csv",
  stringsAsFactors = FALSE
)

# Match tree accessions with metadata
tree_accessions <- sub("\\..*", "", tree$tip.label)
metadata$Accession <- sub("\\..*", "", metadata$Accession)

lineage <- metadata$Pangolin[
  match(tree_accessions, metadata$Accession)
]

# Assign colours to lineages
lineage_colours <- c(
  "B.1.1.7" = "red",
  "B.1.617.2" = "blue",
  "BA.1" = "darkgreen",
  "BA.2" = "purple"
)

tip_colours <- unname(lineage_colours[lineage])

# Make the reference genome black
tip_colours[tree_accessions == "NC_045512"] <- "black"

# Plot the coloured tree
png(
  "Figures/phylogenetic_tree_lineages.png",
  width = 1800,
  height = 1200
)

plot(
  tree,
  main = "Phylogenetic Tree of SARS-CoV-2 Genomes by Lineage",
  tip.color = tip_colours,
  cex = 0.8
)

legend(
  "topleft",
  legend = c(
    "Reference",
    "Alpha (B.1.1.7)",
    "Delta (B.1.617.2)",
    "Omicron BA.1",
    "Omicron BA.2"
  ),
  col = c("black", "red", "blue", "darkgreen", "purple"),
  pch = 19,
  cex = 0.8,
  bty = "n"
)

dev.off()
# Save the tree
write.tree(
  tree,
  file = "results/phylogenetic_tree.nwk"
)

# Plot the tree
png(
  "Figures/phylogenetic_tree.png",
  width = 1600,
  height = 1200
)

plot(
  tree,
  main = "Phylogenetic Tree of SARS-CoV-2 Genomes",
  cex = 0.7
)

dev.off()