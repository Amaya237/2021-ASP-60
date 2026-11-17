library(ape)

# Load phylogenetic tree and metadata
tree <- read.tree("results/phylogenetic_tree.nwk")

metadata <- read.csv(
  "data/raw/genome-metadata.csv",
  stringsAsFactors = FALSE
)

# Match tree accessions with metadata
lineage <- metadata$Pangolin[
  match(tree$tip.label, metadata$Accession)
]

# Check matching
print(data.frame(
  Accession = tree$tip.label,
  Lineage = lineage
))

# Define lineage colours
lineage_levels <- unique(metadata$Pangolin)

lineage_colors <- setNames(
  rainbow(length(lineage_levels)),
  lineage_levels
)

# Assign colours to tree tips
tip_colors <- rep("black", length(tree$tip.label))

matched <- !is.na(lineage)

tip_colors[matched] <- lineage_colors[lineage[matched]]

# Save lineage-coloured tree
png(
  "Figures/lineage_phylogenetic_tree.png",
  width = 1800,
  height = 1200
)

plot(
  tree,
  main = "Phylogenetic Tree by SARS-CoV-2 Lineage",
  tip.color = tip_colors,
  cex = 0.8
)

legend(
  "topleft",
  legend = names(lineage_colors),
  col = lineage_colors,
  pch = 19,
  cex = 0.8,
  bty = "n"
)

dev.off()


# Calculate pairwise genetic distances
aligned_sequences <- Biostrings::readDNAStringSet(
  "results/aligned_genomes.fasta"
)

sequence_matrix <- do.call(
  rbind,
  strsplit(as.character(aligned_sequences), split = "")
)

rownames(sequence_matrix) <- names(aligned_sequences)

dna_data <- ape::as.DNAbin(sequence_matrix)

distance_matrix <- ape::dist.dna(
  dna_data,
  model = "raw",
  pairwise.deletion = TRUE
)

# Calculate average within-lineage distances
distance_values <- as.matrix(distance_matrix)

lineage_results <- data.frame()

for (group in unique(metadata$Pangolin)) {
  
  accessions <- metadata$Accession[
    metadata$Pangolin == group
  ]
  
  group_distances <- distance_values[
    accessions,
    accessions,
    drop = FALSE
  ]
  
  average_distance <- mean(
    group_distances[upper.tri(group_distances)]
  )
  
  lineage_results <- rbind(
    lineage_results,
    data.frame(
      Lineage = group,
      Average_Genetic_Distance = average_distance
    )
  )
}

print(lineage_results)

write.csv(
  lineage_results,
  "results/lineage_distance_summary.csv",
  row.names = FALSE
)


# Interpretation:
# BA.1 showed the lowest average within-lineage
# genetic distance (0.000266).
#
# BA.2 also showed relatively low genetic variation.
#
# Alpha showed a higher average distance than
# both Omicron groups.
#
# Delta showed the highest average within-lineage
# genetic distance (0.001163).
#
# These findings apply only to the selected samples.

