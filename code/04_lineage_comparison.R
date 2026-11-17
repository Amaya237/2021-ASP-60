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