
# 1. Load package
library(Biostrings)

# 2. Import the saved alignment
aligned_sequences <- readDNAStringSet(
  "results/aligned_genomes.fasta"
)

# 3. Convert aligned sequences into a character matrix
sequence_matrix <- do.call(
  rbind,
  strsplit(as.character(aligned_sequences), split = "")
)

rownames(sequence_matrix) <- names(aligned_sequences)

# 4. Extract the reference sequence
reference_seq <- sequence_matrix[1, ]

# 5. Define standard nucleotides
valid_bases <- c("A", "C", "G", "T")

# 6. Create an empty results table
mutation_results <- data.frame()

# 7. Compare every sample against the reference
for (i in 2:nrow(sequence_matrix)) {
  
  sample_seq <- sequence_matrix[i, ]
  
  # Positions where both sequences have standard bases
  valid_positions <- reference_seq %in% valid_bases &
    sample_seq %in% valid_bases
  
  # Count substitutions
  substitutions <- sum(
    valid_positions & reference_seq != sample_seq
  )
  
  # Count insertions
  insertions <- sum(
    reference_seq == "-" & sample_seq %in% valid_bases
  )
  
  # Count deletions
  deletions <- sum(
    reference_seq %in% valid_bases & sample_seq == "-"
  )
  
  # Count positions involving ambiguous symbols
  ambiguous_positions <- sum(
    !(reference_seq %in% c(valid_bases, "-")) |
      !(sample_seq %in% c(valid_bases, "-"))
  )
  
  # Store results
  mutation_results <- rbind(
    mutation_results,
    data.frame(
      Accession = names(aligned_sequences)[i],
      Substitutions = substitutions,
      Insertions = insertions,
      Deletions = deletions,
      Ambiguous_Positions = ambiguous_positions
    )
  )
}

# 8. Display mutation results
print(mutation_results)

# 9. Save mutation results
write.csv(
  mutation_results,
  "results/mutation_summary.csv",
  row.names = FALSE
)

# 10. Load visualization package
library(ggplot2)

# 11. Convert mutation results into long format
mutation_long <- data.frame(
  Accession = rep(mutation_results$Accession, 3),
  Mutation_Type = rep(
    c("Substitutions", "Insertions", "Deletions"),
    each = nrow(mutation_results)
  ),
  Count = c(
    mutation_results$Substitutions,
    mutation_results$Insertions,
    mutation_results$Deletions
  )
)

# 12. Create mutation comparison plot
mutation_plot <- ggplot(
  mutation_long,
  aes(x = Accession, y = Count, fill = Mutation_Type)
) +
  geom_col(position = "dodge") +
  labs(
    title = "Mutation Comparison of SARS-CoV-2 Genomes",
    x = "Accession",
    y = "Number of differences",
    fill = "Mutation type"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1)
  )

print(mutation_plot)

# 13. Save the plot
ggsave(
  "Figures/mutation_comparison.png",
  plot = mutation_plot,
  width = 12,
  height = 7
)
