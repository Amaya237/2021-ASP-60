# 1. Load the required package
library(Biostrings)

# 2. Import genome sequences
samples <- readDNAStringSet("data/raw/sample_genomes.fasta")
reference <- readDNAStringSet("data/raw/reference.fasta")

# 3. Import metadata
metadata <- read.csv(
  "data/raw/genome-metadata.csv",
  stringsAsFactors = FALSE
)

# 4. Check the number of records
length(samples)
length(reference)
nrow(metadata)

# 5. Check metadata columns
colnames(metadata)

# 6. View sample sequence names
head(names(samples))

# 7. View metadata
head(metadata)

# 8. Check lineage distribution
table(metadata$Pangolin)

# 9. Calculate nucleotide frequencies
base_counts <- letterFrequency(
  samples,
  letters = c("A", "C", "G", "T")
)

# 10. Calculate GC content
gc_content <- (
  base_counts[, "G"] + base_counts[, "C"]
) / rowSums(base_counts) * 100

# 11. Calculate ambiguous bases
ambiguous_bases <- width(samples) - rowSums(base_counts)

# 12. Create a quality-control table
qc_results <- data.frame(
  Accession = sub("\\..*", "", names(samples)),
  Genome_Length = width(samples),
  GC_Content = round(gc_content, 2),
  Ambiguous_Bases = ambiguous_bases
)

# 13. Display the results
print(qc_results)

# 14. Summary statistics
summary(qc_results)

# 15. Inspect all nucleotide symbols
all_counts <- alphabetFrequency(samples, baseOnly = FALSE)

# 16. Display counts of ambiguous symbols
print(all_counts[, c("N", "M", "R", "W", "S",
                     "Y", "K", "V", "H", "D", "B")])

# 17. Calculate percentage of ambiguous bases
qc_results$Ambiguous_Percentage <- round(
  (qc_results$Ambiguous_Bases / qc_results$Genome_Length) * 100,
  2
)

# 18. Classify sequence quality for screening
qc_results$QC_Status <- ifelse(
  qc_results$Ambiguous_Percentage > 5,
  "Review",
  "Pass"
)

# 19. Display updated QC results
print(qc_results)

# 20. Count QC statuses
table(qc_results$QC_Status)

# 21. Save quality-control results
write.csv(
  qc_results,
  "data/raw/genome_quality_control.csv",
  row.names = FALSE
)

# 22. Load DECIPHER
library(DECIPHER)

# 23. Combine reference and sample genomes
all_sequences <- c(reference, samples)

# 24. Assign accession names
names(all_sequences) <- c(
  "NC_045512",
  sub("\\..*", "", names(samples))
)

# 25. Check the combined sequences
length(all_sequences)
head(names(all_sequences))

# 26. Perform multiple sequence alignment
aligned_sequences <- AlignSeqs(
  all_sequences,
  processors = 1,
  verbose = TRUE
)

# 27. Check alignment dimensions
length(aligned_sequences)
width(aligned_sequences)

# 28. Check the first few sequence names
head(names(aligned_sequences))

# 29. Create results folder
dir.create("results", showWarnings = FALSE)

# 30. Save aligned sequences as FASTA
writeXStringSet(
  aligned_sequences,
  "results/aligned_genomes.fasta",
  format = "fasta"
)