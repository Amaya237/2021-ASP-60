
library(ggplot2)

lineage_results <- read.csv(
  "results/lineage_distance_summary.csv"
)

lineage_plot <- ggplot(
  lineage_results,
  aes(
    x = Lineage,
    y = Average_Genetic_Distance,
    fill = Lineage
  )
) +
  geom_col(width = 0.65) +
  labs(
    title = "Genetic Variation Within SARS-CoV-2 Lineages",
    x = "Viral lineage",
    y = "Average genetic distance"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold")
  )

print(lineage_plot)

ggsave(
  "Figures/lineage_genetic_distance.png",
  plot = lineage_plot,
  width = 9,
  height = 6
)
