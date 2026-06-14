
# ============================================================
# Differential Protein Abundance Analysis and Volcano Plot
# ============================================================
#
# Purpose:
# Identify proteins associated with lactogenic treatment and
# visualise differential abundance using a volcano plot.
#
# Input:
#   Sheet 1: Protein abundance matrix
#   Sheet 2: Sample metadata
#
# Output:
#   Limma differential abundance results
#   Publication-ready volcano plot
#
# ============================================================


# ============================================================
# Load Packages
# ============================================================

library(tidyverse)
library(readxl)
library(limma)
library(ggrepel)


# ============================================================
# Import Protein Abundance Data
# ============================================================

df_m <- read_xlsx(
  "mer_proteins.xlsx",
  sheet = 1
)

# Store protein identifiers

pr_col <- df_m[, 1]

# Remove identifier column from abundance matrix

df_m0 <- df_m[, -1]


# ============================================================
# Remove Proteins with Zero Variance
# ============================================================

no_var_rows <- which(
  apply(df_m0, 1, sd) == 0
)

df_m <- df_m0[-no_var_rows, ]

pr_col <- pr_col[-no_var_rows, ]


# ============================================================
# Log Transformation
# ============================================================

df_m <- log(df_m + 1)


# ============================================================
# Import Sample Metadata
# ============================================================

des_data <- read_xlsx(
  "mer_proteins.xlsx",
  sheet = 2
) %>%
  mutate(
    treatment = factor(treatment),
    sample_nr = factor(sample_nr)
  )


# ============================================================
# Create Design Matrix
# ============================================================

d <- model.frame(
  ~ treatment,
  data = des_data
)

design <- model.matrix(
  ~ treatment,
  data = d
)


# ============================================================
# Estimate Correlation Structure
# ============================================================

corfit <- duplicateCorrelation(
  df_m,
  design
)

corfit$consensus


# ============================================================
# Differential Abundance Analysis (limma)
# ============================================================

M <- lmFit(
  df_m,
  design = design,
  correlation = corfit$consensus
)

M <- eBayes(M)

v_data <- topTable(
  M,
  n = 1188,
  genelist = pr_col,
  adjust.method = "BH"
)


# ============================================================
# Volcano Plot Classification
# ============================================================

v_data$diffexpressed <- "NO"

v_data$diffexpressed[
  v_data$logFC > 0.6 &
    v_data$P.Value < 0.05
] <- "UP"

v_data$diffexpressed[
  v_data$logFC < -0.6 &
    v_data$P.Value < 0.05
] <- "DOWN"


# ============================================================
# Protein Labels
# ============================================================

v_data$delabel <- ifelse(
  v_data$Accession %in%
    head(
      v_data[order(v_data$P.Value), "Accession"],
      1188
    ),
  v_data$Accession,
  NA
)

v_data <- v_data %>%
  mutate(
    delabel = fct_recode(
      delabel,
      "FASN" = "Q71SP7",
      "SFN" = "Q0VC36",
      "TGFB1" = "P18341",
      "Cofilin-1" = "Q5E9F7",
      "Calponin-2" = "Q3SYU6",
      "α-actinin-1" = "Q3B7N2",
      "Dystroglycan1" = "O18738",
      "B2M" = "P01888",
      "CXCL6" = "P80221",
      "Lactotransferrin" = "P24627",
      "β-Lactoglobulin" = "P02754",
      "αS1-Casein" = "P02662",
      "αS2-Casein" = "P02663",
      "κ-Casein" = "P02668",
      "β-Casein" = "P02666"
    )
  )


# ============================================================
# Proteins Selected for Annotation
# ============================================================

specific_labels <- c(
  "FASN",
  "SFN",
  "TGFB1",
  "Cofilin-1",
  "Calponin-2",
  "α-actinin-1",
  "Dystroglycan1",
  "B2M",
  "CXCL6",
  "Lactotransferrin",
  "β-Lactoglobulin",
  "αS1-Casein",
  "αS2-Casein",
  "κ-Casein",
  "β-Casein"
)

labeled_data <- v_data[
  v_data$delabel %in% specific_labels,
]


# ============================================================
# Plot Theme
# ============================================================

theme_set(
  theme_classic(base_size = 20) +
    theme(
      axis.title.y = element_text(
        face = "bold",
        margin = margin(0, 20, 0, 0),
        size = rel(1.1),
        color = "black"
      ),
      axis.title.x = element_text(
        hjust = 0.5,
        face = "bold",
        margin = margin(20, 0, 0, 0),
        size = rel(1.1),
        color = "black"
      ),
      plot.title = element_text(
        hjust = 0.5
      )
    )
)


# ============================================================
# Volcano Plot
# ============================================================

p <- ggplot(
  data = v_data,
  aes(
    x = logFC,
    y = -log10(P.Value),
    col = diffexpressed
  )
) +
  geom_vline(
    xintercept = c(-0.6, 0.6),
    linetype = "dashed"
  ) +
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed"
  ) +
  geom_point(size = 2) +
  scale_color_manual(
    values = c(
      "#00AFBB",
      "black",
      "#bb0c00"
    ),
    labels = c(
      "Downregulated",
      "Not significant",
      "Upregulated"
    )
  ) +
  coord_cartesian(
    ylim = c(0, 3),
    xlim = c(-7, 7)
  ) +
  labs(
    color = "",
    x = expression("log"[2] * "FC"),
    y = expression("-log"[10] * "p-value")
  ) +
  scale_x_continuous(
    breaks = seq(-7, 7, 3)
  ) +
  geom_label_repel(
    data = labeled_data,
    aes(label = delabel),
    hjust = 0.8,
    vjust = 0.8,
    box.padding = 0.2,
    point.padding = 0.2,
    max.overlaps = Inf
  )

p


# ============================================================
# Identify Highly Ranked Milk Proteins
# ============================================================

v_data$MP <- ifelse(
  v_data$Accession %in%
    head(
      v_data[order(v_data$P.Value), "Accession"],
      60
    ),
  v_data$Accession,
  NA
)

