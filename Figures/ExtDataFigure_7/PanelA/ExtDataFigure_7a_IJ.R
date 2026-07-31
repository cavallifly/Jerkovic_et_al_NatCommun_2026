rm(list=ls())
library(paleomorph)
library(FactoMineR)
library(factoextra)
library(ggplot2)
library(reshape2)
library(ggpubr)

# Import and select the coordinates of triangles made by Prom/EnhA/EnhB
setwd("/Users/ivana.jerkovic/Downloads/Extended Data Figure 7/Raw Data/")

# Load ESC WT
ESC_WT <- rbind(
  read.csv2("ESC_WT_Rep1_merge_25Fields.csv", header = TRUE, stringsAsFactors = FALSE),
  read.csv2("ESC_WT_Rep2_merge_38Fields.csv", header = TRUE, stringsAsFactors = FALSE)
)
ESC_WT[,-1] <- apply(ESC_WT[,-1], 2, as.numeric)

ESC_WT$F1F2 <- sqrt(((ESC_WT$xFC2 - ESC_WT$xFC1) * 40)^2 + ((ESC_WT$yFC2 - ESC_WT$yFC1) * 40)^2 + ((ESC_WT$zFC2 - ESC_WT$zFC1) * 140)^2)
ESC_WT$F2F3 <- sqrt(((ESC_WT$xFC3 - ESC_WT$xFC2) * 40)^2 + ((ESC_WT$yFC3 - ESC_WT$yFC2) * 40)^2 + ((ESC_WT$zFC3 - ESC_WT$zFC2) * 140)^2)
ESC_WT$F1F3 <- sqrt(((ESC_WT$xFC3 - ESC_WT$xFC1) * 40)^2 + ((ESC_WT$yFC3 - ESC_WT$yFC1) * 40)^2 + ((ESC_WT$zFC3 - ESC_WT$zFC1) * 140)^2)

# --- Cutoff Function ---
cut_off_mydist <- function(a) {
  a <- a[a$F1F2 %in% subset(a$F1F2, a$F1F2 < 1500), ]
  a <- a[a$F2F3 %in% subset(a$F2F3, a$F2F3 < 1500), ]
  a <- a[a$F1F3 %in% subset(a$F1F3, a$F1F3 < 1500), ]
  return(a)
}

ESC_WT_final <- cut_off_mydist(ESC_WT)
rm(ESC_WT)

# Load NPC WT
NPC_WT <- rbind(
  read.csv2("NPC_WT_Rep1_merge_21Fields.csv", header = TRUE, stringsAsFactors = FALSE),
  read.csv2("NPC_WT_Rep2_merge_34Fields.csv", header = TRUE, stringsAsFactors = FALSE)
)
NPC_WT[,-1] <- apply(NPC_WT[,-1], 2, as.numeric)

# --- Calculate distances for NPC WT (FIXED DATA MIX-UP HERE) ---
NPC_WT$F1F2 <- sqrt(((NPC_WT$xFC2 - NPC_WT$xFC1) * 40)^2 + ((NPC_WT$yFC2 - NPC_WT$yFC1) * 40)^2 + ((NPC_WT$zFC2 - NPC_WT$zFC1) * 140)^2)
NPC_WT$F2F3 <- sqrt(((NPC_WT$xFC3 - NPC_WT$xFC2) * 40)^2 + ((NPC_WT$yFC3 - NPC_WT$yFC2) * 40)^2 + ((NPC_WT$zFC3 - NPC_WT$zFC2) * 140)^2)
NPC_WT$F1F3 <- sqrt(((NPC_WT$xFC3 - NPC_WT$xFC1) * 40)^2 + ((NPC_WT$yFC3 - NPC_WT$yFC1) * 40)^2 + ((NPC_WT$zFC3 - NPC_WT$zFC1) * 140)^2)

NPC_WT_final <- cut_off_mydist(NPC_WT)
rm(NPC_WT)

ESC_WT_final$condition <- "ESC WT"
NPC_WT_final$condition <- "NPC WT"

df <- rbind(ESC_WT_final, NPC_WT_final)

df <- reshape(
  df,
  varying = c("F1F2", "F2F3", "F1F3"),
  v.names = "valeur",
  timevar = "distance",
  times = c("F1F2", "F2F3", "F1F3"),
  direction = "long"
)

df$condition <- factor(df$condition, levels = c("ESC WT", "NPC WT"))
df$distance  <- factor(df$distance, levels = c("F1F2", "F1F3", "F2F3"), labels = c("P-A", "P-B", "A-B"))

# ==============================================================================
# --- NEW: STATISTICAL COMPARISONS (P-A vs A-B within each condition) ---
# ==============================================================================

cat("\n--- Within-Sample Wilcoxon Test: P-A vs A-B ---\n")
# 1. ESC WT: P-A vs A-B (Un paired test on same locus triplets)
p_val_esc_unpaired <- wilcox.test(ESC_WT_final$F1F2, ESC_WT_final$F2F3, paired = FALSE)$p.value
cat(sprintf("ESC WT (Unpaired Wilcoxon) : p-value = %g\n", p_val_esc_unpaired))

# 2. NPC WT: P-A vs A-B (Unpaired test on same locus triplets)
p_val_npc_unpaired <- wilcox.test(NPC_WT_final$F1F2, NPC_WT_final$F2F3, paired = FALSE)$p.value
cat(sprintf("NPC WT (Unpaired Wilcoxon) : p-value = %g\n", p_val_npc_unpaired))
cat("--------------------------------------------------\n\n")

# Calculate sample sizes
n_wt <- sum(df$condition == "ESC WT") / length(unique(df$distance))
n_da <- sum(df$condition == "NPC WT") / length(unique(df$distance))

legend_labels <- c(
  "ESC WT" = paste0("ESC WT (n=", n_wt, ")"),
  "NPC WT" = paste0("NPC WT (n=", n_da, ")")
)

# --- PERFECT SEAMLESS SPLIT VIOLIN GEOM ---
GeomSplitViolin <- ggproto(
  "GeomSplitViolin", GeomViolin,
  draw_group = function(self, data, ..., draw_quantiles = NULL) {
    data <- transform(data, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
    grp <- data$group[1]
    
    # Split precisely on x without offsetting
    if (data$fill[1] == "#504A50") {
      newdata <- transform(data, x = xminv)
    } else {
      newdata <- transform(data, x = xmaxv)
      newdata <- newdata[rev(seq_len(nrow(newdata))), ]
    }
    
    # Close polygon flat along central x line
    newdata <- rbind(
      newdata[1, ], 
      newdata, 
      newdata[nrow(newdata), ]
    )
    newdata$x[c(1, nrow(newdata))] <- data$x[1]
    
    GeomPolygon$draw_panel(newdata, ...)
  }
)

geom_split_violin <- function(mapping = NULL, data = NULL, stat = "ydensity", 
                              position = "identity", ..., draw_quantiles = NULL, 
                              trim = TRUE, scale = "area", na.rm = FALSE, 
                              show.legend = NA, inherit.aes = TRUE) {
  layer(
    data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin, 
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles, na.rm = na.rm, ...)
  )
}

# --- PLOT GENERATION ---

p <- ggplot(df, aes(x = distance, y = valeur, fill = condition)) +
  # 1. Split Violins
  geom_split_violin(
    aes(color = condition),
    width = 0.375,
    alpha = 0.35,
    trim = TRUE,
    na.rm = TRUE
  ) +
  
  # 2. Boxplots
  geom_boxplot(
    width = 0.05,
    color = "#333333",
    position = position_dodge(width = 0.05),
    outlier.size = 0.4,
    outlier.alpha = 0.4,
    fatten = 1.2,
    show.legend = FALSE,
    na.rm = TRUE
  ) +
  
  # 3. Scale definitions
  scale_fill_manual(
    values = c("ESC WT" = "#504A50", "NPC WT" = "#81B441"),
    labels = legend_labels,
    name = "Condition"
  ) +
  scale_color_manual(
    values = c("ESC WT" = "#5B535C", "NPC WT" = "#5B535C"),
    labels = legend_labels,
    name = "Condition"
  ) +
  
  # 4. P-Values
  stat_compare_means(
    aes(group = condition),
    method = "wilcox.test",
    label = "p.format",
    label.y = 1550,
    vjust = -0.5,
    size = 4,
    na.rm = TRUE
  ) +
  
  # 5. Expand y-axis limits to fit p-values cleanly
  scale_y_continuous(limits = c(0, 1700), expand = expansion(mult = c(0, 0.05))) +
  
  # 6. Styling and Theme
  labs(
    title = "Pairwise Loci Distance Distributions",
    subtitle = "ESC WT vs NPC WT",
    x = "",
    y = "Distance (nm)"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 11, color = "grey30"),
    legend.position = "top",
    axis.text.x = element_text(color = "black", size = 11, face = "bold"),
    axis.text.y = element_text(color = "black")
  )

# Display Plot
p

# Save PDF directly to target folder
ggsave(
  filename = "/Users/ivana.jerkovic/Downloads/Extended Data Figure 7/ExtDataFigure7a.pdf", 
  plot = p, 
  width = 8, 
  height = 6, 
  units = "in",
  dpi = 300
)