rm(list=ls())
library(paleomorph)
library(FactoMineR)
library(factoextra)
library(ggplot2)
library(reshape2)


rm(list=ls())


# Import and select the coordinates of triangles made by Prom/EnhA/EnhB for NPC WT
setwd("/Users/ivana.jerkovic/Downloads/Figure 4/Raw Data/")

# Load NPC WT
NPC_WT <- rbind(
  read.csv2("NPC_wt_rep1_merge.csv", header = TRUE, stringsAsFactors = FALSE),
  read.csv2("NPC_wt_rep3_merge.csv", header = TRUE, stringsAsFactors = FALSE)
)
NPC_WT[,-1] <- apply(NPC_WT[,-1], 2, as.numeric)

NPC_WT$F1F2 <- sqrt(((NPC_WT$xFC2 - NPC_WT$xFC1) * 40)^2 + ((NPC_WT$yFC2 - NPC_WT$yFC1) * 40)^2 + ((NPC_WT$zFC2 - NPC_WT$zFC1) * 140)^2)
NPC_WT$F2F3 <- sqrt(((NPC_WT$xFC3 - NPC_WT$xFC2) * 40)^2 + ((NPC_WT$yFC3 - NPC_WT$yFC2) * 40)^2 + ((NPC_WT$zFC3 - NPC_WT$zFC2) * 140)^2)
NPC_WT$F1F3 <- sqrt(((NPC_WT$xFC3 - NPC_WT$xFC1) * 40)^2 + ((NPC_WT$yFC3 - NPC_WT$yFC1) * 40)^2 + ((NPC_WT$zFC3 - NPC_WT$zFC1) * 140)^2)

# --- Updated Cutoff Function (FIX 2) ---
cut_off_mydist <- function(a) {
  a <- a[a$F1F2 %in% subset(a$F1F2, a$F1F2 < 1500), ]
  a <- a[a$F2F3 %in% subset(a$F2F3, a$F2F3 < 1500), ]
  a <- a[a$F1F3 %in% subset(a$F1F3, a$F1F3 < 1500), ]
  return(a)
}

NPC_WT_final <- cut_off_mydist(NPC_WT)

# Load NPC dA
NPC_WT <- rbind(
  read.csv2("NPC_dA_rep1_merge.csv", header = TRUE, stringsAsFactors = FALSE),
  read.csv2("NPC_dA_rep2_merge.csv", header = TRUE, stringsAsFactors = FALSE),
  read.csv2("NPC_dA_rep3_merge.csv", header = TRUE, stringsAsFactors = FALSE)
)
NPC_WT[,-1] <- apply(NPC_WT[,-1], 2, as.numeric)

# --- Corrected Formula (FIX 1: Changed zFC3 - zFC2 to zFC2 - zFC1) ---
NPC_WT$F1F2 <- sqrt(((NPC_WT$xFC2 - NPC_WT$xFC1) * 40)^2 + ((NPC_WT$yFC2 - NPC_WT$yFC1) * 40)^2 + ((NPC_WT$zFC2 - NPC_WT$zFC1) * 140)^2)
NPC_WT$F2F3 <- sqrt(((NPC_WT$xFC3 - NPC_WT$xFC2) * 40)^2 + ((NPC_WT$yFC3 - NPC_WT$yFC2) * 40)^2 + ((NPC_WT$zFC3 - NPC_WT$zFC2) * 140)^2)
NPC_WT$F1F3 <- sqrt(((NPC_WT$xFC3 - NPC_WT$xFC1) * 40)^2 + ((NPC_WT$yFC3 - NPC_WT$yFC1) * 40)^2 + ((NPC_WT$zFC3 - NPC_WT$zFC1) * 140)^2)

NPC_dA_final <- cut_off_mydist(NPC_WT)
rm(NPC_WT)

NPC_WT_final$condition <- "NPC WT"
NPC_dA_final$condition <- "NPC dA"

df <- rbind(NPC_dA_final, NPC_WT_final)

df <- reshape(
  df,
  varying = c("F1F2", "F2F3", "F1F3"),
  v.names = "valeur",
  timevar = "distance",
  times = c("F1F2", "F2F3", "F1F3"),
  direction = "long"
)

df$condition <- factor(df$condition, levels = c("NPC WT", "NPC dA"))
df$distance  <- factor(df$distance, levels = c("F1F2", "F1F3", "F2F3"), labels = c("P-A", "P-B", "A-B"))

# Calculate sample sizes
n_wt <- sum(df$condition == "NPC WT") / length(unique(df$distance))
n_da <- sum(df$condition == "NPC dA") / length(unique(df$distance))

legend_labels <- c(
  "NPC WT" = paste0("NPC WT (n=", n_wt, ")"),
  "NPC dA" = paste0("NPC dA (n=", n_da, ")")
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
  # 1. Split Violins (Width set to exactly 50%: 0.375)
  geom_split_violin(
    aes(color = condition),
    width = 0.375,
    alpha = 0.35,
    trim = TRUE,
    na.rm = TRUE
  ) +
  
  # 2. Boxplots (Width and Dodge set to exactly 50%: 0.05)
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
    values = c("NPC WT" = "#504A50", "NPC dA" = "#81B441"),
    labels = legend_labels,
    name = "Condition"
  ) +
  scale_color_manual(
    values = c("NPC WT" = "#5B535C", "NPC dA" = "#5B535C"),
    labels = legend_labels,
    name = "Condition"
  ) +
  
  # 4. P-Values positioned directly above each split violin
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
    subtitle = "NPC wt vs NPC dA",
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
  filename = "/Users/ivana.jerkovic/Downloads/Figure 4/Raw Data/Figure4b.pdf", 
  plot = p, 
  width = 8, 
  height = 6, 
  units = "in",
  dpi = 300
)

