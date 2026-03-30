#NMDS definitius, pero les dades de abd rel mitjanes no estan bé!
library(tidyverse)
library(dplyr)
library(tibble)
library(tidyr)
library(ggplot2)
library(permute)
library(vegan)
library(geometry)

###AMB COOOI
metadata<- read.xlsx("metadata_repliques_f.xlsx")
arel_mitjCOI<- read.xlsx("Clean_Data_COI_Reps_AbRel.xlsx")
rownames(arel_mitjCOI) <- arel_mitjCOI$X1

arel_mitjCOI1 <- arel_mitjCOI[,c(2:127)]

arel_mitjCOI1_t<- as.data.frame(t(arel_mitjCOI1))

arel_mitjCOI1_t$Sample<- rownames(arel_mitjCOI1_t)

arel_mitjCOI1_t[, 1:(ncol(arel_mitjCOI1_t)-1)] <- lapply(arel_mitjCOI1_t[, 1:(ncol(arel_mitjCOI1_t)-1)], as.numeric)


any(is.na(arel_mitjCOI1_t))     
which(is.na(arel_mitjCOI1_t))
      
dist_COI <- vegdist(arel_mitjCOI1_t[, -ncol(arel_mitjCOI1_t)], method = "bray")

any(is.na(dist_COI))     
which(is.na(dist_COI))   


#k=2 dimensions
nmds2_COI_result<- metaMDS(dist_COI, k = 2, trymax = 100)


#stress nmds 2 dim
nmds2_COI_result$stress

nmds2_COI_points <- as.data.frame(nmds2_COI_result$points)
nmds2_COI_points$Sample<- rownames(nmds2_COI_points)
nmds2_COI_points<- merge(nmds2_COI_points, metadata, by = "Sample")
nmds2_COI_points$Depth<- as.factor(nmds2_COI_points$Depth) 
nmds2_COI_points$TechRep<- as.factor(nmds2_COI_points$TechRep) 
nmds2_COI_points$Core_Rep <- as.factor(nmds2_COI_points$Core_Rep)

nmds2_COI_points$agegroup <- ifelse(
  nmds2_COI_points$Depth %in% c(4, 5) | nmds2_COI_points$Depth %in% c("4", "5"),
  "recent",
  ifelse(
    nmds2_COI_points$Depth %in% c(45, 50) | nmds2_COI_points$Depth %in% c("45", "50"),
    "old",
    NA
  )
)

nmds2_COI_points$agegroup <- factor(nmds2_COI_points$agegroup, levels = c("recent", "old"))

colors <- c("CCO" = "#8B3A3A", "CLR" = "cyan4", "STO" = "#CDAD00")

stress2_label_COI <- paste0("Stress = ", round(nmds2_COI_result$stress, 3))

shape_values <- c(21, 22, 24, 25, 23)
names(shape_values) <- levels(nmds2_COI_points$Core_Rep)

ggplot(nmds2_COI_points, aes(x = MDS1, y = MDS2)) +
  stat_ellipse(
    aes(fill = Site, group = Site),
    geom = "polygon",
    alpha = 0.2,
    color = NA,
    show.legend = FALSE
  ) +
  geom_point(
    data = subset(nmds2_COI_points, agegroup == "recent"),
    aes(fill = Site, shape = Core_Rep),
    size = 4,
    alpha = 0.9,
    color = "grey30",
    stroke = 0.8
  ) +
  geom_point(
    data = subset(nmds2_COI_points, agegroup == "old"),
    aes(fill = Site, shape = Core_Rep),
    size = 4,
    alpha = 0.9,
    color = "black",
    stroke = 1.2
  ) +
  scale_fill_manual(values = colors) +
  scale_shape_manual(values = shape_values) +
  theme_classic(base_size = 14) +
  labs(
    title = "NMDS by Sample (COI)",
    fill = "Location",
    shape = "CoreRep"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.title = element_text(size = 13, face = "bold"),
    legend.text = element_text(size = 10)
  ) +
  annotate(
    "text",
    x = 0, y = 4,
    label = stress2_label_COI,
    hjust = 1.1, vjust = -0.5,
    size = 4.5
  )


###AMB 18SSSS
arel_mitj18S<- read.xlsx("Clean_Data_18S_Reps_AbRel.xlsx")

arel_mitj18S1 <- column_to_rownames(arel_mitj18S, var = colnames(arel_mitj18S)[1])
arel_mitj18S1 <- arel_mitj18S1[,c(2:140)]

arel_mitj18S1_t<- as.data.frame(t(arel_mitj18S1))

arel_mitj18S1_t$Sample<- rownames(arel_mitj18S1_t)

arel_mitj18S1_t[, 1:(ncol(arel_mitj18S1_t)-1)] <- lapply(arel_mitj18S1_t[, 1:(ncol(arel_mitj18S1_t)-1)], as.numeric)

#Calcular NMDS
dist_18S <- vegdist(arel_mitj18S1_t[, -ncol(arel_mitj18S1_t)], method = "bray")


#k=2 dimensions
nmds2_18S_result<- metaMDS(dist_18S, k = 2, trymax = 100)


#stress nmds 2 dim
nmds2_18S_result$stress


nmds2_18S_points <- as.data.frame(nmds2_18S_result$points)
nmds2_18S_points$Sample<- rownames(nmds2_18S_points)
nmds2_18S_points<- merge(nmds2_18S_points, metadata, by = "Sample")
nmds2_18S_points$Depth<- as.factor(nmds2_18S_points$Depth) 
nmds2_18S_points$TechRep<- as.factor(nmds2_18S_points$TechRep) 
nmds2_18S_points$Core_Rep <- as.factor(nmds2_18S_points$Core_Rep)

nmds2_18S_points$agegroup <- ifelse(
  nmds2_18S_points$Depth %in% c(4, 5) | nmds2_18S_points$Depth %in% c("4", "5"),
  "recent",
  ifelse(
    nmds2_18S_points$Depth %in% c(45, 50) | nmds2_18S_points$Depth %in% c("45", "50"),
    "old",
    NA
  )
)

nmds2_18S_points$agegroup <- factor(nmds2_18S_points$agegroup, levels = c("recent", "old"))

colors <- c("CCO" = "#8B3A3A", "CLR" = "cyan4", "STO" = "#CDAD00")

stress2_label_18 <- paste0("Stress = ", round(nmds2_18S_result$stress, 3))


shape_values <- c(21, 22, 24, 25, 23)
names(shape_values) <- levels(nmds2_18S_points$Core_Rep)

ggplot(nmds2_18S_points, aes(x = MDS1, y = MDS2)) +
  stat_ellipse(
    aes(fill = Site, group = Site),
    geom = "polygon",
    alpha = 0.2,
    color = NA,
    show.legend = FALSE
  ) +
  geom_point(
    aes(fill = Site, shape = Core_Rep, color = agegroup),
    size = 6,
    alpha = 0.9,
    stroke = 1.1
  ) +
  scale_fill_manual(values = colors) +
  scale_shape_manual(values = shape_values) +
  scale_color_manual(
    values = c("recent" = "grey40", "old" = "black"),
    name = "Age group"
  ) +
  theme_classic(base_size = 14) +
  labs(
    title = "NMDS by Sample (18S)",
    fill = "Location",
    shape = "CoreRep"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.title = element_text(size = 13, face = "bold"),
    legend.text = element_text(size = 10)
  ) +
  annotate(
    "text",
    x = -2, y = 4,
    label = stress2_label_18,
    hjust = 1.1, vjust = -0.5,
    size = 4.5
  )




### PERMANOVA Stadistics

run_dispersion_tests <- function(dist_obj, meta, factors, nperm = 999) {
  out <- list()
  for (f in factors) {
    g <- meta[[f]]
    g <- droplevels(as.factor(g))
    
    tab <- table(g)
    if (length(tab) < 2 || any(tab < 2)) {
      message("Skipping betadisper for ", f, " (needs >=2 groups and >=2 obs per group).")
      next
    }
    
    bd <- betadisper(dist_obj, group = g, type = "centroid")
    pt <- permutest(bd, permutations = nperm)
    
    out[[f]] <- list(betadisper = bd, permutest = pt)
  }
  out
}

# -------------------------
# 1) PERMANOVA (location)
# -------------------------

samples_in_dist <- rownames(as.matrix(dist_COI))

metadata_COI <- metadata[match(samples_in_dist, metadata$Sample), ]

stopifnot(nrow(metadata_COI) == attr(dist_COI, "Size"))
stopifnot(all(metadata_COI$Sample == samples_in_dist))

metadata_COI$SampleID  <- interaction(metadata_COI$Site, metadata_COI$Core_Rep, metadata_COI$Depth)
metadata_COI$SiteDepth <- interaction(metadata_COI$Site, metadata_COI$Depth)
metadata_COI$CoreID    <- interaction(metadata_COI$Site, metadata_COI$Core_Rep)

metadata$SampleID <- interaction(metadata$Site, metadata$Core_Rep, metadata$Depth)

disp_tech <- betadisper(dist_COI, metadata_COI$SampleID)
disp_bio  <- betadisper(dist_COI, metadata_COI$SiteDepth)

df_tech <- data.frame(
  value = disp_tech$distances,
  Level = "Technical"
)

comm_COI_df <- as.data.frame(arel_mitjCOI1_t)
comm_COI_df$Sample <- rownames(comm_COI_df)

comm_COI_avg <- comm_COI_df |>
  left_join(metadata_COI[, c("Sample", "SampleID", "Site", "Core_Rep", "Depth")], by = "Sample") |>
  group_by(SampleID, Site, Core_Rep, Depth) |>
  summarise(across(where(is.numeric), mean), .groups = "drop")

comm_COI_avg_mat <- as.matrix(comm_COI_avg[, !(names(comm_COI_avg) %in% c("SampleID", "Site", "Core_Rep", "Depth"))])
rownames(comm_COI_avg_mat) <- comm_COI_avg$SampleID

metadata_COI_avg <- comm_COI_avg[, c("SampleID", "Site", "Core_Rep", "Depth")]
metadata_COI_avg$SiteDepth <- interaction(metadata_COI_avg$Site, metadata_COI_avg$Depth)
metadata_COI_avg$CoreID <- interaction(metadata_COI_avg$Site, metadata_COI_avg$Core_Rep)


dist_COI_avg <- vegdist(comm_COI_avg_mat, method = "bray")
dist_mat_avg <- as.matrix(dist_COI_avg)

bio_vals <- list()

for (sd in unique(metadata_COI_avg$SiteDepth)) {
  idx <- which(metadata_COI_avg$SiteDepth == sd)
  if (length(idx) > 1) {
    vals <- dist_mat_avg[idx, idx][upper.tri(dist_mat_avg[idx, idx])]
    bio_vals[[sd]] <- vals
  }
}

df_bio <- data.frame(
  value = unlist(bio_vals),
  Level = "Biological"
)

df_compare <- rbind(df_tech, df_bio)

medians <- aggregate(value ~ Level, data = df_compare, median)

ggplot(df_compare, aes(x = Level, y = value)) +
  geom_boxplot(outlier.shape = NA, width = 0.5, size = 0.8, alpha = 0.6) +
  geom_jitter(width = 0.12, alpha = 0.6, size = 4) +
  
  geom_point(data = medians, aes(x = Level, y = value),
             color = "black", size = 4, shape = 18) +
  
  geom_text(data = medians,
            aes(x = Level, y = value, label = round(value, 3)),
            vjust = -1.2, hjust= 4, size = 5, fontface = "bold") +
  
  theme_classic(base_size = 16) +
  
  labs(
    y = "Bray–Curtis dissimilarity",
    x = ""
  ) +
  
  theme(
    legend.position = "none",
    axis.title.y = element_text(size = 18, face = "bold"),
    axis.text.x  = element_text(size = 16, face = "bold"),
    axis.text.y  = element_text(size = 14)
  )


wilcox.test(value ~ Level, data = df_compare)
aggregate(value ~ Level, data = df_compare, median)




#### 18S

samples_in_dist <- rownames(as.matrix(dist_18S))
metadata<- read.xlsx("metadata_repliques_f.xlsx")

metadata_COI <- metadata[match(samples_in_dist, metadata$Sample), ]

stopifnot(nrow(metadata_COI) == attr(dist_18S, "Size"))
stopifnot(all(metadata_COI$Sample == samples_in_dist))

metadata_COI$SampleID  <- interaction(metadata_COI$Site, metadata_COI$Core_Rep, metadata_COI$Depth)
metadata_COI$SiteDepth <- interaction(metadata_COI$Site, metadata_COI$Depth)
metadata_COI$CoreID    <- interaction(metadata_COI$Site, metadata_COI$Core_Rep)

metadata$SampleID <- interaction(metadata$Site, metadata$Core_Rep, metadata$Depth)

disp_tech_18 <- betadisper(dist_18S, metadata_COI$SampleID)
disp_bio_18  <- betadisper(dist_18S, metadata_COI$SiteDepth)

df_tech_18 <- data.frame(
  value = disp_tech_18$distances,
  Level = "Technical"
)

comm_COI_df <- as.data.frame(arel_mitj18S1_t)
comm_COI_df$Sample <- rownames(comm_COI_df)

comm_COI_avg <- comm_COI_df |>
  left_join(metadata_COI[, c("Sample", "SampleID", "Site", "Core_Rep", "Depth")], by = "Sample") |>
  group_by(SampleID, Site, Core_Rep, Depth) |>
  summarise(across(where(is.numeric), mean), .groups = "drop")

comm_COI_avg_mat <- as.matrix(comm_COI_avg[, !(names(comm_COI_avg) %in% c("SampleID", "Site", "Core_Rep", "Depth"))])
rownames(comm_COI_avg_mat) <- comm_COI_avg$SampleID

metadata_COI_avg <- comm_COI_avg[, c("SampleID", "Site", "Core_Rep", "Depth")]
metadata_COI_avg$SiteDepth <- interaction(metadata_COI_avg$Site, metadata_COI_avg$Depth)
metadata_COI_avg$CoreID <- interaction(metadata_COI_avg$Site, metadata_COI_avg$Core_Rep)


dist_COI_avg <- vegdist(comm_COI_avg_mat, method = "bray")
dist_mat_avg <- as.matrix(dist_COI_avg)

bio_vals <- list()

for (sd in unique(metadata_COI_avg$SiteDepth)) {
  idx <- which(metadata_COI_avg$SiteDepth == sd)
  if (length(idx) > 1) {
    vals <- dist_mat_avg[idx, idx][upper.tri(dist_mat_avg[idx, idx])]
    bio_vals[[sd]] <- vals
  }
}

df_bio_18 <- data.frame(
  value = unlist(bio_vals),
  Level = "Biological"
)

df_compare_18 <- rbind(df_tech_18, df_bio_18)

medians <- aggregate(value ~ Level, data = df_compare_18, median)

ggplot(df_compare_18, aes(x = Level, y = value)) +
  geom_boxplot(outlier.shape = NA, width = 0.5, size = 0.8, alpha = 0.6) +
  geom_jitter(width = 0.12, alpha = 0.6, size = 4) +
  
  geom_point(data = medians, aes(x = Level, y = value),
             color = "black", size = 4, shape = 18) +
  
  geom_text(data = medians,
            aes(x = Level, y = value, label = round(value, 3)),
            vjust = -1.2, hjust= 4, size = 5, fontface = "bold") +
  
  theme_classic(base_size = 16) +

  labs(
    y = "Bray–Curtis dissimilarity",
    x = ""
  ) +
  
  theme(
    legend.position = "none",
    axis.title.y = element_text(size = 18, face = "bold"),
    axis.text.x  = element_text(size = 16, face = "bold"),
    axis.text.y  = element_text(size = 14)
  )


wilcox.test(value ~ Level, data = df_compare_18)
aggregate(value ~ Level, data = df_compare_18, median)

