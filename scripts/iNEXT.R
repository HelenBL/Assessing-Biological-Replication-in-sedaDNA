install.packages("iNEXT")

library(iNEXT)
library(ggplot2)

################ COI

tab <- read.xlsx("Clean_Data_COI_Reps_AbTot.xlsx")
metadata <- read.xlsx("metadata_repliques_f.xlsx")

rownames(tab) <- tab$X1
tab$X1 <- NULL
tabt <- t(tab)

inext_input <- apply(tabt, 1, function(x) x[x > 0])

out <- iNEXT(
  inext_input,
  q = 0,
  datatype = "abundance",
  nboot = 50
)

df_plot <- out$iNextEst$coverage_based %>%
  left_join(metadata, by = c("Assemblage" = "Sample"))


ggplot(df_plot, aes(x = SC, y = qD, color = Site, group = Assemblage)) +
  geom_line(alpha = 0.3) +
  theme_bw()

df_summary <- df_plot %>%
  group_by(Site, SC) %>%
  summarise(qD = mean(qD), .groups = "drop")

ggplot(df_summary, aes(x = SC, y = qD, color = Site)) +
  geom_line(size = 1.2) +
  theme_bw()

ggiNEXT(out, type = 1) +
  theme_bw()

ggiNEXT(out, type = 3) +
  theme_bw()




#################### 18S
tab18 <- read.xlsx("Clean_Data_18S_Reps_AbTot.xlsx")
rownames(tab18) <- tab18$X1
tab18$X1 <- NULL
tab18t <- t(tab18)

inext_input18 <- apply(tab18t, 1, function(x) x[x > 0])

out18 <- iNEXT(
  inext_input18,
  q = 0,
  datatype = "abundance",
  nboot = 50
)

ggiNEXT(out18, type = 1) +
  theme_bw()

ggiNEXT(out18, type = 3) +
  theme_bw()