library(pheatmap)

# if not installed
install.packages("pheatmap")



pheatmap(t(expr), show_colnames = F, 
         color = hcl.colors(100, palette = "Roma"),
         scale = "row", border_color = NA,
         annotation_col = meta.[,1:2], cluster_cols = F, 
         gaps_col = 41, gaps_row = 6, cluster_rows = F,
         main = "Heatmap", display_numbers = F, labels_col = c(1:60))

