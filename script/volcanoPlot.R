library(EnhancedVolcano)


args <- commandArgs(trailingOnly = TRUE)

config_file <- args[1]

# Load volcano config
source(config_file, local=TRUE)

deseqTable <- args[2]
print(paste0("Imported from : ", deseqTable))

# Read data
res <- read.table(deseqTable, header=TRUE, sep ="\t")


# Set rownames based on Id or old_locus_tag
rownames(res) <- ifelse(grepl("^CDR20291_RS", res$Id),
                        ifelse(!is.na(res$old_locus_tag),
                               gsub("CDR20291_", "", res$old_locus_tag),
                               gsub("CDR20291_", "", res$Id)),
                        gsub("CDR20291_", "", res$Id))

# Define custom colors for points
keyvals <- ifelse(res$padj < padjSeuil & res$log2FoldChange > lg2Seuil, 'red',
                  ifelse(res$padj < padjSeuil & res$log2FoldChange < -lg2Seuil, 'blue', 'grey'))
keyvals[is.na(keyvals)] <- 'grey'
names(keyvals)[keyvals == 'red'] <- 'up'
names(keyvals)[keyvals == 'grey'] <- 'NS'
names(keyvals)[keyvals == 'blue'] <- 'down'

# Define labels based on rcd22type
if (rcd22type == "long") {
    spec_label <- c("1558", "2768", "1829", "3357", "0538", "1462", "RCd22", "RS16617", "1951", "RCd2")
} else if (rcd22type == "short") {
    spec_label <- c("1558", "2768", "1829", "3357", "0538", "1462", "RCd22", "RS16617", "RCd21", "RS10482", "RS19465")
}

# Generate personalized filename with timestamp
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
file_name <- paste0("volcano_plot.pdf")

# Save the plot as PNG
pdf(file_name, width=10, height=6)

#png(filename = file_name, width = 1000, height = 600)

# Plot volcano
EnhancedVolcano(
  res,
  title = titleVal,
  subtitle = subtitleVal,
  lab = rownames(res),
  x = xVal,
  xlab = x_labelVal,
  y = yVal,
  ylab = y_labelVal,
  colCustom = keyvals,
  legendPosition = legendPositionVal,
  pCutoff = padjSeuil,
  FCcutoff = lg2Seuil,
  caption = captionVal,
  labSize = 3,
  selectLab = spec_label
) + coord_flip()

# Close PNG device
dev.off()

# Notify user of saved file
print(paste0("Volcano plot saved as: ", file_name))
