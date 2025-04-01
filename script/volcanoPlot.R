library(EnhancedVolcano)


source("/home/joe.ueda/Bureau/cDiff/maps_deseq2/rcd22/results/volcano_config.R", local=TRUE)
print(paste0("Imported from : ",path2DESeq2))


res <- read.table(path2DESeq2, header=TRUE, sep ="\t")

rownames(res) <- ifelse(grepl("^CDR20291_RS", res$Id), # true if value in column "id" starts with "CDR20291_RS", 
                        ifelse(!is.na(res$old_locus_tag), gsub("CDR20291_","",res$old_locus_tag), gsub("CDR20291_","",res$Id)), # check if old_locus_tag is not NA
                        gsub("CDR20291_","",res$Id))

keyvals <- keyvals <- ifelse(res$padj < padjSeuil & res$log2FoldChange > lg2Seuil, 'red', ifelse(res$padj < padjSeuil & res$log2FoldChange < -lg2Seuil, 'blue', 'grey'))
keyvals[is.na(keyvals)] <- 'grey'
names(keyvals)[keyvals == 'red'] <- 'up'
names(keyvals)[keyvals == 'grey'] <- 'NS'
names(keyvals)[keyvals == 'blue'] <- 'down'

dev.new(width=10, height=6)

# label for RCd22 short
spec_label <- c("1558", "2768", "1829", "3357", "0538", "1462", "RCd22","RS16617","RCd21","RS10482", "RS19465")
# label for RCd22 long
spec_label2 <- c("1558", "2768", "1829", "3357", "0538", "1462", "RCd22","RS16617","1951", "RCd2")

if(rcd22type=="long"){
    # label for RCd22 long
    spec_label <- c("1558", "2768", "1829", "3357", "0538", "1462", "RCd22","RS16617","1951", "RCd2")
    
}else if (rcd22type=="short"){
    # label for RCd22 short
    spec_label <- c("1558", "2768", "1829", "3357", "0538", "1462", "RCd22","RS16617","RCd21","RS10482", "RS19465")
}

EnhancedVolcano(res, title = titleVal, subtitle = subtitleVal, lab = rownames(res), x = xVal, xlab = x_labelVal, y = yVal, ylab = y_labelVal, colCustom = keyvals, legendPosition = legendPositionVal, pCutoff = padjSeuil, FCcutoff = lg2Seuil, caption = captionVal, labSize = 3, selectLab = spec_label) + coord_flip()

