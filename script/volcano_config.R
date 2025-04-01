# conda env:
# conda activate R-EnhancedVolcano_ce
# DESeq2 table
path2DESeq2 <- "/home/joe.ueda/Bureau/cDiff/maps_deseq2/rcd22/fichier_envoi/rcd22_short_MAPS.txt"

lg2Seuil <- 2                                                                                                                                                          
padjSeuil <- 0.05

titleVal <- 'MAPS RCd22 short(a&b) wt:a&b'
subtitleVal <- 'RCd22vsWt.complete.txt, WT (r1,r2), RCd22 (short)'


xVal <- 'log2FoldChange'
x_labelVal <- bquote(~Log[2]~ 'fold change')
yVal <- 'padj'
y_labelVal <- bquote(~Log[10]~ 'padj')

legendPositionVal <- 'bottom'

captionVal <- paste("Log2FC cutoff:",lg2Seuil,", padj cutoff:",padjSeuil) 

rcd22type <- "short" # choose between "short" or "long"
