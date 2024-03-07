#setwd('documents/1研究/k_distance_limited_MFVS/ASD_oncogene')
library(clusterProfiler)
library(org.Hs.eg.db)

gene_list <- c("AHR", "AR", "ATF4", "BCL6", "BRCA1", "CIITA", "DDIT3", "DNMT1", "E2F1", "E2F4",
                 "EP300", "ERG", "ESR1", "FLI1", "FOS", "GATA1", "GATA3", "GLI1", "GLI2", "HDAC7",
                 "HIF1A", "HNF4A", "HOXD9", "JUN", "KLF4", "MEF2A", "MYB", "MYBL2", "MYC", "MYCN",
                 "NF1", "NFKB1", "NR0B1", "NR1I2", "NR3C1", "NR5A2", "PML", "POU1F1", "POU5F1",
                 "RARB", "RB1", "RELA", "SALL4", "SCD5", "SIRT1", "SNAI1", "SOX9", "SP1", "STAT1",
                 "STAT3", "TFAP2A", "TP53", "TWIST1", "USF2", "VDR", "WT1")

gene_list_MFVS <- c('ARNTL', 'BBC3', 'GATA3', 'GLI1', 'GLI2', 'GLI3', 'RARA', 'SOX2')
gene_list_k2MFVS <- c('ARNTL', 'CEBPA', 'CTNNB1', 'ELK1', 'FOS', 'FOSB', 'GATA3', 'GLI1', 'GLI2', 'GLI3', 'HIF1A', 'IFNG', 'JUN', 'KLF2', 'KLF4', 'MYC', 'NANOG', 'NFATC1', 'NFKB2', 'NFX1', 'NR1H3', 'ONECUT1', 'PDX1', 'PGR', 'PML', 'PPARD', 'PSEN1', 'RARA', 'TBX21', 'THRA', 'TNF', 'TP53')
print(length(gene_list_k2MFVS))

entrez_ids <- bitr(gene_list, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
entrez_ids_MFVS <- bitr(gene_list_MFVS, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
entrez_ids_k2MFVS <- bitr(gene_list_k2MFVS, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
#print(entrez_ids$ENTREZID)

ego_result <- enrichGO(gene=gene_list,
                OrgDb=org.Hs.eg.db,
                keyType='SYMBOL',
                ont='MF',
                pAdjustMethod="BH",
                qvalueCutoff=0.05,
                readable=TRUE)

ego_result_MFVS <- enrichGO(gene=gene_list_MFVS,
                OrgDb=org.Hs.eg.db,
                keyType='SYMBOL',
                ont='MF',
                pAdjustMethod="BH",
                qvalueCutoff=0.05,
                readable=TRUE)
ego_result_k2MFVS <- enrichGO(gene=gene_list_k2MFVS,
                OrgDb=org.Hs.eg.db,
                keyType='SYMBOL',
                ont='MF',
                pAdjustMethod="BH",
                qvalueCutoff=0.05,
                readable=TRUE)

png("GO_enrichment_plot_in_RegNetwork_MF.png",width=800,height=600,units="px",res=100)
plot(barplot(ego_result,x='GeneRatio',color='qvalue',showCategory=20,font.size=13,label_format=90,title='RegNetwork; GO; MF; \nGround truth'))
dev.off()

png("GO_enrichment_plot_in_RegNetwork_MFVS_MF.png",width=800,height=600,units="px",res=100)
plot(barplot(ego_result_MFVS,x='GeneRatio',color='qvalue',showCategory=20,font.size=13,label_format=90,title='RegNetwork; GO; MF; \nMFVS'))
dev.off()

png("GO_enrichment_plot_in_RegNetwork_k2MFVS_MF.png",width=800,height=600,units="px",res=100)
plot(barplot(ego_result_k2MFVS,x='GeneRatio',color='qvalue',showCategory=20,font.size=13,label_format=90,title='RegNetwork; GO; MF; \nk2MFVS'))
dev.off()

ego_result <- enrichGO(gene=gene_list,
                OrgDb=org.Hs.eg.db,
                keyType='SYMBOL',
                ont='CC',
                pAdjustMethod="BH",
                qvalueCutoff=0.05,
                readable=TRUE)

ego_result_MFVS <- enrichGO(gene=gene_list_MFVS,
                OrgDb=org.Hs.eg.db,
                keyType='SYMBOL',
                ont='CC',
                pAdjustMethod="BH",
                qvalueCutoff=0.05,
                readable=TRUE)
ego_result_k2MFVS <- enrichGO(gene=gene_list_k2MFVS,
                OrgDb=org.Hs.eg.db,
                keyType='SYMBOL',
                ont='CC',
                pAdjustMethod="BH",
                qvalueCutoff=0.05,
                readable=TRUE)

png("GO_enrichment_plot_in_RegNetwork_CC.png",width=800,height=600,units="px",res=100)
plot(barplot(ego_result,x='GeneRatio',color='qvalue',showCategory=20,font.size=13,label_format=90,title='RegNetwork; GO; CC; \nGround truth'))
dev.off()

png("GO_enrichment_plot_in_RegNetwork_MFVS_CC.png",width=800,height=600,units="px",res=100)
plot(barplot(ego_result_MFVS,x='GeneRatio',color='qvalue',showCategory=20,font.size=13,label_format=90,title='RegNetwork; GO; CC; \nMFVS'))
dev.off()

png("GO_enrichment_plot_in_RegNetwork_k2MFVS_CC.png",width=800,height=600,units="px",res=100)
plot(barplot(ego_result_k2MFVS,x='GeneRatio',color='qvalue',showCategory=20,font.size=13,label_format=90,title='RegNetwork; GO; CC; \nk2MFVS'))
dev.off()

ego_result <- enrichGO(gene=gene_list,
                OrgDb=org.Hs.eg.db,
                keyType='SYMBOL',
                ont='BP',
                pAdjustMethod="BH",
                qvalueCutoff=0.05,
                readable=TRUE)

ego_result_MFVS <- enrichGO(gene=gene_list_MFVS,
                OrgDb=org.Hs.eg.db,
                keyType='SYMBOL',
                ont='BP',
                pAdjustMethod="BH",
                qvalueCutoff=0.05,
                readable=TRUE)
ego_result_k2MFVS <- enrichGO(gene=gene_list_k2MFVS,
                OrgDb=org.Hs.eg.db,
                keyType='SYMBOL',
                ont='BP',
                pAdjustMethod="BH",
                qvalueCutoff=0.05,
                readable=TRUE)

png("GO_enrichment_plot_in_RegNetwork_BP.png",width=800,height=600,units="px",res=100)
plot(barplot(ego_result,x='GeneRatio',color='qvalue',showCategory=20,font.size=13,label_format=90,title='RegNetwork; GO; MP; \nGround truth'))
dev.off()

png("GO_enrichment_plot_in_RegNetwork_MFVS_BP.png",width=800,height=600,units="px",res=100)
plot(barplot(ego_result_MFVS,x='GeneRatio',color='qvalue',showCategory=20,font.size=13,label_format=90,title='RegNetwork; GO; MP; \nMFVS'))
dev.off()

png("GO_enrichment_plot_in_RegNetwork_k2MFVS_BP.png",width=800,height=600,units="px",res=100)
plot(barplot(ego_result_k2MFVS,x='GeneRatio',color='qvalue',showCategory=20,font.size=13,label_format=90,title='RegNetwork; GO; MP; \nk2MFVS'))
dev.off()

kegg_result <- enrichKEGG(gene = entrez_ids$ENTREZID,
                          organism='hsa',
                          keyType='kegg',
                          pvalueCutoff=0.05,
                          pAdjustMethod='BH',
                          qvalueCutoff=0.1)

kegg_result_MFVS <- enrichKEGG(gene = entrez_ids_MFVS$ENTREZID,
                          organism='hsa',
                          keyType='kegg',
                          pvalueCutoff = 0.05,
                          pAdjustMethod= 'BH',
                          qvalueCutoff = 0.1)
kegg_result_k2MFVS <- enrichKEGG(gene = entrez_ids_k2MFVS$ENTREZID,
                          organism='hsa',
                          keyType='kegg',
                          pvalueCutoff = 0.05,
                          pAdjustMethod= 'BH',
                          qvalueCutoff = 0.1)

png("KEGG_enrichment_plot_in_RegNetwork.png", width = 800, height = 600, units = "px", res = 100)
plot(barplot(kegg_result,x='GeneRatio',color='qvalue',showCategory=20,font.size=13,label_format=90,title='RegNetwork; KEGG; \nGround truth'))
dev.off()

png("KEGG_enrichment_plot_in_RegNetwork_MFVS.png", width = 800, height = 600, units = "px", res = 100)
plot(barplot(kegg_result_MFVS,x='GeneRatio',color='qvalue',showCategory=20,font.size=13,label_format=90,title='RegNetwork; KEGG; \nMFVS'))
dev.off()

png("KEGG_enrichment_plot_in_RegNetwork_k2MFVS.png", width = 800, height = 600, units = "px", res = 100)
plot(barplot(kegg_result_k2MFVS,x='GeneRatio',color='qvalue',showCategory=20,font.size=13,label_format=90,title='RegNetwork; KEGG; \nk2MFVS'))
dev.off()
