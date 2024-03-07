#setwd('documents/1研究/k_distance_limited_MFVS/ASD_oncogene')
library(clusterProfiler)
library(org.Hs.eg.db)

#df <- read.csv('')
gene_list <- c("AHR", "AR", "ATF4", "BCL6", "BRCA1", "CIITA", "DDIT3", "DNMT1", "E2F1", "E2F4",
                 "EP300", "ERG", "ESR1", "FLI1", "FOS", "GATA1", "GATA3", "GLI1", "GLI2", "HDAC7",
                 "HIF1A", "HNF4A", "HOXD9", "JUN", "KLF4", "MEF2A", "MYB", "MYBL2", "MYC", "MYCN",
                 "NF1", "NFKB1", "NR0B1", "NR1I2", "NR3C1", "NR5A2", "PML", "POU1F1", "POU5F1",
                 "RARB", "RB1", "RELA", "SALL4", "SCD5", "SIRT1", "SNAI1", "SOX9", "SP1", "STAT1",
                 "STAT3", "TFAP2A", "TP53", "TWIST1", "USF2", "VDR", "WT1")

gene_list_MFVS <- c('AR', 'ARNT', 'ATF4', 'BCL6', 'BRCA1', 'CIITA', 'DDIT3', 'DNMT1', 'E2F1', 'E2F4', 'EGR2', 'EP300', 'ERG', 'ESR1', 'FLI1', 'FOS', 'GATA1', 'GATA3', 'GLI1', 'GLI2', 'HDAC7', 'HIF1A', 'HNF4A', 'HOXD9', 'JUN', 'KLF4', 'MEF2A', 'MYB', 'MYBL2', 'MYC', 'MYCN', 'NFKB1', 'NR0B1', 'NR3C1', 'NR5A2', 'PML', 'POU1F1', 'POU2AF1', 'POU5F1', 'RARB', 'RB1', 'RELA', 'SALL4', 'SIRT1', 'SNAI1', 'SOX9', 'SP1', 'SP3', 'STAT1', 'STAT3', 'TFAP2A', 'TP53', 'TWIST1', 'USF2', 'VDR', 'WT1')
gene_list_k2MFVS <- c('AR', 'ARNT', 'ARNTL', 'ASCL1', 'ATF1', 'ATF4', 'ATM', 'BCL6', 'BRCA1', 'CIITA', 'CREB1', 'CREBBP', 'CREM', 'CTCF', 'CTNNB1', 'DDIT3', 'DNMT1', 'DR1', 'E2F1', 'E2F3', 'E2F4', 'EGR1', 'EGR3', 'ELK1', 'ELL', 'EP300', 'EPAS1', 'ERG', 'ESR1', 'ESR2', 'ETS1', 'FHL2', 'FLI1', 'FOS', 'FOSL1', 'FOXA2', 'FOXM1', 'FOXO1', 'FOXO3', 'GATA1', 'GATA3', 'GLI1', 'GLI2', 'HDAC1', 'HDAC3', 'HDAC9', 'HEY1', 'HIF1A', 'HINFP', 'HNF4A', 'HOXA10', 'HOXA9', 'HOXD9', 'HSF1', 'ID1', 'ID3', 'IRF1', 'IRF4', 'IRF7', 'ISL1', 'JUN', 'KAT2B', 'KLF4', 'LEF1', 'LMO2', 'MDM2', 'MEF2A', 'MEF2C', 'MEF2D', 'MITF', 'MSX1', 'MTA1', 'MYB', 'MYBL2', 'MYC', 'MYCN', 'NANOG', 'NF1', 'NFKB1', 'NFKB2', 'NKX2-1', 'NR0B1', 'NR1H4', 'NR1I2', 'NR2F2', 'NR3C1', 'NR5A1', 'NR5A2', 'PAX6', 'PGR', 'PML', 'POU1F1', 'POU2F1', 'POU5F1', 'PPARG', 'RAD51', 'RARA', 'RARB', 'RB1', 'RBL1', 'RELA', 'REST', 'RUNX2', 'RUNX3', 'SALL4', 'SATB1', 'SCD5', 'SIRT1', 'SMAD3', 'SMAD4', 'SNAI1', 'SNAI2', 'SOX2', 'SOX4', 'SOX9', 'SP1', 'SP3', 'SPI1', 'SRY', 'STAT1', 'STAT3', 'STAT6', 'TBP', 'TBX21', 'TCF4', 'TFAP2A', 'TP53', 'TP63', 'TWIST1', 'VDR', 'WT1', 'XBP1', 'YY1')
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

png("GO_enrichment_plot_in_TRRUSTv2_MF.png",width=800,height=600,units="px",res=100)
plot_obj <-barplot(ego_result,x='GeneRatio',color='qvalue',showCategory=20,font.size=13,label_format=90,title='TRRUST v2; GO; MF; \nGround truth')
plot_obj$width <- 0.5
plot(plot_obj)
dev.off()

png("GO_enrichment_plot_in_TRRUSTv2_MFVS_MF.png",width=800,height=600,units="px",res=100)
plot(barplot(ego_result_MFVS,x='GeneRatio',color='qvalue',showCategory=20,font.size=13,label_format=90,title='TRRUST v2; GO; MF; \nMFVS'))
dev.off()

png("GO_enrichment_plot_in_TRRUSTv2_k2MFVS_MF.png",width=800,height=600,units="px",res=100)
plot(barplot(ego_result_k2MFVS,x='GeneRatio',color='qvalue',showCategory=20,font.size=13,label_format=90,title='TRRUST v2; GO; MF; \nk2MFVS'))
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

png("GO_enrichment_plot_in_TRRUSTv2_CC.png",width=800,height=600,units="px",res=100)
plot(barplot(ego_result,x='GeneRatio',color='qvalue',showCategory=20,font.size=13,label_format=90,title='TRRUST v2; GO; CC; \nGround truth'))
dev.off()

png("GO_enrichment_plot_in_TRRUSTv2_MFVS_CC.png",width=800,height=600,units="px",res=100)
plot(barplot(ego_result_MFVS,x='GeneRatio',color='qvalue',showCategory=20,font.size=13,label_format=90,title='TRRUST v2; GO; CC; \nMFVS'))
dev.off()

png("GO_enrichment_plot_in_TRRUSTv2_k2MFVS_CC.png",width=800,height=600,units="px",res=100)
plot(barplot(ego_result_k2MFVS,x='GeneRatio',color='qvalue',showCategory=20,font.size=13,label_format=90,title='TRRUST v2; GO; CC; \nk2MFVS'))
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

png("GO_enrichment_plot_in_TRRUSTv2_BP.png",width=800,height=600,units="px",res=100)
plot(barplot(ego_result,x='GeneRatio',color='qvalue',showCategory=20,font.size=13,label_format=90,title='TRRUST v2; GO; BP; \nGround truth'))
dev.off()

png("GO_enrichment_plot_in_TRRUSTv2_MFVS_BP.png",width=800,height=600,units="px",res=100)
plot(barplot(ego_result_MFVS,x='GeneRatio',color='qvalue',showCategory=20,font.size=13,label_format=90,title='TRRUST v2; GO; BP; \nMFVS'))
dev.off()

png("GO_enrichment_plot_in_TRRUSTv2_k2MFVS_BP.png",width=800,height=600,units="px",res=100)
plot(barplot(ego_result_k2MFVS,x='GeneRatio',color='qvalue',showCategory=20,font.size=13,label_format=90,title='TRRUST v2; GO; BP; \nk2MFVS'))
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

png("KEGG_enrichment_plot_in_TRRUSTv2.png", width = 800, height = 600, units = "px", res = 100)
plot(barplot(kegg_result,x='GeneRatio',color='qvalue',showCategory=20,font.size=13,label_format=90,title='TRRUST v2; KEGG; \nGround truth'))
dev.off()

png("KEGG_enrichment_plot_in_TRRUSTv2_MFVS.png", width = 800, height = 600, units = "px", res = 100)
plot(barplot(kegg_result_MFVS,x='GeneRatio',color='qvalue',showCategory=20,font.size=13,label_format=90,title='TRRUST v2; KEGG; \nMFVS'))
dev.off()

png("KEGG_enrichment_plot_in_TRRUSTv2_k2MFVS.png", width = 800, height = 600, units = "px", res = 100)
plot(barplot(kegg_result_k2MFVS,x='GeneRatio',color='qvalue',showCategory=20,font.size=13,label_format=90,title='TRRUST v2; KEGG; \nk2MFVS'))
dev.off()
