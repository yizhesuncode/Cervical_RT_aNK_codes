
#Read the GSE297041 data and visualize the UMAP and 
cervical_RT_tumor=readRDS("~/Wanmeng/cervical_scRNA_RT/GSE297041_CESC_18_scRNA_rmdoublet_0.2_cluster_ident_without_anchor.rds")
DimPlot(cervical_RT_tumor,reduction = "umap",label = T)

#Annotate the cells using SingleR
library(SingleR)
library(celldex)
library(ExperimentHub)
library(ggsci)
cors=pal_igv()(35)
cors=c(cors,"purple","orange")

sc_data=as.matrix(cervical_RT_tumor@assays$RNA@counts)
immune.se=DatabaseImmuneCellExpressionData()
hpca.se=HumanPrimaryCellAtlasData()
monaimmune.se=MonacoImmuneData()

hpca.se$label.main
pred.hesc = SingleR(test = sc_data, ref = hpca.se, assay.type.test=1,
                    labels = hpca.se$label.main)
pred.hesc$labels
table(pred.hesc$labels,cervical_RT_tumor$seurat_clusters)
cervical_RT_tumor$cell_annotation_singleR=pred.hesc$labels
DimPlot(cervical_RT_tumor,group.by = "cell_annotation_singleR",reduction = "umap",label = T,cols=cors)

tab <- table(pred.hesc$labels, cervical_RT_tumor$seurat_clusters)
cluster_annotation <- apply(tab, 2, function(x) {
  names(which.max(x))
})

cluster_annotation
Idents(cervical_RT_tumor)="seurat_clusters"

cluster_annot_per_cell <- cluster_annotation[
  as.character(cervical_RT_tumor$seurat_clusters)
]

cervical_RT_tumor =RenameIdents(cervical_RT_tumor,
                                "0"="T_cells",
                                "1"="Neutrophils",
                                "2"="T_cells",
                                "3"="Monocyte",
                                "4"="Epithelial/Cancer_cells",
                                "5"="Macrophage/Monocyte",
                                "6"="Epithelial/Cancer_cells",
                                "7"="Epithelial/Cancer_cells",
                                "8"="Epithelial/Cancer_cells",
                                "9"="Epithelial/Cancer_cells",
                                "10"="Neutrophils",
                                "11"="Epithelial/Cancer_cells",
                                "12"="Fibroblasts/Tissue_stem_cells",
                                "13"="NK_cells",
                                "14"="Epithelial/Cancer_cells",
                                "15"="Endothelial_cells",
                                "16"="B_cells",
                                "17"="Tissue_stem_cells",
                                "18"="Epithelial/Cancer_cells",
                                "19"="Epithelial/Cancer_cells",
                                "20"="Epithelial/Cancer_cells",
                                "21"="Epithelial/Cancer_cells",
                                "22"="B_cells",
                                "23"="Not_defined",
                                "24"="CMP"
)


cervical_RT_tumor$new_idents=Idents(cervical_RT_tumor)
cervical_RT_tumor@meta.data[1:5,]

Idents(cervical_RT_tumor)="new_idents"
DimPlot(cervical_RT_tumor,reduction = "umap",label = T,cols=cors)

base::saveRDS(cervical_RT_tumor,"~/Wanmeng/cervical_scRNA_RT/total_data_processed.rds")
cervical_RT_tumor=readRDS("~/Wanmeng/cervical_scRNA_RT/total_data_processed.rds")

#Visualize the UMAP with NK and other cells
cervical_RT_tumor@meta.data <- cervical_RT_tumor@meta.data %>%
  mutate(idents_for_NK = case_when(
    new_idents == "NK_cells" ~ "NK_cells",
    TRUE ~ "Others"
  ))

DimPlot(cervical_RT_tumor, group.by = "idents_for_NK", reduction = "umap", label = TRUE)
DimPlot(
  cervical_RT_tumor,
  group.by = "idents_for_NK",
  reduction = "umap",
  label = TRUE,
  repel = TRUE,
  label.size = 10,
  pt.size = 1.2,
  cols = c("#fc8d62", "#a6cee3")  
) +
  theme_classic(base_size = 15) +
  theme(
    axis.title = element_text(size = 25),
    axis.text = element_text(size = 25),
    legend.text = element_text(size = 25),
    legend.position = "right",
    legend.title = element_text(size = 25),
    plot.title = element_text(hjust = 0.5, size = 25)
  ) +
  labs(x = "UMAP_1", y = "UMAP_2") +
  ggtitle("NK cells in scRNA data")



###################Extract the NK cells#################################
NK_subset=subset(cervical_RT_tumor,new.idents="NK_cells")

# Mito-genes
mito_genes <- grep("^MT-", rownames(NK_subset), value = TRUE)

#Workflow for the scRNA process
NK_subset[["percent.mt"]] =PercentageFeatureSet(NK_subset, pattern = "^MT-")
NK_subset <- NormalizeData(NK_subset)
NK_subset <- FindVariableFeatures(NK_subset, selection.method = "vst", nfeatures = 3000)
NK_subset <- ScaleData(NK_subset,vars.to.regress = c("percent.mt"))
NK_subset <- RunPCA(NK_subset)


print(NK_subset[["pca"]], dims = 1:5, nfeatures = 10)
VizDimLoadings(NK_subset, dims = 1:5, reduction = "pca")
ElbowPlot(NK_subset,ndims = 50)

NK_subset <- FindNeighbors(NK_subset, dims = 1:17)
NK_subset <- FindClusters(NK_subset, resolution = 0.5)
NK_subset <- RunUMAP(NK_subset, dims = 1:17)


DimPlot(NK_subset,reduction = "umap",label = T,pt.size = 1)


#Annotation for adaptive NK cells using our NK annotation dataset

lapply(c("dplyr","Seurat","HGNChelper","openxlsx"), library, character.only = T)

source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")

NKdb6_ = "~/Single_cell_Parse_data_6_conditions/New annotations/New_annotations.xlsx"
tissue = "Immune system"
gs_list = gene_sets_prepare(NKdb6_, tissue)


es.max = sctype_score(scRNAseqData = NK_subset[["SCT"]]@scale.data, scaled = TRUE, 
                      gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)


cL_resutls = do.call("rbind", lapply(unique(NK_subset@meta.data$seurat_clusters), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(NK_subset@meta.data[NK_subset@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(NK_subset@meta.data$seurat_clusters==cl)), 10)
}))

NKsctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  
NKsctype_scores
NKsctype_scores$type[as.numeric(as.character(NKsctype_scores$scores)) < NKsctype_scores$ncells/4] = "Unknown"
NKsctype_scores
NKsctype_scores=as.data.frame(NKsctype_scores)


# Add annotation information into the rna_NK object
NK_subset$NK_cell_type <- case_when(
  NK_subset$seurat_clusters == 2 ~ "Adaptive_NK_cells",
  NK_subset$seurat_clusters == 5 ~ "Adaptive_NK_cells",
  TRUE ~ paste0("NK_cluster_", NK_subset$seurat_clusters)
)

Idents(NK_subset)="NK_cell_type"

#Visualize the NK clusters UMAP
cluster_order=c("Adaptive_NK_cells","NK_cluster_0","NK_cluster_1","NK_cluster_3","NK_cluster_4","NK_cluster_6","NK_cluster_7","NK_cluster_8","NK_cluster_9")

NK_subset@meta.data$NK_cell_type=factor(NK_subset@meta.data$NK_cell_type,levels = cluster_order)
Idents(NK_subset)="NK_cell_type"
NK_subset@meta.data[1:5,]
colors = c("#1F77B4", "#AEB6E8", "#FF7F0E", "#FFBB78", "#2CA02C", 
           "#98DF8A", "#D62728", "#FF9896", "#9467BD", "#C5B0D5", 
           "#8C564B", "#C49C94")
DimPlot(
  NK_subset,
  reduction = "umap",
  label = TRUE,
  cols = colors,
  label.size = 6,      
  pt.size = 1.1        
) + theme(
  legend.text = element_text(size = 22), 
  legend.title = element_text(size = 22),
  axis.text.x = element_text(size = 22),     
  axis.text.y = element_text(size = 22),
  axis.title = element_text(size = 22)
)+labs(x = "UMAP_1", y = "UMAP_2")


#Use addmodulescore to score the adaptive NK cluster (Cluster 1)
#Read our annotation geneset
library(readxl)
NK_annotation=read_xlsx("~/Single_cell_Parse_data_6_conditions/New annotations/New_annotations.xlsx")
aNK_gene <- strsplit(NK_annotation$geneSymbolmore1[4], ",")[[1]]
aNK_gene<- trimws(aNK_gene)  
aNK_gene<- unique(aNK_gene)

NK_subset <- AddModuleScore(
  object = NK_subset,
  features = list(aNK_gene),
  name = "aNK_Score"
)
NK_subset@meta.data[1:5,]

#Set adaptive NK cell cluster as the reference group 

ref_group <- "Adaptive_NK_cells"
df <- NK_subset@meta.data %>%
  dplyr::select(aNK_Score1, NK_cell_type) %>%
  dplyr::rename(score = aNK_Score1, group = NK_cell_type) %>%
  dplyr::mutate(group = as.factor(group))

mat <- pairwise.wilcox.test(df$score, df$group, p.adjust.method = "BH")$p.value

# Initialize the list
pval_list <- list()

# One situation：ref_group  in rownames
if (ref_group %in% rownames(mat)) {
  idx <- which(!is.na(mat[ref_group, ]))
  if (length(idx) > 0) {
    pval_list[[length(pval_list) + 1]] <- data.frame(
      group1 = ref_group,
      group2 = names(idx),
      p.adj = unname(mat[ref_group, idx])
    )
  }
}

# One situation:ref_group in colnames
if (ref_group %in% colnames(mat)) {
  idx <- which(!is.na(mat[, ref_group]))
  if (length(idx) > 0) {
    pval_list[[length(pval_list) + 1]] <- data.frame(
      group1 = rownames(mat)[idx],
      group2 = ref_group,
      p.adj = unname(mat[idx, ref_group])
    )
  }
}

# Combine all the results
pval_df <- do.call(rbind, pval_list)

# Add y.position
pval_df$y.position <- seq(max(df$score) * 1.05, by = 0.2, length.out = nrow(pval_df))

# Transder p vaule to asterisk marks
pval_df$p.signif <- cut(pval_df$p.adj,
                        breaks = c(-Inf, 0.0001, 0.001, 0.01, 0.05, Inf),
                        labels = c("****", "***", "**", "*", "ns"))
library(ggpubr)
#Visualize the violin map showing the comparsion of aNK score across all the clusters
ggviolin(
  df, x = "group", y = "score",
  fill = "group",
  palette = cors,
  trim = FALSE,           
  draw_quantiles = NULL   
) +
  geom_jitter(
    width = 0.15,
    size = 1,
    alpha = 0.8,
    shape = 16,            
    color = "black"
  ) +
  stat_pvalue_manual(
    pval_df,
    label = "p.signif",
    y.position = "y.position",
    tip.length = 0.01,
    size = 4.5,
  ) +
  scale_y_continuous(
    limits = c(0, 2.2),
    breaks = seq(0, 2.2, by = 0.2) 
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid = element_blank(),        
    axis.line = element_line(color = "black", size = 0.4),
    axis.ticks = element_line(size = 0.3),
    axis.text.x = element_text(angle = 60, hjust = 1, size = 25),
    axis.text.y = element_text(size = 25),
    axis.title = element_text(size = 25),
    legend.position = "none"
  ) +
  labs(
    x = NULL,
    y = "Adaptive NK Module Score"
  )

########Visualize the distribution of NK clusters across all the conditions 

library(patchwork)
bar.df=NK_subset@meta.data
text.df=as.data.frame(table(bar.df$NK_cell_type,bar.df$category))

colnames(text.df)=c("NK_Cluster","RT_status","Freq")
library(dplyr)
library(ggplot2)

bar.df = NK_subset@meta.data
bar.df$category <- factor(bar.df$category,
                          levels = c("preRT", "onRT1", "onRT2"))

count.df <- as.data.frame(table(bar.df$NK_cell_type, bar.df$category))
colnames(count.df) = c("NK_Cluster", "CRT_status", "Count")

prop.df <- count.df %>%
  group_by(CRT_status) %>%
  mutate(Freq = Count / sum(Count))

prop.df$RT_status <- factor(prop.df$CRT_status,
                            levels = c("preRT", "onRT1", "onRT2"))

ggplot(prop.df, aes(x = CRT_status, y = Freq, fill = NK_Cluster)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_y_continuous("Frequency", expand = c(0.02, 0), limits = c(0, 1)) +
  scale_fill_manual(values = cors) +
  geom_text(aes(label = sprintf("%.2f", Freq)),
            position = position_stack(vjust = 0.5), size = 5.5) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 22),
    axis.text.y = element_text(size = 22),
    axis.title.x = element_text(size = 22),
    axis.title.y = element_text(size = 22),
    legend.text = element_text(size = 22),
    legend.title = element_text(size = 22)
  )


table(NK_subset$category)
table(NK_subset$category, NK_subset$NK_cell_type)

base::saveRDS(NK_subset, "~/Wanmeng/cervical_scRNA_RT/NK_subset_RT.rds")
NK_subset=readRDS("~/Wanmeng/cervical_scRNA_RT/NK_subset_RT.rds")



#######################################################################
#################DEGs analysis #######################################
###################################################################
NK_subset=readRDS("~/Wanmeng/cervical_scRNA_RT/NK_subset_RT.rds")
library(ggrepel)
NK_subset@meta.data[1:5,]

##Extract the adaptive NK cells

aNK <- subset(NK_subset, subset = NK_cell_type == "Adaptive_NK_cells")
table(aNK$category) 

Idents(aNK) <- "category"  
unique(aNK@meta.data$sampleID_v2)
NK_subset@meta.data$patientID <- sub("_.*$", "", NK_subset@meta.data$sampleID_v2)
aNK@meta.data$patientID <- sub("_.*$", "", aNK@meta.data$sampleID_v2)


deg_pre_vs_on1 <- FindMarkers(
  aNK,
  ident.1 = "onRT1",
  ident.2 = "preRT",
  test.use = "wilcox",
  logfc.threshold = 0,
  min.pct = 0.1
)

deg_pre_vs_on1["PRDM1", ]
write.csv(deg_pre_vs_on1,file="~/Wanmeng/cervical_scRNA_RT/DE_files/OnRT1_vs_preRT.markers.csv")

deg_pre_vs_on1$group = 'ns'
deg_pre_vs_on1$group[which((deg_pre_vs_on1$p_val_adj<=0.1)&(deg_pre_vs_on1$avg_log2FC>=0.5)&(deg_pre_vs_on1$pct.1-deg_pre_vs_on1$pct.2>0))]='up'
deg_pre_vs_on1$group[which((deg_pre_vs_on1$p_val_adj<=0.1)&(deg_pre_vs_on1$avg_log2FC<=-0.5)&(deg_pre_vs_on1$pct.1-deg_pre_vs_on1$pct.2<0))]='down'
deg_pre_vs_on1$group=factor(deg_pre_vs_on1$group,levels=c("down","ns","up"))

table(deg_pre_vs_on1$group)

volcano.DEs=ggplot(deg_pre_vs_on1,aes(x=avg_log2FC,y=-log10(p_val_adj),color=group))+geom_point(alpha=0.4, size=2.5)+
  theme_bw(base_size = 12)+
  xlab("Log2FC") + # x轴名字
  ylab("-Log10Padj")+
  scale_colour_manual(values = c('blue','gray','brown')) + # 各自的颜色# y轴名字
  geom_vline(xintercept=c(-0.25,0.25),lty=4,colour="grey11")+
  geom_hline(yintercept=1.3,lty=4,colour="grey11") + 
  labs(title = "OnRT1.vs.PreRT in aNK") +
  theme(
    legend.title = element_blank(),
    legend.text = element_text(size=22),
    axis.text.x = element_text(colour="black", size=22), # X轴字体大小
    axis.text.y = element_text(colour="black", size=22), # Y轴字体大小
    plot.title = element_text(hjust=0.5, size=22), # 标题字体大小
    axis.title.x = element_text(size=22), # X轴标题字体大小
    axis.title.y = element_text(size=22)  # Y轴标题字体大小
  )+guides(color = guide_legend(override.aes = list(shape = 16, alpha = 1, size = 4)))

volcano.DEs

#add the gene names 

subsetDEs.down=subset(deg_pre_vs_on1,group=="down")
subsetDEs.up=subset(deg_pre_vs_on1,group=="up")

subsetDEs=rbind(subsetDEs.down,subsetDEs.up)
subsetDEs$gene=rownames(subsetDEs)

prdm1_df <- deg_pre_vs_on1["PRDM1", , drop = FALSE]
prdm1_df$gene <- "PRDM1"


volcano.DEs_2=volcano.DEs+geom_text_repel(data=subsetDEs,aes(label=gene),size=4.5)+
  geom_text_repel(data = prdm1_df, aes(label = gene),
                  size = 4.5, color = 'brown')
volcano.DEs_2


deg_pre_vs_on2 <- FindMarkers(
  aNK,
  ident.1 = "onRT2",
  ident.2 = "preRT",
  test.use = "wilcox",
  logfc.threshold = 0,
  min.pct = 0.1
)

write.csv(deg_pre_vs_on2,file="~/Wanmeng/cervical_scRNA_RT/DE_files/OnRT2_vs_preRT.markers.csv")

deg_pre_vs_on2$group = 'ns'
deg_pre_vs_on2$group[which((deg_pre_vs_on2$p_val_adj<=0.1)&(deg_pre_vs_on2$avg_log2FC>=0.5)&(deg_pre_vs_on2$pct.1-deg_pre_vs_on2$pct.2>0))]='up'
deg_pre_vs_on2$group[which((deg_pre_vs_on2$p_val_adj<=0.1)&(deg_pre_vs_on2$avg_log2FC<=-0.5)&(deg_pre_vs_on2$pct.1-deg_pre_vs_on2$pct.2<0))]='down'
deg_pre_vs_on2$group=factor(deg_pre_vs_on2$group,levels=c("down","ns","up"))

table(deg_pre_vs_on2$group)

volcano.DEs=ggplot(deg_pre_vs_on2,aes(x=avg_log2FC,y=-log10(p_val_adj),color=group))+geom_point(alpha=0.4, size=2.5)+
  theme_bw(base_size = 12)+
  xlab("Log2FC") + 
  ylab("-Log10Padj")+
  scale_colour_manual(values = c('blue','gray','brown')) +
  geom_vline(xintercept=c(-0.25,0.25),lty=4,colour="grey11")+
  geom_hline(yintercept=1.3,lty=4,colour="grey11") + 
  labs(title = "OnRT2.vs.PreRT in aNK") +
  theme(
    legend.title = element_blank(),
    legend.text = element_text(size=22),
    axis.text.x = element_text(colour="black", size=22), 
    axis.text.y = element_text(colour="black", size=22), #
    plot.title = element_text(hjust=0.5, size=22), 
    axis.title.x = element_text(size=22), 
    axis.title.y = element_text(size=22)  
  )+guides(color = guide_legend(override.aes = list(shape = 16, alpha = 1, size = 4)))

volcano.DEs

#add the gene names
subsetDEs.down=subset(deg_pre_vs_on2,group=="down")
subsetDEs.up=subset(deg_pre_vs_on2,group=="up")

subsetDEs=rbind(subsetDEs.down,subsetDEs.up)
subsetDEs$gene=rownames(subsetDEs)

volcano.DEs_2=volcano.DEs+geom_text_repel(data=subsetDEs,aes(label=gene),size=4.5)
volcano.DEs_2


subsetDEs.up1=subset(deg_pre_vs_on1,group=="up")
subsetDEs.up2=subset(deg_pre_vs_on2,group=="up")
common_genes=intersect(rownames(subsetDEs.up1),rownames(subsetDEs.up2))

write.csv(common_genes,"~/Wanmeng/cervical_scRNA_RT/DE_files/common_DE_genes_up.csv")

##Visualize the overlapped genes 
library(UpSetR)
library(ggvenn)

gene_list=list(onRT1_vs_preRT=rownames(subsetDEs.up1),
               onRT2_vs_preRT=rownames(subsetDEs.up2))


gene_list=list(onRT1_vs_preRT=rownames(subsetDEs.up1),
               onRT2_vs_preRT=rownames(subsetDEs.up2))
data_veen = list_to_data_frame(gene_list)

ggvenn(data_veen,
       show_percentage = T,
       show_elements = F,
       text_size=10,
       digits = 1,
       set_name_size=10,
       stroke_color = "grey30",
       fill_color = c("#FF8C00","#4DAF4A","#B64E89"),
       set_name_color =c("#FF8C00","#4DAF4A","#B64E89"))




######GO enrichment for the overlapped genes########################

library(clusterProfiler)
library(org.Hs.eg.db)

cellID=bitr(common_genes,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = "org.Hs.eg.db")
cellID
Common_genes_GO=enrichGO(gene=cellID$SYMBOL,OrgDb = "org.Hs.eg.db",keyType="SYMBOL",ont = "BP"
                         ,pAdjustMethod = "BH",pvalueCutoff = 0.01,qvalueCutoff = 0.05)


Common_genes_GO=as.data.frame(Common_genes_GO@result)
Common_genes_GO=Common_genes_GO[order(Common_genes_GO$p.adjust),]
Common_genes_GO_30=Common_genes_GO[c(1:30),]
Common_genes_GO_30=Common_genes_GO_30[order(Common_genes_GO_30$Count,decreasing = F),]
Common_genes_GO_30$Description=factor(Common_genes_GO_30$Description,levels = Common_genes_GO_30$Description)


ggplot(Common_genes_GO_30, aes(x = Description, y = Count)) +
  geom_bar(aes(fill = p.adjust), stat = "identity", width = 0.8) +
  coord_flip() +
  scale_fill_gradient(limits = c(0, max(Common_genes_GO_30$p.adjust)), low = "red", high = "blue") +
  xlab("GO_BP_term") +
  ylab("Num_of_genes") +
  theme_bw() +
  theme(
    axis.text = element_text(size = 20),        
    axis.title = element_text(size = 20),      
    legend.title = element_text(size = 20),     
    legend.text = element_text(size = 18),      
    plot.title = element_text(size = 20, hjust = 0.5),
    axis.title.x = element_text(size = 20),     
    axis.title.y = element_text(size = 20)
  )

#Rename the pathways 
short_name <-c(
  "Response to virus",
  "Antiviral defense",
  "Tissue remodeling",
  "p53 signaling",
  "Viral process regulation",
  "Viral genome replication",
  "Host process modulation",
  "Viral life cycle regulation",
  "p53-mediated apoptosis",
  "Viral genome repl. regulation",
  "Negative viral regulation",
  "Negative viral genome repl.",
  "γ-radiation response",
  "Vascular remodeling",
  "Transposition",
  "Retrotransposition",
  "TE silencing",
  "Placenta development",
  "Ca2+ transport (cytosol)",
  "DNA modification",
  "Viral RNA replication",
  "Ca2+ transport regulation",
  "Base editing",
  "ssRNA repl. via dsDNA",
  "Regulation of ssRNA repl.",
  "Upregulation of Ca2+ transport",
  "DNA deamination",
  "C-to-U RNA editing",
  "Granzyme-mediated apoptosis",
  "Negative ssRNA repl."
)
Common_genes_GO_30$short_name <- rev(short_name)

Common_genes_GO_30=Common_genes_GO_30[order(Common_genes_GO_30$Count,decreasing = F),]
Common_genes_GO_30$short_name=factor(Common_genes_GO_30$short_name,levels = Common_genes_GO_30$short_name)

ggplot(Common_genes_GO_30, aes(x = short_name, y = Count)) +
  geom_bar(aes(fill = p.adjust), stat = "identity", width = 0.8) +
  coord_flip() +
  scale_fill_gradient(limits = c(0, max(Common_genes_GO_30$p.adjust)), low = "red", high = "blue") +
  xlab("GO_BP_term") +
  ylab("Num_of_genes") +
  theme_bw() +
  theme(
    axis.text = element_text(size = 20),        
    axis.title = element_text(size = 20),      
    legend.title = element_text(size = 20),    
    legend.text = element_text(size = 18),     
    plot.title = element_text(size = 20, hjust = 0.5),
    axis.title.x = element_text(size = 20),    
    axis.title.y = element_text(size = 20)
  )


############################Psedotime analysis##########################
#############################Monocle3#################################
library(leidenbase)
library(Seurat)
library(monocle3)
library(Matrix)
library(ggplot2)


NK_subset=readRDS("~/Wanmeng/cervical_scRNA_RT/NK_subset_RT.rds")
expression_matrix = as(as.matrix(NK_subset@assays$SCT@counts), 'sparseMatrix')
cell_metadata = NK_subset@meta.data
gene_annotation = data.frame(gene_short_name = rownames(expression_matrix))
rownames(gene_annotation) = rownames(expression_matrix)

cds = new_cell_data_set(expression_matrix,
                        cell_metadata = cell_metadata,
                        gene_metadata = gene_annotation)


cds = preprocess_cds(cds, num_dim = 100)


saveRDS(cds,"~/Wanmeng/cervical_scRNA_RT/monocle3/NK_subset_mono3.rds")
cds=readRDS("~/Wanmeng/cervical_scRNA_RT/monocle3/NK_subset_mono3.rds")
plot_pc_variance_explained(cds)
cds = reduce_dimension(cds,reduction_method='UMAP',
                       preprocess_method = 'PCA')


cds = cluster_cells(cds)
plot_cells(cds, reduction_method="UMAP",color_cells_by="NK_cell_type",cell_size=0.5,group_label_size=5)
plot_cells(cds,cell_size=0.5,group_label_size=5)

##replace the UMAP information with the previous oneee
cds.embed = cds@int_colData$reducedDims$UMAP
int.embed = Embeddings(NK_subset, reduction = "umap")
int.embed = int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP = int.embed 
#cds@int_colData$reducedDims$tSNE = int.embed
plot_cells(cds, color_cells_by="NK_cell_type",reduction_method = "UMAP",
           cell_size=0.5,group_label_size=4) 


saveRDS(cds,"~/Wanmeng/cervical_scRNA_RT/monocle3/NK_subset_mono3.rds")
cds=readRDS("~/Wanmeng/cervical_scRNA_RT/monocle3/NK_subset_mono3.rds")

mycds = cds
mycds = learn_graph(mycds,use_partition = F,
                    verbose=T,
                    learn_graph_control=list(minimal_branch_len=5,
                                             euclidean_distance_ratio=0.8, prune_graph=T))

plot_cells(mycds, 
           color_cells_by = "NK_cell_type",
           label_groups_by_cluster=F,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           cell_size=0.5,group_label_size=4)

mycds1 = mycds
mycds1 = order_cells(mycds1)

#pseudotime
plot_cells(mycds1, 
           label_cell_groups = FALSE, 
           color_cells_by = "pseudotime", 
           label_leaves = FALSE, 
           label_branch_points = T, 
           graph_label_size = 4.5, 
           cell_size = 0.9, 
           show_trajectory_graph = TRUE,
           trajectory_graph_segment_size = 1, 
           trajectory_graph_color = "gray10") + 
  xlab("UMAP_1") + 
  ylab("UMAP_2") + 
  theme(legend.title = element_text(size = 25),
        legend.text = element_text(size =  25), 
        axis.title.x = element_text(size = 25), 
        axis.title.y = element_text(size = 25),
        axis.text.x = element_text(size = 25),  
        axis.text.y = element_text(size = 25)
  )

#RT status
plot_cells(mycds1, 
           label_cell_groups = FALSE, 
           color_cells_by = "category", 
           label_leaves = FALSE, 
           label_branch_points = T, 
           graph_label_size = 4.5, 
           cell_size = 0.9, 
           show_trajectory_graph = TRUE,
           trajectory_graph_segment_size = 1, 
           trajectory_graph_color = "gray10") + 
  scale_color_manual(values = c(
    "preRT" = "#ECECEC",  
    "onRT1" = "#D55E00",  
    "onRT2" = "#0072B2" 
  ))+
  xlab("UMAP_1") + 
  ylab("UMAP_2") + 
  theme(legend.title = element_text(size = 25),
        legend.text = element_text(size =  25),  
        axis.title.x = element_text(size = 25),  
        axis.title.y = element_text(size = 25),
        axis.text.x = element_text(size = 25),   
        axis.text.y = element_text(size = 25)
  )







#或者函数推断
mycds2 = mycds
get_earliest_principal_node <- function(cds){
  cell_ids <- dim(cds)[2]
  closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <- igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names(which.max(table(closest_vertex[cell_ids,]))))]
  root_pr_nodes
}
mycds2 <- mycds
mycds2 <- order_cells(mycds2, root_pr_nodes=get_earliest_principal_node(mycds2))



plot_cells(mycds2, 
           label_cell_groups = FALSE, 
           color_cells_by = "pseudotime", 
           label_leaves = FALSE, 
           label_branch_points = T, 
           graph_label_size = 4.5, 
           cell_size = 0.9, 
           show_trajectory_graph = TRUE,
           trajectory_graph_segment_size = 1.5, 
           trajectory_graph_color = "gray10") + 
  xlab("UMAP_1") + 
  ylab("UMAP_2") + 
  theme(legend.title = element_text(size = 18),
        legend.text = element_text(size =  15),  # 调整图例字体大小
        axis.title.x = element_text(size = 20),  # 调整 X 轴标题字体大小
        axis.title.y = element_text(size = 20),
        axis.text.x = element_text(size = 15),    # 调整 X 轴刻度字体大小
        axis.text.y = element_text(size = 15)
  )


#################Visualize the PRDM1 expression in pseudotime UMAP
pd=pseudotime(mycds1, reduction_method = 'UMAP')
NK_subset = AddMetaData(NK_subset,metadata = pd,col.name = 'pseudotime')
NK_subset@meta.data[1:5,]

mycds1_res = graph_test(mycds1, 
                        neighbor_graph="principal_graph", cores=4)
write.csv(mycds1_res, file = "~/Wanmeng/cervical_scRNA_RT/monocle3/NK_subset_Monocle3_mycds1_res.csv")
mycds1_res=read.csv("~/Wanmeng/cervical_scRNA_RT/monocle3/NK_subset_Monocle3_mycds1_res.csv",row.names = 1)


res_ids = row.names(subset(mycds1_res, q_value < 0.05))
gene_module_df = find_gene_modules(mycds1[res_ids,], 
                                   resolution=c(10^seq(-6,-1)))

cell_group_df = tibble::tibble(cell=row.names(colData(mycds1)), 
                               cell_group=colData(mycds1)$NK_cell_type)

agg_mat = aggregate_gene_expression(mycds1, 
                                    gene_module_df, 
                                    cell_group_df)
row.names(agg_mat) = stringr::str_c("Module", row.names(agg_mat))
pheatmap::pheatmap(agg_mat,
                   scale="column", clustering_method="ward.D2",cluster_rows = T,cluster_cols = F)


library(dplyr)
genes_sig = mycds1_res %>% top_n(n=10, morans_I) %>%
  pull(gene_short_name) %>% as.character()


plot_genes_in_pseudotime(mycds1[genes_sig,], color_cells_by="NK_cell_type", 
                         min_expr=0.5, ncol = 2)+ xlab("UMAP_1") + ylab("UMAP_2")

plot_genes_in_pseudotime(mycds1[genes_sig,], color_cells_by="pseudotime", 
                         min_expr=0.5, ncol = 2)+ xlab("UMAP_1") + ylab("UMAP_2")

plot_cells(mycds1, genes=genes_sig, show_trajectory_graph=T,
           label_cell_groups=T,  label_leaves=FALSE)+ xlab("UMAP_1") + ylab("UMAP_2")

library(viridis)
#Visualize the PRDM1 expression in pseudotime UMAP
p <- plot_cells(mycds1,
                genes=c("PRDM1"),           
                label_cell_groups=TRUE,
                #color_cells_by = "NK_cell_type",
                show_trajectory_graph=TRUE, 
                cell_size=1, trajectory_graph_color="black", 
                label_branch_points=TRUE, 
                label_roots=FALSE, label_leaves=FALSE) +
  scale_color_viridis(option="viridis",direction = -1) +
  xlab("UMAP_1") + ylab("UMAP_2") +
  ggtitle("PRDM1 expression") + 
  theme(
    plot.title = element_text(size=20, hjust=0.5),
    axis.title.x = element_text(size=18),
    axis.title.y = element_text(size=18),
    axis.text.x  = element_text(size=18),
    axis.text.y  = element_text(size=18),
    legend.title = element_text(size=18),
    legend.text  = element_text(size=18),
    strip.text = element_text(size=18),
    text = element_text(size=15)   
  )

p



#########Pseudotime module heatmap 
mycds1_res=read.csv("~/Wanmeng/cervical_scRNA_RT/monocle3/NK_subset_Monocle3_mycds1_res.csv",row.names = 1)
genes = row.names(subset(mycds1_res, q_value< 0.05 & morans_I > 0.2))

plot_matrix = exprs(mycds1)[match(genes,
                                  rownames(rowData(mycds1))),
                            order(pseudotime(mycds1))]
plot_matrix = t(apply(plot_matrix,1,function(x){smooth.spline(x,df=3)$y}))
plot_matrix = t(apply(plot_matrix,1,function(x){(x-mean(x))/sd(x)}))
rownames(plot_matrix) = genes;
dim(plot_matrix)

############Construct a function##
cutColumn_Means <- function(data_exp, 
                            cut       
) {
  plot_matrix_combin <- list() # 
  nums <- ncol(data_exp) / cut #
  if (nums - round(nums, 0) == 0) {
    splits <- seq(1, ncol(data_exp), by = cut) 
  } else {
    splits <- seq(1, ncol(data_exp) - cut + 1, by = cut)
  }
  
  for (i in seq_along(splits)) {
    col_range <- splits[i]:(min(splits[i] + cut - 1, ncol(data_exp)))
    A <- as.data.frame(rowMeans(data_exp[, col_range]))
    plot_matrix_combin[[i]] <- A
  }
  plot_matrix_combin <- do.call(cbind, plot_matrix_combin)
  rownames(plot_matrix_combin) <- rownames(data_exp)
  colnames(plot_matrix_combin) <- seq_len(ncol(plot_matrix_combin))
  
  return(plot_matrix_combin)
}

plot_test = cutColumn_Means(plot_matrix,cut = 25)

callback = function(hc, mat){
  sv = svd(t(mat))$v[,1]
  dend = reorder(as.dendrogram(hc), wts = sv)
  as.hclust(dend)
}
library(pheatmap)
p2 = pheatmap::pheatmap(plot_test,
                        useRaster = T,
                        cluster_cols=FALSE,
                        cluster_rows=T,
                        show_rownames=F,
                        show_colnames=F,
                        clustering_method = "ward.D2",
                        cutree_rows=4,
                        filename=NA,
                        border_color = NA,
                        fontsize_row = 8,
                        color=colorRampPalette(c('#1A5592','white',"#B83D3D"))(100),
                        clustering_callback = callback)

annotation_row = data.frame(Cluster=factor(cutree(p2$tree_row, 4)))
row.names(annotation_row) = rownames(plot_test)
rowcolor = c("#85B22E","#E29827","#922927",'#57C3F3')
names(rowcolor) = c("1","2","3","4") 
ann_colors = list(Cluster=rowcolor)

p3 = pheatmap::pheatmap(plot_test,
                        cluster_cols=FALSE,
                        cluster_rows=T,
                        show_rownames=F,
                        show_colnames=F,
                        clustering_method = "ward.D2",
                        filename=NA,
                        cutree_rows=4,
                        border_color = NA,
                        fontsize_row = 8,
                        color=colorRampPalette(c('#1A5592','white',"#B83D3D"))(100),
                        annotation_colors=ann_colors,
                        annotation_row = annotation_row,
                        clustering_callback = callback,
                        annotation_names_col = F,
                        annotation_names_row = F,
                        main="Pseudotime")

celltype_marker=c("CD69","KLRC1","KLRC2","KIR2DL1","KIR2DL3","IFNG","PDCD1","HAVCR2","LAG3","NCR2","SYK","FCER1G","ZNF384","CD7",
                  "NCAM1","ZBTB16","TIGIT","XCL2","GZMA","GZMB","SELE","SELL","MKI67","KLRD1","B3GAT1","FCGR3A","GPRC5A","ITGA1","ITGA2","ITGAE","ENTPD1")
genes=union(celltype_marker,common_genes)
gene=genes
p4 = pheatmap::pheatmap(plot_test,
                        cluster_cols=FALSE,
                        cluster_rows=T,
                        show_rownames=T,
                        show_colnames=F,
                        clustering_method = "ward.D2",
                        filename=NA,
                        border_color = NA,
                        fontsize_row = 8,
                        color=colorRampPalette(c('#1A5592','white',"#B83D3D"))(100),
                        annotation_colors=ann_colors,
                        annotation_row = annotation_row,
                        clustering_callback = callback,
                        annotation_names_col = F,
                        annotation_names_row = F,
                        main="Pseudotime")

library(grid)
source('~/data_my_single_cell/R_code_file/add.flag.R')
add.flag(p4,kept.labels = gene,repel.degree = 0.4,label.size=13)



##############GO enrichment for the pseudotime modules 
library(clusterProfiler)
library(org.Hs.eg.db)
###Extract the module genes 
module_gene <- as.data.frame(cutree(p3$tree_row, k=4))
colnames(module_gene) <- "Module"
module_gene$gene <- rownames(module_gene)

Module_GO=data.frame()
for (i in unique(module_gene$Module)) {
  data=filter(module_gene,module_gene$Module==i)
  df=bitr(data$gene,
          fromType="SYMBOL",
          toType=c("ENTREZID"),
          OrgDb="org.Hs.eg.db")
  go <- enrichGO(gene= unique(df$ENTREZID),
                 OrgDb= org.Hs.eg.db,
                 keyType= 'ENTREZID',
                 ont= "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff= 0.05,
                 qvalueCutoff= 0.05,
                 readable= TRUE)
  go_res=go@result
  if (dim(go_res)[1] != 0) {
    go_res$cluster=i
    Module_GO=rbind(Module_GO,go_res)
  }
}

###PRDM1 expression across the pseudotime and in three CRT conditions 

NK_subset$Cell_barcode=rownames(NK_subset@meta.data)
peuGene_exp = normalized_counts(mycds1, norm_method = "log")
peuGene_exp = as.data.frame(peuGene_exp)
dim(peuGene_exp)
pseudotime = as.data.frame(pseudotime(mycds1))
colnames(pseudotime) = "pseudotime"
pseudotime$Cell_barcode=rownames(pseudotime)
peuGene_exp = peuGene_exp[,colnames(peuGene_exp) %in% rownames(pseudotime)]
dim(peuGene_exp)

df = data.frame(
  Status = NK_subset@meta.data$category,
  Cell_barcode = NK_subset@meta.data$Cell_barcode
)

pseudotime = merge(df, pseudotime, by = "Cell_barcode")

##PreRT 
target_cells_PreRT = pseudotime$Cell_barcode[pseudotime$Status == "preRT"]
peuGene_exp_PreRT = peuGene_exp[, colnames(peuGene_exp) %in% target_cells_PreRT]
pseudotime_PreRT= pseudotime[pseudotime$Status=="preRT",]
##onRT1
target_cells_onRT1 = pseudotime$Cell_barcode[pseudotime$Status == "onRT1"]
peuGene_exp_onRT1 = peuGene_exp[, colnames(peuGene_exp) %in%target_cells_onRT1]
pseudotime_onRT1= pseudotime[pseudotime$Status=="onRT1",]
##onRT2
target_cells_onRT2 = pseudotime$Cell_barcode[pseudotime$Status == "onRT2"]
peuGene_exp_onRT2 = peuGene_exp[, colnames(peuGene_exp) %in% target_cells_onRT2]
pseudotime_onRT2= pseudotime[pseudotime$Status=="onRT2",]



PreRT_PRDM1 = peuGene_exp_PreRT['PRDM1',]
PreRT_PRDM1 = as.data.frame(t(PreRT_PRDM1))
PreRT_PRDM1$pseudotime = pseudotime_PreRT$pseudotime
PreRT_PRDM1$group = 'preRT'


OnRT1_PRDM1 = peuGene_exp_onRT1['PRDM1',]
OnRT1_PRDM1 = as.data.frame(t(OnRT1_PRDM1))
OnRT1_PRDM1$pseudotime = pseudotime_onRT1$pseudotime
OnRT1_PRDM1$group = 'onRT1'

OnRT2_PRDM1 = peuGene_exp_onRT2['PRDM1',]
OnRT2_PRDM1  = as.data.frame(t(OnRT2_PRDM1))
OnRT2_PRDM1 $pseudotime = pseudotime_onRT2$pseudotime
OnRT2_PRDM1 $group = 'onRT2'


peu_trend = rbind(PreRT_PRDM1,OnRT1_PRDM1,OnRT2_PRDM1)

peu_trend_cleaned <- peu_trend[is.finite(peu_trend$pseudotime), ]

cols <- c(preRT="#8A2BE2", onRT1="#FF6347", onRT2="#32CD32")
base_font <- 18  

ggplot(peu_trend_cleaned, aes(x=pseudotime, y=PRDM1, color=group)) +
  geom_point(size=1, alpha=0.25) +
  geom_smooth(method="loess", span=0.5, se=FALSE, linewidth=1.2) +
  scale_color_manual(values=cols, name=NULL) +
  labs(x="Pseudotime", y="Relative expression", title="PRDM1") +
  theme_classic(base_size=base_font) +
  theme(
    text = element_text(size = base_font),
    plot.title = element_text(size = base_font + 4, hjust=0.5),
    axis.title = element_text(size = base_font + 2),
    axis.text  = element_text(size = base_font),
    legend.text = element_text(size = 20),
    legend.position = "top"
  )


##############################################################################
###############scTenifoldKnk virtual KO#######################################

#devtools::install_github("cailab-tamu/scTenifoldKnk")

library(scTenifoldKnk)
NK_subset=readRDS("~/Wanmeng/cervical_scRNA_RT/NK_subset_RT.rds")
NK_subset@meta.data[1:5,]
NK_matrix <- GetAssayData(NK_subset, assay = "SCT", slot = "counts")
aNK <- subset(NK_subset, subset = NK_cell_type == "Adaptive_NK_cells")
table(aNK$category) 
aNK_matrix <- GetAssayData(aNK, assay = "SCT", slot = "counts")

rm_pattern <- "^(RPL|RPS|MRPL|MRPS|MT-|MTRNR|MTRNR2L|HIST)"
keep <- !grepl(rm_pattern, rownames(aNK_matrix))
aNK_matrix_noRibo <- aNK_matrix[keep, ]

saveRDS(aNK_matrix_noRibo,"~/Wanmeng/cervical_scRNA_RT/scTenifoldKnk/aNK_matrix_noRibo.rds")

###

aNK@meta.data
median(aNK$nFeature_RNA)
mean(aNK$nFeature_RNA)
hist(aNK$nCount_RNA, breaks = 50)

res <- scTenifoldKnk(
  countMatrix = aNK_matrix,
  qc = FALSE,
  gKO = "PRDM1",
  nc_nCells = 500,    
  nc_nNet   = 10,       
  nc_nComp  = 30,
  nCores=10
)

str(res$diffRegulation)
head(res$diffRegulation)
diff_reg <- res$diffRegulation
diff_reg <- diff_reg[order(diff_reg$p.value), ]

head(diff_reg)
gene_affected=diff_reg$gene[diff_reg$Z>2]
gene_affected <- gene_affected[gene_affected != "PRDM1"]
df=as.data.frame(gene_affected)
write.csv(df,"~/Wanmeng/cervical_scRNA_RT/scTenifoldKnk/genes.affected.csv")

library(ggplot2)
library(ggrepel)
library(scales)
diff_reg <- diff_reg[diff_reg$gene != "PRDM1", ]
diff_reg$group <- ifelse(diff_reg$Z > 2, "Significant", "Not significant")
diff_reg$group <- factor(diff_reg$group, levels = c("Not significant","Significant"))
label_df <- diff_reg[diff_reg$gene %in% gene_affected, ]


y_cap <- 1.0
diff_reg$group <- ifelse(diff_reg$Z > 2, "Significant", "Not significant")
diff_reg$group <- factor(diff_reg$group, levels = c("Not significant","Significant"))

diff_reg$y <- -log10(diff_reg$p.value)
diff_reg$capped <- diff_reg$y > y_cap
diff_reg$y_plot <- pmin(diff_reg$y, y_cap)
label_df <- diff_reg[diff_reg$gene %in% gene_affected, ]

##Visualize the affected genes under PRDM1 KO
ggplot(diff_reg, aes(x = Z, y = y_plot, color = group)) +
  geom_point(size = 1.6, alpha = 0.7) +
  geom_point(
    data = subset(diff_reg, capped),
    aes(x = Z, y = y_cap),
    shape = 24, size = 2.2, stroke = 0.5
  ) +
  geom_vline(xintercept = c(-2, 2), linetype = "dashed", linewidth = 0.6) +
  geom_text_repel(
    data = label_df,
    aes(label = gene),
    size = 4.3,
    max.overlaps = Inf,
    box.padding = 0.4,
    point.padding = 0.25,
    min.segment.length = 0
  ) +
  scale_color_manual(values = c("Not significant" = "grey75",
                                "Significant"     = "#D73027")) +
  scale_y_continuous(
    limits = c(0, y_cap),
    breaks = seq(0, y_cap, 0.2),
    expand = expansion(mult = c(0.02, 0.18))
  ) +
  coord_cartesian(clip = "off") +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.title = element_blank(),
    plot.margin = margin(10, 30, 10, 10)  
  ) + theme(
    plot.title = element_text(size = 20, hjust = 0.5),
    axis.title = element_text(size = 20),
    axis.text  = element_text(size = 20),
    legend.text  = element_text(size = 20),  
    legend.title = element_blank(),
    #legend.position = c(0.5, 0.5),
    legend.background = element_rect(fill = "white", color = "grey85", linewidth = 0.3),
    axis.line = element_line(linewidth = 0.6)
  ) +
  labs(title = "PRDM1 virtual KO",
       x = "Z-score",
       y = "-log10(p value) (capped)")


##GO enrichment for the affected genes 
library(clusterProfiler)
cellID=bitr(gene_affected,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = "org.Hs.eg.db")
cellID
gene_affected_GO=enrichGO(gene=cellID$SYMBOL,OrgDb = "org.Hs.eg.db",keyType="SYMBOL",ont = "BP"
                          ,pAdjustMethod = "BH",pvalueCutoff = 0.01,qvalueCutoff = 0.05)
gene_affected_GO=as.data.frame(gene_affected_GO@result)
gene_affected_GO=gene_affected_GO[order(gene_affected_GO$p.adjust),]
gene_affected_GO_30=gene_affected_GO[c(1:30),]

gene_affected_GO_30=gene_affected_GO_30[order(gene_affected_GO_30$Count,decreasing = F),]
gene_affected_GO_30$Description=factor(gene_affected_GO_30$Description,levels = gene_affected_GO_30$Description)

ggplot(gene_affected_GO_30, aes(x = Description, y = Count)) +
  geom_bar(aes(fill = p.adjust), stat = "identity", width = 0.8) +
  coord_flip() +
  scale_fill_gradient(limits = c(0, max(gene_affected_GO_30$p.adjust)), low = "red", high = "blue") +
  xlab("GO_BP_term") +
  ylab("Num_of_genes") +
  theme_bw() +
  theme(
    axis.text = element_text(size = 20),        
    axis.title = element_text(size = 20),       
    legend.title = element_text(size = 20),     
    legend.text = element_text(size = 18),     
    plot.title = element_text(size = 20, hjust = 0.5),
    axis.title.x = element_text(size = 20),   
    axis.title.y = element_text(size = 20)
  )

##GSEA GO analysis for the affected genes

library(dplyr)
diff_reg2 <- diff_reg[diff_reg$gene != "PRDM1", ]
gene_rank <- diff_reg2$Z
names(gene_rank) <- diff_reg2$gene
gene_rank <- sort(gene_rank, decreasing = TRUE)

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(patchwork)
gsea_go <- gseGO(
  geneList     = gene_rank,
  OrgDb        = org.Hs.eg.db,
  keyType      = "SYMBOL",
  ont          = "BP",
  minGSSize    = 10,
  maxGSSize    = 500,
  pvalueCutoff = 0.05,
  verbose      = FALSE
)
dotplot(gsea_go, showCategory = 20) +
  ggtitle("GSEA-GO (ranked by network Z-score)") +
  theme(
    plot.title = element_text(size = 18),  # 标题
    axis.title = element_text(size = 18),                
    axis.text.x =  element_text(size = 18), 
    axis.text.y =  element_text(size = 16),# 坐标轴刻度
    legend.title = element_text(size = 18),               
    legend.text  = element_text(size = 18)              
  )
gsea_go@result$Description

p <- gseaplot2(
  gsea_go,
  geneSetID = 2,
  title = gsea_go@result$Description[2],
  base_size = 25,
  color = "#D73027",                
  rel_heights = c(1.2, 0.35, 0.6)    
)
p
p <- gseaplot2(
  gsea_go,
  geneSetID = 3,
  title = gsea_go@result$Description[3],
  base_size = 25,
  color = "#D73027",                
  rel_heights = c(1.2, 0.35, 0.6)     
)
p

##GSEA KEGG
gene_df <- bitr(
  names(gene_rank),
  fromType = "SYMBOL",
  toType   = "ENTREZID",
  OrgDb    = org.Hs.eg.db
)
gene_rank_kegg <- gene_rank[gene_df$SYMBOL]
names(gene_rank_kegg) <- gene_df$ENTREZID
gene_rank_kegg <- sort(gene_rank_kegg, decreasing = TRUE)

gsea_kegg <- gseKEGG(
  geneList     = gene_rank_kegg,
  organism     = "hsa",
  minGSSize    = 10,
  maxGSSize    = 500,
  pvalueCutoff = 0.05
)
library(ggplot2)
dotplot(gsea_kegg, showCategory = 20) +
  ggtitle("GSEA-KEGG (ranked by network Z-score)") +
  theme(
    plot.title = element_text(size = 18),  
    axis.title = element_text(size = 18),                
    axis.text.x =  element_text(size = 18), 
    axis.text.y =  element_text(size = 16),
    legend.title = element_text(size = 18),               
    legend.text  = element_text(size = 18)                
  )


#######################################################
###########Construct the data set for Cytoscape#########
#######################################################
res=readRDS("~/Wanmeng/cervical_scRNA_RT/scTenifoldKnk/res_PRDM1_scTenifoldKnk_2.rds")
names(res)
str(res, max.level = 1)

names(res$tensorNetworks)
str(res$tensorNetworks, max.level = 2)

W_WT <- res$tensorNetworks$WT
W_KO <- res$tensorNetworks$KO
Wdiff <- W_KO - W_WT

pr <- "PRDM1"
stopifnot(pr %in% rownames(Wdiff), pr %in% colnames(Wdiff))
out <- Wdiff[pr, ]      # 1 x genes
inn <- Wdiff[, pr]      # genes x 1

out_v <- as.numeric(out); names(out_v) <- colnames(Wdiff)
inn_v <- as.numeric(inn); names(inn_v) <- rownames(Wdiff)

N <- 60 
nei <- unique(c(
  names(sort(abs(out_v), decreasing = TRUE))[1:N],
  names(sort(abs(inn_v), decreasing = TRUE))[1:N]
))
nei <- setdiff(nei, pr)
genes_sub <- c(pr, nei)
Wsub <- as.matrix(Wdiff[genes_sub, genes_sub])
idx <- which(Wsub != 0, arr.ind = TRUE)
edge_sub <- data.frame(
  source = rownames(Wsub)[idx[,1]],
  target = colnames(Wsub)[idx[,2]],
  weight = Wsub[idx],
  abs_weight = abs(Wsub[idx]),
  direction = ifelse(Wsub[idx] > 0, "Up_in_KO", "Down_in_KO"),
  stringsAsFactors = FALSE
)
edge_sub <- edge_sub[edge_sub$source != edge_sub$target, ]  
M <- 400  
edge_sub <- edge_sub[order(edge_sub$abs_weight, decreasing = TRUE), ]
edge_sub <- head(edge_sub, M)

##Node data
diff_reg <- res$diffRegulation
diff_reg <- diff_reg[diff_reg$gene != "PRDM1", ]               
gene_affected <- gene_affected[gene_affected != "PRDM1"]         

node_df <- data.frame(
  id = diff_reg$gene,
  Z  = diff_reg$Z,
  pvalue = diff_reg$p.value,
  significant = diff_reg$Z > 2,
  is_affected = diff_reg$gene %in% gene_affected,
  stringsAsFactors = FALSE
)

node_df_sub <- node_df[node_df$id %in% unique(c(edge_sub$source, edge_sub$target)), ]
node_df_sub <- rbind(
  data.frame(id="PRDM1", Z=NA, pvalue=NA, significant=NA, is_affected=FALSE),
  node_df_sub
)
node_df_sub$Z <- as.numeric(node_df_sub$Z)
write.csv(edge_sub,
          "~/Wanmeng/cervical_scRNA_RT/scTenifoldKnk/PRDM1_subnetwork_edges_clean.csv",
          row.names = FALSE,
          quote = FALSE)

write.csv(node_df_sub,
          "~/Wanmeng/cervical_scRNA_RT/scTenifoldKnk/PRDM1_subnetwork_nodes_clean.csv",
          row.names = FALSE,
          quote = FALSE,
          na = "")





##########Cell oracle preparation 
library(Seurat)
library(SeuratDisk)
library(purrr)
library(stringr)
NK_subset_for_celloracle=readRDS("~/Wanmeng/cervical_scRNA_RT/Cell_oracle/NK_subset_RT.rds")

###First file: Replace the scale data with the counts matrix
obj_2 <- NK_subset_for_celloracle
obj_2$NK_cell_type=as.character(obj_2$NK_cell_type)
genes_use <- rownames(obj_2@assays$SCT@scale.data)
tail(genes_use)
counts_3000 <- obj_2@assays$SCT@counts[genes_use, ]
obj_2@assays$SCT@scale.data <- as.matrix(counts_3000)

source("~/Wanmeng/cervical_scRNA_RT/function_seurat_janitor.R") 
obj_2 <- fix_seurat_SCT(obj_2)    
SaveH5Seurat(obj_2, "~/Wanmeng/cervical_scRNA_RT/Cell_oracle/NK_subset_2.h5Seurat", overwrite = TRUE)
Convert("NK_subset_2.h5Seurat", dest = "h5ad", overwrite = TRUE)


######Second file
obj <- NK_subset_for_celloracle
obj$NK_cell_type=as.character(obj$NK_cell_type)

source("~/Wanmeng/cervical_scRNA_RT/function_seurat_janitor.R") 
obj <- fix_seurat_SCT(obj)                                      
SaveH5Seurat(obj, "~/Wanmeng/cervical_scRNA_RT/Cell_oracle/NK_subset.h5Seurat", overwrite = TRUE)
Convert("NK_subset.h5Seurat", dest = "h5ad", overwrite = TRUE)










