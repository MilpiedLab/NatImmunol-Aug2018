# #######################################
# all figures for Fig3, Fig5 and Fig3 Sup
# #######################################
library(dplyr)

# ##################
# Load other scripts
# ##################

# Load the global constants
SCRIPT_DIR = "~/scripts/"
CONSTANT_SOURCE = paste( SCRIPT_DIR, "loading_data_scrnaseq.R", sep="")
source( CONSTANT_SOURCE)

# Load functions
ALL_FUNCTIONS_SOURCE = paste( SCRIPT_DIR, "own_functions.R", sep="")
source( ALL_FUNCTIONS_SOURCE)

# #########################
# Start working on the data
# #########################

# ##########
# ##########
# FigSup3_c
# ##########
# ##########

physio_seurat_object = CreateSeuratObject(raw.data = as.matrix(t(physio_df[,all_physio_genes_v])),
                                                  min.cells = 1, min.genes = 1, project = "all")
physio_seurat_object <- FindVariableGenes(object = physio_seurat_object, mean.function = ExpMean, dispersion.function = LogVMR, 
                                                  x.low.cutoff = 0.1, x.high.cutoff = 6, y.cutoff = 0.5, do.plot=FALSE, display.progress = FALSE)
physio_variable_genes_v = physio_seurat_object@var.genes

# we remove the Ig genes from the variable genes
noIG_physio_variable_genes_v = physio_variable_genes_v[!grepl("^IGH|^IGL|^IGK",physio_variable_genes_v)]

pca_results = pca_function(physio_df[,noIG_physio_variable_genes_v], 
                                PlotHist = FALSE,
                                VectorToColorPCA = physio_df$BCL6, 
                                NameOfVectorToColor = "BCL6",dim_all=FALSE,
                                title=paste0("physio_df cells: ",nrow(physio_df),
                                             "\n noIG_physio_variable_genes_v: ",
                                             length(noIG_physio_variable_genes_v)))

physio_df$PhenoCluster = "GC"
physio_df[pca_results$li$Axis1>0,]$PhenoCluster = "non-GC"

allData_df$PhenoCluster = NA
allData_df[rownames(physio_df),]$PhenoCluster = physio_df$PhenoCluster

pca_results = pca_function(physio_df[,noIG_physio_variable_genes_v], 
                           PlotHist = FALSE,
                           VectorToColorPCA = physio_df$PhenoCluster, 
                           NameOfVectorToColor = "PhenoCluster",dim_all=FALSE,
                           title=paste0("physio_df cells: ",nrow(physio_df),
                                        "\n noIG_physio_variable_genes_v: ",
                                        length(noIG_physio_variable_genes_v)))


# ##########
# ##########
# FigSup3_c
# ##########
# ##########


FL_seurat_object = CreateSeuratObject(raw.data = as.matrix(t(FL_df[,all_FL_genes_v])),
                                       min.cells = 1, min.genes = 1, project = "all")
FL_seurat_object <- FindVariableGenes(object = FL_seurat_object, mean.function = ExpMean, dispersion.function = LogVMR, 
                                       x.low.cutoff = 0.1, x.high.cutoff = 6, y.cutoff = 0.5, do.plot=FALSE, display.progress = FALSE)
FL_variable_genes_v = FL_seurat_object@var.genes

# we remove the Ig genes from the variable genes
noIG_FL_variable_genes_v = FL_variable_genes_v[!grepl("^IGH|^IGL|^IGK",FL_variable_genes_v)]

for(gene in c("CD3D","MS4A1")){
  pca_results = pca_function(FL_df[,noIG_FL_variable_genes_v], 
                             PlotHist = FALSE,
                             VectorToColorPCA = FL_df[,gene], 
                             NameOfVectorToColor = gene,dim_all=FALSE,
                             title=paste0("FL_df cells: ",nrow(FL_df),
                                          "\n noIG_FL_variable_genes_v: ",
                                          length(noIG_physio_variable_genes_v)))
}

PC1_coord = pca_results$li$Axis1
PC2_coord = pca_results$li$Axis2

FL_df$PhenoCluster = "Other Cells"
FL_df[PC1_coord+0.4*PC2_coord>(-1),]$PhenoCluster = "FL5"

allData_df[rownames(FL_df),]$PhenoCluster = FL_df$PhenoCluster

pca_results = pca_function(FL_df[,noIG_FL_variable_genes_v], 
                           PlotHist = FALSE,
                           VectorToColorPCA = FL_df$PhenoCluster, 
                           NameOfVectorToColor = "PhenoCluster",dim_all=FALSE,
                           title=paste0("FL_df cells: ",nrow(FL_df),
                                        "\n noIG_FL_variable_genes_v: ",
                                        length(noIG_FL_variable_genes_v)))


# ##########
# ##########
# Fig5_e
# ##########
# ##########

GC_df = allData_df[allData_df$PhenoCluster %in% c("GC"),]
GC_FL_df = allData_df[allData_df$PhenoCluster %in% c("GC","FL5"),]
FL_df = allData_df[allData_df$PhenoCluster %in% c("FL5"),]


GC_seurat = CreateSeuratObject(raw.data = as.matrix(t(GC_df[,all_genes_v])),
                                      min.cells = 1, min.genes = 1, project = "all")
GC_seurat <- FindVariableGenes(object = GC_seurat, mean.function = ExpMean, dispersion.function = LogVMR, 
                                      x.low.cutoff = 0.1, x.high.cutoff = 6, y.cutoff = 0.5, do.plot=FALSE, display.progress = FALSE)
GC_variable_genes_v = GC_seurat@var.genes

# we remove the Ig genes from the variable genes
noIG_GC_variable_genes_v = GC_variable_genes_v[!grepl("^IGH|^IGL|^IGK",GC_variable_genes_v)]


tsne_coord_GC_FL_df = run_tsne_function(GC_FL_df[,noIG_GC_variable_genes_v])

plot_tsne_function(GC_FL_df[,noIG_GC_variable_genes_v], tsne_coord_GC_FL_df,
                   VectorToColor = GC_FL_df$PhenoCluster, 
                   NameOfVectorToColor = "PhenoCluster",
                   title="GC_FL_df on\n noIG_GC_variable_genes_v ")


pca_results = pca_function(GC_FL_df[,noIG_GC_variable_genes_v], 
                           PlotHist = FALSE,
                           VectorToColorPCA = GC_FL_df$PhenoCluster, 
                           NameOfVectorToColor = "PhenoCluster",dim_all=FALSE,
                           title=paste0("GC_FL_df cells: ",nrow(GC_FL_df),
                                        "\n noIG_GC_variable_genes_v: ",
                                        length(noIG_GC_variable_genes_v)))


# ##########
# ##########
# Fig5_e
# ##########
# ##########


pca_results = pca_function(allData_df[allData_df$PhenoCluster!="Other Cells",noIG_GC_variable_genes_v], 
                           PlotHist = FALSE,
                           VectorToColorPCA = allData_df[allData_df$PhenoCluster!="Other Cells","PhenoCluster"], 
                           NameOfVectorToColor = "PhenoCluster",dim_all=FALSE,
                           title="GC non-GC FL cells")


# ###############################################
# ###############################################
# GC
# ###############################################
# ###############################################


# list of cell cycle markers, from Tirosh et al, 2015
cc.genes <- readLines(con = paste0(DATA_FOLDER,"regev_lab_cell_cycle_genes.txt"))
s.genes <- cc.genes[1:43]
g2m.genes <- cc.genes[44:98]
GC_seurat <- CellCycleScoring(object = GC_seurat, s.genes = s.genes, g2m.genes = g2m.genes, set.ident = FALSE)

GC_df$CellCyclePhase = as.vector(GC_seurat@meta.data$Phase)
GC_df$S.Score = as.vector(GC_seurat@meta.data$S.Score)
GC_df$G2M.Score = as.vector(GC_seurat@meta.data$G2M.Score)

gg_plot = ggplot(GC_df, aes(S.Score,G2M.Score, color= CellCyclePhase))+ geom_point()+
  ggtitle("list of cell cycle markers, from Tirosh et al, 2015\n S.genes=42\n G2M.genes=55")
print(gg_plot)

# ##########
# ##########
# Fig3_a
# ##########
# ##########
pca_results = pca_function(GC_df[,noIG_GC_variable_genes_v], 
                           PlotHist = FALSE,
                           VectorToColorPCA = GC_df$CellCyclePhase, 
                           NameOfVectorToColor = "CellCyclePhase",dim_all=FALSE,
                           title="GC_df")


# ##########
# ##########
# FigSup3_i
# ##########
# ##########
pca_results_GC = pca_function(GC_df[,noIG_GC_variable_genes_v], 
                           PlotHist = FALSE,
                           VectorToColorPCA = GC_df$CellCyclePhase, 
                           NameOfVectorToColor = "CellCyclePhase",dim_all=FALSE,
                           dim_x="PC2", dim_y="PC3",
                           title="GC_df")

GC_df$PC2 = pca_results_GC$li$Axis2


# ############################################################################
# DZ/LZ: scoring the GC cells with the list of genes from Victora Blood 2012
# ############################################################################

victora_genes_df = read.table(paste0(DATA_FOLDER,"victora_blood_2012.csv"), 
                              stringsAsFactors = FALSE, header = TRUE, sep=",")
victora_genes_df = victora_genes_df[victora_genes_df$GeneSymbol!="---",]

victora_genes_df = victora_genes_df[!duplicated(victora_genes_df$GeneSymbol),]
victora_genes_LZ_v = victora_genes_df$GeneSymbol[victora_genes_df$EnrichedIn=="LZ"]
victora_genes_DZ_v = victora_genes_df$GeneSymbol[victora_genes_df$EnrichedIn=="DZ"]

victora_genes_LZ_v = victora_genes_LZ_v[victora_genes_LZ_v %in% all_genes_v]
victora_genes_DZ_v = victora_genes_DZ_v[victora_genes_DZ_v %in% all_genes_v]
#victora_genes_DZ_v_withoutCycleGenes = victora_genes_DZ_v[!victora_genes_DZ_v %in% c(cc.genes)]

GC_seurat <- CellCycleScoring(object = GC_seurat, s.genes = victora_genes_LZ_v, g2m.genes = victora_genes_DZ_v, set.ident = FALSE)
GC_df$LZDZPhase = as.vector(GC_seurat@meta.data$Phase)
GC_df$LZDZPhase[GC_df$LZDZPhase=="S"]="LZ"
GC_df$LZDZPhase[GC_df$LZDZPhase=="G2M"]="DZ"
GC_df$LZ.Score = as.vector(GC_seurat@meta.data$S.Score)
GC_df$DZ.Score = as.vector(GC_seurat@meta.data$G2M.Score)

# ##########
# ##########
# FigSup3_j
# ##########
# ##########
gg_plot = ggplot(GC_df, aes(LZ.Score,DZ.Score, color= CellCyclePhase))+ geom_point()+
  ggtitle("list of DZLZ markers")
print(gg_plot)


# ############################################################################
# scoring the GC cells with the list of genes from Victora Blood 2012 without cell cycle genes
# ############################################################################
GO_CellCycle_df = read.csv(paste0(DATA_FOLDER,"goa_human_merged_cellcycle_onlyGeneNames.gaf"),
                            stringsAsFactors = FALSE, header = FALSE, sep=",",check.names=FALSE)
GO_CellCycle_v = as.vector(GO_CellCycle_df[,1])
GO_CellCycle_v = unique(GO_CellCycle_v)

victora_genes_LZ_v_withoutCycleGenes = victora_genes_LZ_v[!victora_genes_LZ_v %in% c(cc.genes,GO_CellCycle_v)]
victora_genes_DZ_v_withoutCycleGenes = victora_genes_DZ_v[!victora_genes_DZ_v %in% c(cc.genes,GO_CellCycle_v)]

GC_seurat <- AddModuleScore(GC_seurat, genes.list = list(victora_genes_LZ_v_withoutCycleGenes),  enrich.name = "LZ.Score")
GC_df$LZ.Score_noCellCycleGenes = as.vector(GC_seurat@meta.data$LZ.Score)
GC_seurat <- AddModuleScore(GC_seurat, genes.list = list(victora_genes_DZ_v_withoutCycleGenes),  enrich.name = "DZ.Score")
GC_df$DZ.Score_noCellCycleGenes = as.vector(GC_seurat@meta.data$DZ.Score)
GC_df$LZDZPhase_noCellCycleGenes = sapply(1:nrow(GC_df), function(x)if(GC_df[x,"LZ.Score_noCellCycleGenes"]>GC_df[x,"DZ.Score_noCellCycleGenes"]){"LZ"}else{"DZ"})



# ##########
# ##########
# Fig3_c
# ##########
# ##########

reg_lin = lm(GC_df[,"DZ.Score_noCellCycleGenes"] ~ GC_df[,"LZ.Score_noCellCycleGenes"])
geom_abline(aes(intercept=reg_lin$coefficients[1], slope=reg_lin$coefficients[2]), colour="blue")

seq_v = seq(from = min(GC_df$LZ.Score_noCellCycleGenes)-0.18, to = max(GC_df$LZ.Score_noCellCycleGenes)+0.05, by = 0.0001)
lines_v = reg_lin$coefficients[2]*seq_v + reg_lin$coefficients[1]
line_df = data.frame(LZ.Score_noCellCycleGenes=seq_v,
                     DZ.Score_noCellCycleGenes=lines_v, 
                     Index=1:length(seq_v))


DZLZ.Score_v = c()
for(i in 1:nrow(GC_df)){
  eucl_dist_v = sqrt( (GC_df[i,"LZ.Score_noCellCycleGenes"]-line_df[,"LZ.Score_noCellCycleGenes"])^2 +
                        (GC_df[i,"DZ.Score_noCellCycleGenes"]-line_df[,"DZ.Score_noCellCycleGenes"])^2 )
  DZLZ.Score_v = c(DZLZ.Score_v, which.min(eucl_dist_v) )
}

GC_df$DZLZ.Score = (DZLZ.Score_v-min(DZLZ.Score_v))/max(DZLZ.Score_v)
gg_plot = ggplot(GC_df, aes(LZ.Score_noCellCycleGenes,DZ.Score_noCellCycleGenes))+
  geom_point(aes( color=CellCyclePhase))+
  geom_line(data=line_df, aes(LZ.Score_noCellCycleGenes,DZ.Score_noCellCycleGenes))+
  ggtitle("list of DZLZ markers, from Victora et al, 2012\n LZ.genes_noCellCycleGenes=180\n DZ.genes_noCellCycleGenes=110")
print(gg_plot)



# ##########
# ##########
# Fig3_d
# ##########
# ##########

G1_spearman_cor = signif( cor( GC_df[GC_df$CellCyclePhase=="G1",]$PC2, 
                               GC_df[GC_df$CellCyclePhase=="G1",]$DZLZ.Score, 
                               method="spearman", use="complete"), digits=3)
S_G2M_spearman_cor = signif( cor( GC_df[GC_df$CellCyclePhase %in% c("G2M","S"),]$PC2, 
                                  GC_df[GC_df$CellCyclePhase %in% c("G2M","S"),]$DZLZ.Score, 
                                  method="spearman", use="complete"), digits=3)
gg_plot = ggplot(GC_df, aes(DZLZ.Score,PC2, color= CellCyclePhase))+ geom_point()+
  ggtitle(paste("G1_spearman_cor:",G1_spearman_cor,"\n S_G2M_spearman_cor",S_G2M_spearman_cor))
print(gg_plot)

# ##########
# ##########
# Fig3_e
# ##########
# ##########


startGreyZoneScore = 0.42
endGreyZoneScore = 0.56
GC_df$LZDZ.score.cutoff2 = sapply(GC_df$DZLZ.Score, 
                                             function(x) if(x<startGreyZoneScore){"DZ"}else if(x<endGreyZoneScore){"GreyZone"}else{"LZ"})
#GC_df$LZDZ.score.cutoff2


gg_plot = ggplot(GC_df, aes(LZ.Score_noCellCycleGenes, DZ.Score_noCellCycleGenes, color=LZDZ.score.cutoff2))+geom_point()+
  ggtitle("list of DZLZ markers, from Victora et al, 2012\n LZ.genes_noCellCycleGenes=180\n DZ.genes_noCellCycleGenes=110")
print(gg_plot)


# ########################################
# Find DE between the 3 zones (DZ/grey/LZ)
# ########################################

tpm_v = as.factor(GC_df$LZDZ.score.cutoff2)
names(tpm_v) = rownames(GC_df)
GC_seurat@ident = tpm_v
GC_seurat.markers <- FindAllMarkers(object = GC_seurat, only.pos = TRUE, min.pct = 0.25, 
                                    thresh.use = 0.25)
GC_seurat.markers$diff.pct = GC_seurat.markers$pct.1 - GC_seurat.markers$pct.2
top20_DiffPerc_LZDZ <- GC_seurat.markers %>% group_by(cluster) %>% top_n(20, diff.pct) # library(dplyr)


# ##########
# ##########
# Fig3_f
# ##########
# ##########
order_GC_df = GC_df[order(GC_df$DZLZ.Score),]
heatmap.3(order_GC_df,dendrogram = "none",Rowv = FALSE, Colv = FALSE,
          genes_v=top20_DiffPerc_LZDZ$gene,
          col=c("grey50",TWO_COLORS_GRADIANT[-1]),
          margins = c(5, 10),
          lhei=c(1.2,4,0.5),lwid=c(1.2,4),
          ColSide_name_v=c("LZDZ.score.cutoff2","CellCyclePhase"),
          ColSideColorsSize = 2.5,trace="none",
          main="order_GC_df\n top20 based on PerDiff")
legend_for_heatmap.3(order_GC_df,ColSide_name_v=c("LZDZ.score.cutoff2","CellCyclePhase"))

# ##########
# ##########
# FigSup3_k
# ##########
# ##########

order_GC_df = GC_df[order(GC_df$DZLZ.Score),]
ordered_GENES_VICTORA = c(victora_genes_DZ_v[order(colSums(GC_df[,victora_genes_DZ_v]), decreasing=TRUE)],
                          victora_genes_LZ_v[order(colSums(GC_df[,victora_genes_LZ_v]))])
heatmap.3(order_GC_df,dendrogram = "none",Rowv = FALSE, Colv = FALSE,
          genes_v=ordered_GENES_VICTORA,
          col=c("grey50",TWO_COLORS_GRADIANT[-1]),
          margins = c(5, 14),
          lhei=c(1.2,4,0.5),lwid=c(1.2,4),
          ColSide_name_v=c("LZDZ.score.cutoff2","CellCyclePhase"),
          ColSideColorsSize = 2.5,trace="none",
          main=" genes are ordered based on their sum of expression\n 
          Victora genes with cell cycle genes")
legend_for_heatmap.3(order_GC_df,ColSide_name_v=c("LZDZ.score.cutoff2","CellCyclePhase"))


# ##########
# ##########
# Fig3_g
# ##########
# ##########

# The Fig3_g is in the script Fig4_scrnaseq.R
