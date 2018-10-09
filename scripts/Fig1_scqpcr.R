# #######################################
# all figures for Fig1 and Fig1 Sup
# #######################################


# ##################
# Load other scripts
# ##################

# Load the global constants
SCRIPT_DIR = "~/scripts/"
CONSTANT_SOURCE = paste( SCRIPT_DIR, "loading_data_scqpcr.R", sep="")
source( CONSTANT_SOURCE)

# Load functions
ALL_FUNCTIONS_SOURCE = paste( SCRIPT_DIR, "own_functions.R", sep="")
source( ALL_FUNCTIONS_SOURCE)

# #########################
# Start working on the data
# #########################


# ##########
# Fig1_a
# ##########

tsne_coord_df = run_tsne_function(raw_data_GC_PB_Mem[,gene_names_v])

plot_tsne_function(raw_data_GC_PB_Mem[,gene_names_v], tsne_coord_df,
                   VectorToColor = raw_data_GC_PB_Mem$IndexedPheno, 
                   NameOfVectorToColor = "IndexedPheno",
                   title="raw_data_GC_PB_Mem")

# ##########
# Fig1_b
# ##########

genes_v = c('NOTCH2',"PRDM1","BCL6","JAM3","SDC1","AICDA")
for(gene in genes_v){
  plot_tsne_function(raw_data_GC_PB_Mem[,gene_names_v], tsne_coord_df,
                     VectorToColor = raw_data_GC_PB_Mem[,gene], 
                     NameOfVectorToColor = gene,
                     title="raw_data_GC_PB_Mem")
}


# ##########
# FigSup1_f
# ##########

hc <- hclust(as.dist(1-cor(t(raw_data_GC_PB_Mem[,gene_names_v]), method="spearman")), method="complete")
hr_eucl <- hclust(dist(t(raw_data_GC_PB_Mem[,gene_names_v])), method="complete")

heatmap.3(raw_data_GC_PB_Mem,
          genes_v=gene_names_v,
          Rowv=as.dendrogram(hr_eucl), 
          Colv=as.dendrogram(hc),
          col=TWO_COLORS_GRADIANT,
          margins = c(5, 14),
          lhei=c(1.2,4,0.5),lwid=c(1.2,4),
          ColSide_name_v=c("IndexedPheno",'Sample'),
          ColSideColorsSize = 2.5,trace="none",
          titleColorGradient='Et',
          main='raw_data_GC_PB_Mem')
legend_for_heatmap.3(raw_data_GC_PB_Mem,ColSide_name_v=c("IndexedPheno",'Sample'))

# ##########
# FigSup1_g
# ##########

pca_result = pca_function(raw_data_GC_PB_Mem[,gene_names_v], PlotHist=FALSE, 
                               VectorToColorPCA=raw_data_GC_PB_Mem$IndexedPheno, 
                               NameOfVectorToColor="IndexedPheno",dim_all=FALSE,
                               title='raw_data_GC_PB_Mem')

# ##########
# FigSup1_h
# ##########

gg_plots = pca_gene_information_function(pca_result, dim_all=FALSE, threshold_norm=0.6,
                                         title='raw_data_GC_PB_Mem')


# ##########
# Fig1_c
# ##########

pca = pca_function(raw_data_PB[,gene_names_v], PlotHist=TRUE , bool_plot=FALSE)

# ##########
# Fig1_d
# ##########

# get the PC1 coordinates from the PCA
raw_data_PB$PC1=pca$li$Axis1
# order the data based on coordinate PC1
raw_data_PB=raw_data_PB[order(raw_data_PB$PC1),]
# take the top gene which drive the PC1
threshold_norm = 0.75
current_vector = pca$c1[[1]]
decreasing_order = order( abs(current_vector), decreasing = TRUE)
genes_vector = c()
for( gene_index in seq(1, length( current_vector))){
  current_norm = sqrt (sum( (pca$c1[[ 1]][decreasing_order][1: gene_index])^2))
  genes_vector=c(genes_vector,row.names( pca$c1)[decreasing_order][gene_index])
  if( current_norm >= threshold_norm){
    break
  }
}
genes_to_study = c("PRDM1", "IRF4",genes_vector) # two important genes
# reordering the genes (for the heatmap)
genes_to_study = c("PRDM1","IRF4","BCL6","MS4A1","KMT2C","POLH","PAX5","HLA.DRA","P2RY8","GNA13","SDC1")

heatmap.3(raw_data_PB,
          genes_v=genes_to_study,
          dendrogram = 'none',
          Rowv=FALSE, 
          Colv=FALSE,
          col=TWO_COLORS_GRADIANT,
          margins = c(5, 14),
          lhei=c(1.2,4,0.5),lwid=c(1.2,4),
          ColSide_name_v=c("IndexedPheno"),
          ColSideColorsSize = 2.5,trace="none",
          titleColorGradient='Et',
          main='raw_data_PB')
legend_for_heatmap.3(raw_data_PB,ColSide_name_v=c("IndexedPheno"))

# ##########
# Fig1_e
# ##########
# monocle DDRTree (plots to view what happens)
results_DDRTree_PB = pseudotimeDDRTree_function(raw_data_PB, GENES_TO_ANALYSE=gene_names_v,
                                                cell_attribute_colorBy=c("Pseudotime","IndexedPheno"),
                                                DDRTree_max_components=3)

raw_data_PB$DDRTreePseudotime = as.vector(results_DDRTree_PB$DDRTree_pseudotime)

DDRTree_pseudotime = raw_data_PB$DDRTreePseudotime
names(DDRTree_pseudotime)=raw_data_PB$UniqueCellID

PC1_pca_PB = raw_data_PB$PC1
names(PC1_pca_PB)=raw_data_PB$UniqueCellID

# Create the article figure 1 e
# ############################
gg_plot = correlation_pseudotime_function(PC1_pca_PB, DDRTree_pseudotime, "PC1 PCA PB","Pseudotime DDRTree PB" )
print(gg_plot)
