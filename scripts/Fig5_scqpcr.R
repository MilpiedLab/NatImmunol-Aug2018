# ###################################
# all figures for Fig5 (scqpcr part) 
# ###################################

# ##################
# Load other scripts
# ##################

# Load the global constants
SCRIPT_DIR = "~/workspace/ciml-scqpcr/articles/article1/public_scripts/scripts/"
CONSTANT_SOURCE = paste( SCRIPT_DIR, "loading_data_scqpcr.R", sep="")
source( CONSTANT_SOURCE)

# Load functions
ALL_FUNCTIONS_SOURCE = paste( SCRIPT_DIR, "own_functions.R", sep="")
source( ALL_FUNCTIONS_SOURCE)

# #########################
# Start working on the data
# #########################


# ##########
# ##########
# Fig5_a
# ##########
# ##########

hc <- hclust(as.dist(1-cor(t(raw_data_bulk10[,gene_names_v]), method="spearman")), method="complete")
hr_eucl <- hclust(dist(t(raw_data_bulk10[,gene_names_v])), method="complete")


heatmap.3(raw_data_bulk10,
          genes_v=gene_names_v,
          Rowv=as.dendrogram(hr_eucl), 
          Colv=as.dendrogram(hc),
          col=TWO_COLORS_GRADIANT,
          margins = c(5, 14),
          lhei=c(1.2,4,0.5),lwid=c(1.2,4),
          ColSide_name_v=c("SimpleSortPheno",'Sample'),
          ColSideColorsSize = 2.5,trace="none",
          titleColorGradient='Et',
          main='raw_data')
legend_for_heatmap.3(raw_data_bulk10,ColSide_name_v=c("SimpleSortPheno",'Sample'))



# ##########
# ##########
# Fig5_b
# ##########
# ##########

# This plot was not done in R.


# ##########
# ##########
# Fig5_c
# ##########
# ##########


tsne_coord_df = run_tsne_function(raw_data[,gene_names_v])

plot_tsne_function(raw_data[,gene_names_v], tsne_coord_df,
                   VectorToColor = raw_data$SimpleSortPheno, 
                   NameOfVectorToColor = "SimpleSortPheno",
                   title="raw_data")


# ##########
# ##########
# Fig5_d
# ##########
# ##########

for(gene in c('LMO2','KMT2D','CD74', 'BCL2', 'GNAI2', 'HLA.DOA')){
  
  plot_tsne_function(raw_data[,gene_names_v], tsne_coord_df,
                     VectorToColor = raw_data[,gene], 
                     NameOfVectorToColor = gene,
                     title="raw_data")
}
