# ###################################
# all figures for Fig7 (scqpcr part) 
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


number_K=5
set.seed(42) 
Kmeans_cluster=kmeans(raw_data_FL[,gene_names_v],number_K, nstart=20)

raw_data_FL$Kmeans=paste("clusterK",Kmeans_cluster$cluster,sep="")
cell_attribute_toSplit="Kmeans"



# ##########
# ##########
# Fig7_a
# ##########
# ##########


tsne_coord_df = run_tsne_function(raw_data_FL[,gene_names_v])

plot_tsne_function(raw_data_FL[,gene_names_v], tsne_coord_df,
                   VectorToColor = raw_data_FL$Kmeans, 
                   NameOfVectorToColor = "Kmeans",
                   title="raw_data_FL")



# ##########################
# volcano plots Combined LRT
# ##########################

results_list=volcano_plot_LRT_function(raw_data_FL[,gene_names_v],raw_data_FL$Kmeans)


cat("\n# ################################################### #\n")
cat("\n# ################################################### #\n")
cat("\n# ################################################### #\n")
print(paste("number of cell in the clusterK1 group: ",nrow(raw_data_FL[which(raw_data_FL$Kmeans %in% c("clusterK1")),]) ))
print(paste("number of cell in the clusterK2 group: ",nrow(raw_data_FL[which(raw_data_FL$Kmeans %in% c("clusterK2")),]) ))
print(paste("number of cell in the clusterK3 group: ",nrow(raw_data_FL[which(raw_data_FL$Kmeans %in% c("clusterK3")),]) ))
print(paste("number of cell in the clusterK4 group: ",nrow(raw_data_FL[which(raw_data_FL$Kmeans %in% c("clusterK4")),]) ))
print(paste("number of cell in the clusterK5 group: ",nrow(raw_data_FL[which(raw_data_FL$Kmeans %in% c("clusterK5")),]) ))
cat("\n# ################################################### #\n")
cat("\n# ################################################### #\n")
cat("\n# ################################################### #\n")


# besoin de trier les genes en fonction de si ils sont positifs ou negatifs dans le volcano plot 
# car il faut prendre seulement les genes toujours positifs par rapport aux autres grps 
# ou les genes toujours negatifs par rapport au autres grps

# for clusterK1:
# genes which are in "clusterK1 vs clusterK2" and "clusterK1 vs clusterK3"
#genes_specific_clusterK1=results_list[["clusterK1 vs clusterK2"]][which (results_list[["clusterK1 vs clusterK2"]] %in% results_list[["clusterK1 vs clusterK3"]])]

# clusterK1 vs clusterK2
genes_specific_clusterK1_pos1 = results_list[["clusterK1 vs clusterK2"]][["clusterK1"]]
genes_specific_clusterK1_neg1 = results_list[["clusterK1 vs clusterK2"]][["clusterK2"]]
genes_specific_clusterK2_neg1 = genes_specific_clusterK1_pos1
genes_specific_clusterK2_pos1 = genes_specific_clusterK1_neg1
# clusterK1 vs clusterK3
genes_specific_clusterK1_pos2 = results_list[["clusterK1 vs clusterK3"]][["clusterK1"]]
genes_specific_clusterK1_neg2 = results_list[["clusterK1 vs clusterK3"]][["clusterK3"]]
genes_specific_clusterK3_neg1 = genes_specific_clusterK1_pos2
genes_specific_clusterK3_pos1 = genes_specific_clusterK1_neg2
# clusterK1 vs clusterK4
genes_specific_clusterK1_pos3 = results_list[["clusterK1 vs clusterK4"]][["clusterK1"]]
genes_specific_clusterK1_neg3 = results_list[["clusterK1 vs clusterK4"]][["clusterK4"]]
genes_specific_clusterK4_neg1 = genes_specific_clusterK1_pos3
genes_specific_clusterK4_pos1 = genes_specific_clusterK1_neg3
# clusterK1 vs clusterK5
genes_specific_clusterK1_pos4 = results_list[["clusterK1 vs clusterK5"]][["clusterK1"]]
genes_specific_clusterK1_neg4 = results_list[["clusterK1 vs clusterK5"]][["clusterK5"]]
genes_specific_clusterK5_neg1 = genes_specific_clusterK1_pos4
genes_specific_clusterK5_pos1 = genes_specific_clusterK1_neg4

# clusterK2 vs clusterK3
genes_specific_clusterK2_pos2 = results_list[["clusterK2 vs clusterK3"]][["clusterK2"]]
genes_specific_clusterK2_neg2 = results_list[["clusterK2 vs clusterK3"]][["clusterK3"]]
genes_specific_clusterK3_neg2 = genes_specific_clusterK2_pos2
genes_specific_clusterK3_pos2 = genes_specific_clusterK2_neg2
# clusterK2 vs clusterK4
genes_specific_clusterK2_pos3 = results_list[["clusterK2 vs clusterK4"]][["clusterK2"]]
genes_specific_clusterK2_neg3 = results_list[["clusterK2 vs clusterK4"]][["clusterK4"]]
genes_specific_clusterK4_neg2 = genes_specific_clusterK2_pos3
genes_specific_clusterK4_pos2 = genes_specific_clusterK2_neg3
# clusterK2 vs clusterK4
genes_specific_clusterK2_pos4 = results_list[["clusterK2 vs clusterK5"]][["clusterK2"]]
genes_specific_clusterK2_neg4 = results_list[["clusterK2 vs clusterK5"]][["clusterK5"]]
genes_specific_clusterK5_neg2 = genes_specific_clusterK2_pos4
genes_specific_clusterK5_pos2 = genes_specific_clusterK2_neg4

# clusterK3 vs clusterK4
genes_specific_clusterK3_pos3 = results_list[["clusterK3 vs clusterK4"]][["clusterK3"]]
genes_specific_clusterK3_neg3 = results_list[["clusterK3 vs clusterK4"]][["clusterK4"]]
genes_specific_clusterK4_neg3 = genes_specific_clusterK3_pos3
genes_specific_clusterK4_pos3 = genes_specific_clusterK3_neg3
# clusterK3 vs clusterK4
genes_specific_clusterK3_pos4 = results_list[["clusterK3 vs clusterK5"]][["clusterK3"]]
genes_specific_clusterK3_neg4 = results_list[["clusterK3 vs clusterK5"]][["clusterK5"]]
genes_specific_clusterK5_neg3 = genes_specific_clusterK3_pos4
genes_specific_clusterK5_pos3 = genes_specific_clusterK3_neg4

# clusterK4 vs clusterK5
genes_specific_clusterK4_pos4 = results_list[["clusterK4 vs clusterK5"]][["clusterK4"]]
genes_specific_clusterK4_neg4 = results_list[["clusterK4 vs clusterK5"]][["clusterK5"]]
genes_specific_clusterK5_neg4 = genes_specific_clusterK4_pos4
genes_specific_clusterK5_pos4 = genes_specific_clusterK4_neg4


genes_specific_clusterK1_pos=genes_specific_clusterK1_pos1[which(genes_specific_clusterK1_pos1 %in% genes_specific_clusterK1_pos2)]
genes_specific_clusterK1_pos=genes_specific_clusterK1_pos[which(genes_specific_clusterK1_pos %in% genes_specific_clusterK1_pos3)]
genes_specific_clusterK1_pos=genes_specific_clusterK1_pos[which(genes_specific_clusterK1_pos %in% genes_specific_clusterK1_pos4)]
genes_specific_clusterK1_neg=genes_specific_clusterK1_neg1[which(genes_specific_clusterK1_neg1 %in% genes_specific_clusterK1_neg2)]
genes_specific_clusterK1_neg=genes_specific_clusterK1_neg[which(genes_specific_clusterK1_neg %in% genes_specific_clusterK1_neg3)]
genes_specific_clusterK1_neg=genes_specific_clusterK1_neg[which(genes_specific_clusterK1_neg %in% genes_specific_clusterK1_neg4)]

genes_specific_clusterK2_pos=genes_specific_clusterK2_pos1[which(genes_specific_clusterK2_pos1 %in% genes_specific_clusterK2_pos2)]
genes_specific_clusterK2_pos=genes_specific_clusterK2_pos[which(genes_specific_clusterK2_pos %in% genes_specific_clusterK2_pos3)]
genes_specific_clusterK2_pos=genes_specific_clusterK2_pos[which(genes_specific_clusterK2_pos %in% genes_specific_clusterK2_pos4)]
genes_specific_clusterK2_neg=genes_specific_clusterK2_neg1[which(genes_specific_clusterK2_neg1 %in% genes_specific_clusterK2_neg2)]
genes_specific_clusterK2_neg=genes_specific_clusterK2_neg[which(genes_specific_clusterK2_neg %in% genes_specific_clusterK2_neg3)]
genes_specific_clusterK2_neg=genes_specific_clusterK2_neg[which(genes_specific_clusterK2_neg %in% genes_specific_clusterK2_neg4)]

genes_specific_clusterK3_pos=genes_specific_clusterK3_pos1[which(genes_specific_clusterK3_pos1 %in% genes_specific_clusterK3_pos2)]
genes_specific_clusterK3_pos=genes_specific_clusterK3_pos[which(genes_specific_clusterK3_pos %in% genes_specific_clusterK3_pos3)]
genes_specific_clusterK3_pos=genes_specific_clusterK3_pos[which(genes_specific_clusterK3_pos %in% genes_specific_clusterK3_pos4)]
genes_specific_clusterK3_neg=genes_specific_clusterK3_neg1[which(genes_specific_clusterK3_neg1 %in% genes_specific_clusterK3_neg2)]
genes_specific_clusterK3_neg=genes_specific_clusterK3_neg[which(genes_specific_clusterK3_neg %in% genes_specific_clusterK3_neg3)]
genes_specific_clusterK3_neg=genes_specific_clusterK3_neg[which(genes_specific_clusterK3_neg %in% genes_specific_clusterK3_neg4)]

genes_specific_clusterK4_pos=genes_specific_clusterK4_pos1[which(genes_specific_clusterK4_pos1 %in% genes_specific_clusterK4_pos2)]
genes_specific_clusterK4_pos=genes_specific_clusterK4_pos[which(genes_specific_clusterK4_pos %in% genes_specific_clusterK4_pos3)]
genes_specific_clusterK4_pos=genes_specific_clusterK4_pos[which(genes_specific_clusterK4_pos %in% genes_specific_clusterK4_pos4)]
genes_specific_clusterK4_neg=genes_specific_clusterK4_neg1[which(genes_specific_clusterK4_neg1 %in% genes_specific_clusterK4_neg2)]
genes_specific_clusterK4_neg=genes_specific_clusterK4_neg[which(genes_specific_clusterK4_neg %in% genes_specific_clusterK4_neg3)]
genes_specific_clusterK4_neg=genes_specific_clusterK4_neg[which(genes_specific_clusterK4_neg %in% genes_specific_clusterK4_neg4)]

genes_specific_clusterK5_pos=genes_specific_clusterK5_pos1[which(genes_specific_clusterK5_pos1 %in% genes_specific_clusterK5_pos2)]
genes_specific_clusterK5_pos=genes_specific_clusterK5_pos[which(genes_specific_clusterK5_pos %in% genes_specific_clusterK5_pos3)]
genes_specific_clusterK5_pos=genes_specific_clusterK5_pos[which(genes_specific_clusterK5_pos %in% genes_specific_clusterK5_pos4)]
genes_specific_clusterK5_neg=genes_specific_clusterK5_neg1[which(genes_specific_clusterK5_neg1 %in% genes_specific_clusterK5_neg2)]
genes_specific_clusterK5_neg=genes_specific_clusterK5_neg[which(genes_specific_clusterK5_neg %in% genes_specific_clusterK5_neg3)]
genes_specific_clusterK5_neg=genes_specific_clusterK5_neg[which(genes_specific_clusterK5_neg %in% genes_specific_clusterK5_neg4)]

genes_specific_clusterK1=c(genes_specific_clusterK1_pos,genes_specific_clusterK1_neg)
genes_specific_clusterK2=c(genes_specific_clusterK2_pos,genes_specific_clusterK2_neg)
genes_specific_clusterK3=c(genes_specific_clusterK3_pos,genes_specific_clusterK3_neg)
genes_specific_clusterK4=c(genes_specific_clusterK4_pos,genes_specific_clusterK4_neg)
genes_specific_clusterK5=c(genes_specific_clusterK5_pos,genes_specific_clusterK5_neg)


cat("\nDifferentially Expressed Genes specific to clusterK1 Cells (OVEREXPRESSED GENES):\n")
print(genes_specific_clusterK1_pos, quote=FALSE)
cat("Differentially Expressed Genes specific to clusterK1 Cells (UNDER EXPRESSED GENES):\n")
print(genes_specific_clusterK1_neg, quote=FALSE)

cat("\nDifferentially Expressed Genes specific to clusterK2 Cells (OVEREXPRESSED GENES):\n")
print(genes_specific_clusterK2_pos, quote=FALSE)
cat("Differentially Expressed Genes specific to clusterK2 Cells (UNDER EXPRESSED GENES):\n")
print(genes_specific_clusterK2_neg, quote=FALSE)

cat("\nDifferentially Expressed Genes specific to clusterK3 Cells (OVEREXPRESSED GENES):\n")
print(genes_specific_clusterK3_pos, quote=FALSE)
cat("Differentially Expressed Genes specific to clusterK3 Cells (UNDER EXPRESSED GENES):\n")
print(genes_specific_clusterK3_neg, quote=FALSE)

cat("\nDifferentially Expressed Genes specific to clusterK4 Cells (OVEREXPRESSED GENES):\n")
print(genes_specific_clusterK4_pos, quote=FALSE)
cat("Differentially Expressed Genes specific to clusterK4 Cells (UNDER EXPRESSED GENES):\n")
print(genes_specific_clusterK4_neg, quote=FALSE)

cat("\nDifferentially Expressed Genes specific to clusterK5 Cells (OVEREXPRESSED GENES):\n")
print(genes_specific_clusterK5_pos, quote=FALSE)
cat("Differentially Expressed Genes specific to clusterK5 Cells (UNDER EXPRESSED GENES):\n")
print(genes_specific_clusterK5_neg, quote=FALSE)

differentially_expr_genes_list=list()
differentially_expr_genes_list[["clusterK1"]]=genes_specific_clusterK1
differentially_expr_genes_list[["clusterK2"]]=genes_specific_clusterK2
differentially_expr_genes_list[["clusterK3"]]=genes_specific_clusterK3
differentially_expr_genes_list[["clusterK4"]]=genes_specific_clusterK4
differentially_expr_genes_list[["clusterK5"]]=genes_specific_clusterK5

differentially_expr_genes_vector = union(genes_specific_clusterK1, genes_specific_clusterK2)
differentially_expr_genes_vector = union(differentially_expr_genes_vector, genes_specific_clusterK3)
differentially_expr_genes_vector = union(differentially_expr_genes_vector, genes_specific_clusterK4)
differentially_expr_genes_vector = union(differentially_expr_genes_vector, genes_specific_clusterK5)


# ##########
# ##########
# Fig7_b
# ##########
# ##########

for(gene in differentially_expr_genes_vector){
  
  plot_tsne_function(raw_data_FL[,gene_names_v], tsne_coord_df,
                     VectorToColor = raw_data_FL[,gene], 
                     NameOfVectorToColor = gene,
                     title="raw_data_FL")
}

# ##########
# ##########
# Fig7_c (not the same order)
# ##########
# ##########

differentially_expr_genes_vector_original = differentially_expr_genes_vector
# we reorder the vector to have the same order of genes as in the paper
differentially_expr_genes_vector=c("CD79B","BACH2","CD38","GNA13","DNMT1","CD86","EZH2","GPR183","CD83","CXCR4","BCL2A1","PRDM1","SEMA7A")



perc_expr_Kcluster_data_frame = data.frame(GeneNames = differentially_expr_genes_vector)
for(clusterK in sort(unique(raw_data_FL$Kmeans))){
  vector_clusterK = c()
  raw_data_for_analysis_K=raw_data_FL[which(raw_data_FL$Kmeans==clusterK),]
  for(gene in differentially_expr_genes_vector){
    
    vector_clusterK = c(vector_clusterK,length(which(raw_data_for_analysis_K[,gene]>0))/length(raw_data_for_analysis_K[,gene])*100)
  }
  perc_expr_Kcluster_data_frame[[clusterK]] = vector_clusterK
}
rownames(perc_expr_Kcluster_data_frame) = perc_expr_Kcluster_data_frame$GeneNames
perc_expr_Kcluster_data_frame$GeneNames<-NULL

heatmap_plot=heatmap.2(as.matrix(perc_expr_Kcluster_data_frame), col=THREE_COLORS_HEATMAP, 
                       tracecol=NA, Colv=FALSE, Rowv=FALSE, dendrogram="none",
                       main="",
                       cexRow=0.55, cexCol=0.55,
                       notecex=0.55,notecol="black")



# ##########
# ##########
# Fig7_d
# ##########
# ##########

plot_tsne_function(raw_data_FL[,gene_names_v], tsne_coord_df,
                   VectorToColor = raw_data_FL$Sample, 
                   NameOfVectorToColor = "Sample",
                   VectorFacetWrap = raw_data_FL$Sample, 
                   title="raw_data_FL")


# ##########
# ##########
# Fig7_e
# ##########
# ##########

for(sample in sort(unique(raw_data_FL$Sample))){
  sample_df = raw_data_FL[raw_data_FL$Sample==sample,]
  sample_df = sample_df[order(sample_df$Kmeans),]
  
  heatmap.3(sample_df,
            genes_v=differentially_expr_genes_vector,
            Rowv=FALSE, 
            Colv=FALSE,
            dendrogram="none",
            col=TWO_COLORS_GRADIANT,
            margins = c(5, 14),
            lhei=c(1.2,4,0.5),lwid=c(1.2,4),
            ColSide_name_v=c("Kmeans"),
            ColSideColorsSize = 2.5,trace="none",
            titleColorGradient='Et',
            main=paste(sample))
  legend_for_heatmap.3(sample_df,ColSide_name_v=c("Kmeans"))
}





