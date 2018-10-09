# #############################################################
# all scqpcr figures for Fig4 and Fig6  and FigSup4 and FigSup6
# #############################################################

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
# ##########
# Fig4_a
# ##########
# ##########


# Calcul the derive of the kernel smoother curve
# ##############################################
x=kernel_smoother_coord_dataframe$PseudotimeTheta
derived_kernel_smoother_coord_df = data.frame(Pseudotime = kernel_smoother_coord_dataframe$PseudotimeTheta[1:(nrow(kernel_smoother_coord_dataframe)-1)])
for( gene in gene_names_v){
  y=kernel_smoother_coord_dataframe[,gene]
  derivedGene_vector=c()
  for(i in 1:(nrow(kernel_smoother_coord_dataframe)-1)){
    
    derivedGene = (y[i+1]-y[i])/(x[i+1]-x[i])
    
    derivedGene_vector=c(derivedGene_vector,derivedGene) 
  }
  derived_kernel_smoother_coord_df[,gene] = derivedGene_vector
}


# list of genes (selected with a threshold on the percentage of expression)
# #########################################################################
percentage_expression_vector = c()
for(gene in gene_names_v){
  percentage_expression_vector = c(percentage_expression_vector,length(which(!raw_data_GC[,gene]==0))/nrow(raw_data_GC)*100)
}
percentage_expression_df=data.frame(Gene=gene_names_v, PercExpr=percentage_expression_vector)
percentage_expression_df= percentage_expression_df[order(percentage_expression_df$PercExpr, decreasing = TRUE),]
percentage_expression_df$NumberOfGenes=1:length(gene_names_v)

list_genes=as.vector(percentage_expression_df$Gene[which(percentage_expression_df$PercExpr>5)])


#firstly we work on all the GC
raw_data_for_genes_heatmap=raw_data_GC
spearman_matrix=cor(raw_data_for_genes_heatmap[,list_genes], method="spearman")
my_palette= colorRampPalette(c("purple", "black", "yellow"))(n = 299)
#hclust_GC <- hclust(as.dist(1-spearman_matrix), method="average")

spearman_matrix_dendr=cor(derived_kernel_smoother_coord_df[,list_genes], method="spearman")
hclust_GC <- hclust(as.dist(1-spearman_matrix_dendr), method="average")

# for a triangle heatmap
for(j in 1:(ncol(spearman_matrix)-1)){
  for(k in (nrow(spearman_matrix)-j+1):nrow(spearman_matrix)){
    spearman_matrix[rev(hclust_GC$order)[k],rev(hclust_GC$order)[ncol(spearman_matrix)-j]]=NA
  }
}

# for the gene labels
cellnote_gene_names = matrix(data=NA, nrow= nrow(spearman_matrix), ncol=ncol(spearman_matrix))
for(i in 1:ncol(cellnote_gene_names)){
  cellnote_gene_names[i,i]=colnames(spearman_matrix)[i]
}

# Create the article Fig4_a
# ############################

Physio_BB=heatmap.2(spearman_matrix, col=THREE_COLORS_HEATMAP, tracecol=NA, 
                    Colv=as.dendrogram(hclust_GC), Rowv=as.dendrogram(hclust_GC),
                    cexRow=0.55, cexCol=0.55,
                    main="hclust: Spearman on Kernel",
                    cellnote=cellnote_gene_names,notecex=0.55,notecol="black",
                    labRow = FALSE,labCol = FALSE)


# ####################
# ####################
# Fig4_b and FigSup4
# ####################
# ####################

cluster_of_genes_list=list()
cluster_of_genes_list[['Cluster 1']] = c("TNFRSF17","GCSAM", "CD38","CD79B","BACH2","IL4R","BCL6")
cluster_of_genes_list[['Cluster 2']] = c("AICDA","MME","ATM","SLAMF6")
cluster_of_genes_list[['Cluster 3']] = c("AURKA","TCF3","MKI67","POLH","RAD51","CCNB1","EZH2","GAPDH","DNMT1")
cluster_of_genes_list[['Cluster 4']] = c("CCR7","CXCR5","IRF4","PRDM1","IL10RA")
cluster_of_genes_list[['Cluster 5']] = c("GPR183","CD83","BCL2A1","SLAMF1")
cluster_of_genes_list[['Cluster 6']] = c("FAS", "HLA.DRA", "CD80", "EGR2","B2M","IL21R","MYC")
cluster_of_genes_list[['Surface B cell markers']] = c("CD79B","MS4A1","MME","CD27")
cluster_of_genes_list[['Somatic Hypermutation']] = c("AICDA","ATM","POLH","TP53")
cluster_of_genes_list[['Antigen Presentation']] = c("CD74", "CD83", "CIITA","HLA.DOA")
cluster_of_genes_list[['Co-stimulation 1']] = c( "CD40","CD72", "CD80", "CD86")
cluster_of_genes_list[['Co-stimulation 2']] = c("SEMA4B", "SEMA7A","SLAMF1", "SLAMF6")
cluster_of_genes_list[['Activation']] = c("BCL2A1","EGR2","IRF4","MYC","NFKB2","TNFAIP3")
cluster_of_genes_list[['Migration']] = c("CXCR4","CXCR5","GNA13","GNAI2","GPR183")
cluster_of_genes_list[['Cell Cycle']] = c("AURKA","CCNB1","MKI67")
cluster_of_genes_list[['Epigenetics Modifiers']] = c("CREBBP","DNMT1","EZH2","KMT2C","KMT2D")
cluster_of_genes_list[['Transcription Factors']] = c("BCL6","FOXO1","LMO2","PAX5","TCF3")
cluster_of_genes_list[['Cytokine Receptors']] = c("IL10RA","IL21R","IL4R")

for(i in 1:length(cluster_of_genes_list)){
  genes_plot_kernel=cluster_of_genes_list[[i]]
  
  genes_kernel_smoother_data_frame_all=data.frame()
  for (gene in genes_plot_kernel){
    genes_kernel_smoother_data_frame=data.frame(PseudotimeAsString=kernel_smoother_coord_dataframe$PseudotimeTheta, 
                                                GeneExpr=kernel_smoother_coord_dataframe[,gene],
                                                Gene = gene,
                                                Gene2 = gene)
    genes_kernel_smoother_data_frame_all=rbind(genes_kernel_smoother_data_frame_all,genes_kernel_smoother_data_frame)
  }
  # Create the article Fig4_b
  # ############################
  ggplot_gene_expression_Kernel=ggplot() +
    geom_segment(data=zones_segment_data_frame, mapping=aes(x=xbegin, xend=xend, y=ybegin, yend=yend), linetype="dashed", color="grey") +
    geom_text(data=zones_data_frame, aes(x=x1+(x2-x1)/2, y=y1+3*(y2-y1)/4, label=r), size=3, angle=90, color="grey") +
    expand_limits(y=limits_allFigs) +
    geom_line(data=genes_kernel_smoother_data_frame_all, aes(x=PseudotimeAsString, y=GeneExpr, colour=Gene), size=1.5) +
    xlab("Theta GC") +
    ylab("Et (Gene Expression)") +
    ggtitle(names(cluster_of_genes_list)[i])+
    coord_fixed(1/4) +
    theme(legend.text = element_text(face="italic"))
  print(ggplot_gene_expression_Kernel)
}



# ##########
# ##########
# Fig6_a
# ##########
# ##########


# we work on each sample of the GC
# #################################
for(sample in sort(unique(raw_data_GC$Sample))){
  raw_data_for_genes_heatmap_sample=raw_data_GC[raw_data_GC$Sample==sample,]
  raw_data_for_genes_heatmap_sample = raw_data_for_genes_heatmap_sample[,list_genes]
  spearman_matrix_sample=cor(raw_data_for_genes_heatmap_sample[,hclust_GC$order], method="spearman")
  spearman_matrix_sample=spearman_matrix_sample[nrow(spearman_matrix_sample):1,]
  for(j in 1:(ncol(spearman_matrix_sample)-1)){
    for(k in (nrow(spearman_matrix_sample)-j+1):nrow(spearman_matrix_sample)){
      spearman_matrix_sample[k,j+1]=NA
    }
  }
  # for the gene labels
  cellnote_gene_names = matrix(data=NA, nrow= nrow(spearman_matrix_sample), ncol=ncol(spearman_matrix_sample))
  for(k in 1:ncol(cellnote_gene_names)){
    cellnote_gene_names[nrow(spearman_matrix_sample)-k+1,k]=colnames(spearman_matrix_sample)[k]
  }
  heatmap_plot=heatmap.2(spearman_matrix_sample, col=THREE_COLORS_HEATMAP, 
                         tracecol=NA, Colv=FALSE, Rowv=FALSE, dendrogram="none",
                         main=paste("Spear Corr for: ",sample),
                         cexRow=0.55, cexCol=0.55,
                         cellnote=cellnote_gene_names,notecex=0.55,notecol="black",
                         labRow = FALSE,labCol = FALSE)
}


# we work on each sample of the FL
# #################################
for(sample in sort(unique(raw_data_FL$Sample))){
  raw_data_for_genes_heatmap_sample=raw_data_FL[raw_data_FL$Sample==sample,]
  raw_data_for_genes_heatmap_sample = raw_data_for_genes_heatmap_sample[,list_genes]
  spearman_matrix_sample=cor(raw_data_for_genes_heatmap_sample[,hclust_GC$order], method="spearman")
  spearman_matrix_sample=spearman_matrix_sample[nrow(spearman_matrix_sample):1,]
  for(j in 1:(ncol(spearman_matrix_sample)-1)){
    for(k in (nrow(spearman_matrix_sample)-j+1):nrow(spearman_matrix_sample)){
      spearman_matrix_sample[k,j+1]=NA
    }
  }
  # for the gene labels
  cellnote_gene_names = matrix(data=NA, nrow= nrow(spearman_matrix_sample), ncol=ncol(spearman_matrix_sample))
  for(k in 1:ncol(cellnote_gene_names)){
    cellnote_gene_names[nrow(spearman_matrix_sample)-k+1,k]=colnames(spearman_matrix_sample)[k]
  }
  heatmap_plot=heatmap.2(spearman_matrix_sample, col=THREE_COLORS_HEATMAP, 
                         tracecol=NA, Colv=FALSE, Rowv=FALSE, dendrogram="none",
                         main=paste("Spear Corr for: ",sample),
                         cexRow=0.55, cexCol=0.55,
                         cellnote=cellnote_gene_names,notecex=0.55,notecol="black",
                         labRow = FALSE,labCol = FALSE)
}


# ##########
# ##########
# Fig6_b
# ##########
# ##########

# There is a script dedicated for this figure (Fig6_b_scqpcr.R)

# ##########
# ##########
# Fig6_c
# ##########
# ##########
raw_data_for_plotGC = raw_data_GC
raw_data_for_plotFL = raw_data_FL

min_CXCR4 = min(raw_data$asinh_CXCR4)
min_CD83 = min(raw_data$asinh_CD83)
max_CXCR4 = max(raw_data$asinh_CXCR4)
max_CD83 = max(raw_data$asinh_CD83)

# modify the genes with zero expression to NA value (to fit the ggplot colors)
for(gene in c('MKI67','AICDA')){
  raw_data_for_plotGC[,gene]=sapply(raw_data_for_plotGC[,gene], function(x) if(x==0){x=NA}else{x})
  raw_data_for_plotFL[,gene]=sapply(raw_data_for_plotFL[,gene], function(x) if(x==0){x=NA}else{x})
}

gg_plot = ggplot(raw_data_for_plotGC, aes(asinh_CXCR4, asinh_CD83)) +
  geom_point(aes(color=MKI67), size=size_geom_point) +
  scale_colour_gradient(low =TWO_COLORS_VECTOR[1], high = TWO_COLORS_VECTOR[2], na.value = "grey", limits=limits_allFigs) +
  theme(aspect.ratio=1,
        legend.title = element_text(face = "italic")) +
  xlab("CXCR4 (Index)") +
  ylab("CD83 (Index)") +
  xlim(min_CXCR4,max_CXCR4) +
  ylim(min_CD83,max_CD83)+
  ggtitle("GC")
print(gg_plot)

gg_plot = ggplot(raw_data_for_plotFL, aes(asinh_CXCR4, asinh_CD83)) +
  geom_point(aes(color=MKI67), size=size_geom_point) +
  scale_colour_gradient(low =TWO_COLORS_VECTOR[1], high = TWO_COLORS_VECTOR[2], na.value = "grey", limits=limits_allFigs) +
  theme(aspect.ratio=1,
        legend.title = element_text(face = "italic")) +
  xlab("CXCR4 (Index)") +
  ylab("CD83 (Index)") +
  xlim(min_CXCR4,max_CXCR4) +
  ylim(min_CD83,max_CD83)+
  ggtitle("FL")
print(gg_plot)


gg_plot = ggplot(raw_data_for_plotGC, aes(asinh_CXCR4, asinh_CD83)) +
  geom_point(aes(color=AICDA), size=size_geom_point) +
  scale_colour_gradient(low =TWO_COLORS_VECTOR[1], high = TWO_COLORS_VECTOR[2], na.value = "grey", limits=limits_allFigs) +
  theme(aspect.ratio=1,
        legend.title = element_text(face = "italic")) +
  xlab("CXCR4 (Index)") +
  ylab("CD83 (Index)") +
  xlim(min_CXCR4,max_CXCR4) +
  ylim(min_CD83,max_CD83)+
  ggtitle("GC")
print(gg_plot)

gg_plot = ggplot(raw_data_for_plotFL, aes(asinh_CXCR4, asinh_CD83)) +
  geom_point(aes(color=AICDA), size=size_geom_point) +
  scale_colour_gradient(low =TWO_COLORS_VECTOR[1], high = TWO_COLORS_VECTOR[2], na.value = "grey", limits=limits_allFigs) +
  theme(aspect.ratio=1,
        legend.title = element_text(face = "italic")) +
  xlab("CXCR4 (Index)") +
  ylab("CD83 (Index)") +
  xlim(min_CXCR4,max_CXCR4) +
  ylim(min_CD83,max_CD83)+
  ggtitle("FL")
print(gg_plot)



# ##########
# ##########
# Fig6_d
# ##########
# ##########

raw_data_for_analysis_FL_AICDApos = raw_data_FL
raw_data_for_analysis_GC_AICDApos = raw_data_GC
raw_data_for_analysis_FL_AICDApos = raw_data_for_analysis_FL_AICDApos[which(raw_data_for_analysis_FL_AICDApos$AICDA>0),]
raw_data_for_analysis_GC_AICDApos = raw_data_for_analysis_GC_AICDApos[which(raw_data_for_analysis_GC_AICDApos$AICDA>0),]
raw_data_for_analysis_AICDApos = rbind(raw_data_for_analysis_FL_AICDApos[,c("SimpleSortPheno","asinh_CXCR4")],
                                       raw_data_for_analysis_GC_AICDApos[,c("SimpleSortPheno","asinh_CXCR4")])

wilcoxon_test_AICDA = wilcox.test(x=raw_data_for_analysis_FL_AICDApos$asinh_CXCR4,
                                  y=raw_data_for_analysis_GC_AICDApos$asinh_CXCR4,
                                  paired=FALSE, conf.int=TRUE)

gg_plot = ggplot(raw_data_for_analysis_AICDApos, aes(SimpleSortPheno,asinh_CXCR4)) +
  geom_violin(aes(fill=SimpleSortPheno,colour=SimpleSortPheno)) +
  geom_jitter() +
  ggtitle("AICDA⁺") +
  theme(plot.title = element_text(face = "italic"),
        legend.position="none") +
  xlab("Indexed Phenotype")+
  ylab("CXCR4 (Index)")
print(gg_plot)


raw_data_for_analysis_FL_MKI67pos = raw_data_FL
raw_data_for_analysis_GC_MKI67pos = raw_data_GC
raw_data_for_analysis_FL_MKI67pos = raw_data_for_analysis_FL_MKI67pos[which(raw_data_for_analysis_FL_MKI67pos$MKI67>0),]
raw_data_for_analysis_GC_MKI67pos = raw_data_for_analysis_GC_MKI67pos[which(raw_data_for_analysis_GC_MKI67pos$MKI67>0),]
raw_data_for_analysis_MKI67pos = rbind(raw_data_for_analysis_FL_MKI67pos[,c("SimpleSortPheno","asinh_CXCR4")],
                                       raw_data_for_analysis_GC_MKI67pos[,c("SimpleSortPheno","asinh_CXCR4")])

gg_plot = ggplot(raw_data_for_analysis_MKI67pos, aes(SimpleSortPheno,asinh_CXCR4)) +
  geom_violin(aes(fill=SimpleSortPheno,colour=SimpleSortPheno)) +
  geom_jitter() +
  ggtitle("MKI67⁺") +
  theme(plot.title = element_text(face = "italic"),
        legend.position="none") +
  xlab("Indexed Phenotype")+
  ylab("CXCR4 (Index)")
print(gg_plot)




# ##########
# ##########
# FigSup6_c
# ##########
# ##########

heatmap.3(raw_data_FL,
          genes_v=gene_names_v,
          col=TWO_COLORS_GRADIANT,
          margins = c(5, 14),
          lhei=c(1.2,4,0.5),lwid=c(1.2,4),
          ColSide_name_v=c("Sample"),
          ColSideColorsSize = 2.5,trace="none",
          titleColorGradient='Et',
          main='raw_data_FL')
legend_for_heatmap.3(raw_data_FL,ColSide_name_v=c("Sample"))


# ##########
# ##########
# FigSup6_d
# ##########
# ##########
pca_result = pca_function(raw_data_FL[,gene_names_v], PlotHist=FALSE, 
                          VectorToColorPCA=raw_data_FL$Sample, 
                          NameOfVectorToColor="Sample",dim_all=TRUE,
                          title='raw_data_FL')

gg_plots = pca_gene_information_function(pca_result, dim_all=FALSE, threshold_norm=0.6,
                                         title='raw_data_FL')


