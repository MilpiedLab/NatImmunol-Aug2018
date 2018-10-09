# #################################################
# all scrnaseq figures for Fig4 and Fig6 
# (You need to run Fig3and5_scrnaseq.R before)
# #################################################

# ##################
# Load other scripts
# ##################

# Load the script Fig3and5_scrnaseq.R
SCRIPT_DIR = "~/scripts/"
SCRIPT_SOURCE = paste( SCRIPT_DIR, "Fig3and5_scrnaseq.R", sep="")
source( SCRIPT_SOURCE)

# #########################
# Start working on the data
# #########################


raw_data_for_analysis=GC_df

intersect_genes_v=all_genes_v 


GENES_TO_ANALYSE=intersect_genes_v
if(TRUE){
  # New list of genes (selected with a threshold)
  # #############################################
  percentage_expression_vector = apply(raw_data_for_analysis[,GENES_TO_ANALYSE],2,function(x) sum(x>0)/length(x))
  
  percentage_expression_df=data.frame(Gene=GENES_TO_ANALYSE, PercExpr=percentage_expression_vector)
  percentage_expression_df= percentage_expression_df[order(percentage_expression_df$PercExpr, decreasing = TRUE),]
  percentage_expression_df$NumberOfGenes=1:length(GENES_TO_ANALYSE)
  
  list_genes=as.vector(percentage_expression_df$Gene[which(percentage_expression_df$PercExpr>0.25)])
}
GENES_TO_ANALYSE=list_genes

# ############################################################################################################################################
# calcul the Kernel Smoother of each gene
# ############################################################################################################################################

pseudotime_used_AsString="DZLZ.Score" 
max_pseudotime = max(raw_data_for_analysis[,pseudotime_used_AsString])
min_pseudotime = min(raw_data_for_analysis[,pseudotime_used_AsString])

kernel_smoother_bandwith=max_pseudotime/8 

raw_data_to_kernel=raw_data_for_analysis


genes_kernel_smoother_data_frame_all=data.frame()

kernel_smoother_coord_dataframe = data.frame()

for (j in 1:length(GENES_TO_ANALYSE)){
  genes_kernel_smoother_list=locpoly(x=raw_data_to_kernel[,pseudotime_used_AsString], y=raw_data_to_kernel[,GENES_TO_ANALYSE[j]], bandwidth = kernel_smoother_bandwith)
  genes_kernel_smoother_data_frame=data.frame(PseudotimeAsString=genes_kernel_smoother_list[[1]], GeneExpr=genes_kernel_smoother_list[[2]])
  genes_kernel_smoother_data_frame$Gene=GENES_TO_ANALYSE[j]
  genes_kernel_smoother_data_frame$Gene2=GENES_TO_ANALYSE[j]
  #kernel_smoother_data_frame_list[[i]]=genes_kernel_smoother_data_frame
  genes_kernel_smoother_data_frame_all=rbind(genes_kernel_smoother_data_frame_all,genes_kernel_smoother_data_frame)
  
  kernel_smoother_coord_dataframe = rbind(kernel_smoother_coord_dataframe,genes_kernel_smoother_data_frame$GeneExpr)
}
kernel_smoother_coord_dataframe = as.data.frame(t(kernel_smoother_coord_dataframe))
colnames(kernel_smoother_coord_dataframe) = GENES_TO_ANALYSE
rownames(kernel_smoother_coord_dataframe) =1:nrow(kernel_smoother_coord_dataframe)
kernel_smoother_coord_dataframe$PseudotimeAsString = genes_kernel_smoother_list[[1]]

# ##############################################
# Calcul the derive of the kernel smoother curve
# ##############################################
x=kernel_smoother_coord_dataframe$PseudotimeAsString
derived_kernel_smoother_coord_df = data.frame(Pseudotime = kernel_smoother_coord_dataframe$PseudotimeAsString[1:(nrow(kernel_smoother_coord_dataframe)-1)])
for( gene in GENES_TO_ANALYSE){
  y=kernel_smoother_coord_dataframe[,gene]
  derivedGene_vector=c()
  for(i in 1:(nrow(kernel_smoother_coord_dataframe)-1)){
    
    derivedGene = (y[i+1]-y[i])/(x[i+1]-x[i])
    
    derivedGene_vector=c(derivedGene_vector,derivedGene) 
  }
  derived_kernel_smoother_coord_df[,gene] = derivedGene_vector
}


#firstly we work on all the GC
raw_data_for_genes_heatmap=raw_data_for_analysis

my_palette= colorRampPalette(c("purple", "black", "yellow"))(n = 299)


spearman_matrix_dendr=cor(kernel_smoother_coord_dataframe[,list_genes], method="spearman")

hclust_GC <- hclust(as.dist(1-spearman_matrix_dendr), method="complete")


spearman_matrix=cor(raw_data_for_genes_heatmap[,list_genes], method="spearman")
spearman_matrix_rearrange = spearman_matrix[rev(hclust_GC$order),hclust_GC$order]
spearman_matrix_withNA = t(sapply(1:ncol(spearman_matrix_rearrange), function(i) x=c(spearman_matrix_rearrange[i,1:(ncol(spearman_matrix_rearrange)-i+1)],rep(NA,i-1 )) ) )
rownames(spearman_matrix_withNA) = rownames(spearman_matrix_rearrange)
spearman_matrix = spearman_matrix_withNA[rownames(spearman_matrix),colnames(spearman_matrix)]

# for the gene labels
cellnote_gene_names = matrix(data=NA, nrow= nrow(spearman_matrix), ncol=ncol(spearman_matrix))
for(i in 1:ncol(cellnote_gene_names)){
  cellnote_gene_names[i,i]=colnames(spearman_matrix)[i]
}

THREE_COLORS_HEATMAP = colorRampPalette(c("blue", "white", "red"))(n = 299)


raw_data_for_analysis$PseudotimeAsString = raw_data_for_analysis$DZLZ.Score


# ##########
# ##########
# Fig3_g
# ##########
# ##########

genes_test_v = c("CXCR4","CD83")
for(gene in (genes_test_v)){
  gg_plot = ggplot(raw_data_for_analysis, aes_string("PseudotimeAsString",gene))+
    geom_point( )+
    geom_vline(xintercept=startGreyZoneScore)+
    geom_vline(xintercept=endGreyZoneScore)+
    geom_line(data=kernel_smoother_coord_dataframe, aes_string("PseudotimeAsString",gene), size=2)+
    labs(x = "DZ_LZ_score")
  print(gg_plot)
}

raw_data_for_analysis$CXCR4_CD83_DblExpr = 0
raw_data_for_analysis[raw_data_for_analysis$CXCR4!=0 & raw_data_for_analysis$CD83!=0,]$CXCR4_CD83_DblExpr = 1

DblExpr_kernel_smoother_list=locpoly(x=raw_data_for_analysis[,pseudotime_used_AsString], y=raw_data_for_analysis[,"CXCR4_CD83_DblExpr"], bandwidth = 0.2)
DblExpr_kernel_smoother_data_frame=data.frame(PseudotimeAsString=DblExpr_kernel_smoother_list[[1]], CXCR4_CD83_DblExpr=DblExpr_kernel_smoother_list[[2]])

gg_plot = ggplot(raw_data_for_analysis, aes(PseudotimeAsString, CXCR4_CD83_DblExpr))+
  geom_point( )+
  geom_vline(xintercept=startGreyZoneScore)+
  geom_vline(xintercept=endGreyZoneScore)+
  geom_line(data=DblExpr_kernel_smoother_data_frame, aes_string("PseudotimeAsString","CXCR4_CD83_DblExpr"), size=2)+
  labs(x = "DZ_LZ_score")
print(gg_plot)

# ############################################################
# ############################################################
# ############################################################

raw_data_for_analysis$DZLZ.Score.bin = cut(raw_data_for_analysis$DZLZ.Score, 3)
raw_data_for_analysis$DZLZ.Score.bin2 = sapply(raw_data_for_analysis$DZLZ.Score.bin, function(x) which(unique(raw_data_for_analysis$DZLZ.Score.bin)==x)  )


spearman_df = as.data.frame(t(spearman_matrix))

cellcycle_genes_v = unique(c(cc.genes))#,GO_CellCycle_v))

spearman_df$GeneList = NA
spearman_df[rownames(spearman_df) %in% victora_genes_LZ_v,]$GeneList = "Victora_LZ"
spearman_df[rownames(spearman_df) %in% victora_genes_DZ_v,]$GeneList = "Victora_DZ"
gene_top20_GreyZone = top20_DiffPerc_LZDZ$gene[21:40]
spearman_df[rownames(spearman_df) %in% gene_top20_GreyZone,]$GeneList = "top20_GreyZone"

cluster_hclust = cutree(hclust_GC, k=8)

spearman_df$cluster_hclust = cluster_hclust


# ##########
# ##########
# Fig4_c
# ##########
# ##########

heatmap.3(spearman_df,
          genes_v=rownames(spearman_df), #c("AICDA","STMN1"),#
          col=THREE_COLORS_HEATMAP,tracecol=NA,
          Colv=as.dendrogram(hclust_GC), Rowv=as.dendrogram(hclust_GC),
          #cellnote=cellnote_gene_names,notecex=0.25,
          #margins = c(5, 14),
          lhei=c(1.2,4,0.5),lwid=c(1.2,4),
          ColSide_name_v=c("GeneList","cluster_hclust"),
          ColSideColorsSize = 2.5,trace="none",
          notecol="black",
          labRow = FALSE,labCol = FALSE, 
          titleColorGradient="SpearCor") 
legend_for_heatmap.3(spearman_df,ColSide_name_v=c("GeneList","cluster_hclust"))

# #################################################
# #################################################
# #################################################
# #################################################

melt_kernel_smoother_coord_dataframe_df =  melt(kernel_smoother_coord_dataframe[which(1:nrow(kernel_smoother_coord_dataframe) %% 5 == 1),
                                                                                c(rownames(spearman_df),"PseudotimeAsString")], 
                                                measure.vars=rownames(spearman_df))
melt_kernel_smoother_coord_dataframe_df$hcluster = spearman_df[melt_kernel_smoother_coord_dataframe_df$variable,]$cluster_hclust



melt_kernel_smoother_coord_dataframe_df = melt_kernel_smoother_coord_dataframe_df[, c("PseudotimeAsString","value","hcluster")] %>%
  group_by(PseudotimeAsString, hcluster) %>%
  summarise_all( funs(mean(., na.rm=TRUE)))

for( i in unique(melt_kernel_smoother_coord_dataframe_df$hcluster)){
  hclust_df = melt_kernel_smoother_coord_dataframe_df[melt_kernel_smoother_coord_dataframe_df$hcluster==i,]
  hclust_df$value = (hclust_df$value-min(hclust_df$value))/max(hclust_df$value-min(hclust_df$value))
  if(i==1){
    norm_melt_kernel_smoother_coord_dataframe_df = hclust_df
  }else{
    norm_melt_kernel_smoother_coord_dataframe_df = rbind(norm_melt_kernel_smoother_coord_dataframe_df, hclust_df)
  }
}


norm_kernel_smoother_coord_dataframe = kernel_smoother_coord_dataframe
norm_kernel_smoother_coord_dataframe[,GENES_TO_ANALYSE] =(apply(norm_kernel_smoother_coord_dataframe[,GENES_TO_ANALYSE],2,function(x) (x-min(x))/max(x-min(x))) )


norm_melt_annotated_df =  melt(norm_kernel_smoother_coord_dataframe[which(1:nrow(norm_kernel_smoother_coord_dataframe) %% 10 == 1),c(GENES_TO_ANALYSE,"PseudotimeAsString")], 
                               measure.vars=GENES_TO_ANALYSE)
norm_melt_annotated_df$variable = as.vector(norm_melt_annotated_df$variable)
norm_melt_annotated_df$GeneList = spearman_df[norm_melt_annotated_df$variable,]$GeneList
norm_melt_annotated_df$hcluster = spearman_df[norm_melt_annotated_df$variable,]$cluster_hclust

# ##########
# ##########
# Fig4_d
# ##########
# ##########

gg_plot = ggplot()+
  geom_line(data=norm_melt_annotated_df, aes(PseudotimeAsString,value, group=variable, color=GeneList), size=0.1)+
  ggtitle("normalized curves of annotated genes\n black line is the mean for each cluster")+
  facet_wrap(~hcluster)+
  geom_line(data=norm_melt_kernel_smoother_coord_dataframe_df, aes(PseudotimeAsString,value))
print(gg_plot)


# ##########################################
# Spearman FL triangle heatmap
# ##########################################

raw_data_FL_onlyGenes = FL_df[,list_genes]
raw_data_FL_onlyGenes_ordered = raw_data_FL_onlyGenes[,hclust_GC$order]
spearman_matrix_FL=cor(raw_data_FL_onlyGenes_ordered, method="spearman")

spearman_matrix_FL=spearman_matrix_FL[nrow(spearman_matrix_FL):1,]
if(FALSE){
  for(j in 1:(ncol(spearman_matrix_FL)-1)){
    for(k in (nrow(spearman_matrix_FL)-j+1):nrow(spearman_matrix_FL)){
      spearman_matrix_FL[k,j+1]=NA
    }
  }
}
spearman_matrix_FL_withNA = t(sapply(1:ncol(spearman_matrix_FL), function(i) x=c(spearman_matrix_FL[i,1:(ncol(spearman_matrix_FL)-i+1)],rep(NA,i-1 )) ) )
rownames(spearman_matrix_FL_withNA) = rownames(spearman_matrix_FL)
spearman_matrix_FL = spearman_matrix_FL_withNA[rownames(spearman_matrix_FL),colnames(spearman_matrix_FL)]

spearman_df_FL = as.data.frame(t(spearman_matrix_FL))
spearman_df_FL$GeneList = NA
spearman_df_FL[rownames(spearman_df_FL) %in% victora_genes_LZ_v,]$GeneList = "Victora_LZ"
spearman_df_FL[rownames(spearman_df_FL) %in% victora_genes_DZ_v,]$GeneList = "Victora_DZ"
#spearman_df_FL[rownames(spearman_df_FL) %in% cellcycle_genes_v,]$GeneList = "cellcycle_list"
spearman_df_FL[rownames(spearman_df_FL) %in% gene_top20_GreyZone,]$GeneList = "top20_GreyZone"
spearman_df_FL$cluster_hclust = cluster_hclust[rownames(spearman_df_FL)]

# ##########
# ##########
# Fig6_e
# ##########
# ##########

heatmap.3(spearman_df_FL,
          genes_v=rownames(spearman_matrix_FL),
          col=THREE_COLORS_HEATMAP,tracecol=NA,
          Colv=FALSE, Rowv=FALSE, dendrogram="none",
          #cellnote=cellnote_gene_names,notecex=0.25,
          #margins = c(5, 14),
          lhei=c(1.2,4,0.5),lwid=c(1.2,4),
          ColSide_name_v=c("GeneList","cluster_hclust"),
          ColSideColorsSize = 2.5,trace="none",
          notecol="black",
          labRow = FALSE,labCol = FALSE, 
          titleColorGradient="SpearCor",
          main="Spearman FL") 
legend_for_heatmap.3(spearman_df,ColSide_name_v=c("GeneList","cluster_hclust"))
