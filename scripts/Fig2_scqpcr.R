# #######################################
# all figures for Fig2 and Fig2 Sup
# #######################################
require(scales) # for FigSup2_e (labels=percent)

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

# #####################################
# #####################################
# Weighted voting predictor analysis
# #####################################
# #####################################

# Get the subset of data linked to cells in light and dark zones
class1_zone_data = raw_data_GC[ which( raw_data_GC$IndexedPheno=="GC DZ"), ]
class2_zone_data = raw_data_GC[ which( raw_data_GC$IndexedPheno=="GC LZ"), ]

# Look at the number of cells in light and dark zones
class1_size = nrow( class1_zone_data)
class2_size = nrow( class2_zone_data)

# Define a dataset with both cell types
both_zone_data = rbind( class1_zone_data, class2_zone_data)

# Define the names of genes to take into account
WEIGHTED_VOTING_MODEL_GENES = gene_names_v

# Apply a Leave One Out Cross-validation method to evaluate the precision of the model
false = 0
bad_indexes_set = vector()
bad_indexes_strength = vector()
all_indexes_strength = vector()
for( turn in seq(1, class1_size + class2_size,1)){
  
  # Prepare the set of indexes to train the model
  # by excluding only the current index
  class1_training_indexes = seq(1, class1_size,1)
  class2_training_indexes = seq(class1_size + 1, class1_size + class2_size,1)
  
  if( turn <= class1_size){
    class1_training_indexes = class1_training_indexes[ -turn]
  }else{
    class2_training_indexes = class2_training_indexes[ -turn]
  }
  
  # For each gene, compute the parameters of the model (P(g,c) and b values) (see article)
  b = list()
  pgc = list()
  for( gene_name in WEIGHTED_VOTING_MODEL_GENES){
    
    ml = mean( both_zone_data[ class1_training_indexes, gene_name])
    sdl = sd( both_zone_data[ class1_training_indexes, gene_name])
    md = mean( both_zone_data[ class2_training_indexes, gene_name])
    sdd = sd( both_zone_data[ class2_training_indexes, gene_name])
    
    b[[ gene_name]] = (ml + md) / 2
    pgc[[ gene_name]] = (ml - md) / (sdl + sdd)
  }
  
  # Compute the prediction of the excluded cell with the model built before
  if( turn <= class1_size){
    class1_result = test_model( c(turn), both_zone_data, pgc, b, TRUE)
    false = false + class1_result[[ "bad_number"]]
    bad_indexes_set = append( bad_indexes_set, class1_result[[ "bad_indexes"]])
    bad_indexes_strength = append( bad_indexes_strength, class1_result[[ "bad_strengths"]])
    all_indexes_strength = append( all_indexes_strength, class1_result[[ "vote_strengths"]])
  }else{
    class2_result = test_model( c(turn), both_zone_data, pgc, b, FALSE)
    false = false + class2_result[[ "bad_number"]]
    bad_indexes_set = append( bad_indexes_set, class2_result[[ "bad_indexes"]])
    bad_indexes_strength = append( bad_indexes_strength, class2_result[[ "bad_strengths"]])
    all_indexes_strength = append( all_indexes_strength, class2_result[[ "vote_strengths"]])
  }
}

# Compute the FDR of the cross-validation computation
cat("\nWeighted voting predictor model FDR (LOOC)= ", false / (class1_size + class2_size))

# Export the complete data with predictor strentgh
both_zone_data$predictorStrength = all_indexes_strength

# Export the cells that were False Positive during the LOOC to file
bad_vote_df = both_zone_data[ bad_indexes_set, ]
bad_vote_df$predictorStrength = bad_indexes_strength

cat("\nNumber of False Positive cells during LOOC :", length( bad_indexes_set))

# Normalize the strength to have them between -100 and 100
# factor = max( abs(max( all_indexes_strength)), abs(min(all_indexes_strength)))
# normalized_strength = 100* all_indexes_strength / factor
factor_pos = abs(max( all_indexes_strength))
factor_neg = abs(min(all_indexes_strength))
normalized_strength = all_indexes_strength
normalized_strength[which( normalized_strength <0 )] = 100* normalized_strength[which( normalized_strength <0 )] / factor_neg
normalized_strength[which( normalized_strength >0 )] = 100* normalized_strength[which( normalized_strength >0 )] / factor_pos



both_zone_data$NormalizedStrength=normalized_strength
rownames(both_zone_data)=both_zone_data$UniqueCellID
normalized_strength_ordered=c()
for(uniquecellid in raw_data_GC$UniqueCellID){
  normalized_strength_ordered = c(normalized_strength_ordered,both_zone_data[uniquecellid,"NormalizedStrength"])
}
raw_data_GC$NormalizedStrength = normalized_strength_ordered

raw_data_GC$FalsePositive=rep(NA,nrow(raw_data_GC))
raw_data_GC$FalsePositive[which(raw_data_GC$UniqueCellID %in% bad_vote_df$UniqueCellID)]="FalsePos"


# ##########
# ##########
# Fig2_a
# ##########
# ##########

gg_plot = ggplot(raw_data_GC, aes(asinh_CXCR4, asinh_CD83)) +
  geom_path(data = gate1, aes(x,y), color="gray50") +
  geom_path(data = gate2, aes(x,y), color="gray50") +
  geom_point(aes( color=NormalizedStrength ) , size=2) +
  geom_point(data=raw_data_GC[which(!is.na(raw_data_GC$FalsePositive)),], aes(fill="black"), size=2) +
  scale_colour_gradient2(low="red", mid="white",high="blue", name="Prediction Strength") +
  coord_fixed() +
  theme(aspect.ratio=1,panel.background = element_rect(fill = 'gray80') ) +
  xlab("CXCR4 (Index)") +
  ylab("CD83 (Index)") +
  scale_fill_manual(name="",values="black",breaks="black", labels="False Positive Cells")
print(gg_plot)  


# ##########
# ##########
# FigSup2a
# ##########
# ##########
raw_data_for_analysis=raw_data_GC
# we remove the ones which are not labeled as LZ or DZ
raw_data_for_analysis=raw_data_for_analysis[which(!raw_data_for_analysis$IndexedPheno=='GC other'),]

result_volcano_plot=volcano_plot_LRT_function(raw_data_for_analysis[,gene_names_v],raw_data_for_analysis$IndexedPheno)


# ##########
# ##########
# Fig2_b
# ##########
# ##########

raw_data_GC_DZLZ = raw_data_GC[raw_data_GC$IndexedPheno!="GC other",]
pca_result = pca_function(raw_data_GC_DZLZ[,gene_names_v], PlotHist=FALSE, 
                          VectorToColorPCA=raw_data_GC_DZLZ$MKI67, 
                          NameOfVectorToColor="MKI67",dim_all=FALSE,
                          title='raw_data_GC_DZLZ')

gg_plots = pca_gene_information_function(pca_result, dim_all=FALSE, threshold_norm=0.6,
                                         title='raw_data_GC_DZLZ')



# ##########
# ##########
# FigSup2_b
# ##########
# ##########
pca_result = pca_function(raw_data_GC_DZLZ[,gene_names_v], PlotHist=FALSE, 
                          VectorToColorPCA=raw_data_GC_DZLZ$Sample, 
                          NameOfVectorToColor="Sample",dim_all=FALSE,
                          title='raw_data_GC_DZLZ')
# ##########
# ##########
# Fig2_c
# ##########
# ##########

raw_data_GC_DZLZ$PC1=pca_result$li$Axis1
raw_data_GC_DZLZ$PC2=pca_result$li$Axis2
raw_data_GC_DZLZ$PC3=pca_result$li$Axis3

gg_plot = ggplot(raw_data_GC_DZLZ, aes(PC1, PC2)) +
  coord_fixed() +
  geom_point(aes( color=NormalizedStrength ) , size=2) +
  geom_point(data=raw_data_GC_DZLZ[which(!is.na(raw_data_GC_DZLZ$FalsePositive)),], aes(fill="black"), size=2) +
  scale_colour_gradient2(low="red", mid="white",high="blue", name="Prediction Strength") +
  #scale_colour_gradient(low="white", high="orange") +
  theme(aspect.ratio=1,panel.background = element_rect(fill = 'gray80') ) +
  scale_fill_manual(name="",values="black",breaks="black", labels="False Positive Cells")+
  ggtitle("Representation of the Weight Voting Predictor\n onto the PCA")
print(gg_plot)


# ##########
# ##########
# Fig2_d
# ##########
# ##########

# take the Principal Composants Coordinates
# #########################################
pca = pca_function(raw_data_GC[,gene_names_v], bool_plot = FALSE)
raw_data_GC$PC1=pca$li$Axis1
raw_data_GC$PC2=pca$li$Axis2
raw_data_GC$PC3=pca$li$Axis3
pca_percentageInfoByAxis = as.integer(100*pca$eig/sum( pca$eig)+0.5)
PC1_mean=mean(raw_data_GC$PC1)
PC2_mean=mean(raw_data_GC$PC2)

# K-means 
# ##############
set.seed(42) # Sets seed for reproducibility
number_K=5
kmeans_cluster = kmeans(raw_data_GC[,gene_names_v],number_K)
raw_data_GC$KmeansCluster = factor(kmeans_cluster$cluster)

PCX="PC1"
PCY="PC2"
# we change the scale of the plot's axis to have a perfect pca plot 
distPCX= max(raw_data_GC[,PCX])-min(raw_data_GC[,PCX])
distPCY= max(raw_data_GC[,PCY])-min(raw_data_GC[,PCY])
addIntervalPCX=0
addIntervalPCY=0
if(distPCX>distPCY){
  addIntervalPCY=(distPCX-distPCY)/2
} else {
  addIntervalPCX=(distPCY-distPCX)/2
}

xend_v=c(-31.5,-24,26,39,6)
rayon=40
yend_v=sqrt(rayon^2-xend_v^2)

Physio_AS <- ggplot(raw_data_GC, aes_string(PCX, PCY, color="KmeansCluster")) +
  geom_segment(aes(x=PC1_mean,y=PC2_mean,xend=xend_v[1],yend=-yend_v[1]), color=color_lines) +
  geom_segment(aes(x=PC1_mean,y=PC2_mean,xend=xend_v[2],yend=yend_v[2]), color=color_lines) +
  geom_segment(aes(x=PC1_mean,y=PC2_mean,xend=xend_v[3],yend=yend_v[3]), color=color_lines) +
  geom_segment(aes(x=PC1_mean,y=PC2_mean,xend=xend_v[4],yend=-yend_v[4]), color=color_lines) +
  geom_segment(aes(x=PC1_mean,y=PC2_mean,xend=xend_v[5],yend=-yend_v[5]), color=color_lines) +
  geom_point(size=size_geom_point) +
  coord_fixed() +
  xlab(paste0(PCX," (",pca_percentageInfoByAxis[1],"%)")) +
  ylab(paste0(PCY," (",pca_percentageInfoByAxis[2],"%)")) +
  xlim(min(raw_data_GC[,PCX])-addIntervalPCX, max(raw_data_GC[,PCX])+addIntervalPCX) +
  ylim(min(raw_data_GC[,PCY])-addIntervalPCY, max(raw_data_GC[,PCY])+addIntervalPCY) +
  #theme(legend.position="none") 
  labs(color="k-means clusters")
print(Physio_AS)

# ##########
# ##########
# Fig2_e
# ##########
# ##########


pseudotime_theta_vector=c()
radius_vector=c()
for(i in 1:nrow(raw_data_GC)){
  x_new=raw_data_GC$PC1[i]-PC1_mean
  y_new=raw_data_GC$PC2[i]-PC2_mean
  radius=sqrt(x_new^2+y_new^2)
  if(x_new<0 && y_new<0){
    theta=atan(y_new/x_new)
  } else if(x_new>0 && y_new<0){
    theta=pi+atan(y_new/x_new)
  } else if(x_new>0 && y_new>0){
    theta=pi+atan(y_new/x_new)
  } else if(x_new<0 && y_new>0){
    theta=2*pi+atan(y_new/x_new)
  }
  pseudotime_theta_vector=c(pseudotime_theta_vector,theta)
  radius_vector=c(radius_vector,radius)
}
names(pseudotime_theta_vector)=raw_data_GC$UniqueCellID

# reorder to fit the pseudotime DDRTree
pseudotime_theta_vector_reordered=pseudotime_theta_vector
for (i in 1:length(pseudotime_theta_vector_reordered)){
  if(pseudotime_theta_vector_reordered[i]>2){
    pseudotime_theta_vector_reordered[i]=pseudotime_theta_vector_reordered[i]-2*pi
  }
}

# get the value from 0 to pi*2
minimum_theta= min(pseudotime_theta_vector_reordered)
if(minimum_theta<0){
  pseudotime_theta_vector_reordered = pseudotime_theta_vector_reordered-minimum_theta
} else {
  pseudotime_theta_vector_reordered = pseudotime_theta_vector_reordered+minimum_theta
}

# inverse the direction of the pseudotimeTheta for our convenience
pseudotime_theta_vector_reordered = -pseudotime_theta_vector_reordered
pseudotime_theta_vector_reordered = pseudotime_theta_vector_reordered - min(pseudotime_theta_vector_reordered)

raw_data_GC$PseudotimeTheta = pseudotime_theta_vector_reordered
raw_data_GC$PseudotimeRadius = radius_vector



radian_bonus = pi -pi/8 # to fit the same position as in the PCA
radius_unitaire_vector=raw_data_GC$PseudotimeRadius/max(raw_data_GC$PseudotimeRadius)
#radius=30+runif(length(pseudotimeTheta),max=20)
radius=13+radius_unitaire_vector*30
x_data=sin(raw_data_GC$PseudotimeTheta+radian_bonus )*radius
y_data=cos(raw_data_GC$PseudotimeTheta+radian_bonus )*radius
raw_data_GC$x_theta = x_data
raw_data_GC$y_theta = y_data

# Create the article Fig2_e
# ############################
Physio_AT = ggplot(raw_data_GC, aes(x_theta, y_theta, color=KmeansCluster)) +
  geom_point(size=size_geom_point) +
  coord_fixed() +
  theme(axis.ticks = element_blank(), 
        axis.text.x = element_blank(), 
        axis.text.y = element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.line=element_blank()) +
  labs(color="k-means clusters")
print(Physio_AT)



# ##########
# ##########
# FigSup2_c
# ##########
# ##########
mean = mean(raw_data_GC$PseudotimeRadius)
sd = sd(raw_data_GC$PseudotimeRadius)
quartils = quantile( raw_data_GC$PseudotimeRadius)
gg_plot = ggplot(raw_data_GC, aes(PseudotimeRadius))+geom_histogram(binwidth=2)+
  ggtitle(paste("Distribution of the radius from the center of gravity",
                "\n mean = ",signif(mean,3),"sd=",signif(sd,3),
                "\n (red=mean and blue=95%CI)"))+
  geom_vline(xintercept = mean, color="red")+
  geom_vline(xintercept = mean-1.96*sd, color="blue")+
  geom_vline(xintercept = mean+1.96*sd, color="blue")
print(gg_plot)

# ##########
# ##########
# FigSup2_d
# ##########
# ##########

raw_data_for_analysis = raw_data_GC
for(i in 3:7){
  set.seed(42) # Sets seed for reproducibility
  number_K=i
  kmeans_cluster = kmeans(raw_data_for_analysis[,gene_names_v],number_K, iter.max=500)
  raw_data_for_analysis$KmeansCluster = paste0("kcluster_",kmeans_cluster$cluster)
  meanKcluster_v = sapply(sort(unique(raw_data_for_analysis$KmeansCluster)), 
                          FUN=function(x) median(raw_data_for_analysis$PseudotimeTheta[raw_data_for_analysis$KmeansCluster==x]))
  names(meanKcluster_v) = sort(unique(raw_data_for_analysis$KmeansCluster))
  meanKcluster_v = sort(meanKcluster_v)
  raw_data_for_analysis$KmeansCluster = factor(paste0("kcluster_",kmeans_cluster$cluster),levels= names(meanKcluster_v))
  
  ggplot_list=list()
  
  ggplot_list[[1]] = ggplot(raw_data_for_analysis, aes(PC1,PC2, color=KmeansCluster ))+
    geom_point()+coord_fixed()+ theme(legend.position="none")+
    ggtitle(paste0("Kclusters= ",number_K))
  ggplot_list[[2]] = ggplot(raw_data_for_analysis, aes(KmeansCluster,PseudotimeTheta, color=KmeansCluster ))+
    geom_violin()+ 
    geom_jitter()+
    theme(legend.position="none")
  
  
  do.call("grid.arrange", c(ggplot_list, nrow=1))
  cat("\n---\n")
}

# ##########
# ##########
# FigSup2_e
# ##########
# ##########


ratio_sample_kcluster_df = data.frame(KmeansCluster=as.vector(sort(unique(raw_data_GC$KmeansCluster))), stringsAsFactors = FALSE)
ratio_sample_kcluster_df$ratioSample_SP1 = sapply(ratio_sample_kcluster_df$KmeansCluster, FUN=function(x) 
  sum(raw_data_GC$KmeansCluster==x & raw_data_GC$Sample=="SP1")/sum(raw_data_GC$KmeansCluster==x) )
ratio_sample_kcluster_df$ratioSample_SP2 = sapply(ratio_sample_kcluster_df$KmeansCluster, FUN=function(x) 
  sum(raw_data_GC$KmeansCluster==x & raw_data_GC$Sample=="SP2")/sum(raw_data_GC$KmeansCluster==x) )
ratio_sample_kcluster_df$ratioSample_SP3 = sapply(ratio_sample_kcluster_df$KmeansCluster, FUN=function(x) 
  sum(raw_data_GC$KmeansCluster==x & raw_data_GC$Sample=="SP3")/sum(raw_data_GC$KmeansCluster==x) )
ratio_sample_kcluster_df$ratioSample_TS1 = sapply(ratio_sample_kcluster_df$KmeansCluster, FUN=function(x) 
  sum(raw_data_GC$KmeansCluster==x & raw_data_GC$Sample=="TS1")/sum(raw_data_GC$KmeansCluster==x) )

ratio_sample_kcluster_df2 = data.frame(KmeansCluster=rep(ratio_sample_kcluster_df$KmeansCluster,4),
                                       ratioSample = c(ratio_sample_kcluster_df$ratioSample_SP1,
                                                       ratio_sample_kcluster_df$ratioSample_SP2,
                                                       ratio_sample_kcluster_df$ratioSample_SP3,
                                                       ratio_sample_kcluster_df$ratioSample_TS1),
                                       SampleName = rep(c("SP1","SP2","SP3","TS1"),each=length(ratio_sample_kcluster_df$KmeansCluster)) )

gg_plot = ggplot(ratio_sample_kcluster_df2, aes(y = ratioSample, x = KmeansCluster, fill = SampleName) ) + geom_bar(
  stat="identity")+scale_y_continuous(labels = percent)+
  ggtitle("The Distribution is not really beautiful...")
print(gg_plot)


# ##########
# ##########
# FigSup2_f
# ##########
# ##########

ratio_LZDZ_kcluster_df = data.frame(KmeansCluster=as.vector(sort(unique(raw_data_GC$KmeansCluster))), stringsAsFactors = FALSE)
ratio_LZDZ_kcluster_df$ratioLZ = sapply(ratio_LZDZ_kcluster_df$KmeansCluster, FUN=function(x) 
  sum(raw_data_GC$KmeansCluster==x & raw_data_GC$IndexedPheno=="GC LZ")/sum(raw_data_GC$KmeansCluster==x) )
ratio_LZDZ_kcluster_df$ratioDZ = sapply(ratio_LZDZ_kcluster_df$KmeansCluster, FUN=function(x) 
  sum(raw_data_GC$KmeansCluster==x & raw_data_GC$IndexedPheno=="GC DZ")/sum(raw_data_GC$KmeansCluster==x) )
ratio_LZDZ_kcluster_df$ratioOther = sapply(ratio_LZDZ_kcluster_df$KmeansCluster, FUN=function(x) 
  sum(raw_data_GC$KmeansCluster==x & raw_data_GC$IndexedPheno=="GC other")/sum(raw_data_GC$KmeansCluster==x) )

ratio_sample_kcluster_df2 = data.frame(KmeansCluster=rep(ratio_LZDZ_kcluster_df$KmeansCluster,3),
                                       ratioLZDZ = c(ratio_LZDZ_kcluster_df$ratioLZ,
                                                     ratio_LZDZ_kcluster_df$ratioDZ,
                                                     ratio_LZDZ_kcluster_df$ratioOther),
                                       LZ_DZ = rep(c("LZ","DZ","other"),each=length(ratio_sample_kcluster_df$KmeansCluster)) )

gg_plot = ggplot(ratio_sample_kcluster_df2, aes(y = ratioLZDZ, x = KmeansCluster, fill = LZ_DZ) ) + geom_bar(
  stat="identity")+scale_y_continuous(labels = percent)+
  ggtitle("The Distribution of LZ DZ")
print(gg_plot)

# ###################
# ###################
# Fig2_f and FigSup2g
# ###################
# ###################

raw_data_for_analysis=raw_data_GC
for(gene in c("MKI67","AICDA","BCL2A1","CD83","CCNB1")){
  raw_data_for_analysis[,gene]=sapply(raw_data_for_analysis[,gene], function(x) if(x==0){x=NA}else{x})
  Physio_AT = ggplot(raw_data_for_analysis, aes_string('x_theta', 'y_theta', color=gene)) +
    geom_point(size=size_geom_point) +
    coord_fixed() +
    theme(axis.ticks = element_blank(), 
          axis.text.x = element_blank(), 
          axis.text.y = element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          axis.line=element_blank()) +
    scale_colour_gradient(low =TWO_COLORS_VECTOR[1], 
                          high = TWO_COLORS_VECTOR[2],
                          limits=limits_allFigs,na.value = color_NA_geom_point, name="Et")+
    ggtitle(paste("GC cells colored by", gene))
  print(Physio_AT) 
}

# ##########
# ##########
# Fig2_g
# ##########
# ##########

# monocle DDRTree (plots to view what happens)
results_DDRTree_GC = pseudotimeDDRTree_function(raw_data_GC, GENES_TO_ANALYSE=gene_names_v,
                                                cell_attribute_colorBy=c("Pseudotime","State","IndexedPheno"),
                                                DDRTree_max_components=3,cellsAnalyse='GC')

raw_data_GC$DDRTreePseudotime = as.vector(results_DDRTree_GC$DDRTree_pseudotime)

DDRTree_pseudotime = raw_data_GC$DDRTreePseudotime
names(DDRTree_pseudotime)=raw_data_GC$UniqueCellID

# rereordering theta (because PseudotimeDDRtreeGC changes between each computer)
pseudotime_theta_vector_rereorder = pseudotime_theta_vector_reordered
pseudotime_theta_vector_rereorder[pseudotime_theta_vector_rereorder>4.75] = pseudotime_theta_vector_rereorder[pseudotime_theta_vector_rereorder>4.75]-pi*2
pseudotime_theta_vector_rereorder = pseudotime_theta_vector_rereorder+pi*2-4.75

# Create the article figure 2 g
# ############################
gg_plot = correlation_pseudotime_function(pseudotime_theta_vector_rereorder, 
                                          DDRTree_pseudotime, "PseudotimeThetaGC","Pseudotime DDRTree GC" )
print(gg_plot)

# ##########
# ##########
# Fig2_h
# ##########
# ##########


vector_genes = c("SEMA7A","CD83","GPR183","AICDA","MKI67","CXCR4","CXCR5","IRF4","MYC")

for (i in 1:length(vector_genes)){
  gene_name=vector_genes[i]
  Physio_BD1=ggplot() +
    geom_rect(data=zones_data_frame, mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2, fill=Zone), alpha=0.2) +
    geom_text(data=zones_data_frame, aes(x=x1+(x2-x1)/2, y=y1+3*(y2-y1)/4, label=r), size=3, angle=90, color="grey30") +
    geom_point(data=raw_data_GC[which(raw_data_GC$PseudotimeThetaGC>=0 & raw_data_GC$PseudotimeThetaGC<2*pi),], 
               aes_string(x="PseudotimeThetaGC", y=gene_name), colour="grey") +
    geom_line(data=genes_kernel_smoother_data_frame_all[which(genes_kernel_smoother_data_frame_all$Gene==gene_name),], aes(x=PseudotimeTheta, y=GeneExpr), size=1.5) +
    #ggtitle("Kernel Smoother on the CD83 ") +
    expand_limits(y=limits_allFigs) +
    xlab("Theta GC") +
    ylab("Et") +
    ggtitle(gene_name)
  print(Physio_BD1)
}




# ##########
# ##########
# Fig2_i
# ##########
# ##########

# ####################################
# Kernel Smoother on CD83 and on CXCR4 (and trigger pulse width)
# ####################################
kernel_smoother_bandwith=0.3
kernel_smoother_gridsize=401*2 # see Arguments in function locpoly() in {KernSmooth}
pseudotime_used_AsString = "PseudotimeThetaGC"

kernel_smoother_list_CD83=locpoly(x=raw_data_to_kernel[,pseudotime_used_AsString], y=raw_data_to_kernel[,"asinh_CD83"], bandwidth = kernel_smoother_bandwith, gridsize=kernel_smoother_gridsize)
kernel_smoother_list_CXCR4=locpoly(x=raw_data_to_kernel[,pseudotime_used_AsString], y=raw_data_to_kernel[,"asinh_CXCR4"], bandwidth = kernel_smoother_bandwith, gridsize=kernel_smoother_gridsize)
kernel_smoother_list_TPW=locpoly(x=raw_data_to_kernel[,pseudotime_used_AsString], y=raw_data_to_kernel[,"TriggerPulseWidth"], bandwidth = kernel_smoother_bandwith, gridsize=kernel_smoother_gridsize)

kernel_smoother_df_metadata=data.frame(PseudotimeThetaGC=kernel_smoother_list_CD83[[1]], 
                                       asinh_CD83=kernel_smoother_list_CD83[[2]], 
                                       asinh_CXCR4=kernel_smoother_list_CXCR4[[2]],
                                       TriggerPulseWidth =kernel_smoother_list_TPW[[2]] )
kernel_smoother_df_metadata=kernel_smoother_df_metadata[which(kernel_smoother_df_metadata$PseudotimeThetaGC>=0 & kernel_smoother_df_metadata$PseudotimeThetaGC<2*pi),]

# Create the article Fig2i
# #############################
gg_plot=ggplot(raw_data_GC, aes_string(x="PseudotimeThetaGC", y="asinh_CD83")) +
  geom_point(colour="grey") +
  geom_line(data=kernel_smoother_df_metadata, size=1.5) +
  xlab("Theta GC") +
  ylab("CD83 (Index)") 
print(gg_plot)

# Create the article Fig2i
# #############################
gg_plot=ggplot(raw_data_GC, aes_string(x="PseudotimeThetaGC", y="asinh_CXCR4")) +
  geom_point(colour="grey") +
  geom_line(data=kernel_smoother_df_metadata, size=1.5) +
  xlab("Theta GC") +
  ylab("CXCR4 (Index)") 
print(gg_plot)


# ##########
# ##########
# FigSup2_h
# ##########
# ##########

# Create the article Fig2i
# #############################
raw_data_for_TPWplot = raw_data_GC
raw_data_for_TPWplot[,"CCNB1"]=sapply(raw_data_for_TPWplot[,"CCNB1"], function(x) if(x==0){x=NA}else{x})

gg_plot = ggplot(raw_data_for_TPWplot, aes(x=PseudotimeThetaGC, y=TriggerPulseWidth )) +
  geom_point(aes(color=CCNB1)) +
  scale_colour_gradient(low =TWO_COLORS_VECTOR[1], high = TWO_COLORS_VECTOR[2], limits=limits_allFigs) +
  geom_line(data=kernel_smoother_df_metadata, aes(x=PseudotimeThetaGC, y=TriggerPulseWidth), size=1.5) +
  xlab("Theta GC") +
  theme(legend.title = element_text(face = "italic"))
print(gg_plot)

# ##########
# ##########
# Fig2_j
# ##########
# ##########

# Count the number of Single-Cells in each Zones (LZ/IZ1/DZ/IZ2)
# ####################################################################
raw_data_GC$Zone = raw_data_GC$PseudotimeThetaGC
raw_data_GC$Zone=sapply(raw_data_GC$Zone, FUN=function(x) if((x<=LZ_end) | (x>IZ2_end) ){x="LZ"}else{x=x})
raw_data_GC$Zone=sapply(raw_data_GC$Zone, FUN=function(x) if((x<=IZ1_end) & (x>LZ_end) ){x="IZ1"}else{x=x})
raw_data_GC$Zone=sapply(raw_data_GC$Zone, FUN=function(x) if((x<=DZ_end) & (x>IZ1_end) ){x="DZ"}else{x=x})
raw_data_GC$Zone=sapply(raw_data_GC$Zone, FUN=function(x) if((x<=IZ2_end) & (x>DZ_end) ){x="IZ2"}else{x=x})
raw_data_GC$Zone=factor(raw_data_GC$Zone)


gg_plot=ggplot(raw_data_GC, aes(asinh_CXCR4, asinh_CD83, color=Zone)) +
  geom_point() +
  xlab("CXCR4 (Index)") +
  ylab("CD83 (Index)") 
print(gg_plot)

