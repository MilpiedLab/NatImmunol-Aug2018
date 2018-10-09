# ###############################################
# This script aims to load the scqpcr data
# ###############################################


DATA_FOLDER = "~/data/"

raw_data_filename = "41590_2018_181_MOESM3_ESM.csv"

# Load the raw data
raw_data = read.table( paste0(DATA_FOLDER, raw_data_filename), 
                       stringsAsFactors = FALSE, header = TRUE, sep=",") # or sep="\t"
rownames(raw_data) = raw_data$UniqueCellID
gene_names_v = colnames(raw_data)[14:104]



raw_data$asinh_CXCR4 = asinh(raw_data$GC_FL_Comp.CXCR4.APC/100)
raw_data$asinh_CD83 = asinh(raw_data$GC_FL_Comp.CD83.PEDazzle594/10)

raw_data_bulk10 = raw_data[raw_data$CellNumber==10,]
# remove bulk wells with no expression
raw_data_bulk10 = raw_data_bulk10[rowSums(raw_data_bulk10[,gene_names_v])!=0,]

# GC only
# ################
raw_data_GC = raw_data[raw_data$SimpleSortPheno=="GC" &
                         raw_data$QCfilters=="pass",]

# PB_PC only
# ################
raw_data_PB = raw_data[raw_data$SimpleSortPheno=="PB_PC" &
                         raw_data$QCfilters=="pass",]

# GC PB Mem (no FL cells)
# ################
raw_data_GC_PB_Mem = raw_data[raw_data$SimpleSortPheno!="FL" &
                         raw_data$QCfilters=="pass",]


# FL only
# ################
raw_data_FL = raw_data[raw_data$SimpleSortPheno=="FL" &
                         raw_data$QCfilters=="pass",]


raw_data = raw_data[raw_data$QCfilters=="pass",]

# #####################
# #####################
# #####################
# some constants
# #####################
# #####################
# #####################

# ggplot parameters
# #################
color_lines = "gray50"
size_geom_point = 2.8 # tSNE colored by genes
color_NA_geom_point = "gray80" # for na.value in scale_colour_gradient()

# for GC plot CXCR4 vs CD83
gate1 = data.frame(x=c(-0.1,1.5,4.2,4.2,-0.1,-0.1), y=c(1.7,1.7,4.5,5.5,5.5,1.7))
gate2 = data.frame(x=c(1.5,1.5,4.2,4.2), y=c(1.7,-1.8,-1.8,4.5))


limits_allFigs = c(0,25.05)


# Zones for Kernel Smoother plots
limits_allFigs = c(0,25.05)
limits_ggplot=c(0,22.9)
radian_sup = 0
LZ_end =0.8+radian_sup#1.18
IZ1_end =1.98+radian_sup
DZ_end =3.9 #3.9 +radian_sup # 3.64 #3.88  
IZ2_end = 5.25+radian_sup #5.28
zones_data_frame=data.frame(x1=c(0,LZ_end,IZ1_end,DZ_end,IZ2_end), x2=c(LZ_end,IZ1_end,DZ_end,IZ2_end,pi*2), y1=c(0,0,0,0,0), y2=c(22.8,22.8,22.8,22.8,22.8), Zone=c("LZ","IZ1","DZ","IZ2","LZ"), r=c("LZ","IZ1","DZ","IZ2","LZ"))
zones_segment_data_frame = data.frame(xbegin=c(LZ_end,IZ1_end,DZ_end,IZ2_end), ybegin=c(0,0,0,0), xend=c(LZ_end,IZ1_end,DZ_end,IZ2_end), yend=c(22.8,22.8,22.8,22.8))

# ############################################################################################################################################
# calcul the Kernel Smoother of each gene (for Fig2 and Fig4)
# ############################################################################################################################################

pseudotime_used_AsString="PseudotimeThetaGC" 
kernel_smoother_bandwith=0.3

# because we are dealing with circular data:
# we copy the values between [max_theta-1; max_theta] and put then at [min_theta-1;min_theta]
copy_max_pseudotime_data_frame=raw_data_GC[which(raw_data_GC[,pseudotime_used_AsString]>=(2*pi-1)),]
copy_max_pseudotime_data_frame[,pseudotime_used_AsString]=copy_max_pseudotime_data_frame[,pseudotime_used_AsString]-2*pi
# we copy the values between [min_theta; min_theta+1] and put then at [max_theta+1;max_theta]
copy_min_pseudotime_data_frame=raw_data_GC[which(raw_data_GC[,pseudotime_used_AsString]<=(0+1)),]
copy_min_pseudotime_data_frame[,pseudotime_used_AsString]=copy_min_pseudotime_data_frame[,pseudotime_used_AsString]+2*pi

raw_data_to_kernel=rbind(raw_data_GC,copy_max_pseudotime_data_frame,copy_min_pseudotime_data_frame)



#kernel_smoother_data_frame_list=list()
genes_kernel_smoother_data_frame_all=data.frame()

kernel_smoother_coord_dataframe = data.frame()

for (j in 1:length(gene_names_v)){
  genes_kernel_smoother_list=locpoly(x=raw_data_to_kernel[,pseudotime_used_AsString], y=raw_data_to_kernel[,gene_names_v[j]], bandwidth = kernel_smoother_bandwith)
  genes_kernel_smoother_data_frame=data.frame(PseudotimeTheta=genes_kernel_smoother_list[[1]], GeneExpr=genes_kernel_smoother_list[[2]])
  genes_kernel_smoother_data_frame$Gene=gene_names_v[j]
  genes_kernel_smoother_data_frame$Gene2=gene_names_v[j]
  # after smoothing, we only keep data between 0 and 2*pi
  genes_kernel_smoother_data_frame=genes_kernel_smoother_data_frame[which(genes_kernel_smoother_data_frame$PseudotimeTheta>=0 & genes_kernel_smoother_data_frame$PseudotimeTheta<=2*pi),]
  genes_kernel_smoother_data_frame_all=rbind(genes_kernel_smoother_data_frame_all,genes_kernel_smoother_data_frame)
  kernel_smoother_coord_dataframe = rbind(kernel_smoother_coord_dataframe,genes_kernel_smoother_data_frame$GeneExpr)
}
kernel_smoother_coord_dataframe = as.data.frame(t(kernel_smoother_coord_dataframe))
colnames(kernel_smoother_coord_dataframe) = gene_names_v
rownames(kernel_smoother_coord_dataframe) =1:nrow(kernel_smoother_coord_dataframe)
kernel_smoother_coord_dataframe$PseudotimeTheta = genes_kernel_smoother_list[[1]][which(genes_kernel_smoother_list[[1]]>=0 & genes_kernel_smoother_list[[1]]<=2*pi)]

