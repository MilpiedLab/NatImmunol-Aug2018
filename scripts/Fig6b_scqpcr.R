# ###################################
# Fig6b (scqpcr)
# ###################################

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
library(Biostrings)


list_genes=c("ATR","ATM","CCND3","EP300","CD27","AICDA","BACH2","MME","CD38","GCSAM","TNFRSF17","IL4R","CREBBP","BCL6",
             "PAX5","IRF8","CXCR4","EZR","GNAI2","AURKA","CCNB1","EZH2","DNMT1","MKI67","RAD51","SMARCA4",
             "POLH","TP53","TCF3","P2RY8","SLAMF6","KMT2C","MYB","PLCG1","PTEN","RUNX1","MEF2B","FOXO1",
             "SEMA4B","CXCR5","GNA13","CCR7","IRF4","IL10RA","PRDM1","CD74","CD86","EGR2","IL21R","CD80","BCL2A1",
             "CD83","GPR183","HLA.DOA","CIITA","HLA.DRA","SLAMF1","SEMA7A","CD40","STAT6","CD72","PTPN6","NFKB2",
             "IFNGR1","FAS","ICAM1","TNFAIP3","KMT2D","LY75","STAT3","LMO2","SELL")


raw_data_to_study = raw_data
rownames(raw_data_to_study)=paste(raw_data_to_study$SimpleSortPheno,rownames(raw_data_to_study),sep="_")

coExpressionGenes_cell_list = function_coexpr_gene(raw_data_to_study[,list_genes])

coExpressionGenes_cell_matrix = matrix(unlist(coExpressionGenes_cell_list), nrow = length(coExpressionGenes_cell_list), byrow = TRUE)
rownames(coExpressionGenes_cell_matrix) = names(coExpressionGenes_cell_list)
coExpressionGenes_cell_matrix=coExpressionGenes_cell_matrix[grep(c("GC|FL"), rownames(coExpressionGenes_cell_matrix)),]


# index of GC cells:
index_GC = grep(c("GC"), rownames(coExpressionGenes_cell_matrix))


# K-FOLD CROSS VALIDATION 
# ########################
k=5
index_GC_v = index_GC
index_GC_remaining_list = list()
nbrForPWM = as.integer(length(index_GC)/k)
for (i in 1:k){
  if(i!=k){
    index_GC_remaining_new = sample(index_GC_v, nbrForPWM)
    index_GC_v = index_GC_v[which(!index_GC_v %in% index_GC_remaining_new)]
  }else{
    index_GC_remaining_new = index_GC_v
  }
  index_GC_remaining_list[[i]] = index_GC_remaining_new
}


score_vector = c()
sample_for_comparison_vector = c()
for(i in 1:k){
  index_GC_remaining=index_GC_remaining_list[[i]]
  index_GC_forPWM = index_GC[which(!index_GC %in% index_GC_remaining)]
  
  Kfold_GC_coExpressionGenes_cell_matrix = coExpressionGenes_cell_matrix[index_GC_forPWM,]
  Kfold_PWM_GC_log = PWM(consensusMatrix( Kfold_GC_coExpressionGenes_cell_matrix ),type="log2probratio") 
  
  index_GC_FL_remaining = 1:nrow(coExpressionGenes_cell_matrix)
  index_GC_FL_remaining = index_GC_FL_remaining[which(! index_GC_FL_remaining %in% index_GC_forPWM)]
  for(j in 1:length(index_GC_FL_remaining) ){
    print(j)
    score = 0
    score=sum(sapply(seq_along(coExpressionGenes_cell_matrix[j,]), function(x){Kfold_PWM_GC_log[coExpressionGenes_cell_matrix[index_GC_FL_remaining[j],x],x]}))
    score_vector = c(score_vector, score)
    sample_for_comparison_vector = c(sample_for_comparison_vector, strsplit(rownames(coExpressionGenes_cell_matrix)[index_GC_FL_remaining[j]],"_")[[1]][2])
  }
  
}

score_df = data.frame(SampleForComparison = sample_for_comparison_vector, Score = score_vector)
density_plot = ggplot(score_df, aes(Score))+
  geom_line(stat="density",aes(colour=SampleForComparison))+
  #scale_colour_manual(values = COLOR_SAMPLES) +
  scale_x_continuous(breaks=seq(0, 1, 0.1), limits = c(0,1))+
  ggtitle(paste("PWM log2probratio GC "))
print(density_plot)
