
DATA_FOLDER = "~/data/"
#DATA_FOLDER = ".../datafolder/"
physio_filename = "GSM3190075_humanBcells_matrixGeneCell.csv"
FL_filename = "GSM3190076_humanFL5_matrixGeneCell.csv"

# Load the raw data
physio_df = read.table( paste0(DATA_FOLDER, physio_filename), 
                       stringsAsFactors = FALSE, header = TRUE, sep=",") # or sep="\t"
rownames(physio_df) = physio_df$UniqueCellID
physio_df$UniqueCellID<-NULL
all_physio_genes_v = colnames(physio_df)

FL_df = read.table( paste0(DATA_FOLDER, FL_filename), 
                        stringsAsFactors = FALSE, header = TRUE, sep=",") # or sep="\t"
rownames(FL_df) = FL_df$UniqueCellID
FL_df$UniqueCellID<-NULL
all_FL_genes_v = colnames(FL_df)


# merging the data
t_physio_df = as.data.frame(as.matrix(t(physio_df)), stringsAsFactors=FALSE)
t_physio_df$Gene = rownames(t_physio_df)
t_FL_df = as.data.frame(t(FL_df), stringsAsFactors=FALSE)
t_FL_df$Gene = rownames(t_FL_df)
t_all_data_df = merge(t_physio_df,t_FL_df, by="Gene",all=TRUE, stringsAsFactors=FALSE)
rownames(t_all_data_df) = t_all_data_df$Gene
t_all_data_df$Gene<-NULL
t_all_data_df[is.na(t_all_data_df)]=0
allData_df = as.data.frame(t(t_all_data_df),stringsAsFactors=FALSE)
all_genes_v = colnames(allData_df)



