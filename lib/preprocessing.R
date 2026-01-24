#EDA
library(SummarizedExperiment)
colData(data_se)
assays(data_se)

count_matrix <- assay(data_se) #just numbers (genes x samples_matrix)
head(count_matrix)

preprocess_cancer_data <- function(cancer_name){

}

