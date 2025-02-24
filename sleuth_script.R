library("sleuth")

stab = read.table("sleuth_input.txt",header=TRUE)

so = sleuth_prep(stab)

#fit a model comparing the two conditions 
so = sleuth_fit(so, ~condition, 'full')

#fit the reduced model to compare in the likelihood ratio test 
so = sleuth_fit(so, ~1, 'reduced')

#perform the likelihood ratio test for differential expression between conditions 
so = sleuth_lrt(so, 'reduced', 'full')

#load the dplyr package for data.frame filtering
library(dplyr)

#extract the test results from the sleuth object 
sleuth_table = sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE) 

#filter most significant results (FDR/qval < 0.05) and sort by pval
sleuth_significant = dplyr::filter(sleuth_table, qval <= 0.05) |> dplyr::arrange(pval) 

#print top 10 transcripts
sleuth_results <- sleuth_significant[,c("target_id", "test_stat", "pval", "qval")]

#write FDR < 0.05 transcripts to file
write.table(sleuth_results, file="PipelineProject.log", sep = "\t", quote = FALSE,row.names = FALSE, append = TRUE)
