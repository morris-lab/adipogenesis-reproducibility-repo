#redirecting the path for commonly shared R packages
.libPaths("/ref/smlab/software/r-envs/common_libs/curr/")

library(Seurat)
library(Capybara)

key <- "23_11_03_super_ref_on_p21"

#loading in new data seurat object and seurat object for reference
query <- readRDS("./p21_processed.rds")
reference <- readRDS("./23_11_03_capy_super_ref_5_datasets_na_rm_subset.rds")

#extracting counts from SO, and converting df to matrix
counts <- as.matrix(query@assays$RNA@counts)

print("starting background qp")
#background qp scoring
single.round.QP.analysis(reference[[3]], reference[[1]], n.cores = 24, save.to.path = "./capy_results/",
                         save.to.filename = paste0(key, "_reference.background.scores"), unix.par = TRUE)


print("starting qp")
#running capybara scoring
single.round.QP.analysis(reference[[3]], counts, n.cores = 24, save.to.path = "./capy_results/",
                         save.to.filename = paste0(key, "_reference.scores"), unix.par = TRUE)
#print("starting background qp")
#background qp scoring
#single.round.QP.analysis(reference[[3]], reference[[1]], n.cores = 24, save.to.path = "./capy_results/",
#                         save.to.filename = paste0(key, "_reference.background.scores"), unix.par = TRUE)

# Read in background and testing identity scores
background.mtx <- read.csv(paste0("/scratch/smlab/ebutka/capy_results/", key, "_reference.background.scores_scale.csv"), header = T, row.names = 1, stringsAsFactors = F)
mtx.test <- read.csv(paste0("/scratch/smlab/ebutka/capy_results/", key, "_reference.scores_scale.csv"), header = T, row.names = 1, stringsAsFactors = F)

col.sub <- ncol(background.mtx) - 2

# Conduct reference randomization to get empirical p-value matrix
ref.perc.list <- percentage.calc(background.mtx[,c(1:col.sub)], background.mtx[,c(1:col.sub)])

# Conduct test randomization to get empirical p-value matrix
perc.list <- percentage.calc(as.matrix(mtx.test[,c(1:col.sub)]), as.matrix(background.mtx[,c(1:col.sub)]))

# Binarization of inference results
bin.count <- binarization.mann.whitney(mtx = mtx.test[,c(1:col.sub)], ref.perc.ls = ref.perc.list, ref.meta = reference[[2]], perc.ls = perc.list)

#classification
classification <- binary.to.classification(bin.count[,c(1:col.sub)])
rownames(classification) <- classification$barcode

#cross-referencing the order and nomenclature of the 10x barcodes between capy output and SO
#output needs to be TRUE
identical(colnames(query), classification$barcode)

#adding in capyresults as a metadata column into the SO
query@meta.data[, paste0(key, "_capy_call")] <- classification$call

# hybrids
bin.count.rowsums <- apply(bin.count, 1, sum)
rownames(bin.count)[which(bin.count.rowsums > 1)]
bin.count.greaterthan1 <- apply(bin.count[which(bin.count.rowsums > 1), ], c(1, 2), function(x) {x > 0})
multi.celltypes <- apply(bin.count.greaterthan1, 1, function(x) {print(colnames(bin.count.greaterthan1)[x])})

hybrids <- c()
for(i in 1:length(multi.celltypes)) {
  #print(multi.celltypes[i])
  hybrids <- c(hybrids, paste(multi.celltypes[[i]], collapse = "-"))
}

hybrids <- gsub("frxn_cell.type_", "", hybrids)

# p-values for classifications
pval_avg_mat <- data.frame()
for(i in 1:length(perc.list)) {
  pval_avg_mat <- rbind(pval_avg_mat, apply(perc.list[[i]][[1]], 2, mean))
  rownames(pval_avg_mat)[i] <- names(perc.list[[i]])
}
colnames(pval_avg_mat) <- names(apply(perc.list[[1]][[1]], 2, mean))

hyb.df <- data.frame("x" = rownames(bin.count)[which(bin.count.rowsums > 1)], "hybrid" = hybrids)


#query@meta.data[, paste0(key, "_hybrid")] <- NA
#for(i in 1:dim(hyb.df)[1]) {
#  query@meta.data[, paste0(key, "_hybrid")][which(colnames(query) == hyb.df$x[i])] <- hyb.df$hybrid[i]
#}


#saveRDS(query, paste0(key, "_updated_query_seurat.rds"))
write.csv(pval_avg_mat, paste0(key, "_p_values.csv"))
write.csv(hyb.df, paste0(key, "_hybrids.csv"))
write.csv(classification, paste0(key, "_classification.csv"))

