f <- "GTEx_ge_brain_hippocampus.purity_filtered.txt"
ndufs2 <- "ENSG00000158864"
tbl <-
    read.table(f, sep="\t", as.is=TRUE, header=TRUE, nrow=-1)
dim(tbl)
grep(ndufs2, tbl$cs_id)
