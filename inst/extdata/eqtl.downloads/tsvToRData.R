suppressPackageStartupMessages({library(EnsDb.Hsapiens.v86)})
db <- EnsDb.Hsapiens.v86
#----------------------------------------------------------------------------------------------------
files <- list.files(pattern="tsv$")

for(file in files){
   tbl <- read.table(file, sep="\t", as.is=TRUE, nrow=-1)[, c(13, 14, 15, 16, 18, 4, 3, 9)]
   colnames(tbl) <- c("chrom", "pos", "ref", "alt", "rsid", "ensg", "pvalue", "beta")
   ensg <- unique(tbl$ensg)
   tbl.idMap <- select(db, keys=ensg, columns=c("SYMBOL"), keytype="GENEID")
   tbl.2 <- merge(tbl, tbl.idMap, by.x="ensg", by.y="GENEID")
   colnames(tbl.2)[9] <- "gene"
   coi <- c("chrom", "pos", "ref", "alt", "rsid", "pvalue", "beta", "gene", "ensg")
   tbl.eqtl <- tbl.2[, coi]
   rdata.filename <- sub(".tsv", ".RData", file, fixed=TRUE)
   save(tbl.eqtl, file=rdata.filename)
   printf("--- writing %d lines to %s", nrow(tbl.eqtl), rdata.filename)
   } # for file

