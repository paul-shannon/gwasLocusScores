library(ADvariantExplorer)
library(RUnit)

targetGene <- NA
  # the region spanned by haploreg's report of variants in LD, R^2 >= 0.2, with tag.snp
loc.chrom <- "chr1"
loc.start <- 161088961
loc.end   <- 161240513
span <- 1 + loc.end - loc.start
printf("span: %5.2fk", span/1000)   # 151k

avx <- ADvariantExplorer$new(NA, loc.chrom, loc.start, loc.end)
tbl.summary <- avx$geteQTLSummary()
colnames(tbl.summary)
grep("cerebellum", tbl.summary$unique_id, ignore.case=TRUE, value=TRUE)
 # [1] "GTEx.brain_hippocampus"    "GTEx.brain_hippocampus"    "GTEx.brain_hippocampus"
 # [4] "GTEx.brain_hippocampus"    "GTEx_V8.Brain_Hippocampus"

studies <- c("GTEx_V8.Brain_Hippocampus", "GTEx_V8.Brain_Cerebellum",
             "GTEx_V8.Brain_Cortex")

tbl.hcCerCtx <- avx$geteQTLsByLocationAndStudyID(loc.chrom, loc.start, loc.end,
                                           studies, method="tabix", simplify=TRUE)
dim(tbl.hcCerCtx) # 53469     6
save(tbl.hcCerCtx, file="tbl.hcCerCtx.151kspan.rs4575098.RData")

colnames(tbl.hc)
dim(tbl.hc)   # 17513 6
save(tbl.hc, file="tbl.hc.151kspan.rs4575098.RData")

     #--------------------------------------------------
     # quick score of NDUFS2 and PPOX
     #--------------------------------------------------

tbl.ndufs2.hc <- subset(tbl.hcCerCtx, id==studies[1] & gene=="NDUFS2")
tbl.ndufs2.hc$score <- -log10(tbl.ndufs2.hc$pvalue)
sum(tbl.ndufs2.hc$score)   # 317
sum(with(tbl.ndufs2.hc, score * abs(beta))) # [1] 64.26891

tbl.ppox.hc <- subset(tbl.hcCerCtx, id==studies[1] & gene=="PPOX")
tbl.ppox.hc$score <- -log10(tbl.ppox.hc$pvalue)
sum(tbl.ppox.hc$score)   # 335
sum(with(tbl.ppox.hc, score * abs(beta))) # [1] 116.2766


tbl.tstd1.hc <- subset(tbl.hcCerCtx, id==studies[1] & gene=="TSTD1")
tbl.tstd1.hc$score <- -log10(tbl.tstd1.hc$pvalue)
sum(tbl.tstd1.hc$score)   # 148
sum(with(tbl.tstd1.hc, score * abs(beta))) # [1] 30

tbl.adamts4.hc <- subset(tbl.hcCerCtx, id==studies[1] & gene=="ADAMTS4")
tbl.adamts4.hc$score <- -log10(tbl.adamts4.hc$pvalue)
sum(tbl.adamts4.hc$score)   # 141
sum(with(tbl.adamts4.hc, score * abs(beta))) # [1] 21

tbl.rest.hc <- subset(tbl.hcCerCtx, id==studies[1] & gene=="REST")
tbl.rest.hc$score <- -log10(tbl.rest.hc$pvalue)
sum(tbl.rest.hc$score)   # 0
sum(with(tbl.rest.hc, score * abs(beta))) # [1] 0

tbl.b4galt3.hc <- subset(tbl.hcCerCtx, id==studies[1] & gene=="B4GALT3")
tbl.b4galt3.hc$score <- -log10(tbl.b4galt3.hc$pvalue)
sum(tbl.b4galt3.hc$score)   # 117
sum(with(tbl.b4galt3.hc, score * abs(beta))) #  16

goi <- c("PPOX","NDUFS2","KLHDC9","NIT1","HSPA7","HSPA6","FCER1G","RP11-312J18.5","TSTD1")
for(GENE in goi){
   tbl.gene.hc <- subset(tbl.hcCerCtx, id==studies[3] & gene==GENE)
   eqtl.count <- length(unique(subset(tbl.gene.hc, pvalue <= cutoff)$rsid))
   tbl.gene.hc$score <- -log10(tbl.gene.hc$pvalue)
   pval.score.sum <- sum(tbl.gene.hc$score)
   pval_beta.score.sum <- sum(with(tbl.gene.hc, score * abs(beta)))
   printf("--- %s: %d %5.2f   %5.2f", GENE, eqtl.count, pval.score.sum, pval_beta.score.sum)
   }



cutoff <- 0.001
as.data.frame(sort(table(subset(tbl.hc, pvalue <= cutoff)$gene), decreasing=TRUE))


tbl.ppox <- subset(tbl.hc,

tbl.all <- avx$getFullGwasTable(trim.columns=FALSE)

tbl.ndufs2 <-
   dim(
    subset(tbl.all,
           MAPPED_GENE=="NDUFS2" &
           P.VALUE <= 1e-3 &
           abs(OR.or.BETA) > 0.0)
    )


    tbl.cat <- avx$geteQTLSummary()
    dim(tbl.cat) # 498 12
    checkTrue(nrow(tbl.cat) > 390)
    checkEquals(ncol(tbl.cat), 12)

       # find study ids like this
    sort(unique(tbl.cat$unique_id))   # 112 of them (24 nov 2021)
    sort(unique(grep("macrophage_naive", tbl.cat$unique_id, value=TRUE)))

    study.1 <- "Alasoo_2018.macrophage_naive"
    study.2 <- "Nedelec_2016.macrophage_naive"
    checkTrue(all(c(study.1, study.2) %in% grep("macrophage_naive", tbl.cat$unique_id, value=TRUE)))
    chrom <- "8"
    start <- 27603335
    end   <- 27608281
    1 + end - start
    tbl.1 <- avx$geteQTLsByLocationAndStudyID(chrom, start, end, study.1, simplify=TRUE)
    dim(tbl.1)
