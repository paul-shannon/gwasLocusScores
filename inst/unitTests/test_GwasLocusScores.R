suppressPackageStartupMessages({
    library(RUnit)
    library(GwasLocusScores)
    })

message(sprintf("--- calling init.snpLocs()"))
init.snpLocs <- function() x <- EndophenotypeExplorer$new(NA, default.genome="hg19",
                                                          vcf.project="AMPAD", initialize.snpLocs=TRUE)

#----------------------------------------------------------------------------------------------------
# these tests are guided by the ADAMTS4 region's tag snp, and haploreg's report
# of variants in LD with it down to  R^2 >= 0.2,
tbl.tagHap <- read.table("~/github/ADvariantExplorer/explore/adamts4-study/haploreg-rs4575098-0.2.tsv",
                          sep="\t", as.is=TRUE, header=TRUE)
region.chrom <- tbl.tagHap$chrom[1]
region.start <- min(tbl.tagHap$hg38) - 1000
region.end <- max(tbl.tagHap$hg38) + 1000
tag.snp.loc <- 161185602
#----------------------------------------------------------------------------------------------------
runTests <- function()
{
    test_ctor()
    test_eqtls()
    test_breakMotifs()
    test_createTmsTable.runTrena()

} # runTests
#----------------------------------------------------------------------------------------------------
test_ctor <- function()
{
    message(sprintf("--- test_ctor"))

    gls <- GwasLocusScores$new("rs4575096", region.chrom, region.start, region.end,
                               "GTEx_V8.Brain_Hippocampus", targetGene="NDUFS2")
    checkTrue(all(c("GwasLocusScores", "R6") %in% class(gls)))

} # test_ctor
#----------------------------------------------------------------------------------------------------
test_eqtls <- function()
{

    checkException(
        gls <- GwasLocusScores$new("rs4575096", region.chrom, region.start, region.end,
                                   tissue.name="intentional.error", targetGene="NDUFS2"),
        silent=TRUE
        )


    gls <- GwasLocusScores$new("rs4575096", region.chrom, region.start, region.end,
                               tissue.name="GTEx_V8.Brain_Hippocampus", targetGene="NDUFS2")

    gls$load.eqtls()
    tbl.eqtl <- gls$get.tbl.eqtl(5e-4)
    checkTrue(nrow(tbl.eqtl) > 100)
    checkTrue(nrow(tbl.eqtl) < 200)
    # checkEquals(colnames(tbl.eqtl), c("rsid", "pvalue", "gene", "total.alleles", "beta", "id"))
    checkEquals(colnames(tbl.eqtl),
                c("chrom", "pos", "ref", "alt", "rsid", "pvalue", "beta", "gene", "ensg"))
    gene.count <- length(unique(tbl.eqtl$gene))
    checkTrue(gene.count > 10 & gene.count < 20)  # 14 on 30 dec 2021)

    tbl.eqtl <- gls$get.tbl.eqtl(1e-5)
    checkTrue(nrow(tbl.eqtl) > 30 & nrow(tbl.eqtl) <= 40) # 39 on 28 dec 2021

} # test_eqtls
#----------------------------------------------------------------------------------------------------
test_breakMotifs <- function()
{
    message(sprintf("--- test_breakMotifs"))

    gls <- GwasLocusScores$new("rs4575096", region.chrom, region.start, region.end,
                               tissue.name="GTEx_V8.Brain_Hippocampus",
                               targetGene="NDUFS2")

    gls$load.eqtls()
    tbl.eqtl <- gls$get.tbl.eqtl(pvalue.cutoff=1)
    dim(tbl.eqtl)
    targetGene <- "NDUFS2"
    pval.cutoff <- 1.22e-5
    tbl.sub <- subset(tbl.eqtl, pvalue <= pval.cutoff & gene==targetGene)
    snps <- unique(tbl.sub$rsid)
    dim(tbl.sub)
    checkTrue(nrow(tbl.sub) > 2)

    tfs.oi.for.speed <- c("SOX21", "ZEB1")  # rather than the full > 500 from MotifDb
    x <- system.time(tbl.breaks <- gls$breakMotifsAtEQTLs(targetGene, pval.cutoff,
                                                          TFs.preselected=tfs.oi.for.speed))
    message(sprintf("%4.1f minutes to break %d snps", round(x[["elapsed"]]/60, digits=1),
                    length(snps)))
    checkEquals(length(unique(tbl.breaks$SNP_id)), length(snps))
    checkTrue(nrow(subset(tbl.breaks, abs(pctDelta) > 0.10)) >= 2)
    # older, longer version is better.  was "save(tbl.breaks, file="tbl.breaks.fromTest.RData")"

} # test_breakMotifs
#----------------------------------------------------------------------------------------------------
test_createTmsTable.runTrena <- function()
{
    message(sprintf("--- test_createTmsTable.runTrena"))

    targetGene <- "NDUFS2"
    tissue <- "GTEx_V8.Brain_Hippocampus"
    gls <- GwasLocusScores$new("rs4575096", region.chrom, region.start, region.end,
                               tissue.name=tissue,
                               targetGene=targetGene)
    require(TrenaProjectAD)
    tpad <- TrenaProjectAD()

    require(EndophenotypeExplorer)
    etx <- EndophenotypeExplorer$new(targetGene,
                                     default.genome="hg38",
                                     vcf.project="ADNI",
                                     verbose=FALSE,
                                     initialize.snpLocs=TRUE)

    mtx.rna <- etx$get.rna.matrix(tissue)
    dim(mtx.rna)

    tbl.fimo.small <- get(load("~/github/gwasLocusScores/inst/extdata/tbl.fimo.ndufs2.small.RData"))
    data.dir <- "~/github/TrenaProjectAD/inst/extdata/genomicRegions"
    filename <- "boca-hg38-consensus-ATAC.RData"
    tbl.boca <- get(load(file.path(data.dir, filename)))

    tbl.tms <- gls$createTrenaMultiScoreTable(tpad,  tbl.fimo, tbl.boca, mtx.rna)
    tbl.tms.strong <- subset(tbl.tms, (chip | oc) & abs(cor.all) > 0.2 & fimo_pvalue < 1e-4)
    tfs <- unique(tbl.tms.strong$tf)
    length(tfs)

    tbl.trena <- gls$runTrena(tfs, mtx.rna, tbl.tms.strong)
    tbl.trena$rfNorm <- tbl.trena$rfScore / max(tbl.trena$rfScore)
    dim(tbl.tms.strong)
    top.tfs <- tbl.trena$gene [1:20]
    tbl.tms.stronger <- subset(tbl.tms.strong, tf %in% top.tfs)
    dim(tbl.tms.stronger)
    tbl.breaks <- get(load("tbl.breaks.fromTest.RData"))

    tbl.eqtl <- gls$get.tbl.eqtl(1e-4)
    if(is.null(tbl.eqtl)){
       gls$load.eqtls()
       tbl.eqtl <- gls$get.tbl.eqtl(1e-4)
       }

    tbl.trena.scored <- gls$scoreEQTLs(tbl.trena, tbl.breaks, tbl.eqtl)
    printf("--- summed breakage: %5.2f", sum(tbl.trena.scored$breakage.score))
    head(tbl.trena.scored, n=10)
    viz <- FALSE
    if(viz){
        igv <- start.igv("NDUFS2", "hg38")
        zoomOut(igv)
        zoomOut(igv)
              # NDUFS2 eqtls
        gls$load.eqtls()
        tbl.eqtl <- gls$get.tbl.eqtl(pvalue.cutoff=1)
        pval.cutoff <- 1e-3
        tbl.eqtl.sub <- subset(tbl.eqtl, pvalue <= pval.cutoff & gene==targetGene)
        dim(tbl.eqtl.sub)
        tbl.rsid.locs <- gls$getETX()$rsidToLoc(tbl.eqtl.sub$rsid)
        tbl.eqtls.filtered <- merge(tbl.eqtl.sub, tbl.rsid.locs[, c("chrom", "hg38", "rsid")], by="rsid")
        tbl.eqtls.filtered$score <- -log10(tbl.eqtls.filtered$pvalue) * tbl.eqtls.filtered$beta
        tbl.eqtls.filtered <- tbl.eqtls.filtered[, c("chrom", "hg38", "hg38", "score", "rsid")]
        colnames(tbl.eqtls.filtered) <- c("chrom", "start", "end", "score", "rsid")
        tbl.eqtls.filtered$start <- tbl.eqtls.filtered$start - 1
        track <- DataFrameQuantitativeTrack("NDUFS2 eQTLs", tbl.eqtls.filtered, autoscale=TRUE, color="brown")
        displayTrack(igv, track)
              # tag snp and LD parnters
        tbl.track <- tbl.tagHap[, c("chrom", "hg38", "hg38", "rSquared")]
        colnames(tbl.track)[2:3] <- c("start", "end")
        tbl.track$start <- tbl.track$start - 1
        track <- DataFrameQuantitativeTrack("gwas+LD", tbl.track, autoscale=TRUE, color="darkred")
        displayTrack(igv, track)
              # broken motifs
        tbl.breaks <- get(load("tbl.breaks.fromTest.RData"))
        tbl.breaks.sub <- subset(tbl.breaks, geneSymbol %in% top.tfs) # & abs(pctDelta) < -0.1)
        dim(tbl.breaks.sub)
        tbl.breaks.sub <- subset(tbl.breaks.sub, abs(pctDelta) > 0.15)
        dim(tbl.breaks.sub)
        tbl.track <- tbl.breaks[, c("seqnames", "start", "end", "pctDelta", "geneSymbol")]
        colnames(tbl.track)[1] <- "chrom"
        tbl.track$chrom <- as.character(tbl.track$chrom)
        tbl.track$start <- tbl.track$start - 1
        for(broken.tfbs in tbl.breaks.sub$geneSymbol){
           tbl.track.tf <- subset(tbl.track, geneSymbol==broken.tfbs)
           track.title <- sprintf("%s.breaks", broken.tfbs)
           track <- DataFrameQuantitativeTrack(track.title, tbl.track.tf, autoscale=FALSE,
                                               color="random", min=-0.3, max=0.3)
           displayTrack(igv, track)
           } # for broken.tfbs
        } # if viz

} # test_createTmsTable.runTrena
#----------------------------------------------------------------------------------------------------
scoreEQTLs = function(tbl.trena, tbl.breaks, tbl.eqtl)
{
    tfs.oi <- tbl.trena$gene
    tbl.eqtl.sub <- subset(tbl.eqtl, gene==targetGene)

    tbl.eqtl.sub$sig.beta <- with(tbl.eqtl.sub, -log10(pvalue) * abs(beta))

       # for each tf:
       #   calculate a trena score, the rfNorm
       #   for every broken tfbs for this tf
       #      get the snps abs(pctDelta)
       #      multiply it by the eqtls sig.beta
       #      add it to the tf's breakage score

    breakage.scores <- list()
    for(tf in tfs.oi){
       tbl.tf.breaks.sub <- subset(tbl.breaks, geneSymbol == tf & SNP_id %in% tbl.eqtl$rsid)
       dups <- which(duplicated(tbl.tf.breaks.sub[, c("SNP_id", "geneSymbol")]))
       if(length(dups) > 0)
           tbl.tf.breaks.sub <- tbl.tf.breaks.sub[-dups,]
       tf.breakage.score <- 0
       for(breaking.snp in tbl.tf.breaks.sub$SNP_id){
          pctDelta <- abs(subset(tbl.tf.breaks.sub, SNP_id==breaking.snp)$pctDelta)
            # multiply this by the tf's rfNorm and the eQTLs sig.beta
          sig.beta <- subset(tbl.eqtl.sub, rsid==breaking.snp)$sig.beta
          tf.rfNorm <- subset(tbl.trena, gene==tf)$rfNorm
          new.increment <- (pctDelta + sig.beta) * tf.rfNorm
          tf.breakage.score <- tf.breakage.score + new.increment
          } # for breaking.snp
       printf("--- %s: %5.2f", tf, tf.breakage.score)
       breakage.scores[[tf]] <- tf.breakage.score
       } # for tf

    tbl.trena$breakage.score <- as.numeric(breakage.scores)
    tbl.trena$rfNorm <- round(tbl.trena$rfNorm, digits=2)
    tbl.trena$breakage.score <- round(tbl.trena$breakage.score, digits=2)

    return(tbl.trena)

} # scoreEQTLs
#----------------------------------------------------------------------------------------------------
if(!interactive())
    runTests()
