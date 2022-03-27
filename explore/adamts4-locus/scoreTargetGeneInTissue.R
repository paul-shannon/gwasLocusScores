suppressPackageStartupMessages({
    library(GwasLocusScores)
    library(TrenaProjectAD)
    library(EndophenotypeExplorer)
    })
#----------------------------------------------------------------------------------------------------
tbl.tagHap <- read.table("~/github/ADvariantExplorer/explore/adamts4-study/haploreg-rs4575098-0.2.tsv",
                          sep="\t", as.is=TRUE, header=TRUE)
region.chrom <- tbl.tagHap$chrom[1]
region.start <- min(tbl.tagHap$hg38) - 1000
region.end <- max(tbl.tagHap$hg38) + 1000
tag.snp.loc <- 161185602
tpad <- TrenaProjectAD()
#----------------------------------------------------------------------------------------------------
scoreTargetGeneInTissue <- function(targetGene, tag.snp, tbl.fimo, tbl.oc, gtex.tissue)
{
   gls <- GwasLocusScores$new("rs4575096", region.chrom, region.start, region.end,
                              tissue.name=gtex.tissue,
                              targetGene=targetGene)

   etx <- EndophenotypeExplorer$new(targetGene, "hg38", vcf.project="ADNI")
   mtx.rna <- etx$get.rna.matrix(gtex.tissue)
   dim(mtx.rna)

   #tbl.fimo.small <- get(load("~/github/gwasLocusScores/inst/extdata/tbl.fimo.ndufs2.small.RData"))
   # data.dir <- "~/github/TrenaProjectAD/inst/extdata/genomicRegions"
   # filename <- "boca-hg38-consensus-ATAC.RData"
   # tbl.boca <- get(load(file.path(data.dir, filename)))

   tbl.tms <- gls$createTrenaMultiScoreTable(tpad,  tbl.fimo, tbl.oc, mtx.rna)
   tbl.tms.strong <- subset(tbl.tms, (chip | oc) & abs(cor.all) > 0.2 & fimo_pvalue < 1e-4)
   tfs <- unique(tbl.tms.strong$tf)
   length(tfs)

   tbl.trena <- gls$runTrena(tfs, mtx.rna, tbl.tms.strong)
   tbl.trena$rfNorm <- tbl.trena$rfScore / max(tbl.trena$rfScore)
   dim(tbl.tms.strong)
   top.tfs <- tbl.trena$gene [1:20]
   tbl.tms.stronger <- subset(tbl.tms.strong, tf %in% top.tfs)
   dim(tbl.tms.stronger)

   tbl.eqtl <- gls$get.tbl.eqtl(1e-4)
   if(is.null(tbl.eqtl)){
      gls$load.eqtls()
      tbl.eqtl <- gls$get.tbl.eqtl(1e-4)
      }

   printf("starting motifbreakR...")
   x <- system.time(tbl.breaks <- gls$breakMotifsAtEQTLs(targetGene, pvalue.cutoff=1e-2))
   print(x[["elapsed"]])
   #tbl.breaks <- get(load("~/github/gwasLocusScores/inst/unitTests/tbl.breaks.fromTest.RData"))


  tbl.trena.scored <- gls$scoreEQTLs(tbl.trena, tbl.breaks, tbl.eqtl)
  filename <- sprintf("scoredInTissue.%s.%s.%s", targetGene, tissue, tag.snp)
  save(tbl.breaks, tbl.tms, tbl.tms.strong, tfs, tbl.trena, tbl.tms.stronger, tbl.trena.scored,
       tbl.fimo, tbl.eqtl, file=filename)
  return(tbl.trena.scored)

} # scoreTargetGeneInTissue
#----------------------------------------------------------------------------------------------------
# <<<<<<< HEAD
#tbl.fimo.ndufs2 <- get(load("~/github/TrenaProjectAD/prep/bigFimo/from-khaleesi/tbl.fimo.NDUFS2.RData"))
##tbl.fimo <- get(load("~/github/TrenaProjectAD/prep/bigFimo/from-khaleesi/tbl.fimo.PPOX.RData"))
#tbl.fimo.ppox <- subset(tbl.fimo, start > region.start & end < region.end)
#data.dir <- "~/github/TrenaProjectAD/inst/extdata/genomicRegions"
#filename <- "boca-hg38-consensus-ATAC.RData"
#tbl.boca <- get(load(file.path(data.dir, filename)))
#tissue <- "GTEx_V8.Brain_Hippocampus"
#targetGene <- "NDUFS2"
#tbl.trena.scored <- scoreTargetGeneInTissue(targetGene, "rs4575098", tbl.fimo.ndufs2, tbl.boca, tissue)
=======
#tbl.fimo <- get(load("~/github/TrenaProjectAD/prep/bigFimo/from-khaleesi/tbl.fimo.NDUFS2.RData"))
#tbl.fimo <- get(load("~/github/TrenaProjectAD/prep/bigFimo/from-khaleesi/tbl.fimo.PPOX.RData"))
targetGenes <- c("B4GALT3", "FCER1G", "HSPA6", "KLHDC9", "NDUFS2", "NIT1", "PPOX", "TSTD1")
tissues <- c(#"GTEx_V8.Brain_Cerebellar_Hemisphere",
             #"GTEx_V8.Brain_Cerebellum",
             "GTEx_V8.Brain_Cortex",
             "GTEx_V8.Brain_Hippocampus")
for(tissue in tissues){
    for(targetGene in targetGenes){
        tbl.fimo <- get(load(sprintf("~/github/TrenaProjectAD/prep/bigFimo/from-khaleesi/tbl.fimo.%s.RData", targetGene)))
        data.dir <- "~/github/TrenaProjectAD/inst/extdata/genomicRegions"
        filename <- "boca-hg38-consensus-ATAC.RData"
        tbl.boca <- get(load(file.path(data.dir, filename)))
        tbl.fimo.gene <- subset(tbl.fimo, start > region.start & end < region.end)
        tbl.trena.scored <- scoreTargetGeneInTissue(targetGene, "rs4575098", tbl.fimo.gene, tbl.boca, tissue)
        save(tbl.trena.scored, file=sprintf("tbl.trena.scored-%s-%s.RData", targetGene, tissue))
        printf("------- %30s, %10s: %5.2f", tissue, targetGene,
               sum(head(tbl.trena.scored$breakage.score, n=20)),na.rm=TRUE)
        print(head(tbl.trena.scored, n=20))
        } # for targetGene
    } # for tissue
