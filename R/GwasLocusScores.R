#' @title GwasLocusScores
#' @description assess RNA-binding protein binding sites, their structure and motif matches
#' @name GwasLocusScores

#' @import R6
#' @import rtracklayer
#' @import TrenaMultiScore
#' @import TrenaProjectAD
#' @import EndophenotypeExplorer
#' @import ADvariantExplorer
#' @import trena
#' @import MotifDb
#' @import motifbreakR
#' @import SNPlocs.Hsapiens.dbSNP151.GRCh38
#' @import BSgenome.Hsapiens.UCSC.hg38
#' @import BiocParallel
#' @import ghdb
#'
#' @field rbp the gene symbol name of the protein
#'
#' @examples
#'   gls <- GwasLocusScores("rs4575096")
#' @export
#----------------------------------------------------------------------------------------------------
GwasLocusScores = R6Class("GwasLocusScores",
 #--------------------------------------------------------------------------------
   private = list(
      tag.snp=NULL,
      targetGene=NULL,
      chrom=NULL,
      start.loc=NULL,
      end.loc=NULL,
      tissue.name=NULL,
      advx=NULL,      # ADvariantExplorer
      etx=NULL,        # EndophenotypeExplorer
      tbl.eqtl=NULL,
      eqtl.catalog=NULL,
      motifbreakR.results=NULL,
      tms=NULL,
      tbl.tms=NULL,
      tbl.trena=NULL
      ),
 #--------------------------------------------------------------------------------
   public = list(

         #--------------------------------------------------------------------------------
         #' @description
         #' Creates a new instance of this class.
         #'
         #' @param tag.snp character, an rsid identifying the locus
         #' @param chrom character, e.g., "chr1"
         #' @param start.loc numeric the start of the genomic region of interest
         #' @param end.loc numeric the end of the genomic region of interest
         #' @param tissue.name character, e.g. "GTEx_V8.Brain_Hippocampus"
         #' @param targetGene character, e.g. "NDUFS2"
         #' @return an object of the GwasLocusScores class

      initialize = function(tag.snp, chrom, start.loc, end.loc, tissue.name, targetGene){
         stopifnot(length(tissue.name) == 1)
         private$targetGene <- targetGene
         private$tag.snp <- tag.snp
         private$chrom <- chrom
         private$start.loc <- start.loc
         private$end.loc <- end.loc
         private$advx <- ADvariantExplorer$new(NULL, chrom, start.loc, end.loc)
         private$etx <- EndophenotypeExplorer$new(targetGene, "hg38")
         private$tissue.name <- tissue.name
         private$eqtl.catalog <- private$advx$geteQTLSummary()
         search.term <- sprintf("^%s$", tissue.name)
         if(length(grep(search.term, private$eqtl.catalog$unique_id)) !=1)
            stop(sprintf("'%s' does not uniquely identify a catalog study", tissue.name))
         },

         #--------------------------------------------------------------------------------
         #' @description
         #' retrieves all eQTLs from the ebi eqtl catalogue
         #'
         #' @returns nothing
      load.eqtls = function(){
          data.dir <- "~/github/gwasLocusScores/inst/extdata/eqtl.downloads"
          files <- list.files(path=data.dir, pattern="GTEx.*RData")
          filename <- grep(private$tissue.name, files, value=TRUE)
          full.path <- file.path(data.dir, filename)
          #tbl.eqtl <- private$advx$geteQTLsByLocationAndStudyID(private$chrom,
          #                                                     private$start.loc,
          #                                                     private$end.loc,
          #                                                     private$tissue.name,
          #                                                     method="tabix", simplify=TRUE)
         private$tbl.eqtl <- get(load(full.path))
         },

         #--------------------------------------------------------------------------------
         #' @description
         #' extract all previously obtained eQTLs at or aboce the specified pvalue threshold
         #'
         #' @param pvalue.cutoff numeric, e.g., 1e-4
         #' @returns a data.frame
      get.tbl.eqtl = function(pvalue.cutoff){
         if(is.null(private$tbl.eqtl))
             return(NULL)
         return(subset(private$tbl.eqtl, pvalue <= pvalue.cutoff))
         },

         #--------------------------------------------------------------------------------
         #' @description
         #' discover which motifs are broken by the eQTL variants, for specified gene and pval
         #'
         #' @param targetGene character, a gene symbol in the eQTL table
         #' @param pvalue.cutoff numeric, e.g., 1e-4
         #' @param TFs.preselected character, e.g., c("SOX21", "ZEB1")
         #' @returns a data.frame
      breakMotifsAtEQTLs = function(targetGene, pvalue.cutoff=NA, TFs.preselected=NA){
         printf("--- starting motifbreakR")
         if(is.null(private$tbl.eqtl))
             self$load.eqtls()
         tbl.sub <- subset(private$tbl.eqtl, gene==targetGene)
         if(!is.na(pvalue.cutoff))
             tbl.sub <- subset(tbl.sub, pvalue <= pvalue.cutoff)
         rsids <- tbl.sub$rsid
         motifs.selected <- query(MotifDb, "sapiens", c("jaspar2018", "hocomoco-core-A"))
         if(!all(is.na(TFs.preselected))){
            motifs.selected <- query(motifs.selected, andStrings="sapiens", orStrings=TFs.preselected)
            }
         if(length(rsids) == 0) browser()
         message(sprintf("--- breakMotifsAtEQTLs, creating snps.gr for %d variants", length(rsids)))
         snps.gr <- snps.from.rsid(rsid = rsids,
                              dbSNP=SNPlocs.Hsapiens.dbSNP151.GRCh38,
                              search.genome=BSgenome.Hsapiens.UCSC.hg38)

         available.workers <- multicoreWorkers()
         printf("--- starting motifbreakR on %d eQTLs with %d workers", length(snps.gr), available.workers)
         bpparam <- MulticoreParam(workers=available.workers)
         private$motifbreakR.results <- motifbreakR(snpList = snps.gr,
                                                    filterp = TRUE,
                                                    pwmList = motifs.selected,
                                                    show.neutral=FALSE,
                                                    method = c("ic", "log", "notrans")[1],
                                                    bkg = c(A=0.25, C=0.25, G=0.25, T=0.25),
                                                    BPPARAM = bpparam,
                                                    verbose=TRUE)
         tbl.breaks <- as.data.frame(private$motifbreakR.results, row.names=NULL)
         tbl.breaks <- subset(tbl.breaks, effect=="strong")
         tbl.breaks$pctDelta <- with(tbl.breaks, pctAlt - pctRef)
         return(invisible(tbl.breaks))
         }, # breakMotifsAtEQTLs

         #--------------------------------------------------------------------------------
         #' @description
         #' use FIMO and open chromatin, and TrenaMultiScore capabilities to create a tbl.tms
         #'
         #' @param trenaProject a concrete TrenaProject object, e.g., TrenaProjectAD()
         #' @param tbl.fimo data.frame
         #' @param tbl.oc data.frame
         #' @param mtx.rna matrix
         #' @returns a data.frame

      createTrenaMultiScoreTable = function(trenaProject, tbl.fimo, tbl.oc, mtx.rna=NA){
         source("~/github/TrenaMultiScore/tools/runner/v2/tmsCore.R")
         tms <- TMS$new(trenaProject, private$targetGene, tbl.fimo, tbl.oc)
         tms$scoreFimoTFBS()   # chip, conservation, genehancer, genic annotations, distance to tss
         tms$add.tf.mrna.correlations(mtx.rna, featureName="cor.all")
         tbl.tms <- tms$getTfTable()
         private$tbl.tms <- tbl.tms
         private$tms <- tms
         return(invisible(tbl.tms))
         },

         #--------------------------------------------------------------------------------
         #' @description
         #' given tfs, implicit target gene and an expression matrix, run trena
         #'
         #' @param tfs character, a vector of transcription factor geneSymbol names
         #' @param mtx.rna matrix
         #' @param tbl.tms.final  data.frame, a TrenaMultiScore table from which the tfs came
         #' @returns a data.frame

      runTrena = function(tfs, mtx.rna, tbl.tms.final=NA){
         solvers <- c("lasso", "Ridge", "Spearman", "Pearson", "RandomForest", "xgboost")
         solver <- EnsembleSolver(mtx.rna,
                                  targetGene=private$targetGene,
                                  candidateRegulators=tfs,
                                  solverNames=solvers,
                                  geneCutoff=1.0)
         tbl.trena <- run(solver)
         tbl.trena <- tbl.trena[order(abs(tbl.trena$rfScore), decreasing=TRUE),]
         tbl.trena$rank <- seq_len(nrow(tbl.trena))
         rownames(tbl.trena) <- NULL
         tbl.trena$target <- private$targetGene
         for(col in 2:7)
             tbl.trena[, col] <- round(tbl.trena[, col], digits=3)
         tbl.trena$tfbs <- unlist(lapply(tbl.trena$gene, function(gene) length(grep(gene, tbl.tms.final$tf))))
         private$tbl.trena <- tbl.trena
         return(tbl.trena)
         },

         #--------------------------------------------------------------------------------
         #' @description
         #' access to the ADvariantExplorer object
         # @returns ADvariantExplorer
      getADVX = function(){
         private$advx
         },

         #--------------------------------------------------------------------------------
         #' @description
         #' access to the EndophenotypeExplorer object
         # @returns EndophenotypeExplorer
      getETX = function(){
         private$etx
         },

         #--------------------------------------------------------------------------------
         #' @description
         #' calculates the aggregate impact of broken motifs for TFs in the trena model
         #'
         #' @param tbl.trena data.frame
         #' @param tbl.breaks data.frame from motifbreakR
         #' @param tbl.eqtl data.frame
         #'
         # @returns data.frame
      scoreEQTLs = function(tbl.trena, tbl.breaks, tbl.eqtl){
         tfs.oi <- tbl.trena$gene
         tbl.eqtl.sub <- subset(tbl.eqtl, gene==private$targetGene)
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
            #printf("--- %s: %5.2f", tf, tf.breakage.score)
            breakage.scores[[tf]] <- tf.breakage.score
            } # for tf

         tbl.trena$breakage.score <- as.numeric(breakage.scores)
         tbl.trena$rfNorm <- round(tbl.trena$rfNorm, digits=2)
         tbl.trena$breakage.score <- round(tbl.trena$breakage.score, digits=2)

         return(tbl.trena)
         } # scoreEQTLs

   ) # public

 ) # class GwasLocusScores
#--------------------------------------------------------------------------------
