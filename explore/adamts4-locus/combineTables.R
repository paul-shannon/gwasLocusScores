filenames <- list.files(pattern="tbl.trena*")
genes <- unlist(lapply(strsplit(filenames, "-"), "[", 2))
tissues <- unlist(lapply(strsplit(filenames, "-"), "[", 3))
tissues.brief <- sub(".RData", "", sub("GTEx_V8.Brain_", "", tissues))
length(genes)
length(tissues.brief)

model.tbls <- list()
breakage.score.tbls <- list()

for(i in seq_len(length(genes))){
    tbl <- head(get(load(filenames[i])), n=20)
    gene <- genes[i]
    tissue <- tissues.brief[i]
    name <- sprintf("%s::%s", tissue, gene)
    score <- sum(tbl$breakage.score, na.rm=TRUE)
    printf("%30s.%s: %5.2f", tissue, gene, score)
    model.tbls$tissue <- tissue
    model.tbls[[name]] <- tbl
    breakage.score.tbls[[name]] <- data.frame(tissue=tissue, gene=gene, score=score, stringsAsFactors=FALSE)
    #browser()
    #xyz <- 99
    }

tbl.bs <- do.call(rbind, breakage.score.tbls)
rownames(tbl.bs) <- NULL

mtx <- matrix(data=0, nrow=4, ncol=8, dimnames=list(unique(tbl.bs$tissue),
                                                    unique(tbl.bs$gene)))
for(r in seq_len(nrow(tbl.bs))){
    tissue <- tbl.bs$tissue[r]
    gene <- tbl.bs$gene[r]
    score <- tbl.bs$score[r]
    mtx[tissue, gene] <- as.integer(score)
    }

mtx
