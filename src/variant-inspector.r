# Analysis of circadian RNA-seq data from Marchantia Polymorpha

# DATA PRE-PROCESSING
# Set WD where all the files are

options(stringsAsFactors = FALSE)

# Hmisc - package for Co-Expression Analysis
library(Hmisc, include.only = c('rcorr', '%nin%'))
library(stringr)
library(plyr, include.only = 'colwise') # colwise
library(dplyr)
library(ggplot2)
library(data.table, include.only = c('data.table', 'melt', 'transpose'))
library(tibble)
library(igraph, include.only = c('graph_from_adjacency_matrix', 'components'))


# Import genes of interest
TF.list   <- read.csv('genes/TF-table.csv')$Gene_ID_New
Meri.list <- read.csv("genes/meristeme-genes.csv", header=T)$id



#=================================================
# Cleaned-up PlantRegMap output (promoter mining)
#==================================================

fimo.counts <- dplyr::rename(read.csv('fimo-clean-counts.csv'), TF_id='id_51')

# My hit-count table, `fimo.counts`, uses two different gene name conventions
# for the TF (X.pattern.name) and the target (sequence_name):
#    TF:      "Mapoly0093s0032"
#    target:  "SCR2"
# neither of which suits us, because we need Tak5.1 notation like "Mp5g1110"

# omit ARF3 which i dropped and GRAS1 which is same as SCR2
# (GRAS1 is SCR2, mistake on my part in pre-prearing data, getting rid of doubles)

fimo.counts <- fimo.counts[
    which(fimo.counts$sequence_name %nin% c('GRAS1', 'ARF3')),
]

# convert the short names used for the target to IDs (e.g. "SCR2" -> "Mp1g20490")

sequence_name_map <- read.csv('genes/meristeme-genes.csv')
rownames(sequence_name_map) <- sequence_name_map$name
fimo.counts$target_id <- sequence_name_map[fimo.counts$sequence_name, 'id']

# some IDs are missing from our name-mapping tables and end up as NA,
# they were deleted as the genome version was updated so i drop them

fimo.counts <- na.omit(fimo.counts)




# ===========================================
# convert IDs to human-readable names (e.g. "Mp1g20490" -> "SCR2")
# ===========================================


names <- read.csv('genes/gene-names.tsv', sep="\t", na.strings = c(''));
names <- names[order(names$id),]
rownames(names) <- names$id

missing.ids <- rownames(names)[is.na(names$name)]
names[missing.ids, 'name'] <- read.csv('genes/gene-names-missing.tsv', row.names='id', sep='\t')[missing.ids, 'name']




#============================================
# Let's inspect some variants!
#============================================

split.variant  <- function (ids) str_match(ids, '(.+?)(?:\\.(\\d+))?$')[,-1] # "Mp123456.1" -> c("Mp123456", "1")
remove.variant <- function (ids) str_remove(ids, '\\.\\d+$') # "Mp123456.1" -> "Mp123456"


transpose.with.names <- function (d) {
    d2 <- transpose(d)
    row.names(d2) <- colnames(d)
    colnames(d2) <- row.names(d)
    d2
}

plot.variants <- function(counts, jitter = NULL, title = NULL) {
    counts2 <- transpose.with.names(counts)
    # if (!is.null(jitter)) {
    #     counts2 <- counts2 %>% (colwise(function (v) jitter(v, jitter)))
    # }
    counts2$time <- seq_along(counts2[,1])
    counts3 <- melt(data.table(counts2), id.vars = c("time"), variable.name = "gene", value.name = "count")
    ggp <- (ggplot(counts3, aes(time, count, col = gene))
        + ggtitle(title)
        + geom_line(position = position_jitter(width = jitter, height = jitter))
    )
    ggp
}

# jitter <- function (v) {
#     max <- max(v)
#     min <- min(v)
#     jitter <-
# }

ts.is.close <- function(m, eps) {
    for (j in 1:ncol(m)) {
        for (k in j:ncol(m)) {
            if (!ts.is.close2(m[,j], m[,k], eps)) {
                return(F)
            }
        }
    }
    return(T)
}

ts.is.close2 <- function(a, b, eps) {
    all(abs(a - b) <= eps)
}

get.variants <- function(df, genes) df[remove.variant(row.names(df)) %in% genes,]

find.variant.groups <- function (variant.counts, eps) {
    similarity.adj.mat <- outer(
        1:nrow(variant.counts),
        1:nrow(variant.counts),
        Vectorize(function(j,k) {
            res <- ts.is.close2(variant.counts[j,], variant.counts[k,], eps = eps)
            # print(paste(length(j), length(k), j, k, res))
            res
        })
    )
    membership <- components(
        graph_from_adjacency_matrix(similarity.adj.mat, mode='undirected')
    )$membership
    
    split(
        row.names(variant.counts),
        f = membership
    )
}


# Import Total Data for Gene Counts
normalized.counts <- read.csv("constant-light-counts.csv", row.names='Genes')

# Get rid of zeros
normalized.counts <- normalized.counts[apply(normalized.counts, 1, function(row) any(row != 0)), ]

# keep only the genes whose non-variant name (w/o ".1", ".2" etc) is in TF.list
normalized.counts <- get.variants(normalized.counts, TF.list)

TRESHOLD <- 0.001

all.variants.close <- do.call(rbind, (
    normalized.counts
    %>% rownames_to_column('id')
    %>% group_by(gene = remove.variant(row.names(normalized.counts)))
    %>% filter(n() > 1)
    %>% group_map(function (counts, group) {
        counts <- transpose.with.names(column_to_rownames(data.frame(counts), 'id'))
        data.frame(id=as.character(group), is.close=ts.is.close(counts, TRESHOLD))
    })
))

differing <- sort(filter(all.variants.close, !is.close)$id)
differing

#=======================
# This is it Luigi...
#=======================

{
    PLOTTED_N <- 10 # change to view a different one
    
    PLOTTED <- differing[PLOTTED_N]
    PLOTTED <- "Mp1g01960"
    
    variant.counts <- get.variants(normalized.counts, PLOTTED)
    title <- paste(PLOTTED, " aka ", names[PLOTTED, 'name'], sep="")
    cat("showing:", title, "\n")
    grp <- find.variant.groups(variant.counts, eps = TRESHOLD)
    # just ignore this monstrous bit here, it's just text formatting
    grp.s <- as.vector(lapply(grp, function (x) do.call(paste, c(as.list(x), list(sep=", ")) )))
    cat("proposed groups:\n",
        do.call(paste, c(
            as.list(paste(seq_along(grp.s), grp.s, sep=": ")),
            list(sep = "\n")
        )),
        sep = ""
    )
    plot.variants(
        variant.counts,
        jitter = 0.5, # turn on if needed to distinguish overlapping
        title = title
    )
}

x <- lapply(differing, function(g) find.variant.groups(get.variants(normalized.counts, g), eps = TRESHOLD))
names(x) <- differing
x
