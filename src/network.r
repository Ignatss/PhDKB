# Analysis of circadian RNA-seq data from Marchantia Polymorpha

# DATA PRE-PROCESSING
# Set WD where all the files are

options(stringsAsFactors = FALSE)

# Hmisc - package for Co-Expression Analysis
library('Hmisc', include.only = c('rcorr', '%nin%'))
library('dplyr')
library('stringr', include.only = c('str_match', 'str_remove'))


# Import genes of interest
TF.list   <- read.csv('data/genes/TF-table.csv')$Gene_ID_New
Meri.list <- read.csv("data/genes/meristeme-genes.csv", header=T)$id

# Cleaned-up PlantRegMap output (promoter mining)
fimo.counts <- rename(read.csv('data/fimo-tf.csv'), TF_id='id_61', target_id='sequence_name')

# https://marchantia.info/nomenclature/nomenlatures.txt
nomenclature.names <- local({
    names <- read.csv('data/genes/mp-nomenclature.tsv', sep="\t")
    # v3.1 and v6.1 ids are glued into one column, separated by "; ".
    # we need to split them up
    ids <- str_match(names$GeneID.Location, '(.*); (.*)')[,3] # keep only the v6.1 ids
    names <- data.frame(id=ids, name = names$gene_symbol)
    (names %>%
            mutate(name = str_remove(name, '^Mp')) %>%
            # some genes have a v3.1 id but were removed in v6.1
            filter(!is.na(id)) %>% 
            # if a gene has multiple names, glue them together with "/"
            group_by(id) %>%
            summarise(name = paste0(name, collapse='/'), .groups = 'drop_last')
    )
})

# Human-friendly gene names
gene.names.custom  <- read.csv('data/genes/gene-names.tsv', sep="\t", na.strings = '') %>% na.omit()
gene.names.missing <- read.csv('data/genes/gene-names-missing.tsv', sep="\t")





# get all rows from df where `by` = `id`
df.lookup <- function (df, ids, by) {
    ids2 <- setNames(data.frame(ids), by)
    left_join(ids2, df, by=by)
}

df.lookup.chained <- function(dfs, ids, by, val) {
    vals <- lapply(dfs, function (df) { df.lookup(df, ids, by)[[val]] })
    Reduce(na.overwrite, vals)
}

na.overwrite <- function(a, b) {
    ifelse(is.na(a), b, a)
}

human.names <- function (ids, name_dfs) {
    names.orig <- df.lookup.chained(name_dfs, ids, val='name', by='id')
    split <- split.variant(ids)
    if (!is.array(split)) {
        split <- array(split, dim = c(1,2))
    }
    ids.novariant <- split[,1]
    variants <- as.numeric(split[,2])
    names.novariant <- df.lookup.chained(name_dfs, ids.novariant, val='name', by='id')

    case_when(
        !is.na(names.orig) ~ names.orig,
        # if we don't have a name for MpXXXXX.Y but have one for MpXXXXXX,
        # use that one and append " [.Y]" at the end
        !is.na(names.novariant) & !is.na(variants) & variants > 1 ~
            paste0(names.novariant, " [.", variants, "]"),
        # if we don't have a name at all, just use the ID
        TRUE ~ ids,
    )
}




split.variant  <- function (ids) str_match(ids, '(.+?)(?:\\.(\\d+))?$')[,-1] # "Mp123456.1" -> c("Mp123456", "1")
remove.variant <- function (ids) str_remove(ids, '\\.\\d+$') # "Mp123456.1" -> "Mp123456"

# Flatten a correlation matrix
# (convert from matrix to dataframe)
flatten.corr.matrix <- function (cormat, pmat) {
    ut <- upper.tri(cormat)
    data.frame(
        row = rownames(cormat)[row(cormat)[ut]],
        column = rownames(cormat)[col(cormat)[ut]],
        cor = cormat[ut],
        p = pmat[ut]
    )
}


#============================================
# Create coexpression network based on normalized edge counts
#============================================

coexpression.network <- function (normalized.counts) {
    # Get rid of zeros
    normalized.counts <- normalized.counts[apply(normalized.counts, 1, function(row) any(row != 0)), ]

    # keep only the genes whose non-variant name (w/o ".1", ".2" etc) is in TF.list
    normalized.counts <- normalized.counts[remove.variant(row.names(normalized.counts)) %in% TF.list,]

    # keep only the lowest splicing variant for each gene (e.g. ".1" out of [".1", ".2", ".3"])
   
    normalized.counts <- local({
        min.variants <- (
            split.variant(row.names(normalized.counts))
            %>% as.data.frame()
            %>% rename(gene='V1', variant='V2')
            %>% mutate_at('variant', as.numeric)
            %>% group_by(gene)
            %>% summarise(variant = min(variant), .groups='drop_last')
        )
        
        normalized.counts[
            paste(min.variants$gene, min.variants$variant, sep="."),
        ]
    })

    #Get Co-expression Data for Genes of Interest.
    #Combine Parametric/Non-Parametric Methods for best results
    normalized.counts.rcorr <- rcorr(t(as.matrix(normalized.counts)))

    normalized.counts.rcorr.flat <- flatten.corr.matrix(
        normalized.counts.rcorr$r,
        normalized.counts.rcorr$P
    )

    rm("normalized.counts.rcorr")

    filtered.network <- filter(
        normalized.counts.rcorr.flat,
        # row %in% TF.list, column %in% TF.list,
        p <= 0.001,
        cor <= -0.8 | cor >= 0.8
    )

    rm("normalized.counts.rcorr.flat")

    filtered.network <- rename(filtered.network, target_id='row', TF_id='column')




    # ===========================================
    # Annotate the graph - does the graph edge appear in fimo.counts?
    # ===========================================

    fimo.edges <- paste(fimo.counts$TF_id, fimo.counts$target_id)
    network.edges <- paste(remove.variant(filtered.network$TF_id), remove.variant(filtered.network$target_id))
    filtered.network$in_fimo <- sapply(
        network.edges,
        function(edge) edge %in% fimo.edges
    )
    
    # add target_name and TF_name
    for (col in c('target', 'TF')) {
        ID_COL   <- paste0(col, '_id')
        NAME_COL <- paste0(col, '_name')
        filtered.network[[NAME_COL]] <- human.names(
            filtered.network[[ID_COL]],
            name_dfs = list(
                gene.names.custom,
                nomenclature.names,
                gene.names.missing
            )
        )
    }
    
    filtered.network$cor_sign <- ifelse(filtered.network$cor > 0, 'pos', 'neg')
    
    filtered.network
}


write.network <- function (...) {
    write.table(..., na="", sep="\t", quote = FALSE, row.names = FALSE)
}




targets <- list(
    sporeling = list(
        input_file  = "data/sporeling-counts.csv",
        output_file = "output/network-sporeling.tsv"
    ),
    constant.light = list(
        input_file  = "data/constant-light-counts.csv",
        output_file = "output/network-constant-light.tsv"
    )
)



networks <- lapply(targets, function (target) {
    data <- read.csv(target$input_file, row.names = 'Genes')
    coexpression.network(data)
})

for (target_name in names(networks)) {
    write.network(
        networks[[target_name]],
        file = targets[[target_name]]$output_file,
    )
}

x <- networks$sporeling
y <- networks$constant.light

combined <- (
    inner_join(x, y, by = c('TF_id', 'target_id'), suffix = c('', '.2')) %>%
    mutate(cor_agreement = ifelse(
        cor_sign == cor_sign.2,
        cor_sign,
        'opposite'
    )) %>%
    select(c('TF_id', 'target_id', 'TF_name', 'target_name', 'in_fimo', 'cor_agreement'))
)

write.network(combined, file='output/combined.tsv')





