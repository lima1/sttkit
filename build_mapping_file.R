library(org.Mm.eg.db)
library(org.Rn.eg.db)
library(biomaRt)

x <- read.delim("~/symbol_mappings_20200626.tsv", as.is=T)

human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
rat   <- useMart("ensembl", dataset = "rnorvegicus_gene_ensembl")
mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

genesV2 <- getLDS(attributes = c("hgnc_symbol"), 
   filters = "hgnc_symbol", 
   values = x$Approved.symbol,
   mart = human,
   attributesL = c("mgi_symbol"),
   martL = mouse, uniqueRows=T)
genesV2r <- getLDS(attributes = c("hgnc_symbol"), 
   filters = "hgnc_symbol", 
   values = x$Approved.symbol,
   mart = human,
   attributesL = c("rgd_symbol"),
   martL = rat, uniqueRows=T)

symbol.mapping <- merge(genesV2, genesV2r, by = 1, all.x = TRUE)

save(symbol.mapping, file = "~/git/sttkit/data/symbol.mapping.rda",
    compress = "xz")

