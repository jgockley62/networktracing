library(dplyr)
library(rtracklayer)

synapseclient <- reticulate::import("synapseclient")
syn_temp <- synapseclient$Synapse()
syn_temp$login(  )

### load Translation
trans <- read.csv( syn_temp$get("syn26876894")$path)

### - laod Pipeline GTF and biomart object
bm <- read.table( synGet('syn26454641')$path, sep = '\t', header = T)

RNA_gtf <- syn_temp$get('syn20692159')
RNA_gtf <- GenomicTools.fileHandler::importGTF(file=RNA_gtf$path,
                                               level="gene",
                                               features=c("gene_id",
                                                          "gene_name",
                                                          "gene_type",
                                                          'hgnc_id'))
RNA_gtf$gene_id <- do.call(rbind,strsplit(RNA_gtf$gene_id,'[.]'))[,1]
parY <- RNA_gtf[duplicated(RNA_gtf$gene_id),]$gene_id
RNA_gtf <- RNA_gtf[ !(RNA_gtf$gene_id %in% parY & RNA_gtf$V1 == 'chrY'),  ]

#sync Names:
RNA_gtf <- as.data.frame(RNA_gtf)
row.names(RNA_gtf) <- RNA_gtf$gene_id
row.names(bm) <- bm$ensembl_gene_id
bm$hgnc_id <- RNA_gtf[row.names(bm),]$hgnc_id
bm$hgnc_id <- gsub("character\\(0\\)",NA,bm$hgnc_id)


# Load Single cell GTF and filtered expression Matrix
sc_allen <- read.csv(syn_temp$get('syn26720729')$path, row.names = 1)
sc_names <- row.names(sc_allen)
poo_0 <- sc_allen %>% filter_at(vars(1:120), any_vars(. > 0))
poo_05 <- sc_allen %>% filter_at(vars(1:120), any_vars(. > 0.5))
poo_1 <- sc_allen %>% filter_at(vars(1:120), any_vars(. > 1))
poo_25 <- sc_allen %>% filter_at(vars(1:120), any_vars(. > 1))

## Total Allen
total_allen <- data.table::fread(syn_temp$get('syn25883506')$path)
allen_genes <- colnames(total_allen)[2:length(colnames(total_allen))]

table(allen_genes %in% bm$hgnc_symbol)
#     FALSE  TRUE 
#     19594 30687
#     30687/(19594+30687)
#     61.03%

################################################################################
#Demarcate the allen translation gene in bm
bm$allen_gene_name <- NA
bm[ bm$hgnc_symbol %in% allen_genes, ]$allen_gene_name <- bm[ bm$hgnc_symbol %in% allen_genes, ]$hgnc_symbol

table(bm$gene_biotype)['protein_coding']

table(bm[bm$hgnc_symbol %in% allen_genes,]$gene_biotype)['protein_coding']
# 92.52894% Protien Coding Coverage
################################################################################
#Find remainder GTF Enrties, locate Alt IDs:
remainder <- allen_genes[!(allen_genes %in% bm$hgnc_symbol)]
table(remainder %in% trans$Approved_symbol)


Allen_gtf <- syn_temp$get('syn26434482')
Allen_gtf <- GenomicTools.fileHandler::importGTF(file=Allen_gtf$path,
                                                 level="gene",
                                                 features=c("gene_id",
                                                            "gene_symbol"))
Allen_gtf <- Allen_gtf[!(Allen_gtf$gene_symbol %in% bm$hgnc_symbol),]
Allen_gtf <- as.data.frame(Allen_gtf)
table(Allen_gtf$gene_id %in% trans$NCBI_Gene_ID)
trans$NCBI_Gene_ID <- as.character(trans$NCBI_Gene_ID)

seeker <- function( x, biom ){
  #biom <- bm
  if(x$HGNC_ID %in% biom$hgnc_id){
    biom[ biom$hgnc_id %in% x$HGNC_ID,]$allen_gene_name <- x$allen_gene_name 
  }else{
  }
  if(x$Ensembl_ID %in% biom$ensembl_gene_id){
    biom[ biom$ensembl_gene_id == x$Ensembl_ID,]$allen_gene_name <- x$allen_gene_name 
  }else{
  }
  return(biom)
}


tran_mini <- trans[ trans$NCBI_Gene_ID %in% Allen_gtf$gene_id, ]
row.names(Allen_gtf) <- Allen_gtf$gene_id

#temp_gtf <- Allen_gtf[row.names(Allen_gtf) %in% trans$NCBI_Gene_ID,  ]
tran_mini$allen_gene_name <- temp_gtf[as.character(tran_mini$NCBI_Gene_ID),]$gene_symbol

#table(tran_mini$HGNC_ID %in% bm$hgnc_id)
#table(tran_mini$Ensembl_ID %in% bm$ensembl_gene_id)

#table(tran_mini$Ensembl_ID %in% bm$ensembl_gene_id | tran_mini$HGNC_ID %in% bm$hgnc_id)
temp_bm <- bm
for( i in 1:dim(tran_mini)[1]){
  temp_bm <- seeker(tran_mini[i,],temp_bm)
}

bm <- temp_bm

table(allen_genes %in% bm$allen_gene_name)
# 69.34%
#table(bm$gene_biotype)
#table(bm[!is.na(bm$allen_gene_name),]$gene_biotype)
################################################################################

remainder <- allen_genes[!(allen_genes %in% bm$allen_gene_name)]
Allen_gtf <- Allen_gtf[!(Allen_gtf$gene_symbol %in% bm$allen_gene_name),]

##### - entrezgene_accession match
bm_temp <- bm[is.na(bm$allen_gene_name),]

mart <- biomaRt::useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
gene_list = bm_temp$ensembl_gene_id
test = biomaRt::getBM(
  attributes = c('ensembl_gene_id',"hgnc_symbol", "entrezgene_accession"), 
  filters = "ensembl_gene_id", 
  values = gene_list, bmHeader = T, mart = mart)

test <- test[ !( test$`NCBI gene (formerly Entrezgene) accession` == ""), ]
test <- test[ test$`NCBI gene (formerly Entrezgene) accession` %in% remainder, ]

test_cust <- test[test$`Gene stable ID` %in% names(table(test$`Gene stable ID`)[table(test$`Gene stable ID`)>1]), ]
test <- test[!(test$`Gene stable ID` %in% names(table(test$`Gene stable ID`)[table(test$`Gene stable ID`)>1])), ]

for(i in 1:dim(test)[1]){
  bm_temp[ bm_temp$ensembl_gene_id %in% test[i,]$`Gene stable ID`, ]$allen_gene_name <- 
    test[i,]$`NCBI gene (formerly Entrezgene) accession`
}
#table(is.na(bm$allen_gene_name))
#table(is.na(bm_temp$allen_gene_name))

bm <- bm_temp

##### - NCBI gene (formerly Entrezgene) NUMBER

remainder <- allen_genes[!(allen_genes %in% bm$allen_gene_name)]
Allen_gtf <- Allen_gtf[!(Allen_gtf$gene_symbol %in% bm$allen_gene_name),]

mart <- biomaRt::useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
#24914
gene_list = bm_temp[is.na(bm_temp$allen_gene_name),]$ensembl_gene_id
test = biomaRt::getBM(
  attributes = c('ensembl_gene_id',"hgnc_symbol", "entrezgene_accession", 'entrezgene_id'), 
  filters = "ensembl_gene_id", 
  values = gene_list, bmHeader = T, mart = mart)

test <- test[ !(is.na(test$`NCBI gene (formerly Entrezgene) ID`)),]

table(as.character(test$`NCBI gene (formerly Entrezgene) ID`) %in% Allen_gtf$gene_id)


test <- test[test$`NCBI gene (formerly Entrezgene) ID` %in% Allen_gtf$gene_id,]

for( i in 1:dim(test)[1]){
  add_name <- Allen_gtf[ Allen_gtf$gene_id %in% as.character(test[i,]$`NCBI gene (formerly Entrezgene) ID`),]$gene_symbol
  bm[ bm$ensembl_gene_id %in% test[i,]$`Gene stable ID`, ]$allen_gene_name <- add_name
  
}
bm_temp <- bm

################################################################################
# Pathway Commons:
pc <- read.csv( syn_temp$get('syn26844067')$path, row.names = 1 )
pc <- pc
# - Name Present Pathway Commons genes:
bm$Path_commons_name <- NA
for( pc_name in pc ){
  if( pc_name %in% bm$hgnc_symbol){
    bm[ bm$hgnc_symbol  == pc_name,]$Path_commons_name <- pc_name
  }
}

# Find ENSG of Pathway Commons Genes:
missing_pc <- pc[ !(pc %in% bm$hgnc_symbol )]

gene_list = missing_pc
test = biomaRt::getBM(
  attributes = c('ensembl_gene_id',"hgnc_symbol"), 
  filters = "hgnc_symbol", 
  values = gene_list, bmHeader = T, mart = mart)

row.names(test) <- test$`Gene stable ID`
test <- test[test$`Gene stable ID` %in% bm$ensembl_gene_id,]
bm[ test$`Gene stable ID`, ]$Path_commons_name <- test$`HGNC symbol`

missing_pc <- missing_pc[!(missing_pc %in%bm$Path_commons_name) ]


trans[ ("DVL1P1" == trans$Prev_symbols) |
         (grepl(", DVL1P1", trans$Prev_symbols)) |
         (grepl("DVL1P1, ", trans$Prev_symbols)),
]

id_names <- NULL
multi <- NULL
multi_id <- NULL
single <- NULL
for( i in missing_pc ){
  if(TRUE %in% grepl(i, trans$Prev_symbols)){
    
    temp <- trans[ (grepl(paste0("\\b",i,"\\b"), trans$Prev_symbols)) |
                     (grepl(paste0(", ",i,"\\b"), trans$Prev_symbols)) |
                     (grepl(paste0("\\b",i,"), "), trans$Prev_symbols)),
    ]
    if( dim(temp)[1] > 1 ){
      multi <- as.data.frame(rbind( multi,temp ))
      multi_id <- c(multi_id,i)
    }else{
      single <- as.data.frame(rbind( single,temp ))
      id_names <- c(id_names, i)
    }
  }
}


## FALSE  TRUE 
## 18995 41563 
#sink <- bm_temp
bm_temp <- bm

# - Replace the single names
fail <- NULL
for( i in id_names ){
  if( i %in% single$Prev_symbols ){
    temp <- single[ (grepl(paste0("\\b",i,"\\b"), single$Prev_symbols)) |
                      (grepl(paste0(", ",i,"\\b"), single$Prev_symbols)) |
                      (grepl(paste0("\\b",i,"), "), single$Prev_symbols)),
    ]
    if( temp$Ensembl_ID %in% bm$ensembl_gene_id ){
      bm[ bm$ensembl_gene_id == temp$Ensembl_ID,]$Path_commons_name <- i
    }else{
      fail <- c(fail,i)
    }
  }else{
    fail <- c(fail,i)
  }
}

table(is.na(bm$Path_commons_name))

# FALSE  TRUE 
# 19019 41539

missing_pc <- c(missing_pc[ !(missing_pc %in% id_names) ], fail)


### Manual Curation:
manual_curation <- rbind(
  c('LRRC75A-AS1',	'ENSG00000175061'),
  c('SSSCA1',	      'ENSG00000173465'),
  c('WISP1',	      'ENSG00000104415'),
  c('CYR61',	      'ENSG00000142871'),
  c('KIAA1551',	    'ENSG00000174718'),
  c('FAM208A',	    'ENSG00000163946'),
  c('C7orf43',	    'ENSG00000146826'),
  c('FAM198B',	    'ENSG00000164125'),
  c('TROVE2',	      'ENSG00000116747'),
  c('FBXL21',	      'ENSG00000164616'),
  c('PLA2G16',	    'ENSG00000176485')
)

bm[ manual_curation[,2],]$Path_commons_name <- manual_curation[,1] 

################################################################################

trans[ trans$Ensembl_ID %in% bm[ (!(is.na(bm$Path_commons_name)) & is.na(bm$allen_gene_name)),]$ensembl_gene_id,]
write.csv(bm, 'Omics_sc_PathwayCommons_Translation_File.csv')
################################################################################






