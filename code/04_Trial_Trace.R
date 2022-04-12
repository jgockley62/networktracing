################################################################################
# Quick Trial Trace for Jesses kinases - 
## Need to add this to the image!
## apt-get install libzmq3-dev
setwd('~/networktracing/code/')
syns_used <- NULL
remotes::install_github("jgockley62/igraphNetworkExpansion")

reticulate::use_python("/usr/bin/python3", required = TRUE)
synapseclient <- reticulate::import("synapseclient")
syn_temp <- synapseclient$Synapse()
syn_temp$login(  )

#Use Co-expression to filter
Use_Cor <- 'NO'
#use directed network
Directed <- 'YES'
#Filter out low edge occurences
Filter_Edges <- 'YES'
#Simple edges ( One Edge per-interaction )
Simple <- 'YES'

syns_used <- c( syns_used, 'syn23283482', 'syn23283475',
                'syn22863899', 'syn22863896', 'syn23283482', 'syn23283475',
                'syn22992753', 'syn22992709')
if( Use_Cor == 'YES' ){
  if( Directed == 'YES' ){
    load( syn_temp$get('syn23283482')$path )
    load( syn_temp$get('syn23283475')$path )
    #Dummy a network (Could choose JS or Non-JS here)
    net <- net_directed
    
  }else{
    #Pull Networks:
    load( syn_temp$get('syn22863899')$path )
    load( syn_temp$get('syn22863896')$path )
    #Dummy a network (Could choose JS or Non-JS here)
    net <- net_undirected
    
  }
}else{
  if( Directed == 'YES' ){
    load( syn_temp$get('syn23283482')$path )
    load( syn_temp$get('syn23283475')$path )
    #Dummy a network (Could choose JS or Non-JS here)
    net <- net_directed
    
  }else{
    #Pull Networks:
    load( syn_temp$get('syn22992753')$path )
    load( syn_temp$get('syn22992709')$path )
    #Dummy a network (Could choose JS or Non-JS here)
    net <- net_undirected
    
  }
}

if( Simple == 'YES' ){
  
}

if( Filter_Edges == 'YES' ){
  #Net is: 13,867 Vertacies and 806,950 Edges
  # A)
  # loose( 246 vertices ) -- A == 13621 Vertacies and 544871 Edges
  test_net <- igraph::subgraph.edges( net,
                                      E(net)[ ( E(net)$EdgeRep == 1  & 
                                                  E(net)$Occurance == 1 & 
                                                  E(net)$Avg_Cortex_CE == 0 ) == F ],
                                      delete.vertices = TRUE
  )
  # loose( 3053 vertices ) -- C 10,814 Vertacies and 179,148 Edges
  test_net <-  igraph::subgraph.edges( test_net,
                                       E(test_net)[ ( E(test_net)$EdgeRep == 1  &
                                                        E(test_net)$Occurance == 2 &
                                                        E(test_net)$Avg_Cortex_CE == 0  ) == F ],
                                       delete.vertices = TRUE
  )
  #loose ( 424 vertacies )
  #B 13,443 vertices and 441,227 Edges
  #test_net <- subgraph.edges( net,
  #                            E(net)[ ( E(net)$EdgeRep == 1  &
  #                                      E(net)$Occurance == 2 & 
  #                                      E(net)$Avg_Cortex_CE == 0  ) == F ],
  #                            delete.vertices = TRUE
  #                          )
  
  # E)
  # loose( 3691 vertices ) -- A == 10176 Vertacies and 151912 Edges
  test_net <-  igraph::subgraph.edges( 
    net,
    igraph::E(net)[ ( igraph::E(net)$Occurance == 1  &
                      igraph::E(net)$Avg_Cortex_CE == 0 ) == F ],
    delete.vertices = TRUE
  )
  net <- test_net
}

syns_used <- c( syns_used, 'syn22758171')
#Annotate vertices on Omics Weights:
OMICS_dep <- read.csv(syn_temp$get('syn22758171')$path)
OMICS <- read.csv( syn_temp$tableQuery( paste0( 'SELECT * FROM syn25575156 WHERE GeneName in (\'',
                                                paste( names(igraph::V(net)), collapse = '\',\'' ),
                                                '\')'),
                                        resultsAs = 'csv' )$filepath)
OMICS <-  OMICS[ , c('ENSG', 'GeneName', 'OmicsScore', 'GeneticsScore', 'Logsdon')]
colnames(OMICS)[ colnames(OMICS) == 'GeneName' ] <- 'GName'
colnames(OMICS)[ colnames(OMICS) == 'Logsdon' ] <- 'Overall'

OMICS$GName <- as.character( OMICS$GName )
OMICS$ENSG <- as.character(OMICS$ENSG)
OMICS_dep$GName <- as.character( OMICS_dep$GName )
OMICS_dep$ENSG <- as.character(OMICS_dep$ENSG)

OMICS_alt <- OMICS[ (OMICS$GName %in% "") == F, ]
OMICS_dep <- OMICS_dep[ (OMICS_dep$GName %in% "") == F, ]

#Pull out pseduo genes and NAs,  also ENSG00000281123 is a diff ENSG for RNA and Protein...:
OMICS_alt <- OMICS_alt[ (OMICS_alt$ENSG %in% c(
  'ENSG00000272655',
  'ENSG00000284770',
  'ENSG00000168255',
  'ENSG00000281123'
)) == F,]
OMICS_dep <- OMICS_dep[ (OMICS_dep$ENSG %in% c(
  'ENSG00000272655',
  'ENSG00000284770',
  'ENSG00000168255',
  'ENSG00000281123'
)) == F,]

OMICS_alt <- OMICS_alt[ is.na(OMICS_alt$GName)==F,]
OMICS_alt <- OMICS_alt[ OMICS_alt$ENSG %in% as.character(OMICS_dep$ENSG), ] 

OMICS_alt  <- OMICS_alt[ !duplicated(OMICS_alt$ENSG), ]
OMICS_alt <- OMICS_alt[!duplicated(OMICS_alt$GName),]
OMICS_alt <- OMICS_alt[ is.na(OMICS_alt$GName) == F, ]
row.names( OMICS_alt ) <- OMICS_alt$GName

OMICS_dep <- OMICS_dep[ is.na(OMICS_dep$GName) == F, ]
row.names( OMICS_dep ) <- OMICS_dep$GName

OMICS_alt$RNA_TE <- OMICS_dep[ row.names(OMICS_alt), ]$RNA_TE
OMICS_alt$Pro_TE <- OMICS_dep[ row.names(OMICS_alt), ]$Pro_TE

#vertex_attr(net, "weight", index = V(net)) <- OMICS_alt[ names( V(net)), ]$Final_Weight 
OMICS_alt[ names( igraph::V(net)), ]$Overall[ is.na(OMICS_alt[ names( igraph::V(net)), ]$Overall) ] <- 0
igraph::vertex_attr(net, "weight", index = igraph::V(net)) <- OMICS_alt[ names( igraph::V(net)), ]$Overall
igraph::vertex_attr(net, "RNA_EffectScore", index = igraph::V(net)) <- OMICS_alt[ names( igraph::V(net)), ]$RNA_TE 
igraph::vertex_attr(net, "Pro_EffectScore", index = igraph::V(net)) <- OMICS_alt[ names( igraph::V(net)), ]$Pro_TE 


# Zero out the TEs
OMICS_alt$PRO_TE_Cor <- OMICS_alt$Pro_TE
OMICS_alt$RNA_TE_Cor <- OMICS_alt$RNA_TE

igraph::vertex_attr(net, "RNA_Cor_EffectScore", index = igraph::V(net)) <- OMICS_alt[ names( igraph::V(net)), ]$RNA_TE_Cor
igraph::vertex_attr(net, "Pro_Cor_EffectScore", index = igraph::V(net)) <- OMICS_alt[ names( igraph::V(net)), ]$PRO_TE_Cor 

################################################################################
################################################################################
#Run the Trace
this_repo <- githubr::getRepo(
  repository = 'jgockley62/networktracing', 
  ref="branch", 
  refName='main'
)
prov <- githubr::getPermlink(
  repository = this_repo,
  repositoryPath = 'code/04_Trial_Trace.R'
  )

syns_used <- c( syns_used, 'syn26126932', 'syn26126931' )


synapseclient <- reticulate::import("synapseclient")
syn_temp <- synapseclient$Synapse()

#synapse <- log_into_synapse(usr='',pass = '')
synap_import <-synapseclient$Synapse()
synap_import$login()
syn_temp$login()

##############################################################################
# Mine the biodomain object
BioDomains <- readRDS(syn_temp$get('syn27368804')$path)
leading_edges <- list()
for(name in names(table(BioDomains$Biodomain))) {
  leading_edges[[name]] <- 'a'
}

for( ent in 1:dim(BioDomains)[1] ){
  if(BioDomains[ent]$padj < 0.05) {
    leading_edges[[BioDomains[ent]$Biodomain]] <- c(
      leading_edges[[BioDomains[ent]$Biodomain]],
      unlist(BioDomains[ent]$leadingEdge)
    )
  }
}

for( name in names(leading_edges) ){
  leading_edges[[name]] <- leading_edges[[name]][ !(leading_edges[[name]] == 'a') ] 
}
### Translate the ENSGs
trans <- read.csv(syn_temp$get('syn26876894')$path, header = T)
trans <- trans[ !(trans$Ensembl_ID %in% c('ENSG00000230417', 'ENSG00000254876', 'ENSG00000276085', '')),]
row.names(trans) <- trans$Ensembl_ID

for( name in names(leading_edges) ){
  leading_edges[[name]] <- trans[ leading_edges[[name]], ]$Approved_symbol
}
for( name in names(leading_edges) ){
  leading_edges[[name]] <- leading_edges[[name]][!duplicated(leading_edges[[name]])]
}
##############################################################################
#Trace Biodomain leading edges:
immune_targets <- leading_edges$`Immune Response`[ leading_edges$`Immune Response` %in% names(V(net))]
trace <- lapply(
  immune_targets,
  igraphNetworkExpansion::short_paths,
  tnet = net,
  targets = immune_targets,
  sentinals = immune_targets,
  cores = 12
)

targ <- NULL
sent <- NULL
for( i in 1:length(trace) ){
  targ <- c(targ,trace[[i]]$Inter)
  sent <- c(sent,trace[[i]]$Sentinal)
}


u_targ <- targ[!duplicated(targ)]
u_sent <- sent[!duplicated(sent)]

# Only genes fount in both target and sentinal traces
opt_1 <- u_targ[u_targ%in%u_sent]

# All genes present more than twice
opt_2 <- c(
  names(table(sent)[table(sent)>1]),
  names(table(targ)[table(targ)>1])
)
opt_2<-opt_2[!duplicated(opt_2)]

# All Genes 
opt_3 <- c(
  names(table(sent)),
  names(table(targ))
)
opt_3<-opt_3[!duplicated(opt_3)]

#trace_filt <- igraphNetworkExpansion::trace_filter(trace)
trace_filt <- opt_3

igraphNetworkExpansion::store_net(
  network = net,
  net_filename = 'Immune_Response',
  net_synname = 'Immune Response Leading Edge Self Trace',
  p_id = 'syn25190666',
  folder = 'Jesses Immune Resoponse Run',
  act_name = 'Jesses Immune Resoponse Subnetwork',
  act_desc = 'Tracces Leading edge genes of significcant Go-Terms of the Immune Response Biodomain to itself',
  synap_import = syn_temp,
  client_import = synapseclient,
  code = NULL,
  repo = NULL,
  syn_used = syns_used,
  subset = trace_filt,
  prov_object = prov
)


##############################################################################
kinases <- igraphNetworkExpansion::list_load(
  'syn26126932',
  network = net,
  is_syn = T, 
  synap_import = syn_temp
)
sentinals <- igraphNetworkExpansion::list_load(
  'syn26126931',
  network = net,
  is_syn = T, 
  synap_import = syn_temp
)

trace <- lapply(
  kinases,
  igraphNetworkExpansion::short_paths,
  tnet = net,
  targets = kinases,
  sentinals = sentinals,
  cores = 12
)

targ <- NULL
sent <- NULL
for( i in 1:length(trace) ){
  targ <- c(targ,trace[[i]]$Inter)
  sent <- c(sent,trace[[i]]$Sentinal)
}

u_targ <- targ[!duplicated(targ)]
u_sent <- sent[!duplicated(sent)]

# Only genes fount in both target and sentinal traces
opt_1 <- u_targ[u_targ%in%u_sent]

# All genes present more than twice
opt_2 <- c(
  names(table(sent)[table(sent)>1]),
  names(table(targ)[table(targ)>1])
)
opt_2<-opt_2[!duplicated(opt_2)]

# All Genes 
opt_3 <- c(
  names(table(sent)),
  names(table(targ))
)
opt_3<-opt_3[!duplicated(opt_3)]

#trace_filt <- igraphNetworkExpansion::trace_filter(trace)
trace_filt <- opt_3

igraphNetworkExpansion::store_net(
  network = net,
  net_filename = 'KinaseRun',
  net_synname = 'Kinase Run',
  p_id = 'syn25190666',
  folder = 'Jesses Kinase Runs',
  act_name = 'Jesses Kinase Subnetwork',
  act_desc = 'Traces Jesses Kinase and Sentinal List',
  synap_import = syn_temp,
  client_import = synapseclient,
  code = NULL,
  repo = NULL,
  syn_used = syns_used,
  subset = trace_filt,
  prov_object = prov
)

net_trace <- igraphNetworkExpansion::network_load(
  syn_id = 'syn26126978', 
  form = 'graphml',
  synap_import = syn_temp
)

igraphNetworkExpansion::store_net(
  network = net,
  net_filename = 'Base Network',
  net_synname = 'Base Network',
  p_id = 'syn25190666',
  folder = 'Jesses Kinase Runs',
  act_name = 'Base Network',
  act_desc = 'Base Network Jesses Kinase and Sentinal List was traced from',
  synap_import = syn_temp,
  client_import = synapseclient,
  code = NULL,
  repo = NULL,
  syn_used = syns_used,
  subset = trace_filt,
  prov_object = prov
)
