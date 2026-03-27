library(phytools)
library(ape)


#################################################
# inserting species in the phylogeny
####################################################

#create df for especies insertion
df_media_spp <- read.csv("spp_mean_data.csv")
df_spp <- as.data.frame(cbind(superfilo_name=df_media_spp$superfilo_name,genus=df_media_spp$genus, section=df_media_spp$section, spp=df_media_spp$spp))

# read tree
tree <- read.tree("mmc_target_common_Oct17_pruned.tre")

#pruuning tree for df species
tips_keep <- df_spp$superfilo_name[!is.na(df_spp$superfilo_name)]
setdiff(tips_keep, tree$tip.label)
tree <- keep.tip(tree, tips_keep)
Ntip(tree)
length(tips_keep)

# separate species present vs missing
present <- df_spp[!is.na(df_spp$superfilo_name), ]
missing <- df_spp[is.na(df_spp$superfilo_name), ]

# function to insert tip randomly within a clade
insert_random_tip <- function(tree, new_tip, ref_tips){
  
  node <- getMRCA(tree, ref_tips)
  
  clade <- extract.clade(tree, node)
  
  edges <- which(tree$edge[,2] %in% match(clade$tip.label, tree$tip.label))
  
  edge <- sample(edges,1)
  
  pos <- runif(1,0,tree$edge.length[edge])
  
  tree <- bind.tip(tree,
                   tip.label = new_tip,
                   where = tree$edge[edge,2],
                   position = pos)
  
  return(tree)
}

# function that inserts all missing species
insert_species <- function(tree, df_spp){
  
  present <- df_spp[!is.na(df_spp$superfilo_name), ]
  missing <- df_spp[is.na(df_spp$superfilo_name), ]
  
  for(i in 1:nrow(missing)){
    
    sp  <- missing$spp[i]
    g   <- missing$genus[i]
    sec <- missing$section[i]
    
    ref <- NULL
    
    # SECTION insertion
    if(!is.na(sec)){
      
      ref <- present$superfilo_name[present$section == sec]
      ref <- ref[!is.na(ref)]
      
      if(length(ref) < 2){
        ref <- NULL
      }
    }
    
    # GENUS fallback
    if(is.null(ref)){
      
      ref <- present$superfilo_name[present$genus == g]
      ref <- ref[!is.na(ref)]
    }
    
    # insertion logic
    if(length(ref) >= 2){
      
      tree <- insert_random_tip(tree, sp, ref)
      
    } else if(length(ref) == 1){
      
      tree <- bind.tip(tree,
                       tip.label = sp,
                       where = which(tree$tip.label == ref),
                       position = runif(1, 0, 0.001))
    }
  }
  
  return(tree)
}

# generate 100 trees with random placements
set.seed(123)

trees <- replicate(
  100,
  insert_species(tree, df_spp),
  simplify = FALSE
)

# check number of tips
sapply(trees, Ntip)

#check inserted spp
expected_tips <- c(
  df_spp$superfilo_name[!is.na(df_spp$superfilo_name)],
  df_spp$spp[is.na(df_spp$superfilo_name)]
)
expected_tips <- unique(expected_tips)
setdiff(expected_tips, trees[[1]]$tip.label)

# save trees
write.tree(trees, file = "trees_with_inserted_species.tre")
