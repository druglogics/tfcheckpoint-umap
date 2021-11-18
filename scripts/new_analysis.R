library(dplyr)
library(tibble)
library(tidyr)
library(readr)
library(uwot)
library(ggplot2)

# Read data in TFch2 (GO data from November 8th) => 3842 proteins
# Format
# Entry\tName\tGO annotations\tInterpro domains\tDeepTF score\tGO annotations no IEA
data_tbl    = readr::read_tsv(file = 'data/data_tbl_16112021.tsv')
protein_ids = data_tbl %>% pull(Entry)
domains     = data_tbl %>% pull(`Interpro domains`)
gene_names  = data_tbl %>% pull(Name)
go_terms    = data_tbl %>% pull(`GO annotations`)
go_terms_no_iea   = data_tbl %>% pull(`GO annotations no IEA`)
go_terms_only_iea = data_tbl %>% pull(`GO annotations only IEA`)
deeptf_scores     = data_tbl %>% pull(`DeepTF score`)

# If matrices have already been stored as compressed R objects, don't rebuild them!
protein_domain_rds_file          = 'data/protein_domain_mat.rds'
protein_go_rds_file              = 'data/protein_go_mat.rds'
protein_go_no_iea_rds_file       = 'data/protein_go_no_iea_mat.rds'
protein_go_only_iea_rds_file     = 'data/protein_go_only_iea_mat.rds'
if ((!file.exists(protein_domain_rds_file)) | (!file.exists(protein_go_rds_file)) | (!file.exists(protein_go_no_iea_rds_file)) | (!file.exists(protein_go_only_iea_rds_file))) {
  # get the unique InterPro domain IDs
  domains_unique = sapply(domains, function(domain_vec) {
    if (!is.na(domain_vec)) {
      vec = stringr::str_split(domain_vec, ";") %>% unlist() %>% head(-1)
      return(vec)
    }
  }) %>% unlist(use.names = FALSE) %>% unique() # 4637 unique domains

  # get the unique GO terms (with IEA)
  go_terms_unique = sapply(go_terms, function(go_term_vec) {
    if (!is.na(go_term_vec)) {
      vec = stringr::str_split(go_term_vec, " ") %>% unlist()
      return(vec)
    }
  }) %>% unlist(use.names = FALSE) %>% unique() # 9676 unique GO terms

  # get the unique GO terms (no IEA terms)
  go_terms_no_iea_unique = sapply(go_terms_no_iea, function(go_term_vec) {
    if (!is.na(go_term_vec)) {
      vec = stringr::str_split(go_term_vec, " ") %>% unlist()
      return(vec)
    }
  }) %>% unlist(use.names = FALSE) %>% unique() # 8096 unique GO terms (subset of the above)
  stopifnot(all(go_terms_no_iea_unique %in% go_terms_unique))

  # get the unique GO terms (only IEA terms)
  go_terms_only_iea_unique = sapply(go_terms_only_iea, function(go_term_vec) {
    if (!is.na(go_term_vec)) {
      vec = stringr::str_split(go_term_vec, " ") %>% unlist() # 4610 unique GO terms (subset of the 9676 above)
      return(vec)
    }
  }) %>% unlist(use.names = FALSE) %>% unique()
  stopifnot(all(go_terms_only_iea_unique %in% go_terms_unique))

  # build matrices of 0's and 1's
  protein_domain_mat      = matrix(0, nrow = length(protein_ids), ncol = length(domains_unique),
    dimnames = list(protein_ids, domains_unique))
  protein_go_mat          = matrix(0, nrow = length(protein_ids), ncol = length(go_terms_unique),
    dimnames = list(protein_ids, go_terms_unique))
  protein_go_no_iea_mat   = matrix(0, nrow = length(protein_ids), ncol = length(go_terms_no_iea_unique),
    dimnames = list(protein_ids, go_terms_no_iea_unique))
  protein_go_only_iea_mat = matrix(0, nrow = length(protein_ids), ncol = length(go_terms_only_iea_unique),
    dimnames = list(protein_ids, go_terms_only_iea_unique))

  # protein - InterPro matrix
  for(id in protein_ids) {
    # get the InterPro domains for that protein ID
    domain_vec = data_tbl %>%
      filter(Entry == id) %>%
      pull(`Interpro domains`) %>%
      unique()

    # ... and set them to 1 in the matrix
    if (!is.na(domain_vec)) {
      vec = stringr::str_split(domain_vec, ";") %>% unlist() %>% head(-1)
      for(domain in vec) {
        protein_domain_mat[id, domain] = 1
      }
    }
  }

  # Protein - GO matrix
  for(id in protein_ids) {
    # get the GO terms for that protein ID
    go_term_vec = data_tbl %>%
      filter(Entry == id) %>%
      pull(`GO annotations`) %>%
      unique()

    # ... and set them to 1 in the matrix
    if (!is.na(go_term_vec)) {
      vec = stringr::str_split(go_term_vec, " ") %>% unlist()
      for(go_term in vec) {
        protein_go_mat[id, go_term] = 1
      }
    }
  }

  # Protein - GO (no IEA) matrix
  for(id in protein_ids) {
    # get the GO terms for that protein ID
    go_term_vec = data_tbl %>%
      filter(Entry == id) %>%
      pull(`GO annotations no IEA`) %>%
      unique()

    # ... and set them to 1 in the matrix
    if (!is.na(go_term_vec)) {
      vec = stringr::str_split(go_term_vec, " ") %>% unlist()
      for(go_term in vec) {
        protein_go_no_iea_mat[id, go_term] = 1
      }
    }
  }

  # Protein - GO (only IEA) matrix
  for(id in protein_ids) {
    # get the GO terms for that protein ID
    go_term_vec = data_tbl %>%
      filter(Entry == id) %>%
      pull(`GO annotations only IEA`) %>%
      unique()

    # ... and set them to 1 in the matrix
    if (!is.na(go_term_vec)) {
      vec = stringr::str_split(go_term_vec, " ") %>% unlist()
      for(go_term in vec) {
        protein_go_only_iea_mat[id, go_term] = 1
      }
    }
  }

  # save the matrices
  saveRDS(object = protein_domain_mat,      file = protein_domain_rds_file)
  saveRDS(object = protein_go_mat,          file = protein_go_rds_file)
  saveRDS(object = protein_go_no_iea_mat,   file = protein_go_no_iea_rds_file)
  saveRDS(object = protein_go_only_iea_mat, file = protein_go_only_iea_rds_file)
} else {
  protein_domain_mat      = readRDS(file = protein_domain_rds_file)
  protein_go_mat          = readRDS(file = protein_go_rds_file)
  protein_go_no_iea_mat   = readRDS(file = protein_go_no_iea_rds_file)
  protein_go_only_iea_mat = readRDS(file = protein_go_only_iea_rds_file)
}

# save the GO matrices as files and zip them afterwards manually :)
if (FALSE) {
  go_mat          = protein_go_mat
  go_mat_no_iea   = protein_go_no_iea_mat
  go_mat_only_iea = protein_go_only_iea_mat

  # change rownames
  rownames(go_mat)          = gene_names
  rownames(go_mat_no_iea)   = gene_names
  rownames(go_mat_only_iea) = gene_names

  # make tibbles and move rownames to a separate column
  go_tbl          = go_mat %>% as_tibble(rownames = 'Name')
  go_no_iea_tbl   = go_mat_no_iea %>% as_tibble(rownames = 'Name')
  go_only_iea_tbl = go_mat_only_iea %>% as_tibble(rownames = 'Name')

  readr::write_csv(x = go_tbl, file = 'data/protein_go_mat.csv')
  readr::write_csv(x = go_no_iea_tbl, file = 'data/protein_go_no_iea_mat.csv')
  readr::write_csv(x = go_only_iea_tbl, file = 'data/protein_go_only_iea_mat.csv')
}

# Check: percentage of 1's (From more to less)
100*sum(protein_go_mat)/(dim(protein_go_mat)[1]*dim(protein_go_mat)[2])
100*sum(protein_go_no_iea_mat)/(dim(protein_go_no_iea_mat)[1]*dim(protein_go_no_iea_mat)[2])
100*sum(protein_go_only_iea_mat)/(dim(protein_go_only_iea_mat)[1]*dim(protein_go_only_iea_mat)[2])
100*sum(protein_domain_mat)/(dim(protein_domain_mat)[1]*dim(protein_domain_mat)[2])

# Check how many GO terms per protein?
rowSums(protein_go_mat) %>% sort(decreasing = T) %>% head(20)
rowSums(protein_go_no_iea_mat) %>% sort(decreasing = T) %>% head(20) # less GO terms per protein

####################################
# Read GREEKC DbTF list            #
# See commit 6cd92c8 for more info #
####################################
greekc_dbtfs_tbl = readr::read_csv(file = 'data/dbTF_online_02112021.csv',
  col_select = c("id", "symbol"))
dbtf_ids = greekc_dbtfs_tbl %>% pull(id) # 1437 DbTF curated proteins

# check: all GREEKC annotated DbTF proteins are in TFch2 (YES)
stopifnot(dbtf_ids %in% protein_ids)

##########################################################
# TFch2 Protein Annotation ('relaxed', No 1)             #
# Annotate which proteins are DbTFs, coTFs, None or Both #
##########################################################

# Rule for coTFs: have at least one of the GO terms most likely associated with coTFs

# Read GO terms most likely associated with coTFs
cotf_go_map    = readr::read_tsv(file = 'data/coTF_GO_terms.tsv', col_types = 'cc')
cotf_go_terms  = cotf_go_map %>% pull(term)
total_go_terms = colnames(protein_go_mat)
total_go_terms_no_iea   = colnames(protein_go_no_iea_mat)
total_go_terms_only_iea = colnames(protein_go_only_iea_mat)

# check: all coTF GO terms are present in the GO matrices
stopifnot(all(cotf_go_terms %in% total_go_terms)) # yes
stopifnot(all(cotf_go_terms %in% total_go_terms_no_iea)) # yes
all(cotf_go_terms %in% total_go_terms_only_iea) # no
cotf_go_terms_only_iea = cotf_go_terms[cotf_go_terms %in% total_go_terms_only_iea] # 11 GO terms are present
all(cotf_go_terms_only_iea %in% total_go_terms_only_iea) # now yes

cotf_ids = rownames(protein_go_mat)[rowSums(protein_go_mat[,cotf_go_terms]) > 0] # 769 proteins
# with the GO matrix with no IEA-supported terms, coTFs are going to be less
cotf_ids_no_iea = rownames(protein_go_no_iea_mat)[rowSums(protein_go_no_iea_mat[,cotf_go_terms]) > 0] # 553 proteins
# and even less with the GO matrix with only IEA-supported terms
cotf_ids_only_iea = rownames(protein_go_only_iea_mat)[rowSums(protein_go_only_iea_mat[,cotf_go_terms_only_iea]) > 0] # 314 proteins

protein_class_1 = sapply(protein_ids, function(id) {
  if (id %in% dbtf_ids & id %in% cotf_ids) {
    return('Both')
  }
  if (id %in% dbtf_ids) {
    return('DbTF')
  }
  if (id %in% cotf_ids) {
    return('coTF')
  }
  return('None')
})
protein_class_1 %>% table()

protein_class_1_no_iea = sapply(protein_ids, function(id) {
  if (id %in% dbtf_ids & id %in% cotf_ids_no_iea) {
    return('Both')
  }
  if (id %in% dbtf_ids) {
    return('DbTF')
  }
  if (id %in% cotf_ids_no_iea) {
    return('coTF')
  }
  return('None')
})
protein_class_1_no_iea %>% table()

protein_class_1_only_iea = sapply(protein_ids, function(id) {
  if (id %in% dbtf_ids & id %in% cotf_ids_only_iea) {
    return('Both')
  }
  if (id %in% dbtf_ids) {
    return('DbTF')
  }
  if (id %in% cotf_ids_only_iea) {
    return('coTF')
  }
  return('None')
})
protein_class_1_only_iea %>% table()

########################################
# Unsupervised UMAP on InterPro matrix #
########################################
neighbors = c(2,6,10,14,20)
for (n_neighbors in neighbors) {
  print(paste0('Number of neighbors: ', n_neighbors))

  data_file = paste0('data/tfc2_umap_', n_neighbors, 'n_interpro.rds')
  if (!file.exists(data_file)) {
    set.seed(42)
    interpro_umap = uwot::umap(X = protein_domain_mat, n_threads = 4,
      n_neighbors = n_neighbors, metric = 'euclidean',
      verbose = TRUE)
    saveRDS(object = interpro_umap, file = data_file)
  }
}

# save plots (coloring with protein class annotation 1)
for (n_neighbors in neighbors) {
  print(paste0('Number of neighbors: ', n_neighbors))

  image_file = paste0('img/tfch2-InterPro/tfc2_umap_', n_neighbors, 'n_interpro.png')

  if (!file.exists(image_file)) {
    umap_file = paste0('data/tfc2_umap_', n_neighbors, 'n_interpro.rds')
    interpro_umap = readRDS(file = umap_file)

    interpro_umap %>%
      `colnames<-` (c("X", "Y")) %>%
      tibble::as_tibble() %>%
      tibble::add_column(protein_class = protein_class_1 %>% unname() %>% as.factor()) %>%
      ggplot(aes(x = X, y = Y, color = protein_class)) +
      geom_point(size = 0.1) +
      scale_colour_brewer(palette = "Set1") +
      guides(colour = guide_legend(title = "Class",
        label.theme = element_text(size = 12),
        override.aes = list(shape = 19, size = 12))) +
      labs(title = paste0("TFcheckpoint2 (InterPro) - UMAP (", n_neighbors," Neighbors)")) +
      theme_classic() +
      theme(plot.title = element_text(hjust = 0.5))
    ggsave(filename = image_file, width = 7, height = 5, dpi = 'print')
  }
}

#########################################
# Unsupervised UMAP on (full) GO matrix #
#########################################
neighbors = c(2,6,10,14,20)
for (n_neighbors in neighbors) {
  print(paste0('Number of neighbors: ', n_neighbors))

  data_file = paste0('data/tfc2_umap_', n_neighbors, 'n_go.rds')
  if (!file.exists(data_file)) {
    set.seed(42)
    go_umap = uwot::umap(X = protein_go_mat, n_threads = 4,
      n_neighbors = n_neighbors, metric = 'euclidean',
      verbose = TRUE)
    saveRDS(object = go_umap, file = data_file)
  }
}

# save plots (coloring with protein class annotation 1)
for (n_neighbors in neighbors) {
  print(paste0('Number of neighbors: ', n_neighbors))

  image_file = paste0('img/tfch2-GO/tfc2_umap_', n_neighbors, 'n_go.png')

  if (!file.exists(image_file)) {
    umap_file = paste0('data/tfc2_umap_', n_neighbors, 'n_go.rds')
    go_umap = readRDS(file = umap_file)

    go_umap %>%
      `colnames<-` (c("X", "Y")) %>%
      tibble::as_tibble() %>%
      tibble::add_column(
        protein_class = protein_class_1 %>% unname() %>% as.factor()) %>%
      ggplot(aes(x = X, y = Y, color = protein_class)) +
      geom_point(size = 0.1) +
      scale_colour_brewer(palette = "Set1") +
      guides(colour = guide_legend(title = "Class",
        label.theme = element_text(size = 12),
        override.aes = list(shape = 19, size = 12))) +
      labs(title = paste0("TFcheckpoint2 (GO) - UMAP (", n_neighbors," Neighbors)")) +
      theme_classic() +
      theme(plot.title = element_text(hjust = 0.5))
    ggsave(filename = image_file, width = 7, height = 5, dpi = 'print')
  }
}

#############################################################
# Unsupervised UMAP on combined GO (full) + InterPro matrix #
#############################################################
combined_mat = cbind(protein_go_mat, protein_domain_mat)

neighbors = c(2,6,10,14,20)
for (n_neighbors in neighbors) {
  print(paste0('Number of neighbors: ', n_neighbors))

  data_file = paste0('data/tfc2_umap_', n_neighbors, 'n_combined.rds')
  if (!file.exists(data_file)) {
    set.seed(42)
    combined_umap = uwot::umap(X = combined_mat, n_threads = 4,
      n_neighbors = n_neighbors, metric = 'euclidean',
      verbose = TRUE)
    saveRDS(object = combined_umap, file = data_file)
  }
}

# save plots (coloring with protein class annotation 1)
for (n_neighbors in neighbors) {
  print(paste0('Number of neighbors: ', n_neighbors))

  image_file = paste0('img/tfch2-combined/tfc2_umap_', n_neighbors, 'n_combined.png')

  if (!file.exists(image_file)) {
    umap_file = paste0('data/tfc2_umap_', n_neighbors, 'n_combined.rds')
    go_umap = readRDS(file = umap_file)

    go_umap %>%
      `colnames<-` (c("X", "Y")) %>%
      tibble::as_tibble() %>%
      tibble::add_column(
        protein_class = protein_class_1 %>% unname() %>% as.factor()) %>%
      ggplot(aes(x = X, y = Y, color = protein_class)) +
      geom_point(size = 0.1) +
      scale_colour_brewer(palette = "Set1") +
      guides(colour = guide_legend(title = "Class",
        label.theme = element_text(size = 12),
        override.aes = list(shape = 19, size = 12))) +
      labs(title = paste0("TFcheckpoint2 (GO full + InterPro) - UMAP (", n_neighbors," Neighbors)")) +
      theme_classic() +
      theme(plot.title = element_text(hjust = 0.5))
    ggsave(filename = image_file, width = 7, height = 5, dpi = 'print')
  }
}

###########################################
# Unsupervised UMAP on GO matrix (no IEA) #
###########################################
neighbors = c(2,6,10,14,20)
for (n_neighbors in neighbors) {
  print(paste0('Number of neighbors: ', n_neighbors))

  data_file = paste0('data/tfc2_umap_', n_neighbors, 'n_go_no_iea.rds')
  if (!file.exists(data_file)) {
    set.seed(42)
    go_umap = uwot::umap(X = protein_go_no_iea_mat, n_threads = 4,
      n_neighbors = n_neighbors, metric = 'euclidean',
      verbose = TRUE)
    saveRDS(object = go_umap, file = data_file)
  }
}

# save plots (coloring with protein class annotation 1)
for (n_neighbors in neighbors) {
  print(paste0('Number of neighbors: ', n_neighbors))

  image_file = paste0('img/tfch2-GO/tfc2_umap_', n_neighbors, 'n_go_no_iea.png')

  if (!file.exists(image_file)) {
    umap_file = paste0('data/tfc2_umap_', n_neighbors, 'n_go_no_iea.rds')
    go_umap = readRDS(file = umap_file)

    go_umap %>%
      `colnames<-` (c("X", "Y")) %>%
      tibble::as_tibble() %>%
      tibble::add_column(
        protein_class = protein_class_1_no_iea %>% unname() %>% as.factor()) %>%
      ggplot(aes(x = X, y = Y, color = protein_class)) +
      geom_point(size = 0.1) +
      scale_colour_brewer(palette = "Set1") +
      guides(colour = guide_legend(title = "Class",
        label.theme = element_text(size = 12),
        override.aes = list(shape = 19, size = 12))) +
      labs(title = paste0("TFcheckpoint2 (GO - no IEA) - UMAP (", n_neighbors," Neighbors)")) +
      theme_classic() +
      theme(plot.title = element_text(hjust = 0.5))
    ggsave(filename = image_file, width = 7, height = 5, dpi = 'print')
  }
}

#############################################
# Unsupervised UMAP on GO matrix (only IEA) #
#############################################
neighbors = c(2,6,10,14,20)
for (n_neighbors in neighbors) {
  print(paste0('Number of neighbors: ', n_neighbors))

  data_file = paste0('data/tfc2_umap_', n_neighbors, 'n_go_only_iea.rds')
  if (!file.exists(data_file)) {
    set.seed(42)
    go_umap = uwot::umap(X = protein_go_only_iea_mat, n_threads = 4,
      n_neighbors = n_neighbors, metric = 'euclidean',
      verbose = TRUE)
    saveRDS(object = go_umap, file = data_file)
  }
}

# save plots (coloring with protein class annotation 1)
for (n_neighbors in neighbors) {
  print(paste0('Number of neighbors: ', n_neighbors))

  image_file = paste0('img/tfch2-GO/tfc2_umap_', n_neighbors, 'n_go_only_iea.png')

  if (!file.exists(image_file)) {
    umap_file = paste0('data/tfc2_umap_', n_neighbors, 'n_go_only_iea.rds')
    go_umap = readRDS(file = umap_file)

    go_umap %>%
      `colnames<-` (c("X", "Y")) %>%
      tibble::as_tibble() %>%
      tibble::add_column(
        protein_class = protein_class_1_only_iea %>% unname() %>% as.factor()) %>%
      ggplot(aes(x = X, y = Y, color = protein_class)) +
      geom_point(size = 0.1) +
      scale_colour_brewer(palette = "Set1") +
      guides(colour = guide_legend(title = "Class",
        label.theme = element_text(size = 12),
        override.aes = list(shape = 19, size = 12))) +
      labs(title = paste0("TFcheckpoint2 (GO - only IEA) - UMAP (", n_neighbors," Neighbors)")) +
      theme_classic() +
      theme(plot.title = element_text(hjust = 0.5))
    ggsave(filename = image_file, width = 7, height = 5, dpi = 'print')
  }
}

##########################################################
# TFch2 Protein Annotation (No 2)                        #
# Annotate which proteins are DbTFs, coTFs, None or Both #
##########################################################

# Rule for coTFs: have the GO:0003712 term or one of its children terms OR be in the Tcof-DB [@Schmeier2017]

schmeier_cotf_ids = readr::read_lines(file = 'data/cotfs_from_tfcof_db_tfc2_master_16112021.tsv') # 952
stopifnot(all(schmeier_cotf_ids %in% protein_ids))

cotf_go_terms = readr::read_lines(file = 'data/coreg_GO_terms_16112021.tsv') # 6 GO co-regulator terms
cotf_go_terms[cotf_go_terms %in% colnames(protein_go_mat)] # no protein has GO:0016989 (6th entry)
cotf_go_terms[cotf_go_terms %in% colnames(protein_go_no_iea_mat)] # same as above

cotf_ids_from_child_terms        = rownames(protein_go_mat)[rowSums(protein_go_mat[,cotf_go_terms[1:5]]) > 0] # 494
cotf_ids_from_child_terms_no_iea = rownames(protein_go_no_iea_mat)[rowSums(protein_go_no_iea_mat[,cotf_go_terms[1:5]]) > 0] # 454

# union => OR logic
cotf_ids        = union(schmeier_cotf_ids, cotf_ids_from_child_terms) # 1185
cotf_ids_no_iea = union(schmeier_cotf_ids, cotf_ids_from_child_terms_no_iea) # 1164

protein_class_2 = sapply(protein_ids, function(id) {
  if (id %in% dbtf_ids & id %in% cotf_ids) {
    return('Both')
  }
  if (id %in% dbtf_ids) {
    return('DbTF')
  }
  if (id %in% cotf_ids) {
    return('coTF')
  }
  return('None')
})
protein_class_2 %>% table()

protein_class_2_no_iea = sapply(protein_ids, function(id) {
  if (id %in% dbtf_ids & id %in% cotf_ids_no_iea) {
    return('Both')
  }
  if (id %in% dbtf_ids) {
    return('DbTF')
  }
  if (id %in% cotf_ids_no_iea) {
    return('coTF')
  }
  return('None')
})
protein_class_2_no_iea %>% table()

##############################################
# Color GO UMAP using 2nd protein annotation #
##############################################

# full GO dataset
umap_file  = paste0('data/tfc2_umap_20n_go.rds')
go_umap    = readRDS(file = umap_file)
image_file = paste0('img/tfch2-GO/tfc2_umap_20n_go_class_2.png')

set1_col = RColorBrewer::brewer.pal(n = 4, 'Set1')

go_umap %>%
  `colnames<-` (c("X", "Y")) %>%
  tibble::as_tibble() %>%
  tibble::add_column(
    protein_class = protein_class_2 %>% unname() %>% as.factor()) %>%
  ggplot(aes(x = X, y = Y, color = protein_class)) +
  geom_point(size = 0.1) +
  scale_colour_brewer(palette = "Set1") +
  guides(colour = guide_legend(title = "Class",
    label.theme = element_text(size = 12),
    override.aes = list(shape = 19, size = 12))) +
  labs(title = paste0("TFcheckpoint2 (GO) - UMAP (20 Neighbors)")) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(filename = image_file, width = 7, height = 5, dpi = 'print')

# full GO dataset + boxify clusters + saving to file the cluster ids
image_file = paste0('img/tfch2-GO/tfc2_umap_20n_go_class_2_with_boxes.png')

go_umap %>%
  `colnames<-` (c("X", "Y")) %>%
  tibble::as_tibble() %>%
  tibble::add_column(
    protein_class = protein_class_2 %>% unname() %>% as.factor()) %>%
  ggplot(aes(x = X, y = Y, color = protein_class)) +
  geom_point(size = 0.1) +
  scale_colour_brewer(palette = "Set1") +
  guides(colour = guide_legend(title = "Class",
    label.theme = element_text(size = 12),
    override.aes = list(shape = 19, size = 12))) +
  labs(title = paste0("TFcheckpoint2 (GO) - UMAP (20 Neighbors)")) +
  # 1st cluster
  annotate("rect", xmin = -21, xmax = -19, ymin = -4, ymax = -2,
    alpha = 0.05, size = 0.3, color = 'black') +
  annotate("text", x = -20, y = 0, size = 8, label = "1") +
  # 2nd cluster
  annotate("rect", xmin = -17, xmax = -9, ymin = -10, ymax = -3,
    alpha = 0.05, size = 0.3, color = 'black') +
  annotate("text", x = -13, y = -1, size = 8, label = "2") +
  # 3rd cluster
  annotate("rect", xmin = 18, xmax = 20, ymin = 11.5, ymax = 13.5,
    alpha = 0.05, size = 0.3, color = 'black') +
  annotate("text", x = 16.5, y = 12.5, size = 8, label = "3") +
  # 4rd cluster
  annotate("rect", xmin = 18.3, xmax = 20.3, ymin = 15.5, ymax = 17.5,
    alpha = 0.05, size = 0.3, color = 'black') +
  annotate("text", x = 17, y = 16.5, size = 8, label = "4") +
  # 5th cluster
  annotate("rect", xmin = 23, xmax = 26, ymin = 16.5, ymax = 18.5,
    alpha = 0.05, size = 0.3, color = 'black') +
  annotate("text", x = 24.5, y = 20, size = 8, label = "5") +
  # 6th cluster
  annotate("rect", xmin = 21, xmax = 27, ymin = 10, ymax = 14.5,
    alpha = 0.05, size = 0.3, color = 'black') +
  annotate("text", x = 28.5, y = 12.5, size = 8, label = "6") +
  # 7th cluster
  annotate("rect", xmin = 26, xmax = 29, ymin = 6.6, ymax = 8.5,
    alpha = 0.05, size = 0.3, color = 'black') +
  annotate("text", x = 28, y = 5, size = 8, label = "7") +
  # 8th cluster
  annotate("rect", xmin = 21, xmax = 24, ymin = 4.5, ymax = 9,
    alpha = 0.05, size = 0.3, color = 'black') +
  annotate("text", x = 19.5, y = 7, size = 8, label = "8") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(filename = image_file, width = 7, height = 5, dpi = 'print')

cluster_data_tbl = go_umap %>%
  `colnames<-` (c("X", "Y")) %>%
  tibble::as_tibble() %>%
  tibble::add_column(
    protein_id    = protein_ids,
    gene_name     = gene_names,
    protein_class = protein_class_2 %>% unname() %>% as.factor()) %>%
  mutate(cluster_id = case_when(
    (X > -21  & X < -19 & Y > -4 & Y < -2) ~ 1,
    (X > -17  & X < -9  & Y > -10 & Y < -3 ) ~ 2,
    (X > 18   & X < 20   & Y > 11.5 & Y < 13.5) ~ 3,
    (X > 18.3 & X < 20.3 & Y > 15.5 & Y < 17.5 ) ~ 4,
    (X > 23   & X < 26  & Y > 16.5  & Y < 18.5 ) ~ 5,
    (X > 21   & X < 27  & Y > 10  & Y < 14.5 ) ~ 6,
    (X > 26   & X < 29  & Y > 6.6 & Y < 8.5) ~ 7,
    (X > 21 & X < 24 & Y > 4.5 & Y < 9 ) ~ 8,
    TRUE ~ 0 # there shouldn't be any point here!
  ))

# check: no '0' cluster ids
cluster_data_tbl %>% count(cluster_id)

# distribution of protein classes in cluster 2
cluster_data_tbl %>% filter(cluster_id == 2) %>% count(protein_class)

# save to CSV
cluster_data_tbl %>%
  dplyr::select(protein_id, gene_name, protein_class, cluster_id) %>%
  readr::write_csv(file = 'data/tfch2_umap_20n_GO_cluster_annot.csv')

# no IEA GO dataset
umap_file = paste0('data/tfc2_umap_20n_go_no_iea.rds')
go_umap = readRDS(file = umap_file)
image_file = paste0('img/tfch2-GO/tfc2_umap_20n_go_no_iea_class_2.png')

go_umap %>%
  `colnames<-` (c("X", "Y")) %>%
  tibble::as_tibble() %>%
  tibble::add_column(
    protein_class = protein_class_2_no_iea %>% unname() %>% as.factor()) %>%
  ggplot(aes(x = X, y = Y, color = protein_class)) +
  geom_point(size = 0.1) +
  scale_colour_brewer(palette = "Set1") +
  guides(colour = guide_legend(title = "Class",
    label.theme = element_text(size = 12),
    override.aes = list(shape = 19, size = 12))) +
  labs(title = paste0("TFcheckpoint2 (GO - no IEA) - UMAP (20 Neighbors)")) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(filename = image_file, width = 7, height = 5, dpi = 'print')

#######################################
# Color GO UMAP using DeepTF scoring  #
#######################################

# full GO dataset
umap_file  = paste0('data/tfc2_umap_20n_go.rds')
go_umap    = readRDS(file = umap_file)
image_file = paste0('img/tfch2-GO/tfc2_umap_20n_go_deeptf.png')

go_umap %>%
  `colnames<-` (c("X", "Y")) %>%
  tibble::as_tibble() %>%
  tibble::add_column(`DeepTF score` = deeptf_scores) %>%
  ggplot(aes(x = X, y = Y, color = `DeepTF score`)) +
  geom_point(size = 0.1) +
  # low => `None`, high = `DbTF`
  scale_colour_gradient2(low = set1_col[4], high = set1_col[3],
    na.value = 'black', midpoint = 0.5, guide = guide_colourbar(barheight = 10, draw.llim = TRUE)) +
  labs(title = paste0("TFcheckpoint2 (GO) - UMAP (20 Neighbors)")) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(filename = image_file, width = 7, height = 5, dpi = 'print')

# no IEA GO dataset
umap_file  = paste0('data/tfc2_umap_20n_go_no_iea.rds')
go_umap    = readRDS(file = umap_file)
image_file = paste0('img/tfch2-GO/tfc2_umap_20n_go_no_iea_deeptf.png')

go_umap %>%
  `colnames<-` (c("X", "Y")) %>%
  tibble::as_tibble() %>%
  tibble::add_column(`DeepTF score` = deeptf_scores) %>%
  ggplot(aes(x = X, y = Y, color = `DeepTF score`)) +
  geom_point(size = 0.1) +
  # low => `None`, high = `DbTF`
  scale_colour_gradient2(low = set1_col[4], high = set1_col[3],
    na.value = 'black', midpoint = 0.5, guide = guide_colourbar(barheight = 10, draw.llim = TRUE)) +
  labs(title = paste0("TFcheckpoint2 (GO - no IEA) - UMAP (20 Neighbors)")) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(filename = image_file, width = 7, height = 5, dpi = 'print')

#######
# END #
#######

if(FALSE) {
# Old code
#########################################################################
# Annotate UMAP result (GO dataset, 20 Neighbors) with 'boxed' clusters #
#########################################################################
umap_file = paste0('data/tfc2_umap_20n_go.rds')
go_umap_20n = readRDS(file = umap_file)

# Trying some clustering methods to find the clusters (DBSCAN seems to be the best!)
# km = kmeans(x = go_umap_20n, centers = 7, iter.max = 1000) # km$centers
# gmm = ClusterR::GMM(data = go_umap_20n, gaussian_comps = 7) # gmm#centroids
# dbs = dbscan::dbscan(x = go_umap_20n, eps = 0.7, minPts = 3) # dbs$cluster (can't get the cluster centers)
# We are still though going to annotate the cluster boxes using the eye-ball technique ;)

go_umap_20n %>%
  `colnames<-` (c("X", "Y")) %>%
  tibble::as_tibble() %>%
  tibble::add_column(
    #cluster       = dbs$cluster %>% as.factor(),
    protein_id    = protein_ids,
    protein_class = protein_class %>% unname() %>% as.factor()) %>%
  ggplot(aes(x = X, y = Y, color = protein_class)) +
  geom_point(size = 0.1) +
  scale_colour_brewer(palette = "Set1") +
  guides(colour = guide_legend(title = "Class",
    label.theme = element_text(size = 12),
    override.aes = list(shape = 19, size = 12))) +
  labs(title = paste0("TFcheckpoint2 (GO) - UMAP (", n_neighbors," Neighbors)")) +
  # 1st cluster
  annotate("rect", xmin = -15, xmax = -13, ymin = -13, ymax = -10,
    alpha = 0.05, size = 0.3, color = 'black') +
  annotate("text", x = -15.8, y = -9.5, size = 8, label = "1") +
  # 2nd cluster
  annotate("rect", xmin = -11, xmax = -3, ymin = -15, ymax = -5,
    alpha = 0.05, size = 0.3, color = 'black') +
  annotate("text", x = -10, y = -9.5, size = 8, label = "2") +
  # 3rd cluster
  annotate("rect", xmin = -0.8, xmax = 0, ymin = -13, ymax = -10,
    alpha = 0.05, size = 0.3, color = 'black') +
  annotate("text", x = -1.5, y = -9.5, size = 8, label = "3") +
  # 4rd cluster
  annotate("rect", xmin = 6, xmax = 13, ymin = 13, ymax = 22,
    alpha = 0.05, size = 0.3, color = 'black') +
  annotate("text", x = 6.8, y = 20.5, size = 8, label = "4") +
  # 5th cluster
  annotate("rect", xmin = 15, xmax = 17, ymin = 21, ymax = 24,
    alpha = 0.05, size = 0.3, color = 'black') +
  annotate("text", x = 14.2, y = 24, size = 8, label = "5") +
  # 6th cluster
  annotate("rect", xmin = 14, xmax = 16, ymin = 14, ymax = 19,
    alpha = 0.05, size = 0.3, color = 'black') +
  annotate("text", x = 17, y = 16.5, size = 8, label = "6") +
  # 7th cluster
  annotate("rect", xmin = 19, xmax = 21, ymin = 2.5, ymax = 4.5,
    alpha = 0.05, size = 0.3, color = 'black') +
  annotate("text", x = 18, y = 5, size = 8, label = "7") +
  # 8th cluster
  annotate("rect", xmin = 17.5, xmax = 19.5, ymin = -9, ymax = -7,
    alpha = 0.05, size = 0.3, color = 'black') +
  annotate("text", x = 18, y = -5, size = 8, label = "8") +
  # geom_point(data = km$centers %>% as.data.frame(), mapping = aes(x = V1, y = V2), size = 5, color = 'black') +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(filename = 'img/tfch2-GO/tfc2_umap_20n_GO_cluster_annot.png', width = 7, height = 5, dpi = 'print')

data_tbl = go_umap_20n %>%
  `colnames<-` (c("X", "Y")) %>%
  tibble::as_tibble() %>%
  tibble::add_column(
    protein_id    = protein_ids,
    gene_name     = gene_names,
    protein_class = protein_class %>% unname() %>% as.factor()) %>%
  mutate(cluster_id = case_when(
    (X > -15  & X < -13 & Y > -13 & Y < -10) ~ 1,
    (X > -11  & X < -3  & Y > -15 & Y < -5 ) ~ 2,
    (X > -0.8 & X < 0   & Y > -13 & Y < -10) ~ 3,
    (X > 6    & X < 13  & Y > 13  & Y < 22 ) ~ 4,
    (X > 15   & X < 17  & Y > 21  & Y < 24 ) ~ 5,
    (X > 14   & X < 16  & Y > 14  & Y < 19 ) ~ 6,
    (X > 19   & X < 21  & Y > 2.5 & Y < 4.5) ~ 7,
    (X > 17.5 & X < 19.5 & Y > -9 & Y < -7 ) ~ 8,
    TRUE ~ 0 # there shouldn't be any point here!
  ))

# check: no '0' cluster ids
data_tbl %>% count(cluster_id)

# distribution of protein classes in cluster 2
data_tbl %>% filter(cluster_id == 2) %>% count(protein_class)

# save to CSV
data_tbl %>%
  select(protein_id, gene_name, protein_class, cluster_id) %>%
  readr::write_csv(file = 'data/tfch2_umap_20n_GO_cluster_annot.csv')

}