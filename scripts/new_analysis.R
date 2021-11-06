library(dplyr)
library(tibble)
library(tidyr)
library(readr)
library(uwot)
library(ggplot2)

##########################################################
# Read Interpro domain + GO annotations for all proteins #
# in TFch.v2 and build expanded 1/0 data matrices        #
##########################################################
tfc2_data = readr::read_tsv(file = 'data/tfc2_all_human_go_interpro_22102021_ALL.tsv',
  col_types = 'ccc') # 3875
protein_ids = tfc2_data %>% pull(Entry)
domains     = tfc2_data %>% pull(`Cross-reference (InterPro)`)
go_terms    = tfc2_data %>% pull(`Gene ontology IDs`)

# If matrices have already been stored as compressed R objects, don't rebuild them!
protein_domain_file = 'data/protein_domain_mat.rds'
protein_go_file     = 'data/protein_go_mat.rds'
if ((!file.exists(protein_domain_file)) | (!file.exists(protein_go_file))) {
  # get the unique InterPro domain IDs
  domains_unique = sapply(domains, function(domain_vec) {
    if (!is.na(domain_vec)) {
      vec = stringr::str_split(domain_vec, ";") %>% unlist() %>% head(-1)
      return(vec)
    }
  }) %>% unlist(use.names = FALSE) %>% unique()

  # get the unique GO terms
  go_terms_unique = sapply(go_terms, function(go_term_vec) {
    if (!is.na(go_term_vec)) {
      vec = stringr::str_split(go_term_vec, "; ") %>% unlist()
      return(vec)
    }
  }) %>% unlist(use.names = FALSE) %>% unique()

  # build matrices of 0's and 1's
  protein_domain_mat = matrix(0, nrow = length(protein_ids), ncol = length(domains_unique),
    dimnames = list(protein_ids, domains_unique))
  protein_go_mat     = matrix(0, nrow = length(protein_ids), ncol = length(go_terms_unique),
    dimnames = list(protein_ids, go_terms_unique))

  for(id in protein_ids) {
    # get the InterPro domains for that protein ID
    domain_vec = tfc2_data %>%
      filter(Entry == id) %>%
      pull(`Cross-reference (InterPro)`)

    # ... and set them to 1 in the matrix
    if (!is.na(domain_vec)) {
      vec = stringr::str_split(domain_vec, ";") %>% unlist() %>% head(-1)
      for(domain in vec) {
          protein_domain_mat[id, domain] = 1
      }
    }
  }

  for(id in protein_ids) {
    # get the GO terms for that protein ID
    go_term_vec = tfc2_data %>%
      filter(Entry == id) %>%
      pull(`Gene ontology IDs`)

    # ... and set them to 1 in the matrix
    if (!is.na(go_term_vec)) {
      vec = stringr::str_split(go_term_vec, "; ") %>% unlist()
      for(go_term in vec) {
        protein_go_mat[id, go_term] = 1
      }
    }
  }

  # save the matrices
  saveRDS(object = protein_domain_mat, file = protein_domain_file)
  saveRDS(object = protein_go_mat, file = protein_go_file)
} else {
  protein_domain_mat = readRDS(file = protein_domain_file)
  protein_go_mat     = readRDS(file = protein_go_file)
}

# Check: matrix 1-sparsity or percentage of 1's (GO matrix has more 1's)
sum(protein_go_mat)/(dim(protein_go_mat)[1]*dim(protein_go_mat)[2])
sum(protein_domain_mat)/(dim(protein_domain_mat)[1]*dim(protein_domain_mat)[2])

####################################
# Read GREEKC DbTF list            #
# See commit 6cd92c8 for more info #
####################################
greekc_dbtfs_tbl = readr::read_csv(file = 'data/dbTF_online_02112021.csv',
  col_select = c("id", "symbol"))
dbtf_ids = greekc_dbtfs_tbl %>% pull(id) # 1437 DbTF curated proteins

# check: all GREEKC annotated DbTF proteins are in TFch2
stopifnot(dbtf_ids %in% protein_ids)

###################################################
# Read GO terms most likely associated with coTFs #
###################################################
cotf_go_map    = readr::read_tsv(file = 'data/coTF_GO_terms.tsv', col_types = 'cc')
cotf_go_terms  = cotf_go_map %>% pull(term)
total_go_terms = colnames(protein_go_mat)

# check if we can find coTF-based GO terms that are not annotated
cotf_go_terms[!cotf_go_terms %in% total_go_terms] # GO:0001098, need to remove it!
cotf_go_terms = cotf_go_terms[cotf_go_terms %in% total_go_terms]

### Some data checks
if(FALSE) {
  # are there dbTFs that are annotated as coTFs?
  test_tbl = protein_go_mat %>%
    as_tibble(rownames = 'id') %>%
    filter(id %in% dbtf_ids) %>% # keep only the DbTFs
    select(c('id', cotf_go_terms)) # keep only coTF GO terms

  # How many DbTFs are annotated with at least one coTF GO term?
  test_tbl %>% select(cotf_go_terms) %>% as.matrix() %>% colSums()

  # Show the amount of co-TF GO terms per DbTF protein in decreasing order
  test_tbl %>% select(cotf_go_terms) %>% as.matrix() %>% rowSums() %>% sort(decreasing = TRUE) %>% head()

  # 240 DbTFs have at least one co-TF GO term
  test_tbl %>% select(cotf_go_terms) %>% as.matrix() %>% rowSums() %>% sort(decreasing = TRUE) %>% head(n = 240)
}

# Annotate which proteins are DbTFs, coTFs, None or Both
cotf_ids = rownames(protein_go_mat)[rowSums(protein_go_mat[,cotf_go_terms]) > 0] # 848 proteins

protein_class = sapply(protein_ids, function(id) {
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

# save plots
for (n_neighbors in neighbors) {
  print(paste0('Number of neighbors: ', n_neighbors))

  image_file = paste0('img/tfch2-InterPro/tfc2_umap_', n_neighbors, 'n_interpro.png')

  if (!file.exists(image_file)) {
    umap_file = paste0('data/tfc2_umap_', n_neighbors, 'n_interpro.rds')
    interpro_umap = readRDS(file = umap_file)

    interpro_umap %>%
      `colnames<-` (c("X", "Y")) %>%
      tibble::as_tibble() %>%
      tibble::add_column(
        protein_id    = protein_ids,
        protein_class = protein_class %>% unname() %>% as.factor()) %>%
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

##################################
# Unsupervised UMAP on GO matrix #
##################################
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

# save plots
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
        protein_id    = protein_ids,
        protein_class = protein_class %>% unname() %>% as.factor()) %>%
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

######################################################
# Unsupervised UMAP on combined GO + InterPro matrix #
######################################################
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

# save plots
for (n_neighbors in neighbors) {
  print(paste0('Number of neighbors: ', n_neighbors))

  image_file = paste0('img/tfch2-combined/tfc2_umap_', n_neighbors, 'n_combined.png')

  if (!file.exists(image_file)) {
    umap_file = paste0('data/tfc2_umap_', n_neighbors, 'n_go.rds')
    go_umap = readRDS(file = umap_file)

    go_umap %>%
      `colnames<-` (c("X", "Y")) %>%
      tibble::as_tibble() %>%
      tibble::add_column(
        protein_id    = protein_ids,
        protein_class = protein_class %>% unname() %>% as.factor()) %>%
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
  select(protein_id, protein_class, cluster_id) %>%
  readr::write_csv(file = 'data/tfc2_umap_20n_GO_cluster_annot.csv')
