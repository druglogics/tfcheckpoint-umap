library(dplyr)
library(tibble)
library(readr)
library(ggplot2)
library(uwot)

# genes-to-GO terms data
gg_data = readr::read_delim(file = 'data/genes2go_result_tfcheckpoint2_data.tsv',
  delim = '\t', skip = 2, progress = TRUE)

# protein list
protein_names = gg_data %>% pull(`Gene_Ids/GO_Terms`)

# remove column with gene names and make it a matrix
gg_mat = gg_data %>% select(-one_of("Gene_Ids/GO_Terms")) %>% as.matrix()

# for testing
# indexes = sample(x = 1:nrow(gg_mat), size = 100)

##########################################################
# Unsupervised UMAP + DNA binding annotation (neighbors) #
##########################################################
neighbors = c(2,4,6,8,10,12,14,15,17,20)
for (n_neighbors in neighbors) {
  print(paste0('Number of neighbors: ', n_neighbors))

  data_file = paste0('data/tf_umap_', n_neighbors, 'n.rds')
  if (!file.exists(data_file)) {
    set.seed(42)
    gg_umap = uwot::umap(X = gg_mat, n_threads = 4,
      n_neighbors = n_neighbors, metric = 'euclidean',
      verbose = TRUE)
    saveRDS(object = gg_umap, file = data_file)
  }
}

# save plots
for (n_neighbors in neighbors) {
  print(paste0('Number of neighbors: ', n_neighbors))

  image_file = paste0('img/tf_umap_', n_neighbors, 'n.png')

  if (!file.exists(image_file)) {
    umap_file = paste0('data/tf_umap_', n_neighbors, 'n.rds')
    gg_umap = readRDS(file = umap_file)
    gg_umap %>%
      `colnames<-` (c("X", "Y")) %>%
      tibble::as_tibble() %>%
      tibble::add_column(
        is_dna_binding = gg_data %>% pull(`DNA binding`) %>% as.factor()) %>%
      ggplot(aes(x = X, y = Y, color = is_dna_binding)) +
      geom_point(shape = '.') +
      scale_colour_brewer(palette = "Set1", labels = c("NO", "YES")) +
      guides(colour = guide_legend(title = "DNA Binding",
        label.theme = element_text(size = 12),
        override.aes = list(shape = 19, size = 12))) +
      labs(title = paste0("TFcheckpoint - UMAP (", n_neighbors," Neighbors)")) +
      theme_classic() +
      theme(plot.title = element_text(hjust = 0.5))
    ggsave(filename = image_file, width = 7, height = 5, dpi = 'print')
  }
}

##################################################################
# AHDC1: a DbTF in the non-DbTF supercluster                     #
# SHF: a set of DbTFs a little bit different then the rest DbTFs #
##################################################################

# UMAP coordinates with 12 neighbors
gg_umap = readRDS(file = 'data/tf_umap_12n.rds')

# tidy up data
gg_umap = gg_umap %>%
  `colnames<-` (c("X", "Y")) %>%
  tibble::as_tibble() %>%
  tibble::add_column(is_dna_binding = gg_data %>% pull(`DNA binding`) %>% as.factor()) %>%
  tibble::add_column(name = protein_names)

# X > 0 clearly defines the non-DbTF supercluster
gg_umap %>%
  filter(gg_umap$X > 0 & gg_umap$is_dna_binding == 1)
# AHDC1 with UMAP coordinates: (17.8, 4.41)

# -5 < X < 0 clearly identifies the extra DbTF that is a little bit more far away
# from the main cluster of DbTFs
gg_umap %>%
  filter(gg_umap$X > -5 & gg_umap$X < 0 & gg_umap$is_dna_binding == 1)
# HSF proteins with UMAP coordinates: (-3.99, 3.22)

print('UMPA, 12 neighbors + annotations of AHDC1 and HSF proteins')
image_file = 'img/12n_extra/tf_umap_12n_annot.png'
if (!file.exists(image_file)) {
  gg_umap %>%
    ggplot(aes(x = X, y = Y, color = is_dna_binding)) +
    geom_point(shape = '.') +
    # add `AHDC1` annotation
    annotate(
      geom = "curve", x = 10, y = 10, xend = 17.7, yend = 4.41,
      curvature = 0.3, arrow = arrow(length = unit(1.5, "mm"))
    ) +
    annotate(geom = "text", x = 10, y = 10.5, label = "AHDC1") +
    # add `HSF` annotation
    annotate(
      geom = "curve", x = 2, y = 5, xend = -3.7, yend = 3.22,
      curvature = -0.3, arrow = arrow(length = unit(1.5, "mm"))
    ) +
    annotate(geom = "text", x = 2, y = 5.6, label = "Heat shock factors (HSF)") +
    scale_colour_brewer(palette = "Set1", labels = c("NO", "YES")) +
    guides(colour = guide_legend(title = "DNA Binding",
      label.theme = element_text(size = 12),
      override.aes = list(shape = 19, size = 12))) +
    labs(title = paste0("TFcheckpoint - UMAP (12 Neighbors)")) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5))
  ggsave(filename = image_file, width = 7, height = 5, dpi = 'print')
  # ggsave(filename = 'img/tf_umap_12n_annot.svg', width = 7, height = 5) # only with library(svglite)
}

#########################################################
# Unsupervised UMAP + DNA binding annotation (min_dist) #
#########################################################
min_dists = c(0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 1)
n_neighbors = 15

for (min_dist in min_dists) {
  print(paste0('15 neighbors and min_dist = ', min_dist))

  data_file = paste0('data/tf_umap_15n_mindist_', min_dist, '.rds')
  if (!file.exists(data_file)) {
    set.seed(42)
    gg_umap = uwot::umap(X = gg_mat, n_threads = 4,
      n_neighbors = n_neighbors, metric = 'euclidean', min_dist = min_dist,
      verbose = TRUE)
    saveRDS(object = gg_umap, file = data_file)
  }
}

# save plots
for (min_dist in min_dists) {
  print(paste0('15 neighbors and min_dist = ', min_dist))

  image_file = paste0('img/tf_umap_15n_mindist_', min_dist, '.png')

  if (!file.exists(image_file)) {
    umap_file = paste0('data/tf_umap_15n_mindist_', min_dist, '.rds')
    gg_umap = readRDS(file = umap_file)
    gg_umap %>%
      `colnames<-` (c("X", "Y")) %>%
      tibble::as_tibble() %>%
      tibble::add_column(
        is_dna_binding = gg_data %>% pull(`DNA binding`) %>% as.factor()) %>%
      ggplot(aes(x = X, y = Y, color = is_dna_binding)) +
      geom_point(shape = '.') +
      scale_colour_brewer(palette = "Set1", labels = c("NO", "YES")) +
      guides(colour = guide_legend(title = "DNA Binding",
        label.theme = element_text(size = 12),
        override.aes = list(shape = 19, size = 12))) +
      labs(title = paste0("TFcheckpoint - UMAP (15 Neighbors, min_dist = ",
        min_dist, ")")) +
      theme_classic() +
      theme(plot.title = element_text(hjust = 0.5))
    ggsave(filename = image_file, width = 7, height = 5, dpi = 'print')
  }
}

#########################################################################
# Unsupervised UMAP + GREEKC curated DbTF list (vs GO DNA Binding term) #
#########################################################################

greekc_dbtfs = readr::read_csv(file = 'data/dbTF_gene_product_set_greekc.csv',
  col_names = c('name', 'id'), col_types = 'cc') %>% pull(name)

# all GREEKC annotated proteins are in the TFcheckpoint list
stopifnot(greekc_dbtfs %in% protein_names)

# UMAP coordinates with 12 neighbors
gg_umap = readRDS(file = 'data/tf_umap_12n.rds')

# tidy up data
gg_umap = gg_umap %>%
  `colnames<-` (c("X", "Y")) %>%
  tibble::as_tibble() %>%
  tibble::add_column(is_dna_binding = gg_data %>% pull(`DNA binding`) %>% as.factor()) %>%
  tibble::add_column(name = protein_names) %>%
  mutate(is_greekc_dbtf = ifelse(name %in% greekc_dbtfs, 1, 0)) %>% # add GREEKC DbTF annotation
  mutate(is_greekc_dbtf = is_greekc_dbtf %>% as.factor()) %>%
  mutate(cluster_id = case_when(X > 10 ~ 3, Y > -1 ~ 2, TRUE ~ 1)) %>% # add cluster id info
  mutate(cluster_id = cluster_id %>% as.factor())

# save data
umap_cluster_res_file = 'data/12n_dbtf_cluster_annot.tsv'
if (!file.exists(umap_cluster_res_file)) {
  readr::write_tsv(x = gg_umap %>% select(name, is_dna_binding, is_greekc_dbtf, cluster_id),
    file = umap_cluster_res_file, col_names = TRUE)
}

# GO-DbTFs vs GREEKC-DbTFs contingency table
dbtf_stats = table(gg_umap$is_dna_binding, gg_umap$is_greekc_dbtf,
  dnn = c('GO-DbTF', 'GREEKC-DbTF')) %>% as_tibble()
stats_file = 'data/dbtfs_go_greekc_stats.rds'
if (!file.exists(stats_file)) {
  saveRDS(object = dbtf_stats, file = stats_file)
}

# save plots
print('12 Neighbors: GREEK DbTFs and comparison with the GO term DNA Binding')
image_file = 'img/12n_extra/tf_umap_12n_greekc_dbtfs.png'
if (!file.exists(image_file)) {
  gg_umap %>%
    ggplot(aes(x = X, y = Y, color = is_greekc_dbtf)) +
    geom_point(shape = '.') +
    scale_colour_brewer(palette = "Set1", labels = c("NO", "YES")) +
    guides(colour = guide_legend(title = "GREEKC (DbTF)",
      label.theme = element_text(size = 12),
      override.aes = list(shape = 19, size = 12))) +
    labs(title = 'TFcheckpoint - UMAP (12 Neighbors)') +
    # 1st cluster
    annotate("rect", xmin = -17, xmax = -8, ymin = -8.5, ymax = -1, alpha = 0.05, size = 0.3, color = 'black') +
    annotate("text", x = -7, y = -5, size = 7, label = "1") +
    # 2nd cluster
    annotate("rect", xmin = -11, xmax = -3, ymin = -0.8, ymax = 4, alpha = 0.05, size = 0.3, color = 'black') +
    annotate("text", x = -2, y = 5, size = 7, label = "2") +
    # 3rd cluster
    annotate("rect", xmin = 12, xmax = 25.2, ymin = 1, ymax = 16, alpha = 0.05, size = 0.3, color = 'black') +
    annotate("text", x = 11, y = 1, size = 7, label = "3") +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5))
  ggsave(filename = image_file, width = 7, height = 5, dpi = 'print')

  data_00 = dbtf_stats %>% filter(`GO-DbTF` == 0, `GREEKC-DbTF` == 0) %>% pull(n)
  data_10 = dbtf_stats %>% filter(`GO-DbTF` == 1, `GREEKC-DbTF` == 0) %>% pull(n)
  data_01 = dbtf_stats %>% filter(`GO-DbTF` == 0, `GREEKC-DbTF` == 1) %>% pull(n)
  data_11 = dbtf_stats %>% filter(`GO-DbTF` == 1, `GREEKC-DbTF` == 1) %>% pull(n)

  # GREEKC DbTFs vs GO DNA-binding annotated proteins
  gg_umap %>%
    ggplot(aes(x = X, y = Y,
      color = interaction(is_dna_binding, is_greekc_dbtf, sep = '-', lex.order = TRUE))) +
    geom_point(shape = '.') +
    scale_colour_brewer(palette = "Set1", labels = c(paste0("0-0 (",data_00,")"),
      paste0("0-1 (", data_01, ")"), paste0("1-0 (", data_10, ")"),
      paste0("1-1 (", data_11, ")"))) +
    guides(colour = guide_legend(title = "DbTFs: GO vs GREEKC",
      label.theme = element_text(size = 12),
      override.aes = list(shape = 19, size = 12))) +
    labs(title = 'TFcheckpoint - UMAP (12 Neighbors)') +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5))
  ggsave(filename = 'img/12n_extra/tf_umap_12n_go_vs_greekc_dbtfs.png', width = 7, height = 5, dpi = 'print')
}

########################################################
# Unsupervised UMAP + GREEKC curated DbTF + coTFs list #
########################################################

greekc_cotfs = readr::read_csv(file = 'data/coTF_list.csv',
  col_names = c('id', 'name'), col_types = 'cc') %>% pull(name)

# not all coTFs are part of the TFcheckpoint dataset, so subset only to the ones that are
greekc_cotfs = greekc_cotfs[greekc_cotfs %in% protein_names]

# tidy up data (we use the UMAP coordinates with 12 neighbors from above)
gg_umap = gg_umap %>%
  mutate(is_greekc_cotf = ifelse(name %in% greekc_cotfs, 1, 0)) %>% # add GREEKC coTF annotation
  mutate(is_greekc_cotf = is_greekc_cotf %>% as.factor())

# GREEKC-DbTFs vs GREEKC-coTFs contingency table
db_cotf_stats = table(gg_umap$is_greekc_dbtf, gg_umap$is_greekc_cotf,
  dnn = c('GREEKC-DbTF', 'GREEKC-coTF')) %>% as_tibble()

stats_file = 'data/db_cotf_stats.rds'
if (!file.exists(stats_file)) {
  saveRDS(object = db_cotf_stats, file = stats_file)
}

data_00 = db_cotf_stats %>% filter(`GREEKC-coTF` == 0, `GREEKC-DbTF` == 0) %>% pull(n)
data_01 = db_cotf_stats %>% filter(`GREEKC-coTF` == 0, `GREEKC-DbTF` == 1) %>% pull(n)
data_10 = db_cotf_stats %>% filter(`GREEKC-coTF` == 1, `GREEKC-DbTF` == 0) %>% pull(n)
# zero
data_11 = db_cotf_stats %>% filter(`GREEKC-coTF` == 1, `GREEKC-DbTF` == 1) %>% pull(n)

print('GREEKC coTFs vs GREEKC DbTFs (12 neighbors)')
gg_umap %>%
  ggplot(aes(x = X, y = Y,
    color = interaction(is_greekc_cotf, is_greekc_dbtf, sep = '-', lex.order = TRUE))) +
  geom_point(shape = '.') +
  scale_colour_brewer(palette = "Set1", labels = c(paste0("Neither (",data_00,")"),
    paste0("DbTFs only (", data_01, ")"), paste0("co-TFs only (", data_10, ")"))) +
  guides(colour = guide_legend(title = "co-TFs vs DbTFs",
    label.theme = element_text(size = 12),
    override.aes = list(shape = 19, size = 12))) +
  labs(title = 'TFcheckpoint - UMAP (12 Neighbors)') +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(filename = 'img/12n_extra/tf_umap_12n_co_vs_dbtfs.png', width = 7, height = 5, dpi = 'print')

print('GREEKC coTFs vs GREEKC DbTFs (15 neighbors, different min_dist)')
for (min_dist in min_dists) {
  print(paste0('min_dist = ', min_dist))

  image_file = paste0('img/cotf_vs_dbtf/tf_umap_15n_mindist_', min_dist, '.png')

  if (!file.exists(image_file)) {
    umap_file = paste0('data/tf_umap_15n_mindist_', min_dist, '.rds')
    gg_umap = readRDS(file = umap_file)
    gg_umap %>%
      `colnames<-` (c("X", "Y")) %>%
      tibble::as_tibble() %>%
      tibble::add_column(is_dna_binding = gg_data %>% pull(`DNA binding`) %>% as.factor()) %>%
      tibble::add_column(name = protein_names) %>%
      mutate(is_greekc_dbtf = ifelse(name %in% greekc_dbtfs, 1, 0)) %>% # add GREEKC DbTF annotation
      mutate(is_greekc_dbtf = is_greekc_dbtf %>% as.factor()) %>%
      mutate(is_greekc_cotf = ifelse(name %in% greekc_cotfs, 1, 0)) %>% # add GREEKC coTF annotation
      mutate(is_greekc_cotf = is_greekc_cotf %>% as.factor()) %>%
      ggplot(aes(x = X, y = Y, color = interaction(is_greekc_cotf, is_greekc_dbtf, sep = '-', lex.order = TRUE))) +
      geom_point(shape = '.') +
      scale_colour_brewer(palette = "Set1", labels = c(paste0("Neither (",data_00,")"),
        paste0("DbTFs only (", data_01, ")"), paste0("co-TFs only (", data_10, ")"))) +
      guides(colour = guide_legend(title = "co-TFs vs DbTFs",
        label.theme = element_text(size = 12),
        override.aes = list(shape = 19, size = 12))) +
      labs(title = paste0('TFcheckpoint - UMAP (15 Neighbors, min_dist = ',
        min_dist, ')')) +
      theme_classic() +
      theme(plot.title = element_text(hjust = 0.5))
    ggsave(filename = image_file, width = 7, height = 5, dpi = 'print')
  }
}

###################
# Supervised UMAP #
###################
neighbors = c(6,8,10,12,14,20)

# target class: 'no-TF: 0, DbTF: 1, co-TF: 2'
target_class = sapply(protein_names, function(name) {
  is_dbtf = name %in% greekc_dbtfs
  is_cotf = name %in% greekc_cotfs
  if (!is_dbtf & !is_cotf) return(0) # none of the 2
  if (is_dbtf) return(1) # DbTF
  else return(2) # co-TF
}, USE.NAMES = FALSE) %>% as.factor()

for (n_neighbors in neighbors) {
  print(paste0('Number of neighbors (supervised UMAP): ', n_neighbors))

  data_file = paste0('data/tf_sumap_', n_neighbors, 'n.rds')
  if (!file.exists(data_file)) {
    set.seed(42)
    gg_umap = uwot::umap(X = gg_mat, n_threads = 4, y = target_class,
      n_neighbors = n_neighbors, metric = 'euclidean',
      verbose = TRUE)
    saveRDS(object = gg_umap, file = data_file)
  }
}

# save plots
for (n_neighbors in neighbors) {
  print(paste0('Number of neighbors: ', n_neighbors))

  image_file = paste0('img/sumap/tf_sumap_', n_neighbors, 'n.png')

  if (!file.exists(image_file)) {
    umap_file = paste0('data/tf_sumap_', n_neighbors, 'n.rds')
    gg_umap = readRDS(file = umap_file)
    gg_umap %>%
      `colnames<-` (c("X", "Y")) %>%
      tibble::as_tibble() %>%
      tibble::add_column(target = target_class) %>%
      ggplot(aes(x = X, y = Y, color = target)) +
      geom_point(shape = '.') +
      scale_colour_brewer(palette = "Set1", labels = c("no-TFs", "DbTFs", "co-TFs")) +
      guides(colour = guide_legend(title = "Target Class",
        label.theme = element_text(size = 12),
        override.aes = list(shape = 19, size = 12))) +
      labs(title = paste0("TFcheckpoint - Supervised UMAP (", n_neighbors," Neighbors)")) +
      theme_classic() +
      theme(plot.title = element_text(hjust = 0.5))
    ggsave(filename = image_file, width = 7, height = 5, dpi = 'print')
  }
}
