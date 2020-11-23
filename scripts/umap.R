library(dplyr)
library(tibble)
library(ggplot2)
library(uwot)

# genes-to-GO terms data
gg_data = readr::read_delim(file = 'data/genes2go_result_tfcheckpoint2_data.txt',
  delim = '\t', skip = 2)

# remove column with gene names and make it a matrix
gg_mat = gg_data %>% select(-one_of("Gene_Ids/GO_Terms")) %>% as.matrix()

indexes = sample(x = 1:nrow(gg_mat), size = 100)

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
      labs(title = paste0("TFchechpoint - UMAP (", n_neighbors," Neighbors)")) +
      theme_classic() +
      theme(plot.title = element_text(hjust = 0.5))
    ggsave(filename = image_file, width = 7, height = 5, dpi = 'print')
  }
}

# find the protein that is a DNA binding TF, but is clustered with the non-binding ones
gg = gg_umap %>%
  `colnames<-` (c("X", "Y")) %>%
  tibble::as_tibble() %>%
  tibble::add_column(is_dna_binding = gg_data %>% pull(`DNA binding`) %>% as.factor())

pr_index = which(gg$Y < 0 & gg$is_dna_binding == 1)
gg_data %>% slice(pr_index) %>% pull(`Gene_Ids/GO_Terms`)
# [1] "AHDC1"

# play with `min_dist`
min_dists = c(0.05, 0.1, 0.3, 0.5, 0.7, 1)
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
      labs(title = paste0("TFchechpoint - UMAP (15 Neighbors, min_dist = ",
        min_dist, ")")) +
      theme_classic() +
      theme(plot.title = element_text(hjust = 0.5))
    ggsave(filename = image_file, width = 7, height = 5, dpi = 'print')
  }
}
