library(dplyr)
library(tibble)
library(ggplot2)
library(uwot)

# genes-to-GO terms data
gg_data = readr::read_delim(file = 'data/genes2go_result_tfcheckpoint2_data.txt',
  delim = '\t', skip = 2)

# remove column with gene names and make it a matrix
gg_mat = gg_data %>% select(-one_of("Gene_Ids/GO_Terms")) %>% as.matrix()

# for testing
# indexes = sample(x = 1:nrow(gg_mat), size = 100)

#############################################
# UMAP + DNA binding annotation (neighbors) #
#############################################
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
  tibble::add_column(name = gg_data %>% pull(`Gene_Ids/GO_Terms`))

# X > 0 clearly defines the non-DbTF supercluster
gg_umap %>%
  filter(gg_umap$X > 0 & gg_umap$is_dna_binding == 1)
# AHDC1 with UMAP coordinates: (17.8, 4.41)

# -5 < X < 0 clearly identifies the extra DbTF that is a little bit more far away
# from the main cluster of DbTFs
gg_umap %>%
  filter(gg_umap$X > -5 & gg_umap$X < 0 & gg_umap$is_dna_binding == 1)
# HSF proteins with UMAP coordinates: (-3.99, 3.22)

image_file = 'img/tf_umap_12n_annot.png'
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

############################################
# UMAP + DNA binding annotation (min_dist) #
############################################
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
      labs(title = paste0("TFchechpoint - UMAP (15 Neighbors, min_dist = ",
        min_dist, ")")) +
      theme_classic() +
      theme(plot.title = element_text(hjust = 0.5))
    ggsave(filename = image_file, width = 7, height = 5, dpi = 'print')
  }
}
