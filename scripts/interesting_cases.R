library(dplyr)
library(tibble)
library(readr)
library(ggplot2)

# genes-to-GO terms data
gg_data = readr::read_delim(file = 'data/genes2go_result_tfcheckpoint2_data.tsv',
  delim = '\t', skip = 2, progress = TRUE)

# protein list
protein_names = gg_data %>% pull(`Gene_Ids/GO_Terms`)

# remove column with gene names and make it a matrix
gg_mat = gg_data %>% select(-one_of("Gene_Ids/GO_Terms")) %>% as.matrix()

# use the supervised UMAP results with 14 neighbors
gg_sumap = readRDS(file = 'data/tf_sumap_14n.rds')

# GREEKC DbTF and co-TF list
greekc_dbtfs = readr::read_csv(file = 'data/dbTF_gene_product_set_greekc.csv',
  col_names = c('name', 'id'), col_types = 'cc') %>% pull(name)
greekc_cotfs = readr::read_csv(file = 'data/coTF_list.csv',
  col_names = c('id', 'name'), col_types = 'cc') %>% pull(name)

# all GREEKC annotated proteins are in the TFcheckpoint list
stopifnot(greekc_dbtfs %in% protein_names)

# not all coTFs are part of the TFcheckpoint dataset, so subset only to the ones that are
greekc_cotfs = greekc_cotfs[greekc_cotfs %in% protein_names]

# make the target class: 'no-TF: 0, DbTF: 1, co-TF: 2'
target_class = sapply(protein_names, function(name) {
  is_dbtf = name %in% greekc_dbtfs
  is_cotf = name %in% greekc_cotfs
  if (!is_dbtf & !is_cotf) return(0) # none of the 2
  if (is_dbtf) return(1) # DbTF
  else return(2) # co-TF
}, USE.NAMES = FALSE) %>% as.factor()

gg_sumap = gg_sumap %>%
  `colnames<-` (c("X", "Y")) %>%
  tibble::as_tibble() %>%
  tibble::add_column(target = target_class) %>%
  tibble::add_column(name = protein_names)

set1_col = RColorBrewer::brewer.pal(n = 9, 'Set1')

gg_sumap %>%
  ggplot(aes(x = X, y = Y, color = target)) +
    geom_point(shape = '.') +
    scale_colour_brewer(palette = "Set1", labels = c("no-TFs", "DbTFs", "co-TFs")) +
    guides(colour = guide_legend(title = "Target Class",
      label.theme = element_text(size = 12),
      override.aes = list(shape = 19, size = 12))) +
    # DbTFs: 2 blocks
    annotate("rect", xmin = -21.5, xmax = -15.5, ymin = 14, ymax = 20, alpha = 0.05, size = 0.3, color = set1_col[2]) +
    annotate("text", x = -14, y = 20, size = 7, label = "1", color = set1_col[2]) +
    annotate("rect", xmin = -22, xmax = -18, ymin = 10, ymax = 13, alpha = 0.05, size = 0.3, color = set1_col[2]) +
    annotate("text", x = -16.5, y = 12, size = 7, label = "2", color = set1_col[2]) +
    # co-TFs: 1 block
    annotate("rect", xmin = -16, xmax = -7, ymin = -3, ymax = 8, alpha = 0.05, size = 0.3, color = set1_col[3]) +
    annotate("text", x = -6, y = 9, size = 7, label = "3", color = set1_col[3]) +
    # no-TFs: 2 blocks
    annotate("rect", xmin = -5, xmax = 6, ymin = -11, ymax = -2, alpha = 0.05, size = 0.3, color = set1_col[1]) +
    annotate("text", x = 7, y = -1, size = 7, label = "4", color = set1_col[1]) +
    annotate("rect", xmin = 18, xmax = 24, ymin = -10, ymax = -5, alpha = 0.05, size = 0.3, color = set1_col[1]) +
    annotate("text", x = 25, y = -4, size = 7, label = "5", color = set1_col[1]) +
    labs(title = paste0("TFcheckpoint - Supervised UMAP (14 Neighbors)")) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5))
ggsave(filename = 'img/sumap/tf_sumap_14n_annot.png', width = 7, height = 5, dpi = 'print')
#ggsave(filename = 'img/sumap/tf_sumap_14n_annot.svg', width = 7, height = 5)

# Find the odd members in each cluster as annotated in the above figure
# These proteins maybe misclassified and are worth looking at

# 1 (DbTF cluster)
## all 3 classes here
stopifnot(gg_sumap %>% filter(X > -21.5, X < -15.5, Y > 14, Y < 20) %>% distinct(target) %>% nrow() == 3)

cotfs_1 = gg_sumap %>%
  filter(X > -21.5, X < -15.5, Y > 14, Y < 20) %>%
  filter(target == 2) %>% pull(name)
readr::write_lines(x = cotfs_1, file = 'data/interesting_cases/cotfs_1.txt')

notfs_1 = gg_sumap %>%
  filter(X > -21.5, X < -15.5, Y > 14, Y < 20) %>%
  filter(target == 0) %>% pull(name)
readr::write_lines(x = notfs_1, file = 'data/interesting_cases/notfs_1.txt')

dbtfs_1 = gg_sumap %>%
  filter(X > -21.5, X < -15.5, Y > 14, Y < 20) %>%
  filter(target == 1) %>% pull(name)
readr::write_lines(x = dbtfs_1, file = 'data/interesting_cases/dbtfs_1.txt')

# 2 (DbTF cluster)
## all 3 classes here
stopifnot(gg_sumap %>% filter(X > -22, X < -18, Y > 10, Y < 13) %>% distinct(target) %>% nrow() == 3)

cotfs_2 = gg_sumap %>%
  filter(X > -22, X < -18, Y > 10, Y < 13) %>%
  filter(target == 2) %>% pull(name)
readr::write_lines(x = cotfs_2, file = 'data/interesting_cases/cotfs_2.txt')

notfs_2 = gg_sumap %>%
  filter(X > -22, X < -18, Y > 10, Y < 13) %>%
  filter(target == 0) %>% pull(name)
readr::write_lines(x = notfs_2, file = 'data/interesting_cases/notfs_2.txt')

dbtfs_2 = gg_sumap %>%
  filter(X > -22, X < -18, Y > 10, Y < 13) %>%
  filter(target == 1) %>% pull(name)
readr::write_lines(x = dbtfs_2, file = 'data/interesting_cases/dbtfs_2.txt')

# 3 (co-TF cluster)
## only co-TFs and some no-TFs here
stopifnot(gg_sumap %>% filter(X > -16, X < -7, Y > -3, Y < 8) %>% distinct(target) %>% nrow() == 2)

notfs_3 = gg_sumap %>%
  filter(X > -16, X < -7, Y > -3, Y < 8) %>%
  filter(target == 0) %>% pull(name)
readr::write_lines(x = notfs_3, file = 'data/interesting_cases/notfs_3.txt')

# 4 (no-TF cluster)
## all 3 classes here
stopifnot(gg_sumap %>% filter(X > -5, X < 6, Y > -11, Y < -2) %>% distinct(target) %>% nrow() == 3)

cotfs_4 = gg_sumap %>%
  filter(X > -5, X < 6, Y > -11, Y < -2) %>%
  filter(target == 2) %>% pull(name)
readr::write_lines(x = cotfs_4, file = 'data/interesting_cases/cotfs_4.txt')

dbtfs_4 = gg_sumap %>%
  filter(X > -5, X < 6, Y > -11, Y < -2) %>%
  filter(target == 1) %>% pull(name)
readr::write_lines(x = dbtfs_4, file = 'data/interesting_cases/dbtfs_4.txt')

# 5 (no-TF cluster)
## all 3 classes here
stopifnot(gg_sumap %>% filter(X > 18, X < 24, Y > -10, Y < -5) %>% distinct(target) %>% nrow() == 3)

cotfs_5 = gg_sumap %>%
  filter(X > 18, X < 24, Y > -10, Y < -5) %>%
  filter(target == 2) %>% pull(name)
readr::write_lines(x = cotfs_5, file = 'data/interesting_cases/cotfs_5.txt')

dbtfs_5 = gg_sumap %>%
  filter(X > 18, X < 24, Y > -10, Y < -5) %>%
  filter(target == 1) %>% pull(name)
readr::write_lines(x = dbtfs_5, file = 'data/interesting_cases/dbtfs_5.txt')
