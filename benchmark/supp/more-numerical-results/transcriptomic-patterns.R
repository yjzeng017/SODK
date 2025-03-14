# Supplementary: transcrptomic examples




# user must specify the correct path
root_path = "C:/Users/Youjie Zeng/Desktop/SODK"

############# load data ########
library('reticulate')


# install spatialDE tool
# if (!require("BiocManager", quietly = TRUE)){
#   install.packages("BiocManager")
#   BiocManager::install("spatialDE")
# }

library('spatialDE') # require miniconda environment


transcrptomic_pattern = function(data, sample_info, num_gene = NULL, C=3, length=1.5) {
  # data = Rep11_MOB_0
  # sample_info = MOB_sample_info
  set.seed(12)
  
  require('spatialDE')
  data = data[rowSums(data) >= 5, ]
  data = data[, row.names(sample_info)]
  sample_info$total_counts = colSums(data)
  X = sample_info[, c("x", "y")]
  
  # stabilize
  norm_expr <- stabilize(data)
  # regression out
  resid_expr <- regress_out(norm_expr, sample_info = sample_info)
  
  #For this example, run spatialDE on the first num_gene genes
  if(is.null(num_gene)) num_gene = dim(data)[1]
  sample_resid_expr = head(resid_expr, num_gene)
  results = spatialDE::run(sample_resid_expr, coordinates = X)
  de_results = results[results$qval < 0.05, ]
  
  # spatial patter
  sp = spatial_patterns(
    sample_resid_expr,
    coordinates = X,
    de_results = de_results,
    n_patterns = C, length = length
  )
  sp$coordinate = X
  sp$patterns_normalized = apply(sp$patterns, MARGIN = 2, FUN = function(x){(x - min(x))/(max(x) - min(x))})
  return(sp)
  
}


# load the BreastCancer data
# BC = readr::read_tsv('./benchmark/supp/more-numerical-results/Layer2_BC_count_matrix-1.tsv',
#                      show_col_types = FALSE)
# BC = t(BC[, 2:dim(BC)[2]])
# BC_sample_info = t(apply(a[, 1], MARGIN = 1, FUN=function(x){as.numeric(strsplit(x, 'x')[[1]])}))
# BC_sample_info = data.frame(BC_sample_info, colSums(BC))
# colnames(BC_sample_info) = c('x', 'y', 'total_counts')
# save(BC, file = './benchmark/supp/more-numerical-results/Layer2_BC_count_matrix.Rdata')
# save(BC_sample_info, file = './benchmark/supp/more-numerical-results/BC_sample_info.Rdata')

data("Rep11_MOB_0")
data("MOB_sample_info")

load(file = paste(root_path, '/benchmark/supp/more-numerical-results/Layer2_BC_count_matrix.Rdata', sep = ''))
load(file = paste(root_path, '/benchmark/supp/more-numerical-results/BC_sample_info.Rdata', sep = ''))

# identify the spatial patterns for MOB and BC
MOB_patterns = transcrptomic_pattern(Rep11_MOB_0, MOB_sample_info, num_gene = NULL, C = 2, length = 1.5)
BC_patterns = transcrptomic_pattern(BC, BC_sample_info, num_gene = NULL, C = 2,length = 1.5)
save(MOB_patterns, file = paste(root_path, '/benchmark/supp/more-numerical-results/MOB_patterns.Rdata', sep = ''))
save(BC_patterns, file = paste(root_path, '/benchmark/supp/more-numerical-results/BC_patterns.Rdata', sep = ''))


# plot the spatial patterns
df.plot = data.frame(
  x = c(rep(unlist(MOB_patterns$coordinate['x']), 2), rep(unlist(BC_patterns$coordinate['x']), 2)),
  y = c(rep(unlist(MOB_patterns$coordinate['y']), 2), rep(unlist(BC_patterns$coordinate['y']), 2)),
  gene_expression = c(MOB_patterns$patterns_normalized[, 1], MOB_patterns$patterns_normalized[, 2],
                      BC_patterns$patterns_normalized[, 1], BC_patterns$patterns_normalized[, 2]),
  
  pattern = c(rep(c('Pattern 1', 'Pattern 2'), each = dim(MOB_patterns$patterns_normalized)[1]),
              rep(c('Pattern 1', 'Pattern 2'), each = dim(BC_patterns$patterns_normalized)[1])),
  example = c(rep('Mouse olfactory bulb', 2*dim(MOB_patterns$patterns_normalized)[1]),
              rep('Breast Cancer', each = 2*dim(BC_patterns$patterns_normalized)[1]))
)


ggplot(data = df.plot, aes(x = x, y = y, color = gene_expression)) +
  geom_point(size = 7) +
  scale_color_viridis_c(option = 'C', direction = -1) +
  facet_grid(~example~pattern) + 
  theme_bw() + 
  theme(
    panel.background = element_rect(fill = "white"),
    plot.background = element_rect(fill = "white"),
    panel.border = element_rect(color = 'black', fill = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 16),
    legend.position = 'bottom',
    legend.title = element_text(size = 16, vjust = 1),
    legend.key.width = unit(1, "cm"), legend.key.height = unit(0.4, "cm")) + 
  labs(color = 'Relative Expression', x='Spatial dim 1', y='Spatial dim 2')

