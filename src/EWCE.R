# Cell-type specific enrichment for LE hits
# Tutorial https://nathanskene.github.io/EWCE/articles/extended.html#create-a-celltypedataset
library(EWCE)
library(MAGMA.Celltyping)

# Tutorial to setup reticulate: https://github.com/ttimbers/intro-to-reticulate/blob/main/setup-instructions/setup-after-installing-python.md
library(reticulate)
library(SingleCellExperiment) 
library(scKirby)

# Number of bootstrap
reps <- 50000

pooltop100 <- read.csv2('~/Desktop/scProteomic/output_data/pool_mean_int_ranked.csv', sep = ',')
negtop100 <- read.csv2('~/Desktop/scProteomic/output_data/pool_mean_int_ranked_neg.csv', sep = ',')
postop100 <- read.csv2('~/Desktop/scProteomic/output_data/pool_mean_int_ranked_pos.csv', sep = ',')

hits <- pooltop100$GeneName[0:100]
poshits <- postop100$GeneName[0:600]
neghits <- negtop100$GeneName[0:600]

## FIRST PASS ANALYSIS USING DEFAULT SINGLE CELL DATA (MOUSE) ## 
# Load mouse hippocampus data
ctd <- ewceData::ctd()

cont_results <- EWCE::bootstrap_enrichment_test(sct_data = ctd,
                                                hits = hits, 
                                                sctSpecies = "mouse",
                                                genelistSpecies = "human",
                                                reps = reps,
                                                annotLevel = 1,
                                                geneSizeControl = TRUE)

try({
  plot_list <- EWCE::ewce_plot(total_res = cont_results$results,
                               mtc_method = "BH")
  print(plot_list$plain)
})

pdf(file='~/Desktop/scProteomic/figures/EWCE/EWCE_defaultMsHipp.pdf', width = 6, height = 4)
print(plot_list$plain)
dev.off()

### OTHER SC DATA - AIBS
ctd <- MAGMA.Celltyping::get_ctd('ctd_AIBS')

cont_results <- EWCE::bootstrap_enrichment_test(sct_data = ctd,
                                                hits = hits, 
                                                sctSpecies = "human",
                                                genelistSpecies = "human",
                                                reps = reps,
                                                annotLevel = 1,
                                                geneSizeControl = TRUE)

try({
  plot_list <- EWCE::ewce_plot(total_res = cont_results$results,
                               mtc_method = "BH")
  print(plot_list$plain)
})

pdf(file='~/Desktop/scProteomic/figures/EWCE/EWCE_AIBS.pdf', width = 6, height = 4)
print(plot_list$plain)
dev.off()


### OTHER SC DATA - KI
ctd <- MAGMA.Celltyping::get_ctd('ctd_allKI')

cont_results <- EWCE::bootstrap_enrichment_test(sct_data = ctd,
                                                hits = hits, 
                                                sctSpecies = "human",
                                                genelistSpecies = "human",
                                                reps = reps,
                                                annotLevel = 1,
                                                geneSizeControl = TRUE)

try({
  plot_list <- EWCE::ewce_plot(total_res = cont_results$results,
                               mtc_method = "BH")
  print(plot_list$plain)
})

pdf(file='~/Desktop/scProteomic/figures/EWCE/EWCE_KI.pdf', width = 6, height = 4)
print(plot_list$plain)
dev.off()

### OTHER SC DATA - KI
ctd <- MAGMA.Celltyping::get_ctd('ctd_BlueLake2018_FrontalCortexOnly')

cont_results <- EWCE::bootstrap_enrichment_test(sct_data = ctd,
                                                hits = hits, 
                                                sctSpecies = "human",
                                                genelistSpecies = "human",
                                                reps = reps,
                                                annotLevel = 1,
                                                geneSizeControl = TRUE)

try({
  plot_list <- EWCE::ewce_plot(total_res = cont_results$results,
                               mtc_method = "BH")
  print(plot_list$plain)
})

pdf(file='~/Desktop/scProteomic/figures/EWCE/EWCE_BlueLake.pdf', width = 6, height = 4)
print(plot_list$plain)
dev.off()

# BlueLake Level2

cont_results <- EWCE::bootstrap_enrichment_test(sct_data = ctd,
                                                hits = hits, 
                                                sctSpecies = "human",
                                                genelistSpecies = "human",
                                                reps = reps,
                                                annotLevel = 2,
                                                geneSizeControl = TRUE)

try({
  plot_list <- EWCE::ewce_plot(total_res = cont_results$results,
                               mtc_method = "BH")
  print(plot_list$plain)
})

pdf(file='~/Desktop/scProteomic/figures/EWCE/EWCE_BlueLake_lvl2.pdf', width = 6, height = 4)
print(plot_list$plain)
dev.off()

write.csv2(cont_results$results, file='~/Desktop/scProteomic/output_data/EWCE_bluelake_lvl2.csv')

# BlueLake Level2 NEG

cont_results <- EWCE::bootstrap_enrichment_test(sct_data = ctd,
                                                hits = neghits, 
                                                sctSpecies = "human",
                                                genelistSpecies = "human",
                                                reps = reps,
                                                annotLevel = 2,
                                                geneSizeControl = TRUE)

try({
  plot_list <- EWCE::ewce_plot(total_res = cont_results$results,
                               mtc_method = "BH")
  print(plot_list$plain)
})

pdf(file='~/Desktop/scProteomic/figures/EWCE/EWCE_BlueLake_lvl2_neghits.pdf', width = 6, height = 4)
print(plot_list$plain)
dev.off()

# BlueLake Level2 OPS

cont_results <- EWCE::bootstrap_enrichment_test(sct_data = ctd,
                                                hits = poshits, 
                                                sctSpecies = "human",
                                                genelistSpecies = "human",
                                                reps = reps,
                                                annotLevel = 2,
                                                geneSizeControl = TRUE)

try({
  plot_list <- EWCE::ewce_plot(total_res = cont_results$results,
                               mtc_method = "BH")
  print(plot_list$plain)
})

pdf(file='~/Desktop/scProteomic/figures/EWCE/EWCE_BlueLake_lvl2_pos.pdf', width = 6, height = 4)
print(plot_list$plain)
dev.off()


