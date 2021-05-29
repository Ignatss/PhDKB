options(stringsAsFactors = FALSE)

library('dplyr')
library('sleuth')
library('ggplot2')
library('splines')

s2c <- read.csv("ConstantLight.csv", header = TRUE, stringsAsFactors = FALSE)


library(splines)
timepoint <- s2c$timepoint
full_design <- model.matrix(formula(~ ns(timepoint, df = 5)))
full_design

so <- sleuth_prep(s2c, full_model = full_design,)
so <- sleuth_fit(so)
so <- sleuth_fit(so, formula = ~ 1, fit_name = "reduced")
so <- sleuth_lrt(so, "reduced", "full")


Normalized_counts <- sleuth_to_matrix(so, 'obs_norm', 'est_counts')

write.csv(Normalized_counts, file = "Normalized.counts.csv")


lrt_results <- sleuth_results(so, 'reduced:full', test_type = 'lrt')




table(lrt_results[,"qval"] < 0.05)

lrt_results %>% head(n = 20) %>% dplyr::select(target_id, qval)

tmp <- so$obs_raw %>% dplyr::filter(target_id == 'Mp1g00050.1')
tmp <- dplyr::full_join(so$sample_to_covariates, tmp, by = 'sample')
tmp

tmp <- transform(tmp, timepoint = as.numeric(timepoint))
(ggplot(tmp, aes(x=timepoint, y=est_counts)) 
  + geom_point(shape=1) 
  + geom_smooth(method = loess))


plot_transcript_heatmap(so, head(lrt_results, n = 20)$target_id, 'est_counts')
plot_qq(so, test = 'reduced:full', test_type = 'lrt', sig_level = 0.05)
