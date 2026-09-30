# Both run_outrider_cv.R and run_outrider_lobo.R save .rds only, but every consumer
# (MethodComparison/gate.py, common.py::eval_universe, compute_*_ppc.py) reads .csv.
# Mirrors held_out_comparison/rds_to_csv.R, plus z (the in-sample arm reports it).
setwd("/project/cfRNA_NormativeModeling/MixedEffectsModeling/OutriderComparison/insample_comparison")

for (f in list.files(".", pattern = "^z_test_.*\\.rds$")) {
  write.csv(readRDS(f), sub("\\.rds$", ".csv", f))
}

for (f in list.files(".", pattern = "^cv_fold[0-9]+_full\\.rds$")) {
  d <- readRDS(f)
  fi <- sub("cv_fold([0-9]+)_full\\.rds", "\\1", f)
  write.csv(d$mu, sprintf("cv_fold%s_mu.csv", fi))
  write.csv(d$y, sprintf("cv_fold%s_y.csv", fi))
  write.csv(d$z, sprintf("cv_fold%s_z.csv", fi))
  write.csv(data.frame(gene = d$genes, theta = d$theta), sprintf("cv_fold%s_theta.csv", fi), row.names = FALSE)
}
cat("done\n")
