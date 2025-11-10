# 5d_pgs_per_person_scatter.R
# Per-person PGS y=x scatter: Baseline vs DP (BMI and T2D).

source(file.path("R", "config.R"))
d5 <- stage_dirs("5_MethodComparison")
safe_png <- function(path, w=850, h=750) png_open(path, w, h)

load_scores <- function() {
  if (!exists("s_base_bmi")) source(file.path("5_MethodComparison", "R Code", "5b_predictive_scores_violin.R"), local = TRUE)
  list(
    BMI = list(base=s_base_bmi, dp=s_dp_bmi),
    T2D = list(base=s_base_t2d, dp=s_dp_t2d)
  )
}

plot_yx <- function(x, y, main, xlab, ylab, out, nmax=6000) {
  ok <- which(is.finite(x) & is.finite(y))
  if (length(ok) > nmax) ok <- sample(ok, nmax)
  safe_png(out, 850, 750)
  par(mar=c(5,5,4,2)+0.1)
  smoothScatter(x[ok], y[ok], xlab=xlab, ylab=ylab, main=main)
  abline(0,1,lty=2,col="gray40")
  dev.off()
}

sc <- load_scores()
plot_yx(sc$BMI$base, sc$BMI$dp,
        "Per-person PGS: Baseline vs DP (BMI)",
        "Baseline score (BMI)", "DP score (BMI)",
        file.path(d5$out_g, "perperson_yx_BMI_base_vs_DP.png"))

plot_yx(sc$T2D$base, sc$T2D$dp,
        "Per-person PGS: Baseline vs DP (T2D)",
        "Baseline score (T2D)", "DP score (T2D)",
        file.path(d5$out_g, "perperson_yx_T2D_base_vs_DP.png"))

writeLines(c("perperson_yx_BMI_base_vs_DP.png","perperson_yx_T2D_base_vs_DP.png"),
           file.path(d5$out_t, "5d_outputs.txt"))
