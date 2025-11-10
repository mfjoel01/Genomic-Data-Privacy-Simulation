# 5c_predictive_scores_pirate.R
# Pirate-style plots (jittered points + mean lines) for the same score sets.

source(file.path("R", "config.R"))
d5 <- stage_dirs("5_MethodComparison")
safe_png <- function(path, w=900, h=700) png_open(path, w, h)

load_scores <- function() {
  # Reuse objects created in 5b without recomputing (read back from environment by sourcing file)
  # If not present, source 5b.
  if (!exists("s_base_bmi")) source(file.path("5_MethodComparison", "R Code", "5b_predictive_scores_violin.R"), local = TRUE)
  list(
    BMI = list(Baseline=s_base_bmi, DP=s_dp_bmi, Federated=s_fed_bmi),
    T2D = list(Baseline=s_base_t2d, DP=s_dp_t2d, Federated=s_fed_t2d)
  )
}

plot_pirate <- function(vecs, title, outfile) {
  nm <- names(vecs); k <- length(vecs)
  safe_png(outfile, 1000, 700)
  par(mar=c(5,6,4,2)+0.1)
  y_rng <- range(unlist(vecs), na.rm = TRUE)
  plot(0,0,type="n", xlim=c(0.5, k+0.5), ylim=y_rng,
       xaxt="n", xlab="", ylab="Score", main=title)
  for (i in seq_len(k)) {
    v <- vecs[[i]]
    xj <- jitter(rep(i, length(v)), amount=0.2)
    points(xj, v, pch=16, cex=0.25, col=adjustcolor("black", 0.2))
    segments(i-0.35, mean(v, na.rm=TRUE), i+0.35, mean(v, na.rm=TRUE), lwd=3)
  }
  axis(1, at=seq_len(k), labels=nm, las=2)
  dev.off()
}

sc <- load_scores()
f1 <- file.path(d5$out_g, "pirate_scores_BMI.png")
plot_pirate(sc$BMI, "Predictive scores (pirate-style) — BMI", f1)

f2 <- file.path(d5$out_g, "pirate_scores_T2D.png")
plot_pirate(sc$T2D, "Predictive scores (pirate-style) — T2D", f2)

writeLines(c(basename(f1), basename(f2)),
           file.path(d5$out_t, "5c_outputs.txt"))
