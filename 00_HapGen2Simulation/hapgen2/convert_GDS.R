library(SNPRelate)

in_gen    <- "chr22_EUR_sim_10k.controls.gen.gz"
in_sample <- "chr22_EUR_sim_10k.controls.sample"
out_gds   <- "chr22_sim_hapgen.gds"

snpgdsGEN2GDS(
  gen.fn         = in_gen,
  sample.fn      = in_sample,
  out.fn         = out_gds,
  chr.code       = "22",
  call.threshold = 0.90,
  snpfirstdim    = TRUE,
  verbose        = TRUE
)

snpgdsSummary(out_gds)
genofile <- snpgdsOpen(out_gds, readonly = TRUE)
on.exit(snpgdsClose(genofile), add = TRUE)
# ... (your downstream reads/tests here as needed)
