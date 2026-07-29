# ============================================================
# Complete 10-replicate benchmark reproduction
# ============================================================
# Study 01 uses simulated genotype data that are publicly accessible from a
# pinned psoerensen/qgdata revision. The setup downloads and checksum-validates
# the five files automatically, then caches them under results/local/, which is
# ignored by Git. Set SBLR_BENCH_GLIST to use an existing compatible Glist.
#
# targets::tar_make() runs only missing or outdated work. Do not normally delete
# _targets/: it contains the cache that prevents expensive recomputation.
# Exact numerical reproduction depends on the pinned files and checksums,
# package versions, seeds, platform, compiler, and numerical libraries.

Sys.setenv(
  SBLR_BENCH_STUDY = "01_finemapping",
  SBLR_BENCH_REPLICATES = "10"
)

targets::tar_make()

# Inspect compact implemented summary targets.
targets::tar_read(marker_metrics)
targets::tar_read(credible_set_summary)
targets::tar_read(credible_set_metrics)
targets::tar_read(computational_summary)
targets::tar_read(replicate_status)

# Optional cache and provenance inspection.
targets::tar_outdated()
targets::tar_meta()
