select_separated_causal_markers <- function(glist, chromosome, marker_ids,
    n_causal, min_distance_bp, seed, min_maf = NULL, max_maf = NULL) {
  marker_ids <- as.character(marker_ids)
  if (!length(marker_ids) || anyNA(marker_ids) || anyDuplicated(marker_ids))
    stop("marker_ids must be unique and non-missing.", call. = FALSE)
  chromosome <- as.integer(chromosome)
  chromosome_ids <- glist$rsids[[chromosome]]
  index <- match(marker_ids, chromosome_ids)
  if (anyNA(index)) stop("Candidate markers are absent from the chromosome.", call. = FALSE)
  position <- glist$pos[[chromosome]][index]
  maf <- glist$maf[[chromosome]][index]
  if (any(!is.finite(position)) || any(!is.finite(maf)))
    stop("Aligned marker positions and MAF values must be finite.", call. = FALSE)
  eligible <- rep(TRUE, length(marker_ids))
  if (!is.null(min_maf)) eligible <- eligible & maf >= min_maf
  if (!is.null(max_maf)) eligible <- eligible & maf <= max_maf
  candidates <- which(eligible)
  set.seed(as.integer(seed))
  starts <- sample(candidates)
  selected <- integer()
  for (candidate in starts) {
    if (!length(selected) ||
        all(abs(position[candidate] - position[selected]) >= min_distance_bp))
      selected <- c(selected, candidate)
    if (length(selected) == as.integer(n_causal)) break
  }
  if (length(selected) != as.integer(n_causal))
    stop("Insufficient separated candidates for the requested causal design.",
      call. = FALSE)
  selected <- sort(selected)
  minimum <- if (length(selected) < 2L) Inf else min(abs(outer(
    position[selected], position[selected], "-")[lower.tri(matrix(0,
      length(selected), length(selected)))]))
  list(marker_ids = marker_ids[selected], marker_indices = index[selected],
    chromosome = chromosome, positions = position[selected], maf = maf[selected],
    pairwise_min_distance_bp = minimum, selection_seed = as.integer(seed),
    settings = list(n_causal = as.integer(n_causal),
      min_distance_bp = min_distance_bp, min_maf = min_maf, max_maf = max_maf))
}

study01_locus_table <- function(selection, spec) {
  data.frame(locus = selection$marker_ids, chromosome = selection$chromosome,
    causal_marker = selection$marker_ids, position_bp = selection$positions,
    left_boundary_bp = selection$positions - spec$locus_design$max_locus_distance,
    right_boundary_bp = selection$positions + spec$locus_design$max_locus_distance,
    credible_set_target = spec$locus_design$credible_set_target,
    stringsAsFactors = FALSE)
}
