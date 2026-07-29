# Frozen configuration for reference benchmark v1.1.
# Do not edit in place; create a new snapshot for changed settings.
list(study = "01_finemapping", architecture = "separated", development_settings = TRUE, 
    replicate_count = 10L, methods = c("st_bed_bayesc", "st_bed_bayesr", 
    "st_csr_sbayesc", "st_csr_sbayesr"), simulation = list(h2 = 0.2, 
        n_causal = 10L, base_seed = 2001L), causal_selection = list(
        min_distance_bp = 1000000L, min_maf = 0.05, max_maf = 0.5), 
    mcmc = list(nit = 500L, nburn = 250L, nthin = 1L, nchains = 1L, 
        ncores = 1L, seed_offset = 10000L, convergence_control = list(
            warn = FALSE)), credible_sets = list(coverage = 0.95, 
        min_r2 = 0.5, pip_cutoff = 0.001, locus_pip_cutoff = 0.01, 
        max_locus_distance = 1000000L), example_data = list(repository = "psoerensen/qgdata", 
        commit = "6cca5819e711d326cfb2614d7e9d9f34942612cd", 
        subdirectory = "simulated_human_data", files = c("human.bed", 
        "human.bim", "human.fam", "human.pheno", "human.covar"
        ), size_bytes = c(human.bed = 62500003, human.bim = 1882359, 
        human.fam = 117786, human.pheno = 92786, human.covar = 641513
        ), md5 = c(human.bed = "e89bea9a6cedd9eeef3fd0a5c807db81", 
        human.bim = "0105119b04c67b7ac7f66cc5e6680963", human.fam = "3c5db3d9eb7f3fc893c75f6f2b89836d", 
        human.pheno = "6a9e7cb1162e43999c170a363863176d", human.covar = "d06002aa2b1b79bdc4c0e92c21f27ae5"
        )))
