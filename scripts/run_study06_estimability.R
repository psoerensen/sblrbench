#!/usr/bin/env Rscript
source(file.path("studies", "06_annotation_models", "estimability-and-contrasts.R"))
run_study06_estimability()
source(file.path("studies", "06_annotation_models", "estimability-finalize.R"))
finalize_study06_estimability()
