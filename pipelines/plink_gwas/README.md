# Case/control GWAS (PLINK2)

This workflow accepts the governed `dap-pack-manifest/v2` produced for a named
`control` and `patient` cohort pair. It follows the Cloudfield GWASTutorial
structure: stage shared pVCF, convert to PLINK triples, pre-GWAS QC, KING/PCA,
Firth logistic association, and post-GWAS summary/plots. Cohort-specific pVCFs
are never written or retained.

The runtime is responsible for resolving and authorizing the two cohort packs.
The EKS workload identity must permit `s3:GetObject` for the release objects:
the staging process uses `aws s3 cp` before PLINK reads local files. Build and
publish `Dockerfile` before execution, then pass its immutable image digest as
`--gwas_image`. The workflow fails for an invalid manifest, overlapping
cohorts, missing pVCF roles, or missing pVCF sample IDs.
