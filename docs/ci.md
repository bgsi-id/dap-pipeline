# CI for dap-pipeline

`.github/workflows/ci.yml` runs the DevSecOps gates for the pipelines in this repository. Nothing is
built or deployed from here: DAP runs the Nextflow workflows pinned to a full Git commit SHA.

## What runs

| Job | What it does | Blocking |
| --- | --- | --- |
| `validate` | `catalog.yaml` contract check (`.github/scripts/validate_catalog.py`) and the tests of the CI scripts | yes |
| `tests` | `pytest` over `pipelines/` | yes, except the known failures below |
| `security` | Central `security-pipeline.yml`: SonarQube SAST, Semgrep SAST, Trivy filesystem SCA (also finds committed secrets). Findings go to SonarQube and DefectDojo | yes (Sonar quality gate, Semgrep ERROR, Trivy CRITICAL/HIGH) |
| `discover-images` | Lists every container image the pipelines reference (`.github/scripts/list_images.py`) | yes (a malformed reference fails) |
| `image-scan` | One Trivy scan per image through the central `image-scan-registry-pipeline.yml`; report artifact, step summary table, DefectDojo import | **no** for now (`non-blocking: true`) |

## Triggers

Pull requests to `dev`, pushes to `dev`, manual runs, and a weekly run (Monday 02:00 UTC) that only
rescans the images, because their CVEs change without a commit here. Pull requests from forks and
Dependabot receive no secrets: SonarQube and the DefectDojo imports are skipped with a warning, while
tests, validation, Semgrep and Trivy still run.

## Image scan

- Each image is resolved once to its `linux/amd64` manifest digest, which is printed in the log and in
  the summary, then scanned with Trivy straight from the registry.
- The gate is `CRITICAL,HIGH` with `gate-ignore-unfixed: true`. The report and the DefectDojo import
  keep every finding, including those without a fix.
- Images that are only tagged produce a warning. Set `ENFORCE_DIGEST_PINS` in `ci.yml` to `"true"`
  (a reviewed change) once they are pinned to digests.
- To make one image blocking, run it from its own entry with `non-blocking: false` after its baseline
  has no fixable CRITICAL/HIGH finding.

## Test baseline (measured on the first run of this work)

`17 passed, 2 skipped, 2 failed`.

- The two failures are `test_load_variants_anti_join_idempotency` and
  `test_load_variants_failure_before_ledger_retry_does_not_duplicate` in
  `pipelines/variant_etl/test_variant_etl_control.py`. They import `variant_ingest` from the
  `dap-variant` package, which lives in another repository. They run in a separate non-blocking step.
- The two skipped tests in `differentially_methylated_regions/tests/test_call_dmrs_container.py` need
  Docker. Their behaviour on a GitHub runner is confirmed on the first run and recorded here.

## Settings

No repository settings are needed: `SONAR_TOKEN`, `DEFECTDOJO_TOKEN` and `SONAR_HOST_URL` are
organization-level. The SonarQube project key and the DefectDojo product are both `dap-pipeline`.
The organization is on the GitHub Free plan, so required checks cannot be enforced yet. When they can,
require `validate`, `tests`, the three `security / ...` jobs and the blocking `image-scan` jobs.

## Rollback

Delete or disable `.github/workflows/ci.yml`. Nothing else depends on it.
