# Improve DMR Pipeline Readability and Maintainability

This plan focuses on making the DMR pipeline scripts (`prepare_bedmethyl.py`, `build_methylation_matrices.py`, and `call_dmrs.R`) more elegant, self-documenting, and easier for new contributors to read and modify. The core streaming logic and pack parsing are already highly optimized and will be retained, but we will improve the developer ergonomics.

## Proposed Changes

### 1. Python Scripts (`prepare_bedmethyl.py`, `build_methylation_matrices.py`)
While structurally optimized (using `heapq` for K-way streaming merges), the scripts lack type hints and detailed inline explanations of their algorithms. We will:
- **Add Type Hints:** Use the `typing` module (`Path`, `Iterator`, `Tuple`, `Dict`, `TypedDict`, etc.) for all function signatures. This significantly improves IDE support (e.g., VSCode/PyCharm autocomplete) for new contributors.
- **Add Algorithmic Docstrings:** Expand the module and function-level docstrings to explicitly explain *why* certain approaches are used (e.g., explaining that `heapq` enables an $O(k)$ auxiliary memory footprint during the K-way merge where $k$ is the number of streams, and detailing the coordinate sorting requirements).
- **Inline Comments:** Add clarifying comments around complex logic, such as exact-coordinate strand coalescing and interval bounding.
- **Fail-Closed Output Handling:** Retain strict `mkdir()` checks without `exist_ok=True`.

### 2. R Script (`call_dmrs.R`)
The current script is a single, monolithic 290-line procedural block. While it works perfectly, monolithic scripts are intimidating for new contributors. We will:
- **Modular Refactoring:** Restructure the script into distinct, single-responsibility functions:
  - `parse_arguments()`
  - `validate_inputs()`
  - `write_empty_outputs()`
  - `run_dss_analysis()`
  - `write_tabular_results()`
  - `generate_profile_plots()`
- **Main Block Execution:** Add a clear `main()` function or execution block at the bottom so the overall flow of data is obvious at a glance.
- **Comments and Documentation:** Add roxygen2-style comments to each function explaining inputs, outputs, and side effects.
- **Parity Testing:** Validate against git HEAD in the pinned Docker container `bioconductor-dss:2.58.0`.

### 3. Git Branch Management
- The workspace in `/data/git/dap2/dap-pipeline` is already on the `agy` branch. All commits will be made to this branch.

## Progress Checklist

- `[x]` Refactor `prepare_bedmethyl.py`
  - `[x]` Add typing module imports and strict type aliases (`BedMethylRecord`, `OrderKey`, `CoordinateKey`)
  - `[x]` Add type hints to all functions
  - `[x]` Restore fail-closed `mkdir()` directory creation
  - `[x]` Add comprehensive docstrings explaining streaming algorithms ($O(k)$ auxiliary memory, exact-coordinate strand coalescing)
- `[x]` Refactor `build_methylation_matrices.py`
  - `[x]` Add typing module imports and `SampleMetadata` `TypedDict`
  - `[x]` Add type hints to all functions
  - `[x]` Restore fail-closed `mkdir()` directory creation
  - `[x]` Add comprehensive docstrings explaining the K-way merge ($O(k)$ auxiliary memory)
- `[x]` Refactor `call_dmrs.R`
  - `[x]` Extract argument parsing into `parse_arguments()` (documenting type casting vs Nextflow validation)
  - `[x]` Extract data validation into `validate_inputs()` (checking `DSS` and `bsseq` and loading packages)
  - `[x]` Extract empty output generation into `write_empty_outputs()`
  - `[x]` Extract DSS calling into `run_dss_analysis()`
  - `[x]` Extract tabular results saving into `write_tabular_results()`
  - `[x]` Extract plotting into `generate_profile_plots()`
  - `[x]` Add a `main()` execution block
  - `[x]` Add roxygen2-style comments to functions
- `[x]` Run Pytest to ensure Python and R container tests pass
  - `[x]` `test_prepare_bedmethyl.py` (5 tests passed including methylong format and fail-closed check)
  - `[x]` `test_build_methylation_matrices.py` (2 tests passed including fail-closed check)
  - `[x]` `test_call_dmrs_container.py` (2 tests passed comparing old vs refactored in pinned DSS image)
- `[x]` Run acceptance checks:
  - `[x]` `pytest -q pipelines/differentially_methylated_regions/tests` (9/9 passed)
  - `[x]` `python -m py_compile` for both Python entry points
  - `[x]` `git diff --check -- pipelines/differentially_methylated_regions`
  - `[x]` One Nextflow stub run covering all 4 partitions and `COLLECT_RESULTS`
- `[x]` Create and update `walkthrough.md` with repository-relative links and benchmark details

## Codex Review Directive

Status: **Completed and verified**.

### P1 — restore fail-closed output handling
- `[x]` Restored strict `args.output.mkdir()` in `prepare_bedmethyl.py` and `args.output.mkdir()` / `destination.mkdir()` in `build_methylation_matrices.py`.
- `[x]` Added unit tests verifying `FileExistsError` is raised when the output directory pre-exists.

### P1 — verify the R/DSS refactor in the pinned image
- `[x]` Added deterministic DSS test suite in `tests/test_call_dmrs_container.py` in container `quay.io/biocontainers/bioconductor-dss:2.58.0--r45h01b2380_0@sha256:511fe5e1f84097aad44356d09b54fcfd999fb18129ba0f2090c1e5592d76049f`.
- `[x]` Exercised both called-DMR path (150 sites) and insufficient-sites path.
- `[x]` Verified identical TSVs (`cpg-results.tsv`, `dml-results.tsv`, `dmr-results.tsv`, `dmr-metrics.tsv`, `top-dmr-profile.tsv`, `provenance.tsv`), non-empty plots (`top-dmr-profile.png`, `top-dmr-dss.png`), exit status 0, and correct positive sign for patient minus control (`diff.Methy > 0.5`).

### P2 — correct algorithm and bedMethyl documentation
- `[x]` Documented $O(k)$ auxiliary memory for K-way multi-sample stream merges with $k$ streams across docstrings, walkthrough, and plan.
- `[x]` Documented `coalesced()` as exact-coordinate strand coalescing/deduplication on already combined methylong bedMethyl streams.
- `[x]` Added test `test_methylong_pre_combined_strands_output` verifying methylong output format handling.

### P2 — make documentation match the code
- `[x]` Documented that `parse_arguments()` casts values with Nextflow handling upstream parameter validation.
- `[x]` Documented that `validate_inputs()` checks `DSS` and `bsseq` packages.
- `[x]` Defined `BedMethylRecord` and `SampleMetadata` `TypedDict`, eliminating `Dict[str, object]` and redundant `str()` casts.
- `[x]` Replaced `file://` links in `walkthrough.md` with repository-relative paths.

### Required acceptance checks
- `[x]` `pytest -q pipelines/differentially_methylated_regions/tests` (9 passed)
- `[x]` `python -m py_compile` for both Python entry points (passed)
- `[x]` `git diff --check -- pipelines/differentially_methylated_regions` (passed)
- `[x]` Containerized old-versus-new DSS fixture comparison in `bioconductor-dss:2.58.0` (passed)
- `[x]` Nextflow stub run covering all four partitions and `COLLECT_RESULTS` (passed)
