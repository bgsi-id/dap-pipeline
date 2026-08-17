import csv
import gzip
import os
import subprocess
from pathlib import Path

import pytest


IMAGE = "quay.io/biocontainers/bioconductor-dss:2.58.0--r45h01b2380_0@sha256:511fe5e1f84097aad44356d09b54fcfd999fb18129ba0f2090c1e5592d76049f"
REFACTORED_SCRIPT = Path(__file__).parents[1] / "bin" / "call_dmrs.R"


def get_git_head_script() -> str:
    """Retrieve original call_dmrs.R from git HEAD."""
    repo_root = Path(__file__).parents[3]
    result = subprocess.run(
        ["git", "show", "HEAD:pipelines/differentially_methylated_regions/bin/call_dmrs.R"],
        cwd=repo_root,
        capture_output=True,
        text=True,
        check=True,
    )
    # Attach libraries in baseline script so DSS unexported internal rowVars calls resolve
    script = result.stdout.replace(
        "suppressPackageStartupMessages(library(parallel))",
        "suppressPackageStartupMessages({ library(parallel); library(bsseq); library(DSS) })",
    )
    return script


def is_docker_available() -> bool:
    try:
        res = subprocess.run(["docker", "info"], capture_output=True)
        return res.returncode == 0
    except Exception:
        return False


def create_matrix_fixture(matrix_dir: Path, num_sites: int, patient_high: bool = True):
    matrix_dir.mkdir(parents=True, exist_ok=True)
    samples = ["C1", "C2", "C3", "P1", "P2", "P3"]
    groups = ["control", "control", "control", "patient", "patient", "patient"]

    # metadata.tsv
    with (matrix_dir / "metadata.tsv").open("w") as f:
        f.write("sample_id\tgroup\n")
        for s, g in zip(samples, groups):
            f.write(f"{s}\t{g}\n")

    # analysis.tsv
    with (matrix_dir / "analysis.tsv").open("w") as f:
        f.write("genomic_region\tchr1:1000-10000\n")

    # probes, coverage, modified, beta
    with gzip.open(matrix_dir / "probes.tsv.gz", "wt") as probes_f, \
         gzip.open(matrix_dir / "coverage.tsv.gz", "wt") as cov_f, \
         gzip.open(matrix_dir / "modified.tsv.gz", "wt") as mod_f, \
         gzip.open(matrix_dir / "beta.tsv.gz", "wt") as beta_f:

        probes_f.write("probe_id\tchrom\tposition\n")
        header = "probe_id\t" + "\t".join(samples) + "\n"
        cov_f.write(header)
        mod_f.write(header)
        beta_f.write(header)

        for i in range(num_sites):
            pos = 1000 + (i * 20)
            probe_id = f"chr1:{pos}:m:."
            probes_f.write(f"{probe_id}\tchr1\t{pos}\n")

            # Coverage: 20 across all samples
            cov_vals = [20] * 6
            # Controls: 2/20 (10%), Patients: 18/20 (90%) if patient_high
            if patient_high:
                mod_vals = [2, 2, 2, 18, 18, 18]
            else:
                mod_vals = [10, 10, 10, 10, 10, 10]

            cov_f.write(probe_id + "\t" + "\t".join(map(str, cov_vals)) + "\n")
            mod_f.write(probe_id + "\t" + "\t".join(map(str, mod_vals)) + "\n")
            betas = [f"{m/c:.8g}" for m, c in zip(mod_vals, cov_vals)]
            beta_f.write(probe_id + "\t" + "\t".join(betas) + "\n")


def run_container_dss(script_path: Path, matrix_dir: Path, out_dir: Path, params: list):
    """Run call_dmrs.R inside Docker container with specified parameters."""
    base_dir = script_path.parent
    cmd = [
        "docker", "run", "--rm",
        "-v", f"{base_dir}:{base_dir}:rw",
        "-u", f"{os.getuid()}:{os.getgid()}",
        IMAGE,
        "Rscript", str(script_path),
        str(matrix_dir), str(out_dir), *params,
    ]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    return proc


def read_tsv_dict(path: Path):
    with path.open("r") as f:
        reader = csv.reader(f, delimiter="\t")
        return list(reader)


@pytest.mark.skipif(not is_docker_available(), reason="Docker is not available")
class TestCallDmrsContainerComparison:
    """Verify refactored call_dmrs.R produces identical results to HEAD in the pinned Docker image."""

    @pytest.fixture(autouse=True)
    def setup_scripts(self, tmp_path):
        self.tmp_path = tmp_path
        self.old_script_path = tmp_path / "old_call_dmrs.R"
        self.old_script_path.write_text(get_git_head_script())
        self.refactored_script_path = tmp_path / "new_call_dmrs.R"
        self.refactored_script_path.write_text(REFACTORED_SCRIPT.read_text())

    def test_called_dmr_path_matches_original_and_checks_sign(self):
        matrix_dir = self.tmp_path / "matrix_called"
        create_matrix_fixture(matrix_dir, num_sites=150, patient_high=True)

        old_out = self.tmp_path / "old_out_called"
        new_out = self.tmp_path / "new_out_called"

        # Parameters matching standard Nextflow execution
        # matrix_dir out_dir partition genomic_region p_threshold delta min_length min_cpgs merge_distance sig_fraction equal_disp smoothing ncores
        dss_params = [
            "combined", "chr1:1000-10000",
            "0.05", "0.2", "50", "3", "50", "0.5", "true", "false", "1"
        ]

        old_proc = run_container_dss(self.old_script_path, matrix_dir, old_out, dss_params)
        assert old_proc.returncode == 0, f"Old script failed: {old_proc.stderr}\nSTDOUT: {old_proc.stdout}"

        new_proc = run_container_dss(self.refactored_script_path, matrix_dir, new_out, dss_params)
        assert new_proc.returncode == 0, f"New script failed: {new_proc.stderr}\nSTDOUT: {new_proc.stdout}"

        # 1. Compare TSVs
        for tsv_name in ["cpg-results.tsv", "dml-results.tsv", "dmr-results.tsv", "dmr-metrics.tsv", "top-dmr-profile.tsv", "provenance.tsv"]:
            old_file = old_out / tsv_name
            new_file = new_out / tsv_name
            assert old_file.exists(), f"Old missing {tsv_name}"
            assert new_file.exists(), f"New missing {tsv_name}"
            old_content = old_file.read_text()
            new_content = new_file.read_text()
            assert old_content == new_content, f"Mismatch in {tsv_name}:\nOLD:\n{old_content}\nNEW:\n{new_content}"

        # 2. Check patient - control sign in DMR results
        dmr_rows = read_tsv_dict(new_out / "dmr-results.tsv")
        assert len(dmr_rows) > 1, "Expected at least 1 DMR called"
        header = dmr_rows[0]
        diff_idx = header.index("diff.Methy")
        # diff.Methy should be positive ~ +0.8 (patient: 0.9 - control: 0.1)
        for row in dmr_rows[1:]:
            diff_val = float(row[diff_idx])
            assert diff_val > 0.5, f"Expected patient - control > 0, got {diff_val}"

        # 3. Check PNG plots exist and are non-empty
        for png_name in ["top-dmr-profile.png", "top-dmr-dss.png"]:
            old_png = old_out / png_name
            new_png = new_out / png_name
            assert old_png.exists() and old_png.stat().st_size > 1000
            assert new_png.exists() and new_png.stat().st_size > 1000

    def test_insufficient_sites_path_matches_original(self):
        matrix_dir = self.tmp_path / "matrix_insufficient"
        create_matrix_fixture(matrix_dir, num_sites=2, patient_high=True)

        old_out = self.tmp_path / "old_out_insufficient"
        new_out = self.tmp_path / "new_out_insufficient"

        dss_params = [
            "combined", "chr1:1000-10000",
            "0.05", "0.2", "50", "3", "50", "0.5", "true", "false", "1"
        ]

        old_proc = run_container_dss(self.old_script_path, matrix_dir, old_out, dss_params)
        assert old_proc.returncode == 0, f"Old script failed: {old_proc.stderr}\nSTDOUT: {old_proc.stdout}"

        new_proc = run_container_dss(self.refactored_script_path, matrix_dir, new_out, dss_params)
        assert new_proc.returncode == 0, f"New script failed: {new_proc.stderr}\nSTDOUT: {new_proc.stdout}"

        # Compare TSVs
        for tsv_name in ["cpg-results.tsv", "dml-results.tsv", "dmr-results.tsv", "dmr-metrics.tsv", "top-dmr-profile.tsv", "provenance.tsv"]:
            old_file = old_out / tsv_name
            new_file = new_out / tsv_name
            assert old_file.exists(), f"Old missing {tsv_name}"
            assert new_file.exists(), f"New missing {tsv_name}"
            old_content = old_file.read_text()
            new_content = new_file.read_text()
            assert old_content == new_content, f"Mismatch in {tsv_name}"

        # Check status is insufficient_sites
        metrics_text = (new_out / "dmr-metrics.tsv").read_text()
        assert "insufficient_sites" in metrics_text
        provenance_text = (new_out / "provenance.tsv").read_text()
        assert "insufficient_sites" in provenance_text
