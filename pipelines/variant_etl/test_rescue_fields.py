import subprocess
from pathlib import Path


SCRIPT = Path(__file__).parent / "bin" / "add_rescue_fields.pl"
PIPELINE = Path(__file__).parent / "main.nf"


def annotate(rows: list[str]) -> list[str]:
    header = [
        "##fileformat=VCFv4.2",
        '##INFO=<ID=CLNSIG,Number=.,Type=String,Description="ClinVar">',
        '##INFO=<ID=CLNSIGCONF,Number=.,Type=String,Description="Conflicts">',
        '##INFO=<ID=SpliceAI,Number=.,Type=String,Description="SpliceAI">',
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    ]
    result = subprocess.run(
        ["perl", str(SCRIPT)],
        input="\n".join([*header, *rows]) + "\n",
        text=True,
        capture_output=True,
        check=True,
    )
    return [line for line in result.stdout.splitlines() if not line.startswith("#")]


def test_rescue_fields_are_numeric_and_clinvar_specific():
    result = annotate(
        [
            "chr1\t1\tv1\tA\tT\t.\tPASS\tCLNSIG=Pathogenic;SpliceAI=T|G|0.01|0.20|0.70|0.02|1|2|3|4",
            "chr1\t2\tv2\tA\tG\t.\tPASS\tCLNSIG=Conflicting_classifications_of_pathogenicity;CLNSIGCONF=Benign(2)|Uncertain_significance(1)",
            "chr1\t3\tv3\tA\tC\t.\tPASS\tCLNSIG=Conflicting_classifications_of_pathogenicity;CLNSIGCONF=Likely_pathogenic(1)|Benign(1)",
            "chr1\t4\tv4\tA\tC\t.\tPASS\tSpliceAI=C|G|0.49|0|0|0|1|2|3|4",
        ]
    )

    assert "DAP_SPLICEAI_DS_MAX=0.7" in result[0]
    assert "DAP_CLINVAR_RESCUE" in result[0]
    assert "DAP_CLINVAR_RESCUE" not in result[1]
    assert "DAP_CLINVAR_RESCUE" in result[2]
    assert "DAP_SPLICEAI_DS_MAX=0.49" in result[3]


def test_rescue_headers_are_declared_once():
    source = (
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "chr1\t1\tv1\tA\tT\t.\tPASS\t.\n"
    )
    result = subprocess.run(
        ["perl", str(SCRIPT)], input=source, text=True, capture_output=True, check=True
    ).stdout
    assert result.count("##INFO=<ID=DAP_SPLICEAI_DS_MAX") == 1
    assert result.count("##INFO=<ID=DAP_CLINVAR_RESCUE") == 1


def test_pipeline_uses_normalized_rescue_fields_and_handles_empty_detail_set():
    source = PIPELINE.read_text()
    assert "process NORMALIZE_AF_FIELDS" in source
    assert "INFO/gnomad_nhomalts INFO/gnomad_nhomalt" in source
    assert "INFO/DAP_CLINVAR_RESCUE=1" in source
    assert "INFO/DAP_SPLICEAI_DS_MAX>=0.5" in source
    assert "INFO/SpliceAI ~" not in source
    assert "detail = INDEX_DETAIL_VCF(vep.vcf, selected.metrics)" in source
    assert 'if [ "\\${selected_count}" -gt 0 ]; then' in source
