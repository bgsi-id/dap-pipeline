import sys
from pathlib import Path

import pytest

# Add bin directory to sys.path
bin_dir = str(Path(__file__).resolve().parent / "bin")
if bin_dir not in sys.path:
    sys.path.insert(0, bin_dir)

import variant_etl_control as control


def test_clinvar_stars():
    assert control.clinvar_stars_from_revstat("practice_guideline") == 4
    assert control.clinvar_stars_from_revstat("reviewed_by_expert_panel") == 3
    assert control.clinvar_stars_from_revstat("criteria_provided,_multiple_submitters,_no_conflicts") == 2
    assert control.clinvar_stars_from_revstat("criteria_provided,_single_submitter") == 1
    assert control.clinvar_stars_from_revstat("no_assertion_criteria_provided") == 0
    assert control.clinvar_stars_from_revstat(None) is None
    assert control.clinvar_stars_from_revstat(".") is None


def test_parse_clinvar_sig():
    assert control.parse_clinvar_sig("Pathogenic/Likely_pathogenic") == ["Pathogenic", "Likely_pathogenic"]
    assert control.parse_clinvar_sig("Benign,Likely_benign") == ["Benign", "Likely_benign"]
    assert control.parse_clinvar_sig(".") == []
    assert control.parse_clinvar_sig(None) == []


def test_parse_clinvar_traits():
    assert control.parse_clinvar_traits("Cardiomyopathy|Long_QT_syndrome") == ["Cardiomyopathy", "Long QT syndrome"]
    assert control.parse_clinvar_traits("not_specified|not_provided") == []


def test_parse_spliceai():
    # SpliceAI format: ALLELE|SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL
    raw = "T|GENE1|0.01|0.85|0.02|0.00|10|-20|3|4,T|GENE2|0.10|0.20|0.05|0.01|1|2|3|4"
    max_ds, max_ag, max_al, max_dg, max_dl = control.parse_spliceai(raw)
    assert max_ds == 0.85
    assert max_ag == 0.01
    assert max_al == 0.85
    assert max_dg == 0.02
    assert max_dl == 0.00


def test_detail_values():
    fields = [
        "Allele", "Consequence", "IMPACT", "SYMBOL", "Gene", "Feature_type", "Feature",
        "BIOTYPE", "EXON", "INTRON", "HGVSc", "HGVSp", "cDNA_position", "CDS_position",
        "Protein_position", "Amino_acids", "Codons", "Existing_variation", "DISTANCE",
        "STRAND", "FLAGS", "PICK", "MANE_SELECT", "CANONICAL", "SIFT", "PolyPhen",
        "HGNC_ID", "REVEL", "CADD_PHRED"
    ]
    csq_entry = (
        "A|missense_variant|MODERATE|BRCA1|ENSG00000012048|Transcript|ENST00000357654|"
        "protein_coding|11/24||ENST00000357654.8:c.5266dupC|ENSP00000350283.3:p.Gln1756ProfsTer74|"
        "5385|5266|1756|Q/P|caa/c-|||1||1|NM_007294.4|1|deleterious(0.01)|probably_damaging(0.98)|"
        "HGNC:1100|0.75|25.4"
    )
    res = control.detail_values(csq_entry, fields)
    assert res["hgnc"] == "BRCA1"
    assert res["hgnc_id"] == "HGNC:1100"
    assert res["vep_consequence"] == ["missense_variant"]
    assert res["vep_impact"] == "MODERATE"
    assert res["vep_transcript_id"] == "ENST00000357654"
    assert res["scores"]["revel"] == 0.75
    assert res["scores"]["cadd_phred"] == 25.4
    assert res["scores"]["sift_score"] == 0.01
    assert res["scores"]["polyphen_score"] == 0.98
    assert res["attributes"]["biotype"] == "protein_coding"
    assert res["attributes"]["mane_select"] == "NM_007294.4"


def test_parse_base_info_aliases():
    # Test primary fields
    info_primary = {
        "gnomad_af": "0.001",
        "gnomad_af_popmax": "0.002",
        "gnomad_nhomalt": "5",
        "CLNSIG": "Pathogenic",
        "CLNREVSTAT": "criteria_provided,_multiple_submitters,_no_conflicts",
        "CLNDN": "Breast-ovarian_cancer",
        "SpliceAI": "T|BRCA1|0.01|0.02|0.03|0.04|1|2|3|4",
    }
    base = control.parse_base_info(info_primary)
    assert base["gnomad_af"] == 0.001
    assert base["gnomad_af_popmax"] == 0.002
    assert base["gnomad_nhomalt"] == 5
    assert base["clinvar_sig"] == ["Pathogenic"]
    assert base["clinvar_stars"] == 2
    assert base["clinvar_traits"] == ["Breast-ovarian cancer"]
    assert base["spliceai_ds_max"] == 0.04

    # Test popmax and nhomalt aliases
    info_alias = {
        "gnomad_popmax_af": "0.003",
        "gnomad_nhomalts": "12",
    }
    base_alias = control.parse_base_info(info_alias)
    assert base_alias["gnomad_af_popmax"] == 0.003
    assert base_alias["gnomad_nhomalt"] == 12

    base_max_only = control.parse_base_info({"gnomad_af_max": "0.004"})
    assert base_max_only["gnomad_af"] == 0.004


def test_load_variants_anti_join_idempotency(monkeypatch, tmp_path):
    import json

    # Create dummy local root with manifest.json
    local_root = tmp_path / "published"
    sample_dir = local_root / "sample-test"
    sample_dir.mkdir(parents=True)
    manifest_file = sample_dir / "manifest.json"
    manifest_file.write_text(json.dumps({
        "sample_id": "sample-test",
        "assembly": "GRCh38",
        "datasets": {"variant_dim": {}, "alt_call": {}}
    }))

    # Mock load_clickhouse
    mock_load_ch_calls = []
    def fake_load_clickhouse(url, db, uri, region):
        mock_load_ch_calls.append((url, db, uri, region))

    import sys
    dap_variant_dir = str(Path(__file__).resolve().parents[3] / "dap-variant")
    if dap_variant_dir not in sys.path:
        sys.path.insert(0, dap_variant_dir)
    import variant_ingest.ingest
    monkeypatch.setattr(variant_ingest.ingest, "load_clickhouse", fake_load_clickhouse)

    queries_run = []
    def fake_post(self, url, **kwargs):
        content = kwargs.get("content", b"")
        if isinstance(content, bytes):
            sql_text = content.decode("utf-8")
        else:
            sql_text = str(content)
        queries_run.append(sql_text)

        class FakeResponse:
            text = "10\n" if "SELECT count() FROM s3" in sql_text else "0\n"
            def raise_for_status(self):
                pass
            def json(self):
                # When checking ledger or calls
                return {"data": []}
        return FakeResponse()

    monkeypatch.setattr(control.httpx.Client, "post", fake_post)

    out_file = tmp_path / "output.json"
    args = type("Args", (), {
        "sample_id": "sample-test",
        "local_root": str(local_root),
        "published_root": "s3://bucket/variants",
        "source_sha256": "abcdef1234567890abcdef1234567890abcdef1234567890abcdef1234567890",
        "region": "ap-southeast-3",
        "clickhouse_url": "http://mock-ch:8123",
        "database": "default",
        "output": str(out_file),
    })()

    # Run load_variants
    control.load_variants(args)

    # Verify single-owner load_clickhouse was invoked
    assert len(mock_load_ch_calls) == 1
    # Verify anti-join was present in call loading query
    call_inserts = [q for q in queries_run if "INSERT INTO default.variant_call" in q]
    assert len(call_inserts) == 1
    assert "LEFT ANTI JOIN" in call_inserts[0]
    assert "FROM default.variant_call" in call_inserts[0]
    assert "WHERE sample_id = 'sample-test'" in call_inserts[0]
    # Verify ledger was written
    ledger_inserts = [q for q in queries_run if "INSERT INTO default.variant_ingestion_ledger" in q]
    assert len(ledger_inserts) == 1


def test_load_variants_failure_before_ledger_retry_does_not_duplicate(monkeypatch, tmp_path):
    import json

    dap_variant_dir = str(Path(__file__).resolve().parents[3] / "dap-variant")
    if dap_variant_dir not in sys.path:
        sys.path.insert(0, dap_variant_dir)
    import variant_ingest.ingest

    local_root = tmp_path / "published"
    sample_dir = local_root / "sample-test"
    sample_dir.mkdir(parents=True)
    (sample_dir / "manifest.json").write_text(json.dumps({"sample_id": "sample-test"}))
    monkeypatch.setattr(variant_ingest.ingest, "load_clickhouse", lambda *args: None)

    state = {"ledger": False, "calls": set(), "fail_ledger_once": True, "call_inserts": 0}

    def fake_sql(client, url, database, statement):
        if "SELECT count() FROM default.variant_ingestion_ledger" in statement:
            return "1" if state["ledger"] else "0"
        if "SELECT count() FROM s3" in statement:
            return "1"
        if "INSERT INTO default.variant_call" in statement:
            state["call_inserts"] += 1
            # Model the LEFT ANTI JOIN against already committed calls.
            state["calls"].add(("release-1", "sample-test", "variant-1"))
            return ""
        if "INSERT INTO default.variant_ingestion_ledger" in statement:
            if state["fail_ledger_once"]:
                state["fail_ledger_once"] = False
                raise RuntimeError("simulated failure before ledger commit")
            state["ledger"] = True
        return ""

    monkeypatch.setattr(control, "sql", fake_sql)
    args = type("Args", (), {
        "sample_id": "sample-test",
        "local_root": str(local_root),
        "published_root": "s3://bucket/variants",
        "source_sha256": "a" * 64,
        "region": "ap-southeast-3",
        "clickhouse_url": "http://mock-ch:8123",
        "database": "default",
        "output": str(tmp_path / "output.json"),
    })()

    with pytest.raises(RuntimeError, match="failure before ledger"):
        control.load_variants(args)
    control.load_variants(args)
    control.load_variants(args)

    assert state["ledger"] is True
    assert state["calls"] == {("release-1", "sample-test", "variant-1")}
    assert state["call_inserts"] == 2


def test_annotation_loader_uses_canonical_hgnc_schema(monkeypatch, tmp_path):
    import gzip
    import json

    base_vcf = tmp_path / "base.vcf.gz"
    detail_vcf = tmp_path / "detail.vcf.gz"
    with gzip.open(base_vcf, "wt") as handle:
        handle.write("##fileformat=VCFv4.2\n")
        handle.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        handle.write("chr1\t1\tv1\tA\tT\t.\tPASS\tgnomad_af=0.001;CLNSIG=Pathogenic\n")
    with gzip.open(detail_vcf, "wt") as handle:
        handle.write("##fileformat=VCFv4.2\n")
        handle.write('##INFO=<ID=CSQ,Number=.,Type=String,Description="Format: Allele|Consequence|IMPACT|SYMBOL|HGNC_ID|PICK|REVEL">\n')
        handle.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        handle.write("chr1\t1\tv1\tA\tT\t.\tPASS\tCSQ=T|missense_variant|MODERATE|GENE1|HGNC:1|1|0.8\n")

    statements = []
    monkeypatch.setattr(control, "sql", lambda client, url, database, statement: statements.append(statement) or "")

    class DummyClient:
        def __init__(self, *args, **kwargs):
            pass
        def __enter__(self):
            return self
        def __exit__(self, *args):
            return False

    monkeypatch.setattr(control.httpx, "Client", DummyClient)
    args = type("Args", (), {
        "base_vcf": str(base_vcf), "detail_vcf": str(detail_vcf),
        "annotation_pack": "grch38-v3", "assembly": "GRCh38",
        "clickhouse_url": "http://mock-ch:8123", "database": "default",
        "completion_table": "variant_annotation_complete",
        "output": str(tmp_path / "receipt.json"),
    })()
    control.load_annotations(args)

    ddl = "\n".join(statement for statement in statements if statement.startswith("CREATE TABLE"))
    assert "hgnc Nullable(String)" in ddl
    assert "vep_consequence Array(String)" in ddl
    assert "gene_symbol" not in ddl
    inserted_rows = [
        json.loads(line)
        for statement in statements if "FORMAT JSONEachRow" in statement
        for line in statement.splitlines()[1:]
    ]
    detail = next(row for row in inserted_rows if row.get("hgnc") == "GENE1")
    assert detail["vep_consequence"] == ["missense_variant"]
    assert detail["scores"]["revel"] == 0.8
