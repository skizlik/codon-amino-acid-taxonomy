"""Reading, repair, filtering, features and groups (the 2022 defect guarded against is noted per test)."""
import numpy as np
import pytest

import codon_utils as cu


def test_comma_shifted_rows_are_repaired(raw):
    """2022: UUU/UUC coerced to 0 in two rows, every other codon shifted."""
    assert raw.shape == (13028, 72)
    assert not raw[cu.CODONS].isna().any().any()
    rep = raw[raw.repaired]
    assert list(rep.index) == [486, 5063]
    assert rep.SpeciesName.str.contains(",").all()
    assert rep.loc[486, "UUU"] == pytest.approx(0.04362)  # was 0 in 2022
    assert rep.loc[486, "UGA"] == pytest.approx(3 / 1238, abs=5e-6)  # recovered from Ncodons
    assert rep.loc[5063, ["UAG", "UGA"]].tolist() == [0.0, 0.0]
    sums = raw[cu.CODONS].sum(axis=1)
    assert sums.between(0.9995, 1.0005).all(), sums.agg(["min", "max"])


def test_repair_row_recovers_one_missing_codon():
    # 4 codons, Ncodons 10, counts 1, 2, 3 present: the 4th must be 4 / 10
    fields = ["vrl", "0", "1", "10", "Virus A", "strain B", "0.10000", "0.20000", "0.30000"]
    out, n_missing = cu.repair_row(fields)
    assert n_missing == 1 and out[4] == "Virus A,strain B" and out[-1] == "0.40000"


def test_repair_row_refuses_ambiguous_rows():
    # two codons missing and 7 counts unaccounted for: cannot split 7 between them
    fields = ["vrl", "0", "1", "10", "Virus", "b", "c", "0.10000", "0.20000"]
    with pytest.raises(ValueError, match="missing codons"):
        cu.repair_row(fields)


def test_dnatype_9_is_a_code_and_filter_is_explicit(raw, clean):
    """2022: code 9 became the string 'NA', which read_csv turned into a missing value."""
    assert (raw.DNAtype == 9).sum() == 2 and (raw.dna_type == "unknown").sum() == 2
    df, log = clean
    assert log.rows.tolist() == [13028, 13010, 12964, 12963]
    assert log.removed.sum() == 65
    assert df.kingdom.nunique() == 10
    assert set(df.dna_type) == {"genomic", "mitochondrial", "chloroplast"}


def test_amino_features_partition_the_codons(clean):
    """2022: a 'start' column duplicated methionine (both AUG): 22 columns, mass > 1."""
    df, _ = clean
    am = cu.amino_features(df)
    assert list(am.columns) == cu.AMINOS and am.shape[1] == 21
    np.testing.assert_allclose(am.sum(axis=1), df[cu.CODONS].sum(axis=1), atol=1e-12)
    np.testing.assert_array_equal(am["methionine"], df["AUG"])


def test_groups(clean):
    df, _ = clean
    assert df.group_species.nunique() == 12331
    assert df.group_genus.nunique() == 4548
    # a genus group always contains whole species groups
    assert (df.groupby("group_species").group_genus.nunique() == 1).all()
    # rows sharing a SpeciesID are different compartments of one organism
    assert not df.duplicated(["SpeciesID", "DNAtype"]).any()


def test_name_token():
    assert cu.name_token("Homo sapiens") == "homo"
    assert cu.name_token("mitochondrion Bos taurus") == "bos"
    assert cu.name_token("x Doritaenopsis sp.") == "doritaenopsis"
    assert cu.name_token("Candidatus Pelagibacter ubique") == "pelagibacter"


ANIMAL = ["vertebrate", "mammal", "rodent", "primate", "invertebrate"]
VERTEBRATE = ["vertebrate", "mammal", "rodent", "primate"]


def test_animal_mitochondria_use_a_different_genetic_code(clean):
    """Not a 2022 defect but a caveat on amino_features: animal mitochondria read
    UGA as tryptophan (not stop), AUA as methionine, and vertebrate mitochondria
    avoid AGA/AGG (stops there). For those rows the standard-code groups are
    codon groups, not amino-acid composition."""
    df, _ = clean
    mito = df[(df.dna_type == "mitochondrial") & df.kingdom.isin(ANIMAL)]
    vert_mito = mito[mito.kingdom.isin(VERTEBRATE)]
    nuclear = df[df.dna_type == "genomic"]
    assert len(mito) == 2747 and len(vert_mito) == 2336
    assert (mito.UGA > mito.UGG).mean() > 0.95      # UGA used as Trp
    assert (nuclear.UGA > nuclear.UGG).mean() < 0.02  # UGA a rare stop
    assert (vert_mito.AGA + vert_mito.AGG).median() < 0.001


def test_raw_file_matches_documented_checksum():
    """data/README.md records this sha256; the raw file must stay unchanged."""
    import hashlib

    digest = hashlib.sha256(cu.DATA.read_bytes()).hexdigest()
    assert digest == "ba8c66381c598ead1c800a66e457e1496c94ecf537db54ca7d89c5a50cd033ca"
