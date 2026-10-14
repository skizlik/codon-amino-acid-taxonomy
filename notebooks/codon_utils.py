"""Shared, tested helpers for the codon-usage notebooks.

Data: read the raw UCI file and repair its two comma-shifted rows; codebooks;
the explicit class filter; amino-acid features; cross-validation groups.
"""
import csv
from pathlib import Path

import numpy as np
import pandas as pd

DATA = Path(__file__).resolve().parent.parent / "data" / "raw" / "codon_usage.csv"

# Column order of the UCI file.
CODONS = [
    "UUU", "UUC", "UUA", "UUG", "CUU", "CUC", "CUA", "CUG", "AUU", "AUC", "AUA", "AUG",
    "GUU", "GUC", "GUA", "GUG", "GCU", "GCC", "GCA", "GCG", "CCU", "CCC", "CCA", "CCG",
    "UGG", "GGU", "GGC", "GGA", "GGG", "UCU", "UCC", "UCA", "UCG", "AGU", "AGC", "ACU",
    "ACC", "ACA", "ACG", "UAU", "UAC", "CAA", "CAG", "AAU", "AAC", "UGU", "UGC", "CAU",
    "CAC", "AAA", "AAG", "CGU", "CGC", "CGA", "CGG", "AGA", "AGG", "GAU", "GAC", "GAA",
    "GAG", "UAA", "UAG", "UGA",
]

# Synonymous-codon groups under the STANDARD genetic code (NCBI table 1).
# AUG appears once (methionine): the 2022 'start' column duplicated it exactly,
# because aggregate frequencies cannot separate initiator from internal AUGs.
# Animal mitochondria (about 21% of rows) use other codes (UGA = Trp, AUA = Met,
# AGA/AGG = stop or Ser), so for those rows these are codon groups, not
# amino-acid composition. Choosing the right code per row would require the
# compartment and lineage -- the targets -- so the standard code is the only
# label-free choice.
AMINO_CODONS = {
    "alanine": ["GCU", "GCC", "GCA", "GCG"],
    "arginine": ["CGU", "CGC", "CGA", "CGG", "AGA", "AGG"],
    "asparagine": ["AAU", "AAC"],
    "aspartic_acid": ["GAU", "GAC"],
    "cysteine": ["UGU", "UGC"],
    "glutamine": ["CAA", "CAG"],
    "glutamic_acid": ["GAA", "GAG"],
    "glycine": ["GGU", "GGC", "GGA", "GGG"],
    "histidine": ["CAU", "CAC"],
    "isoleucine": ["AUU", "AUC", "AUA"],
    "leucine": ["CUU", "CUC", "CUA", "CUG", "UUA", "UUG"],
    "lysine": ["AAA", "AAG"],
    "methionine": ["AUG"],
    "phenylalanine": ["UUU", "UUC"],
    "proline": ["CCU", "CCC", "CCA", "CCG"],
    "serine": ["UCU", "UCC", "UCA", "UCG", "AGU", "AGC"],
    "threonine": ["ACU", "ACC", "ACA", "ACG"],
    "tryptophan": ["UGG"],
    "tyrosine": ["UAU", "UAC"],
    "valine": ["GUU", "GUC", "GUA", "GUG"],
    "stop": ["UAA", "UAG", "UGA"],
}
AMINOS = list(AMINO_CODONS)
assert sorted(c for cs in AMINO_CODONS.values() for c in cs) == sorted(CODONS)

# UCI codebook (archive.ics.uci.edu/dataset/577/codon+usage).
KINGDOM_CODES = {
    "arc": "archaea", "bct": "bacteria", "phg": "bacteriophage", "plm": "plasmid",
    "pln": "plant", "inv": "invertebrate", "vrt": "vertebrate", "mam": "mammal",
    "rod": "rodent", "pri": "primate", "vrl": "virus",
}
DNATYPE_CODES = {
    0: "genomic", 1: "mitochondrial", 2: "chloroplast", 3: "cyanelle", 4: "plastid",
    5: "nucleomorph", 6: "secondary_endosymbiont", 7: "chromoplast", 8: "leucoplast",
    9: "unknown",  # "NA" in the codebook: keep it as a code, never as a missing value
    10: "proplastid", 11: "apicoplast", 12: "kinetoplast",
}


def _is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False


def repair_row(fields):
    """Realign a row whose SpeciesName contains unquoted commas.

    The name spills into the codon columns, so the last codon value(s) are
    absent from the file. Frequencies are counts / Ncodons (5 dp), so the
    missing counts are Ncodons minus the counts present. Returns the repaired
    fields and the number of codons recovered; raises ValueError if the
    missing values cannot be identified.
    """
    k = 5
    while not all(_is_number(x) for x in fields[k:]):
        k += 1
    values = np.array(fields[k:], dtype=float)
    n_missing = len(fields) - 5 - len(values)
    ncodons = int(fields[3])
    counts = values * ncodons
    if np.abs(counts - counts.round()).max() > 0.05:
        raise ValueError(f"row {fields[2]}: frequencies are not counts / Ncodons")
    residual = ncodons - int(counts.round().sum())
    if n_missing == 1:
        recovered = [residual / ncodons]
    elif residual == 0:
        recovered = [0.0] * n_missing  # non-negative counts that sum to 0
    else:
        raise ValueError(f"row {fields[2]}: {n_missing} missing codons, residual {residual}")
    name = ",".join(fields[4:k])
    return fields[:4] + [name] + fields[k:] + [f"{v:.5f}" for v in recovered], n_missing


def read_raw(path=DATA):
    """Read codon_usage.csv, repairing comma-shifted rows.

    Adds `repaired` (bool) and the decoded labels `kingdom` and `dna_type`.
    DNAtype stays an integer code (code 9 is not turned into a missing value).
    """
    with open(path, newline="") as fh:
        header, *body = list(csv.reader(fh))
    rows, repaired = [], []
    for fields in body:
        if len(fields) != len(header):
            raise ValueError(f"{len(fields)} fields, expected {len(header)}")
        is_bad = not all(_is_number(x) for x in fields[5:])
        rows.append(repair_row(fields)[0] if is_bad else fields)
        repaired.append(is_bad)
    df = pd.DataFrame(rows, columns=header)
    df[CODONS] = df[CODONS].astype(float)
    df[["DNAtype", "SpeciesID", "Ncodons"]] = df[["DNAtype", "SpeciesID", "Ncodons"]].astype(int)
    df["repaired"] = repaired
    df["kingdom"] = df.Kingdom.map(KINGDOM_CODES)
    df["dna_type"] = df.DNAtype.map(DNATYPE_CODES)
    if df[["kingdom", "dna_type"]].isna().any().any():
        raise ValueError("a Kingdom or DNAtype code is not in the codebook")
    return df


def filter_classes(df):
    """The project's class filter, step by step. Returns (rows kept, row-count log)."""
    steps = [("raw", len(df))]
    df = df[df.kingdom != "plasmid"]
    steps.append(("drop Kingdom 'plm' (plasmid)", len(df)))
    df = df[df.dna_type.isin(["genomic", "mitochondrial", "chloroplast"])]
    steps.append(("keep DNAtype in {0 genomic, 1 mitochondrial, 2 chloroplast}", len(df)))
    # 'chloroplast Toxoplasma gondii', filed under inv: an apicoplast, not a chloroplast
    df = df[~((df.kingdom == "invertebrate") & (df.dna_type == "chloroplast"))]
    steps.append(("drop inv + chloroplast (Toxoplasma apicoplast) row", len(df)))
    log = pd.DataFrame(steps, columns=["step", "rows"])
    log["removed"] = (-log.rows.diff()).fillna(0).astype(int)
    return df.copy(), log


def amino_features(df):
    """Sum codon frequencies within each amino acid (and stop): 21 columns."""
    return pd.DataFrame({a: df[cs].sum(axis=1) for a, cs in AMINO_CODONS.items()}, index=df.index)


_SKIP = {"mitochondrion", "chloroplast", "x", "candidatus", "uncultured", "unidentified", "[candida]"}


def name_token(name):
    """First informative word of SpeciesName: the genus for most cellular
    organisms, often a host word ('human', 'tomato') for viruses."""
    words = [w for w in str(name).replace("[", " [").split() if w.lower() not in _SKIP]
    return words[0].lower() if words else str(name).lower()


def add_groups(df):
    """group_species = SpeciesID (the NCBI taxid); group_genus = name token,
    made constant within SpeciesID so that genus groups contain whole species."""
    df = df.copy()
    df["group_species"] = df.SpeciesID
    df["group_genus"] = df.SpeciesName.map(name_token).groupby(df.SpeciesID).transform("first")
    return df


def load_project(path=DATA):
    """read_raw -> filter_classes -> add_groups. Returns (modelling frame, row-count log)."""
    df, log = filter_classes(read_raw(path))
    return add_groups(df), log
