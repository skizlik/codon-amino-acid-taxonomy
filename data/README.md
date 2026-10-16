# Data

## Source

`raw/codon_usage.csv` is the **Codon Usage** dataset from the UCI Machine Learning Repository
(dataset 577, <https://archive.ics.uci.edu/dataset/577/codon+usage>), donated by B. B. Khomtchouk in
connection with Khomtchouk (2020), *Codon usage bias levels predict taxonomic identity and genetic
composition*, bioRxiv 10.1101/2020.10.26.356295. The rows are derived from the Codon Usage
Tabulated from GenBank (CUTG) database (Nakamura, Gojobori & Ikemura, *Nucleic Acids Research*
28:292, 2000).

The file is committed unchanged.

| | |
|---|---|
| size | 7,001,200 bytes; 1 header line + 13,028 data rows; 69 fields per line |
| sha256 | `ba8c66381c598ead1c800a66e457e1496c94ecf537db54ca7d89c5a50cd033ca` (checked by `tests/test_data.py`) |
| columns | `Kingdom`, `DNAtype`, `SpeciesID`, `Ncodons`, `SpeciesName`, then 64 codon frequencies (`UUU` … `UGA`) |

- **Citation (VERIFY before publishing):** cite the dataset as the UCI page instructs, including its
  DOI, which is not recorded here.
- **License (VERIFY before publishing):** believed to be CC BY 4.0, as for most UCI datasets; confirm
  on the dataset page. The code in this repository is MIT-licensed; the data keep their own license.

`SpeciesID` is the NCBI Taxonomy identifier (for example 9606 *Homo sapiens*, 10090 *Mus musculus*,
3702 *Arabidopsis thaliana*, 7227 *Drosophila melanogaster*). Codon frequencies are codon counts
divided by `Ncodons`, rounded to 5 decimal places.

## Codebook

`Kingdom` is the CUTG (GenBank) division, not a biological kingdom. Its classes are nested
(primate ⊂ mammal ⊂ vertebrate; each entry is filed in the most specific division) and partly
administrative (the plant division also holds fungi and algae).

| code | label | raw rows | modelled rows |
|---|---|---:|---:|
| `arc` | archaea | 126 | 126 |
| `bct` | bacteria | 2,920 | 2,918 |
| `phg` | bacteriophage | 220 | 220 |
| `plm` | plasmid | 18 | 0 |
| `pln` | plant | 2,523 | 2,490 |
| `inv` | invertebrate | 1,345 | 1,333 |
| `vrt` | vertebrate (other than mammals) | 2,077 | 2,077 |
| `mam` | mammal (other than rodents and primates) | 572 | 572 |
| `rod` | rodent | 215 | 215 |
| `pri` | primate | 180 | 180 |
| `vrl` | virus | 2,832 | 2,832 |

| `DNAtype` | label | raw rows | modelled rows |
|---:|---|---:|---:|
| 0 | genomic | 9,267 | 9,249 |
| 1 | mitochondrial | 2,899 | 2,899 |
| 2 | chloroplast | 816 | 815 |
| 3 | cyanelle | 2 | 0 |
| 4 | plastid | 31 | 0 |
| 5 | nucleomorph | 2 | 0 |
| 6 | secondary endosymbiont | 1 | 0 |
| 7 | chromoplast | 1 | 0 |
| 9 | "NA" in the codebook (decoded as `unknown`, never parsed as missing) | 2 | 0 |
| 11 | apicoplast | 2 | 0 |
| 12 | kinetoplast | 5 | 0 |

Codes 8 (leucoplast) and 10 (proplastid) are in the codebook but do not occur.

## Known defects and repairs

**Two rows with commas in `SpeciesName`.** The names are not quoted, so they spill into the codon
columns. Every line still has 69 fields, which means the *last* codon values are absent from the file
rather than displaced:

| data row | `SpeciesName` | name fields | codons missing | recovered from `Ncodons` |
|---:|---|---:|---:|---|
| 486 | Non-A,non-B hepatitis virus | 2 | 1 (`UGA`) | 3 of 1,238 codons unaccounted for → `UGA` = 0.00242 |
| 5063 | Salmonella enterica subsp. enterica serovar 4,12;I,- | 3 | 2 (`UAG`, `UGA`) | 0 codons unaccounted for → both 0 |

`read_raw` in `notebooks/codon_utils.py` joins the name fields, left-aligns the values (`UUU` first)
and recovers the missing trailing codons from `Ncodons`. Left alignment is the only reading with
plausible stop-codon usage: right-aligned, row 486 would have `UAA` ≈ 2.3% of all codons, and row 5063
`UAA` ≈ 4.5% and `UAG` ≈ 2.5%; left-aligned, stop codons are 0.40% and 0.12% of their codons. After
the repair every row's frequencies sum to 0.99984–1.00017. The 2022 pipeline coerced the stray text
to 0, which left both rows in the modelling set with every codon shifted by one or two positions.
Dropping the two rows would also be defensible and changes nothing material; repair is preferred
because it is exact and tested. The derivation is in `notebooks/01_wrangling.ipynb`.

## Genetic code

The 21 amino-acid features (`amino_features` in `notebooks/codon_utils.py`) sum codon frequencies
by the standard genetic code (NCBI translation table 1). That is right for nuclear genomes,
chloroplasts and most plant mitochondria, but not for animal mitochondria (vertebrate: table 2;
invertebrate: table 5), in which UGA encodes tryptophan instead of stop and AUA encodes methionine
instead of isoleucine; AGA and AGG are stops in vertebrate mitochondria and serine in invertebrate
mitochondria. After the class filter, 2,747 of the 12,963 rows (21%) are animal mitochondria, and UGA
outnumbers UGG in 96.5% of them against 0.9% of nuclear rows. For those rows the features are groups
of synonymous codons under the standard code, not amino-acid composition, and the "stop" group
mostly counts tryptophan.

The features are not corrected row by row because choosing a row's genetic code requires its
compartment and lineage, which are the prediction targets. `notebooks/02_eda.ipynb` shows the code
signatures by compartment; `notebooks/04_modeling.ipynb` checks that the DNA-type models do not
depend on the "stop" group.

The 2022 `start` feature was the AUG frequency again. AUG is both the usual initiation codon and the
methionine codon, and aggregate codon counts cannot separate the first AUG of each gene from the
others, so `start` duplicated `methionine` exactly; it is removed.

## Rows used for modelling

The class filter follows the 2022 project and is applied explicitly by `filter_classes` in
`notebooks/codon_utils.py` (also written to `results/row_counts.csv`):

| step | rows | removed |
|---|---:|---:|
| raw | 13,028 | 0 |
| drop `Kingdom` = `plm` (plasmid) | 13,010 | 18 |
| keep `DNAtype` ∈ {0 genomic, 1 mitochondrial, 2 chloroplast} | 12,964 | 46 |
| drop the `inv` + chloroplast row (*Toxoplasma gondii*, an apicoplast) | 12,963 | 1 |

The 46 rows of the second step include the two `DNAtype` = 9 rows (*x Doritaenopsis* sp. and
*x Tritordeum* sp.), which the 2022 notebooks lost silently as missing values.

## Relatedness and grouping

- 12,331 distinct `SpeciesID`s; 596 of them have 2–3 rows. Those rows are always different
  compartments of one organism (genomic, mitochondrial, chloroplast): no `SpeciesID` has two rows with
  the same `DNAtype`, and their codon vectors are no closer than random pairs of rows
  (`notebooks/03_preprocessing.ipynb`).
- Strains and close relatives often have *different* `SpeciesID`s and near-identical codon vectors
  (23 groups of identical vectors, all spanning several `SpeciesID`s). Cross-validation therefore holds
  out whole groups: by `SpeciesID` (*species*) or by a genus-level group (*genus*: the first
  informative word of `SpeciesName`, made constant within `SpeciesID`; 4,548 groups). For viruses the
  first word is often a host name (for example "Human"), so the genus grouping is a heuristic. Because
  `SpeciesID` is the NCBI taxid, real NCBI ranks could replace it.
