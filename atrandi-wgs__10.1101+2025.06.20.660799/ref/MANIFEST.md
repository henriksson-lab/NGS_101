# Reference material for the Atrandi SPC + PTA page

What this directory holds, what it is for, and how to obtain anything that is **not**
committed. Third-party documents are deliberately absent: they are copyrighted and not
ours to redistribute. Every fact taken from one lives in `tools/atrandi.py` and is pinned
in `tools/selftest.py`, so the derived chemistry is checkable and citable without shipping
the source document.

## Committed here

| file | what it is | depends on it |
|---|---|---|
| `bascet_atrandi_barcodes.tsv` | The 96 round barcodes (24 per round, 8 nt) as Bascet reads them. Verbatim from [henriksson-lab/bascet](https://github.com/henriksson-lab/bascet) `crates/bascet-cli/src/barcode/atrandi_barcodes.tsv` (MIT licensed, so redistributable). | `03_our_protocol.md` §9; the whitelist comparison in `02_library_prep.md` §10 |
| `our_index_primers.tsv` | **Ours.** The i7/P5 index PCR primers as we order them, with the index as it sits in the oligo and as the sample sheet wants it. | `tools/selftest.py` checks `atrandi.OUR_I7`/`OUR_P5` against it base for base; rendered on the page |

## Not committed — obtain from the vendor

Atrandi Biosciences / Droplet Genomics, Research Use Only. Both are free downloads from
the vendor on request; cite by document number **and revision**, since the oligo tables
have changed between revisions.

| document | number / revision | what we took from it |
|---|---|---|
| Single Microbe DNA Barcoding Kit | `DGPM02323198001` **V3** | kit contents, SPC handling, lysis and MDA conditions, barcoding rounds — `01_barcoding_kit.md` |
| Single Microbe DNA Barcoding Kit — Library Prep for Sequencing | `DGPM02323206001` **V3** | the ligation adapter and indexing primer sequences, fragmentation and PCR programs — `02_library_prep.md`, and `LIGADAPT_*` / `ATRANDI_*` in `tools/atrandi.py` |

The adapter and primer strings quoted from `DGPM02323206001` V3 are reproduced once, as a
comment in `tools/atrandi.py`, and the self-test parses that comment and compares it to
the strands the model derives — so the transcription is checked even though the manual
itself is not here.

## Papers

Cited by DOI in the notes, never stored: the Bascet/Zorn preprint
([10.1101/2025.06.20.660799](https://doi.org/10.1101/2025.06.20.660799)) for our protocol,
and the PTA, Atrandi-platform and contaminant references listed in `06_pta.md` and
`03_our_protocol.md` §10.

No check in this directory reads a third-party file, so `tools/selftest.py` runs in full on
a fresh clone.
