# Safe-SeqS — research notes

**Defining paper:** Kinde *et al.*, “Detection and quantification of rare mutations with massively parallel sequencing,” *PNAS* 108:9530–9535 (2011), DOI [10.1073/pnas.1105422108](https://doi.org/10.1073/pnas.1105422108).

## Sources read

- 🟢 Defining article (PMC3111315 full text): Fig. 3 exogenous-UID strategy (extracted text line 1283) and Materials and Methods, Endogenous/Exogenous UIDs (lines 1652–1655).
- 🟢 SI Materials and Methods and SI Table S4, read in the original ten-page PNAS SI facsimile reproduced in a US PTAB exhibit (local `download-documents.pdf`, fetched 2026-10-17; SI text at `.txt` lines 894–1035; Table S4 on exhibit pages 0015–0016). The publisher and PMC SI copies remain gated.
- 🟡 Table S4 is an image: text extraction captures only its caption and legend. The CTNNB1 block (page 10 of 10) and the polymerase-fidelity block (page 9) were read visually from the rendered page at 400 dpi on 2026-10-17; the CTNNB1 N-run was measured as 14 characters. This is one reading; a second independent visual reading is still needed before the transcription counts as checked, and the repository's text-based oligo guard cannot verify it.
- 🔴 Not inspected: the cited Herman *et al.* capture predecessor; the other Table S4 blocks (endogenous capture, inverse PCR, mitochondrial, oligonucleotide composition) beyond orientation.

## Oligos (SI Table S4, CTNNB1 mutations in DNA from normal human cells), 5′→3′

- 🟢 UID assignment amplification, for: `CGACGTAAAACGACGGCCAGT NNNNNNNNNNNNNN GCAGCAACAGTCTTACCTGGACT` (black universal tag · green UID · blue template-complementary arm).
- 🟢 UID assignment amplification, rev: `CACACAGGAAACAGCTATGACCATG TCCACATCCTCTTCCTCAGGATT`.
- 🟢 Library amplification, for: `AATGATACGGCGACCACCGAGATCTACAC CGACGTAAAACGACGGCCA*G*T` (red grafting sequence; `*` phosphorothioate).
- 🟢 Library amplification, rev-1 … rev-10: `CAAGCAGAAGACGGCATACGAGAT [index] CACACAGGAAACAGCTATGACCA*T*G` with indices ATCACG, CGATGT, TGACCA, GCCAAT, CAGATC, ACTTGA, GATCAG, TAGCTT, GGCTAC, CTTGTA.
- 🟢 Sequencing primer (UID and internal sequences): `CGACGTAAAACGACGGCCAGT`; index primer: `CATGGTCATAGCTGTTTCCTGTGTG`.
- 🟡 Derived in the module: the sequencing primer equals the forward tag and the index primer is the reverse complement of the reverse tag; both relations are checked at import.

## Construct-changing steps (SI Materials and Methods, Exogenous UIDs)

1. 🟢 DNA quantified by digital PCR. UID assignment: Phusion Hot Start II, 45 µL, 1× HF buffer, 0.25 mM dNTPs, 0.5 µM each primer, 2 U polymerase; 98 °C 30 s, then two cycles of 98 °C 10 s, 61 °C 120 s, 72 °C 10 s. Final template <1.5 ng/µL (multiple wells if needed).
2. 🟢 Exonuclease I (Enzymatics), 60 U per well, 37 °C 1 h; 98 °C 5 min heat inactivation; no purification.
3. 🟢 Library primers (0.5 µM each; two terminal phosphorothioates; 5′ grafting sequences; index between graft and universal tag) added directly; 25 cycles of 98 °C 10 s, 72 °C 15 s; no intermediate purification. 🟡 The model assumes the carried-over Phusion performs this PCR (no fresh polymerase is mentioned).
4. 🟢 Wells consolidated, QIAquick purified (50 µL EB), expected-size band excised from a 10% TBE polyacrylamide gel (2% agarose Size Select gel for mtDNA libraries), crushed and eluted.
5. 🟢 Illumina GA IIx; read lengths 36–73 bases across experiments (not assigned per assay). Reads kept if the first 25 bases passed chastity, every base Q≥20, ≤3 mismatches to expected sequence, and UID bases Q≥30.
6. 🟢 Families with ≥2 members are analysed; a “supermutant” requires ≥95% of the family to share the mutation.

## Final library and reads

P5 graft – forward tag – 14-nt UID – CTNNB1 forward arm – amplicon interior – CTNNB1 reverse arm (rc) – reverse tag (rc) – 6-nt index (rc) – P7 graft (rc). The custom sequencing primer reads the UID first, then the target; the index primer reads the index. No second genomic read is described.

## Changes in the 2026-10-17 bounded pass

The former representative model used an INFERRED TruSeq-compatible shell. It was replaced by the Table S4 CTNNB1 primers, which install the P5/P7 grafting sequences directly on non-TruSeq universal tags, with custom sequencing and index primers. Reaction conditions, the polyacrylamide size selection and its waste fate were added to the graph.

## Open questions

- Second independent visual reading of Table S4 (bases, N-run length, phosphorothioate positions).
- CTNNB1 amplicon interior and length; acquired read length for this assay.
- Whether fresh polymerase was added for the 25-cycle PCR.
- The endogenous-UID capture/inverse-PCR libraries and the other exogenous assays (index on the P5 side in the polymerase-fidelity and mtDNA blocks) are not drawn.
