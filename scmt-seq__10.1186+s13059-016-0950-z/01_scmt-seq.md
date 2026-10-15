# scMT-seq — source intake

Defining paper: Hu et al., *Genome Biology* (2016), DOI
[10.1186/s13059-016-0950-z](https://doi.org/10.1186/s13059-016-0950-z).

🟢 The physically separated fractions enter transcriptome and bisulfite-methylome
branches. scMT-seq remains distinct from scM&T-seq, a separate named protocol.

## Sources and inherited chemistries

- 🟢 Hu et al. explicitly state that cytoplasmic cDNA is amplified by Smart-seq2 (with SuperScript III substituted for SuperScript II), followed by Nextera tagmentation and library construction.
- 🟢 The nuclear fraction follows the previously published scRRBS workflow: MspI digestion, end repair/dA-tailing, indexed premethylated adapter ligation, bisulfite conversion and two PCR rounds.
- 🟡 Smart-seq2 oligo-dT, TSO and ISPCR bases are taken from the Picelli supplementary protocol through the shared `smartseq2` model. The scRRBS adapter and final-library structures are taken through the shared `rrbs` model.

The defining paper does not introduce a new oligo set for either branch. The earlier page's generic post-bisulfite random-primer branch was therefore incorrect for scMT-seq and has been replaced by the inherited scRRBS chemistry.

🔴 The inherited scRRBS paper names standard premethylated indexed Illumina adapters but
does not print its library-PCR primer oligos. The final arms are assembled from the
Illumina adapter guide; they do not establish exact primer landing on an unmodelled
bisulfite-converted pre-PCR intermediate. The scMT-seq DNA PCR geometry stays unavailable.

## The TSO linker is inherited, not published here — 🟡

🟢 Hu et al. 2016 print **no oligo sequence** for the RNA branch. The whole of
"Single-cell RNA-seq library construction" reads: *"Single-cell complementary DNA was
amplified from the tubes containing cytosol according to the Smart2-seq protocol. Instead
of using Superscript II, we used Superscript III for reverse transcription. After
amplification and purification, 0.1 ng cDNA was used for Nextera Tagmentation and library
construction."* (`s13059-016-0950-z_PMC4858893.pdf.txt`, around line 745.) The reference
list cites both Smart-seq2 papers — ref. 10 Picelli et al., *Nat. Methods* 10:1096–8
(2013) and ref. 18 Picelli et al., *Nat. Protoc.* 9:171–81 (2014) — but neither the main
text nor Additional file 1 (`13059_2016_950_MOESM1_ESM.pdf`) contains an `AAGCAGTGG…`
string.

🟡 That makes the TSO an **inherited inference**, and the inheritance is genuinely
ambiguous, because the two cited papers disagree about the linker:

| Linker | Where | TSO as printed |
|---|---|---|
| `AC` | ref. 10, Picelli et al. 2013 Supplementary Table 2, row `rGrG+G` — the oligo Supplementary Table 4 names as Smart-seq2's own | `AAGCAGTGGTATCAACGCAGAGTAC` `rGrG+G` |
| `ACAT` | ref. 18, the *Nature Protocols* form | `AAGCAGTGGTATCAACGCAGAGTACAT` `rGrG+G` |

In the copy of Supplementary Table 2 fetched into this bundle
(`nmeth.2639_…_MOESM301_ESM.xlsx`, lines 4–16) three of thirteen oligos (`C6 Amino`,
`rGrGrG`, `rG3p`) carry `ACAT` and the other ten, `rGrG+G` among them, carry `AC`.

🟡 **The model assumes the `ACAT` form** — `smart-seq2__10.1038+nmeth.2639/tools/smartseq2.TSO_LINKER`, the constant four
SMART-seq2 derivatives share — because that is the protocol as practised. Citing both
papers does not choose between them, so the schematic now draws the **whole linker**
dotted, with all three required signals: a preamble caveat naming both forms, an
`INFERRED — …` caption on the template-switch step, and `<inf>` around `ACAT` in the oligo
list and the drawn cDNA states. `tools/chem.py:_tso_linker_inferred` applies the mark to
the state itself and raises if the shared segment is renamed, so the uncertainty cannot be
lost by a caller.

🔴 What would settle it: the 2014 *Nature Protocols* Smart-seq2 paper is not in this
bundle, and Hu et al. never say which version of "Smart2-seq" they ran. Only a statement
by these authors would make the linker green for scMT-seq.
