# Large-scale detection of antigen-specific T cells using peptide-MHC-I multimers labeled with DNA barcodes

Defining DOI: 10.1038/nbt.3662. Status: capped source pass; incomplete private model, no finished public schematic.

Sources actually inspected:

- https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fnbt.3662/MediaObjects/41587_2016_BFnbt3662_MOESM2_ESM.pdf: `nbt.3662_41587_2016_BFnbt3662_MOESM2_ESM.pdf.txt` — Published SI FigS1 caption2–26 and Table1 lines119–193; rendered original p11 checked for Biotin/primer/6N/25nt/annealing color boundaries.
- https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fnbt.3662/MediaObjects/41587_2016_BFnbt3662_MOESM3_ESM.xlsx: `nbt.3662_41587_2016_BFnbt3662_MOESM3_ESM.xlsx.txt` — Published Table2 workbook header and B1–B24, exact B1 body and named component roles.
- https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fnbt.3662/MediaObjects/41587_2016_BFnbt3662_MOESM4_ESM.xlsx: `nbt.3662_41587_2016_BFnbt3662_MOESM4_ESM.xlsx.txt` — Published Table3 F1 and R1 whole primer strings and all source named columns; F1 rich-text XML inspected, no unsupported fine linker/sample boundary inferred.

🟢 Published SI FigS1 caption states two identical DNA barcodes are attached per multimerization backbone before pMHC co-attachment. This is an intended source assembly, not computed exact occupancy. Caption reports FACS-sorted cell-associated barcode recovery; qPCR validation is not substituted for the core NGS recipe.

🟢 Table1 rendered original p11 identifies biotin, blue primer body, green6N UMI, black25nt barcode and turquoise annealing arm. Table2 provides OligoB1 with its independently named6N UMI and25nt barcode. These are reporter molecules, not cell barcodes. No individual N6 sequence is instantiated.

- Biotin OligoA1: `BiotinGAAGTTCCAGCCAGCGTCACAGTTTNNNNNNCGAGGGCAATGGTTAACTGACACGTGGTCAGCATCATTTCC`
- OligoB1: `CTGTGACTATGTGAGGCTTTCTCGANNNNNNGCCTGTAGTCCCACGCGATCTAACAGGAAATGATGCTGACC`
- Table3 F1: `CCATCTCATCCCTGCGTGTCTCCGACTCAGCTGGGGGAAGTTCCAGCCAGCGTC`
- Table3 R1: `CCTCTCTATGGGCAGTCGGTGATCTGTGACTATGTGAGGCTTTC`

🟡 The two printed16nt annealing arms are reverse complements by computation. This checks that relation only; full extension, strand synthesis and attachment geometry remain unassembled. Table3 identifies sample-ID/A-Key forward and P1-Key reverse primers, so IonTorrent adapter provenance is positive; exact instrument/cycle/run layout is unreviewed. F1 remains a whole printed primer: the short sample/linker sub-boundaries are not inferred from color prose alone.

The private graph exposes source overview interfaces and known reporter molecules, with explicit opaque main-library continuation. Other sorted populations remain unresolved rather than being silently discarded. It does not infer a TetTCR or current commercial recipe.

Actual blocker: standard source manifest reports no open defining full text at the DOI URL. The UCL repository lead is an accepted manuscript, not automatically a published-method substitute. No additional retrieval was pursued beyond the capped pass. The main-source barrier blocks further preparation ordering, not an accepted claim that details are absent.

Component lineage: potential predecessor relationship to TetTCR must be checked in TetTCR's own source/version before a scoped edge; shared pMHC DNA tags alone do not establish recipe identity.

Deferred targets (not claims of absence):

- Defining published main Online Methods: reporter extension/assembly, staining, recovery, all library PCR/cleanup/acquisition order
- Full barcode duplex extension geometry and exact sample-ID/linker sub-boundaries in primer F1
- Other published SI validation figures and library experiment scope; source-linked manufacturing/predecessor documents

Verification: mandatory standard oligo scraper was run over each fetched source bundle; source pointers and relevant modifications were checked against the originals stated above. Local tests check transcription and consequential fraction/order distinctions; graph validation is separate from molecular/source completeness. No absence audit or absence acceptance.

Source transcription guard: a narrowly scoped `oligo_a` exception records the published Table1 page11 line-wrap extraction issue (row label A1 interposed at lines119–126). Two readers independently verified the whole original biotinylated strand in the rendered PDF. This is an extraction limitation, not a publication error; the guard reports the exemption explicitly.

## Bounded publication follow-up (2026-10-30)

The published supplement positively establishes DNA-barcode attachment before pMHC co-attachment, staining and sorted-cell barcode recovery at overview level (FigS1caption2–26), and its tables identify actual reporter and PCR-primer components. These facts do not specify the complete reporter-extension, recovery, PCR, cleanup or acquisition sequence from the defining Online Methods. The required `unassembled` item remains uncovered. Source primer inventory alone does not establish a DNA-library-producing reaction or a sequencing-ready product. The current opaque continuation and generic co-attachment therefore remain private modeling boundaries, not completed chemical assertions.

This bounded pass reread only three existing published originals: MOESM2 FigS1/Table1 and actualp11 A1/color-role image; Table2 header/B1–B4; Table3 header/F1–F4 and source role key. The actual defining main remains unread under the historical standard-manifest limitation; there was no fresh retrieval attempt or current access-failure claim. No fine F1 sample/linker boundary was inferred from its sequence alone. N6 UMI and25nt reporter identity remain distinct source roles, not cell barcodes.

No module changes are justified merely to pass publication gates. A later source-supported implementation should replace runtime graph-enumerated inventory with a literal reviewed target list and preserve molecular assembly/stocks/fates, but must first satisfy the real required main/assembly scope. Notes status remains a source-scope hold; no absence claim or full-source review.
