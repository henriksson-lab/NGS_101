# scNMT-seq

## 1. What it is

scNMT-seq jointly measures chromatin accessibility, endogenous DNA methylation and RNA expression by GpC marking followed by physical RNA/DNA separation and branch-specific libraries.

## 2. Authoritative source

🟢 Clark et al., *scNMT-seq enables joint profiling of chromatin accessibility, DNA methylation and transcription in single cells*, Nature Communications (2018), DOI `10.1038/s41467-018-03149-4`.

## 3. Protocol boundary

🟢 The RNA branch uses bead-capture Smart-seq2; the DNA branch uses M.CviPI marking followed by scBS-seq. Both branches belong on this protocol page because their pairing is the defining operation.

## 4. Construct path

🟢 Cells are treated with M.CviPI for 15 minutes before lysis. Oligo-dT magnetic beads capture RNA; wash supernatants are transferred to the DNA plate.

🟢 The RNA branch is Smart-seq2 cDNA followed by one-fifth-volume Nextera XT. The DNA branch is bisulfite conversion followed by the published scBS-seq two-round random-priming chemistry and PE1.0/iPCRTag amplification.

🟡 The defining paper names its RNA oligos as the Smart-seq2 oligo-dT, TSO and ISPCR
primer and cites Picelli et al. (2013, 2014). Their sequences and the canonical
Nextera s5/s7 indexing-primer geometry are therefore reused from those authoritative
chemistries rather than left as protocol-local unknowns. The scNMT-seq paper does not
claim a different sequence variant — but it does not print one either, which matters for
the TSO linker (below).

🟢 Both branches were sequenced paired-end. The schematic keeps their distinct Nextera and historical scBS-seq primer layouts separate.

## 5. The TSO linker is inherited, not published here — 🟡

🟢 Clark et al. 2018 print **no oligo sequence**. The Library-preparation paragraph names
each reagent by citation: *"mRNA was captured using Smart-seq2¹⁰,²⁷ oligo-dT pre-annealed
to magnetic beads (MyOne C1, Invitrogen) … 1 μM Template-Switching Oligo¹⁰,²⁷ (Exiqon) …
1 μl of 2 μM ISPCR primer¹⁰,²⁷"* (`s41467-018-03149-4_PMC5823944.pdf.txt`, lines
1092–1104). Refs. 10 and 27 are Picelli et al. 2013 (*Nat. Methods*) and Picelli et al.
2014 (*Nat. Protoc.*). Neither the main text, the three supplementary PDFs
(`41467_2018_3149_MOESM1..3`) nor the supplementary tables contain an `AAGCAGTGG…` string.

🟡 So the TSO is an **inherited inference**, and because the two cited papers are cited
*together* the inheritance is ambiguous: they print different linkers.

| Linker | Where | TSO as printed |
|---|---|---|
| `AC` | ref. 10, Picelli et al. 2013 Supplementary Table 2, row `rGrG+G` — the oligo Supplementary Table 4 names as Smart-seq2's own | `AAGCAGTGGTATCAACGCAGAGTAC` `rGrG+G` |
| `ACAT` | ref. 27, the *Nature Protocols* form | `AAGCAGTGGTATCAACGCAGAGTACAT` `rGrG+G` |

In the copy of Supplementary Table 2 fetched into this bundle
(`nmeth.2639_…_MOESM301_ESM.xlsx`, lines 4–16) three of thirteen oligos (`C6 Amino`,
`rGrGrG`, `rG3p`) carry `ACAT` and the other ten, `rGrG+G` among them, carry `AC`.

🟡 **The model uses the `ACAT` form** — `smart-seq2__10.1038+nmeth.2639/tools/smartseq2.TSO_LINKER`, the constant four
SMART-seq2 derivatives share — because Clark cites the 2014 detailed protocol that prints it. This remains an inferred selection. "Exiqon" names
the supplier of the LNA `+G`, not a linker. The schematic therefore draws the **whole
linker** dotted, with all three required signals: a preamble caveat naming both forms, an
`INFERRED — …` caption on the template-switch step, and `<inf>` around `ACAT` in the oligo
list and the drawn cDNA states. The local graph marks its TSO input linker inferred;
`smartseq2.preamp_template()` marks the same linker inferred in the inherited cDNA state.

🟢 An author-uploaded copy of Picelli et al. 2014 *Nature Protocols* is readable
[online](https://www.researchgate.net/publication/259566573_Full-length_RNA-seq_from_single_cells_using_Smart-seq2):
its Reagent Setup (web lines 763–789) prints the `ACATrGrG+G` TSO, oligo-dT30VN
and ISPCR. It confirms that inherited option, but Clark et al. still cite the
2013 and 2014 protocols together without identifying which TSO they ordered.
The exact scNMT-seq TSO remains inferred.

🟡 The DNA branch inherits the published scBS-seq iPCRTag primer family. Clark et al.
name an iTAG index primer but do not identify the selected sequence member in the
defining Methods or the inspected supplements. The final-library diagram therefore
shows a variable sample-index window, not an experiment-specific iTAG sequence.
