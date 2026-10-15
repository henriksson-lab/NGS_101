# scM&T-seq — research notes

## Sources read

- 🟢 Angermueller et al. 2016, DOI `10.1038/nmeth.3728`, PMC full text and supplementary protocol.
- 🟢 SMART-seq2 and scBS-seq defining papers and supplements for the two downstream branches.

## Construct path

🟢 Biotinylated oligo(dT) beads capture cellular mRNA. The bead-bound fraction is processed by Smart-seq2 and Nextera XT; the supernatant and washes retain genomic DNA for scBS-seq. 🟢 The common cell identity is plate provenance, not a covalent shared barcode, so the schematic does not invent one.

🟡 The schematic reuses the oligos from the branch-defining sources scM&T-seq delegates
to: Smart-seq2 oligo-dT30VN, TSO and ISPCR; canonical Nextera N5xx/N7xx landing
regions; and scBS-seq oligo1/oligo2 handles with PE1/indexed iPCRTag. The scM&T paper's
Methods explicitly change both bisulfite strand-synthesis primers to **3′ random
hexamers**. The original scBS-seq oligo1/oligo2 use N9; those N9 regions must not be
inherited into scM&T. The DNA branch now models N6 on each primer, while retaining the
source-printed scBS-seq handles.

## The TSO linker is inherited, not published here — 🟡

🟢 Angermueller et al. 2016 print **no oligo sequence at all** for the RNA branch. The
Online Methods say only that cells were "processed using the G&T-seq protocol⁷, but
following physical separation of mRNA and genomic DNA from single cells, the DNA was
eluted into 10 μl of H₂O" (ref. 7 = Macaulay et al., *Nature Methods* 12:519–522, 2015).
Searching the whole fetched bundle — `nmeth.3728_PMC4770512.pdf/.xml`, the author-
manuscript supplements `NIHMS66375-supplement-1..6` and the Springer ESM files
`MOESM593…MOESM614` — returns no `AAGCAGTGG…` hit and no mention of a template-switch
oligo. The reference list has no Picelli entry: the Smart-seq2 inheritance is two hops
away, through G&T-seq.

🟡 So the TSO is an **inherited inference**. It would be one even at one hop, because
SMART-seq2 has two published linkers between the 23-nt handle and the G tail:

| Linker | Where | TSO as printed |
|---|---|---|
| `AC` | Picelli et al. 2013 (`10.1038/nmeth.2639`) Supplementary Table 2, row `rGrG+G` — the oligo Supplementary Table 4 names as Smart-seq2's own | `AAGCAGTGGTATCAACGCAGAGTAC` `rGrG+G` |
| `ACAT` | the *Nature Protocols* form of the protocol (Picelli et al. 2014, `10.1038/nprot.2014.006`) | `AAGCAGTGGTATCAACGCAGAGTACAT` `rGrG+G` |

Supplementary Table 2 is internally consistent about this: in the copy fetched into this
bundle (`nmeth.2639_…_MOESM301_ESM.xlsx`, lines 4–16) three of thirteen oligos (`C6 Amino`,
`rGrGrG`, `rG3p`) carry `ACAT` and the other ten, `rGrG+G` among them, carry `AC`.

🟡 **The model assumes the `ACAT` form** — `lib/smartseq2.TSO_LINKER` — because that is the
protocol as practised and the constant four SMART-seq2 derivatives share. Nothing in
scM&T-seq chooses between the two, so the schematic now draws the **whole linker** dotted,
with all three required signals: a preamble caveat naming both forms, an `INFERRED — …`
caption on the template-switch step, and `<inf>` around `ACAT` in the oligo list and the
drawn cDNA states. `tools/chem.py:_tso_linker_inferred` applies the mark to the state
itself and raises if the shared segment is ever renamed, so the uncertainty cannot be lost
by a caller.

🔴 What would settle it: the 2014 *Nature Protocols* Smart-seq2 paper and the G&T-seq
paper (Macaulay 2015) are not in this bundle. If G&T-seq prints a TSO, that is the form
scM&T-seq actually inherited, and this page could go green.
