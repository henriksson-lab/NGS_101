# scM&T-seq — research notes

## Sources read

- 🟢 Angermueller et al. 2016, DOI `10.1038/nmeth.3728`, PMC full text and Supplementary Fig. 3e caption (Springer `MOESM608` and PMC `NIHMS66375-supplement-1`).
- 🟢 SMART-seq2 and Smallwood scBS-seq (`10.1038/nmeth.3035`, PMC4117646) defining articles; the Smallwood binary Word supplement remains unread.
- 🟢 Macaulay G&T-seq (`10.1038/nmeth.3370`) is the cited RNA-separation predecessor. Its original publisher PDF, uploaded by an author, was read at Online Methods p. 4–5. The separate S1 supplement remains unread.

## Construct path

🟢 Biotinylated oligo(dT) beads capture cellular mRNA. The bead-bound fraction is processed by Smart-seq2 and Nextera XT; the supernatant and washes retain genomic DNA for scBS-seq. 🟢 The common cell identity is plate provenance, not a covalent shared barcode, so the schematic does not invent one.

🟡 The schematic reuses the oligos from the branch-defining sources scM&T-seq delegates
to: Smart-seq2 oligo-dT30VN, TSO and ISPCR; canonical Nextera N5xx/N7xx landing
regions; and scBS-seq oligo1/oligo2 handles with PE1/indexed iPCRTag. The scM&T paper's
Methods explicitly change both bisulfite strand-synthesis primers to **3′ random
hexamers**. The original scBS-seq oligo1/oligo2 use N9; those N9 regions must not be
inherited into scM&T. The DNA branch now models N6 on each primer, while retaining the
source-printed scBS-seq handles.

🟢 Angermueller Supplementary Fig. 3e says many **sequenced DNA fragments** are
concatemers of the primer used in first-strand synthesis, reducing genomic alignment.
The graph includes this observed off-product at the DNA PCR output, with an explicit
discard reason referring to genomic alignment. 🟡 Its placement at PCR is a bookkeeping
point for the final sequenced product: the source does not establish the exact stage or
molecular junction at which the concatemers formed, so no bases are drawn for them.

## The TSO linker is inherited from G&T-seq — 🟡

🟢 Angermueller et al. 2016 print **no oligo sequence at all** for the RNA branch. The
Online Methods say only that cells were "processed using the G&T-seq protocol⁷, but
following physical separation of mRNA and genomic DNA from single cells, the DNA was
eluted into 10 μl of H₂O" (ref. 7 = Macaulay et al., *Nature Methods* 12:519–522, 2015).
The original G&T-seq publisher PDF by Macaulay et al. 2015, Online Methods p. 4,
prints the RNA-capture oligo as
`5′-biotin-triethyleneglycol-AAGCAGTGGTATCAACGCAGAGTACT30VN-3′` (lines 569–576),
the TSO as `5′-AAGCAGTGGTATCAACGCAGAGTACrGrG+G-3′` (lines 591–600), and
the cDNA PCR primer as `5′-AAGCAGTGGTATCAACGCAGAGT-3′` (lines 604–611).
Its four bead washes are pooled with the genomic-DNA supernatant (lines 580–590).
The scM&T paper delegates this branch to G&T-seq and does not repeat its sequences.

🟢 The G&T predecessor establishes the **`AC`** linker for its own TSO. The
Smart-seq2 source family has two published linker forms:

| Linker | Where | TSO as printed |
|---|---|---|
| `AC` | Picelli et al. 2013 (`10.1038/nmeth.2639`) Supplementary Table 2, row `rGrG+G` — the oligo Supplementary Table 4 names as Smart-seq2's own | `AAGCAGTGGTATCAACGCAGAGTAC` `rGrG+G` |
| `ACAT` | the *Nature Protocols* form of the protocol (Picelli et al. 2014, `10.1038/nprot.2014.006`) | `AAGCAGTGGTATCAACGCAGAGTACAT` `rGrG+G` |

Supplementary Table 2 is internally consistent about this: in the copy fetched into this
bundle (`nmeth.2639_…_MOESM301_ESM.xlsx`, lines 4–16) three of thirteen oligos (`C6 Amino`,
`rGrGrG`, `rG3p`) carry `ACAT` and the other ten, `rGrG+G` among them, carry `AC`.

🟡 **The model uses G&T-seq's printed `AC` form** for scM&T because Angermueller
explicitly delegates the RNA preparation to that protocol. Since Angermueller does
not itself print the TSO, the model retains the inherited-inference marking: a
preamble caveat with both published forms, an `INFERRED — …` template-switch
caption, and `<inf>` around `AC` in the oligo list and the drawn cDNA state.

🔴 The G&T S1 supplement and historical Nextera XT guide still need review before
the linked-source checklist can be closed. The direct scM&T source does not verify
whether its laboratory used any undocumented TSO change from the cited G&T method.

## Bounded source pass (2026-10-17, round breadth-first-17)

Read: PMC4770512 JATS (`nmeth.3728_PMC4770512.xml.txt`) Results lines 66–68 and Online
Methods lines 80–94 ("Sample collection & Single-cell sequencing", "BS-seq read
alignment"); Springer `MOESM608` Supplementary Fig. 3 caption lines 20–32. The bundle
scraper finds no oligo printed by nmeth.3728 itself; its hits are the Smallwood
scBS-seq and Picelli Smart-seq2 tables already used above.

🟢 Cells are sorted into RLT Plus lysis buffer with 1 U/µl SUPERase-In and processed by
G&T-seq; after physical mRNA/gDNA separation the DNA is eluted into 10 µl H₂O.

🟢 Bisulfite conversion uses EZ Methylation Direct (Zymo) on purified DNA with AMPure XP
beads present; desulphonation and purification are on Zymo magnetic beads, eluting into
the first-strand mastermix. First- and second-strand primers carry 3′ random hexamers,
biotin capture of first-strand products is **omitted**, and an extra 0.8× AMPure XP
purification is done between second-strand synthesis and PCR. The graph now draws that
cleanup and its discarded supernatant before the DNA library PCR.

🟢 Nextera XT is used at one-fifth volumes. Both libraries are sequenced on a HiSeq2000
with 125-bp paired-end reads; the model now declares 125 cycles for Read 1 and Read 2 of
each library. Index-read cycles are not stated and remain unfixed.

🟢 Reads are trimmed by 6 bases (the 6N random-priming portion) and aligned as
non-directional, "due to multiple rounds of random priming with oligo 1".

Open: G&T-seq S1 supplement, Smallwood binary Word supplement and historical Nextera XT
guide (unchanged deferrals); number of oligo1 priming rounds and strand-synthesis
polymerase/conditions (delegated to scBS-seq, not re-read in this pass); index sequences
and index-read cycles.
