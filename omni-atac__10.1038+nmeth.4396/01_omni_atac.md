# Omni-ATAC — research notes

Corces et al., *Nature Methods* (2017),
[doi:10.1038/nmeth.4396](https://doi.org/10.1038/nmeth.4396).

## Sources read

- 🟢 The defining article and Online Methods, available as the PMC author manuscript
  (PMCID: PMC5623106).
- 🟢 Supplementary Protocol 1, *Omni-ATAC protocol*, linked by the article and mirrored
  by the Corces laboratory.
- 🟢 Buenrostro et al., *Current Protocols in Molecular Biology* (2015),
  [doi:10.1002/0471142727.mb2129s109](https://doi.org/10.1002/0471142727.mb2129s109),
  for the inherited ATAC amplification workflow. Its methods refer the complete primer
  list to Buenrostro et al., *Nature Methods* (2013),
  [doi:10.1038/nmeth.2688](https://doi.org/10.1038/nmeth.2688); that paper's
  Supplementary Methods (MOESM200, opening lines) print Ad1_noMX and Ad2 members.

## What Omni changes

1. 🟢 Fifty thousand cells are resuspended in 50 µl ATAC-RSB containing 0.1% NP-40,
   0.1% Tween-20 and 0.01% digitonin and held on ice for 3 min. Multiple detergents make
   permeabilization less cell-type-dependent.
2. 🟢 One millilitre ATAC-RSB with 0.1% Tween-20 washes the permeabilized nuclei before
   they are pelleted. This removes released mitochondria rather than carrying them into
   transposition.
3. 🟢 The 50-µl transposition mix is 25 µl 2× TD buffer, 2.5 µl transposase (100 nM
   final), 16.5 µl PBS, 0.5 µl 1% digitonin, 0.5 µl 10% Tween-20 and 5 µl water. It is
   incubated for 30 min at 37 °C with mixing at 1,000 r.p.m. PBS is an intentional
   Omni-ATAC component (Online Methods, PMC5623106 JATS lines 285–315; Supplementary
   Protocol 1, transposition steps 1–9).
   - 🟢 The transposase is **not** commercial TDE1: the Online Methods cite Picelli et
     al. 2014 (ref. 26) for it, and the Author Contributions state that S.W.C.
     "produced the Tn5 transposase complex that was used in all of the experiments"
     (JATS line 898). 2× TD buffer is home-made: 20 mM Tris-HCl pH 7.6, 10 mM MgCl2,
     20% dimethylformamide (Supplementary Protocol 1, Buffers and Reagents). Model
     correction (2026-10-17): earlier versions of this note and page said TDE1.
   - 🟡 The defining article does not print the Tn5 loading oligos. The S5/S7 plus
     mosaic-end adapter geometry is inferred from the 3' ends of the Ad1/Ad2 primers
     that must anneal to it; Picelli 2014 has not been read in this pass.
4. 🟢 Transposed DNA is purified on a Zymo DNA Clean & Concentrator-5 column (eluted
   in 21 µl). All ~20 µl enters a 50-µl NEBNext 2× Master Mix PCR with 2.5 µl each
   25 µM Ad1 and Ad2: 72 °C 5 min, 98 °C 30 s, then 5 cycles of 98 °C 10 s, 63 °C
   30 s, 72 °C 1 min. A 5-µl aliquot is run as a 15-µl SYBR qPCR for 20 cycles to
   choose the additional cycles; the remaining 45 µl returns to the cycler with no
   added reagents, and the product is purified on a second Zymo column (20 µl water).
   Libraries are quantified with the KAPA Library Quantification Kit (Supplementary
   Protocol 1, Pre-amplification, qPCR and Final amplification sections).

## Shared ATAC library geometry

🟢 Omni-ATAC changes nuclei preparation and the transposition buffer, not the adapter
architecture. 🟡 The in-house Tn5 is taken to transfer the canonical 19-base mosaic end
with S5/S7 adapters, because the Ad1/Ad2 primers end in S5/S7 plus mosaic-end bases. Of the possible
S5/S5, S7/S7 and S5/S7 end combinations, only the heterologous S5/S7 fragment can acquire
both flow-cell arms in PCR.

🟢 The defining protocol delegates amplification to the original ATAC procedure. The
published example pair is:

- `Ad1_noMX`: `AATGATACGGCGACCACCGAGATCTACACTCGTCGGCAGCGTCAGATGTG`
- `Ad2.1_TAAGGCGA`: `CAAGCAGAAGACGGCATACGAGATTCGCCTTAGTCTCGTGGGCTCGGAGATGT`

The eight bases physically present in Ad2.1 are `TCGCCTTA`; the sample-sheet Index 1
value named by the primer is its reverse complement, `TAAGGCGA`. There is no i5 index in
this historical construct. Later dual-index ATAC implementations should not be projected
back onto the defining Omni-ATAC page.

## Sequencing boundary

🟢 The final library is paired-end Nextera geometry: Read 1 starts after the S5-side
mosaic end, Index 1 reports the eight-base Ad2 index, and Read 2 starts after the S7-side
mosaic end. The Online Methods (Sequencing, JATS lines 265–272) state that all deep
sequencing used 2 × 75-bp reads on a HiSeq4000 and all low-depth sequencing 2 × 75 bp
on a MiSeq; pooled libraries for the HiSeq4000 were first PAGE size-selected (>100 bp)
to remove excess primers. Model correction (2026-10-17): the page previously declared
no genomic read length; it now uses 75 + 75. 🟡 The Index 1 cycle count is not printed;
the page shows eight cycles because the Ad2 index is eight bases.

## Not drawn / open

- Frozen-tissue (Supplementary Protocol 2, density gradient) and 500-cell variants.
- Pre-treatment of cultured cells with 200 U/ml DNase and FACS of primary T cells.
- Which Ad2 index members were used per library; Picelli 2014 Tn5 loading oligos and
  Protocol Exchange protex.2017.096 not read in the 2026-10-17 bounded pass.

## Source-owned validation migration, 2026-10-18

Local declaration owns the original Ad1/Ad2 single-index reference, in-house Tn5/DMF/PBS/detergent conditions and source-specific recovery/PAGE/QC branches. Removed shared variant/endpoint/primer recipes and test imports. Full loading-derived template geometry and reference primer assignments retain inferred construction styling. A distinct gap-fill intermediate precedes five-cycle PCR; both PCR geometries use actual predecessor states instead of the final library. Candidate non-library end products remain unresolved rather than presumed physically discarded. Acquired reads are2×75 only: eight molecular index bases no longer imply an acquired index read.

Source scope: Corces2017 published PMC5623106 Sequencing lines261–272 and in-house transposase Methods lines307–313 plus AuthorContributions line898; published MOESM3 SupplementaryProtocol1 cleanup/preamp/qPCR/final amplification lines88–169; MOESM1 inspected for identity/figure context only and not used as oligo authority.

Deferred: Exact in-house Tn5 loading oligos and study-specific index member remain unassigned; original Ad1/Ad2 reference does not verify those reagents; Cited Picelli2014 loading chemistry not re-read within this cap; Transposition/gap filling and full PCR product relations remain unresolved despite primer-site checks. This bounded check does not establish complete molecular validation or an accepted source-absence claim.
