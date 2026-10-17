# Pixel-seq — research notes

Fu et al., *Cell* (2022),
[doi:10.1016/j.cell.2022.10.021](https://doi.org/10.1016/j.cell.2022.10.021).

## Source read

- 🟢 Defining article as the NIH author manuscript (PMC9691594, PMC Cloud JATS,
  licence TDM): STAR Methods "Polony template construction" through "Library
  construction and sequencing" and "Transcript mapping", plus the Key Resources Table
  oligonucleotides (round 2026-10-17-breadth-first-15 re-read).
- The publisher version of record (doi.org/10.1016/j.cell.2022.10.021) and the six
  NIHMS1846275 supplements (Tables S1–S3, Figures) are behind the publisher page or
  the PMC download gate; they remain unread. The author manuscript is the accepted
  text, not a superseded preprint.

## Array replication

1. 🟢 Covalently graft plus and minus primers to a crosslinked polyacrylamide gel, seed
   random templates and bridge-amplify approximately one-micrometre clonal polonies.
2. 🟢 USER-linearize the stamp polonies. Press stamp and copy gels together at 95 °C,
   cool to 60 °C for Taq extension across the gel interface, then re-denature. The template
   remains on the stamp while its complementary copy stays covalently attached to the copy
   gel.
3. 🟢 Bridge-amplify the copied strands for 20 cycles. Copy gels preserve the coordinate
   pattern closely enough that only sibling gels need in-situ sequencing to map the
   24-base random barcode; the other replicas reuse that map.
4. 🟢 Digest copy-gel polonies with TaqI for 1.5 h at 60 °C to expose their 3′ poly(dT)
   RNA-capture ends.

Stamp-gel details (🟢 STAR Methods "Preparation of stamp gels"): templates are first
PCR-amplified with BA(+)/BA(−) (Taq, 15 cycles) and gel-purified; seeded templates are
extended by Taq at 74 °C, then bridge-amplified with Bst for 22 formamide-denaturation
cycles at 60 °C; USER (37 °C, 1 h) linearizes the polonies at the BA(+) uracil.

## Transcript library

1. 🟢 Place a fresh 10-µm frozen section on the dry capture gel and hybridize RNA for
   15 min. Maxima H-minus reverse transcriptase copies captured poly(A) RNA for 1 h at
   42 °C while the modified TSO adds the opposite universal handle.
2. 🟢 Proteinase-K digest the tissue (55 °C, 30 min) and wash with 2× SSC/0.1% SDS and
   0.1× SSC.
3. 🟢 Synthesize the second strand with the biotinylated TSO primer and Bst 2.0
   WarmStart polymerase (65 °C, 15 min), then elute it with 0.08 M KOH and neutralize
   with Tris.
4. 🟢 Extend once with the seven-random-base UMI primer (Q5 Ultra II; 95/65/72 °C),
   remove unused primer with thermolabile Exonuclease I (degraded in the tube, not
   separated), and PCR-amplify with the "TSO" and TruSeq sequencing primers and Q5
   HotStart (4 + 8 cycles for an OB section) to 5–10 ng.
5. 🟢 Tagment about 1 ng cDNA with Nextera XT and amplify with the TruSeq libP5
   primer and Nextera index primers.
6. 🔴 Run read lengths are not printed in the inspected author manuscript; NovaSeq S4
   appears only in a cost estimate (Discussion). The earlier note's "PE150 on NovaSeq"
   could not be traced to an inspected passage and was removed.

## Exact published oligos

🟢 STAR Methods, "Polony template construction" (original JATS lines 183–200), prints
the entire seed template. It joins P5, a T7 promoter, 24 random spatial-barcode bases,
a polony sequencing site, a 19-T capture tract, `TCGA`, a Cy5-probe site, the long
internal segment ending in a second `TCGA`, and the P7-side sequence. The public page
renders the full source string from the protocol module; the random barcode remains
`N` in that ordering sequence.

- 🟢 The STAR Methods resource table lists bridge-amplification primers BA(+)
  `AATGATACGGCGACCACCGAGAUCTACAC` and BA(−)
  `CAAGCAGAAGACGGCATACGAGAT`. The plus primer contains a published U base for
  USER cleavage; it is not a repository typo for T.
- 🟢 The gel-grafted variants are PS-BA(+)
  `/5PS/TTTTTTTTTTAATGATACGGCGACCACCGAGAUCTACAC` and PS-BA(−)
  `/5PS/TTTTTTTTTTCAAGCAGAAGACGGCATACGAGAT`.
- 🟢 Cy5 probe: `/5Cy5/CTTGCGCTAGTCCGGCAACCTCGGTGG`.
- 🟢 UMI primer: `CTACACGACGCTCTTCCGATCTNNNNNNNCCGCTAATACGACTCACTATAGGG`.
- 🟢 Template-switching oligo: `AAGCTGGTATCAACGCAGAGTGAATrG+GrG`.
- 🟢 Biotinylated TSO primer: `/5’BIOSG/AAGCAGTGGTATCAACGCAGAGT`.
- 🟢 TruSeq sequencing primer (used in cDNA PCR): `CCTACACGACGCTCTTCCGATC`.
- 🟢 TruSeq libP5 primer:
  `AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT`.
- 🟢 Nextera index primer: `CAAGCAGAAGACGGCATACGAGAT[INDEX]GTCTCGTGGGCTCGG`.
- 🟢 Polony sequencing primer: `TTCCATCATTGTTCCCGGGTTCCTCATTCTC`.

## Final library and reads

🟡 Derived from the oligos above (computed in `tools/chem.py`): the UMI primer's 3′
anchor `CCGCTAATACGACTCACTATAGGG` matches the template T7-promoter segment, so it
anneals to the eluted second strand and extends across `ATCT`, the 24-base barcode,
the polony sequencing-primer site and poly(T) into antisense cDNA. Its 5′ arm is the
3′ end of TruSeq Read 1, which the TruSeq libP5 primer extends to P5. The library is
therefore a TruSeq/Nextera hybrid with **no i5 window**: P5–TruSeq Read 1–UMI(7)–T7
anchor–barcode(24)–polony site–T19–cDNA–ME–S7–i7–P7. TruSeq Read 1 reads the UMI
first and the barcode after 28 anchor bases; Read 2 reads transcript from the mosaic
end. The earlier repository model (dual-indexed Nextera library, UMI after the barcode)
was a modelling error, corrected in this round. The i7 length is unprinted
(`[INDEX]`); the page shows an inferred eight-position placeholder.

The spatial barcode is a random 24-base polony identity rather than a finite bead
whitelist. An informative final Nextera fragment must retain the barcode, seven-base UMI
and transcript sequence; other random tagmentation products do not contribute mapped
spatial molecules.

## Open questions

- Version of record and supplements (Tables S1–S3) unread: gated publisher/PMC routes.
- Sequencing run read lengths and i7 index length: not in the inspected manuscript;
  GEO GSE186097 run metadata not inspected (named follow-up).
- Whether the cDNA-PCR "TSO" primer is the biotinylated TSO primer (Methods say "TSO
  and TruSeq sequencing primers").
