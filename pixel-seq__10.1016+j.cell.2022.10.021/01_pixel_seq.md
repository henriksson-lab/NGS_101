# Pixel-seq — research notes

Fu et al., *Cell* (2022),
[doi:10.1016/j.cell.2022.10.021](https://doi.org/10.1016/j.cell.2022.10.021).

## Source read

- 🟢 Defining article, STAR Methods and resource table (PMCID: PMC9691594).

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

## Transcript library

1. 🟢 Place a fresh 10-µm frozen section on the dry capture gel and hybridize RNA for
   15 min. Maxima H-minus reverse transcriptase copies captured poly(A) RNA for 1 h at
   42 °C while the modified TSO adds the opposite universal handle.
2. 🟢 Proteinase-K digest the tissue, extend the second strand and elute with KOH.
3. 🟢 Extend once with the seven-random-base UMI primer, remove unused primer with
   thermolabile Exonuclease I, and PCR-amplify with the TSO and TruSeq primers.
4. 🟢 Tagment 1 ng cDNA and complete the sequencing library with Nextera XT. Libraries
   were run PE150 on NovaSeq.

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
- 🟢 TruSeq libP5 primer:
  `AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT`.
- 🟢 Nextera index primer: `CAAGCAGAAGACGGCATACGAGAT[INDEX]GTCTCGTGGGCTCGG`.
- 🟢 Polony sequencing primer: `TTCCATCATTGTTCCCGGGTTCCTCATTCTC`.

The spatial barcode is a random 24-base polony identity rather than a finite bead
whitelist. An informative final Nextera fragment must retain the barcode, seven-base UMI
and transcript sequence; other random tagmentation products do not contribute mapped
spatial molecules.
