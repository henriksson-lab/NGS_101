# Light-Seq — research notes

Kishi et al., *Nature Methods* (2022),
[doi:10.1038/s41592-022-01604-1](https://doi.org/10.1038/s41592-022-01604-1).

## Sources read

- 🟢 Defining open-access article and Methods (PMCID: PMC9636025).
- 🟢 Authors' complete Light-Seq protocol and oligo table,
  [doi:10.17504/protocols.io.x54v9jno4g3e/v1](https://doi.org/10.17504/protocols.io.x54v9jno4g3e/v1).
- 🟢 Authors' deposited analysis code, especially `ExtractTwoBarcodes.py` and
  `ExtractThreeBarcodes.py`, for the actual Read 1 parsing boundary.

## In-situ capture and light-directed indexing

1. 🟢 Fix and permeabilize cells or tissue, then reverse-transcribe total RNA with
   RT.5N.3G, `TTTACACGATTGAGTTATNNNNNGGG`. Five degenerate bases followed by three G
   bases prime without requiring poly(A); the 5′ constant becomes a docking site.
2. 🟢 Terminal transferase adds a poly(dA) tract to the free 3′ cDNA end.
3. 🟢 Hybridize one cnvK-bearing barcode strand throughout the sample, illuminate the
   chosen ROI and wash stringently. The photoactivated cnvK makes a covalent interstrand
   crosslink only in the illuminated area. Repeat with orthogonal barcodes for more ROIs.
4. 🟢 RNase H cleaves the RNA in RNA:cDNA hybrids and releases the selected, crosslinked
   barcode–cDNA complexes without destroying the specimen.

The three published barcode strands contain the same 27-base GATE PCR handle and 12
IUPAC-D UMI positions. Their ROI-identifying prefixes are `TATGGA`, `GTTAGG` and
`AGGGTA`; the authors' parser matches those six bases, followed by fixed `TGAGTTAT`.
The protocol preserves the ordered modifications as `[cnvK]`, `[Cy5-3]`, `[Cy3-3]` and
`[Fl-3]`.

## Crossing a covalent junction

🟢 The crosslink prevents ordinary PCR. GATC.20T primes the poly(dA) tract and Bst 2.0
WarmStart extends through cDNA until it stops at the crosslink. A shared seven-base P
domain—potentially eight after Bst A-tailing—allows the newly synthesized end to invade
the opposite strand and resume synthesis on the barcode template. The uncrosslinked
product therefore carries GATC–transcript–P-domain–ROI–UMI–GATE on one strand. GATC and
GATE primers amplify that product; a small qPCR aliquot determines the endpoint cycle.

## Nextera endpoint and reads

🟢 Nextera XT tagmentation fragments the Light-Seq amplicon. Indexing PCR is asymmetric:
an ordinary Next.N7XX primer amplifies the tagmented end, while S5XX.GATE replaces an
ordinary Nextera S5 primer and selects the original GATE end. Thus Read 1 is not an
ordinary Nextera insert boundary. The required P5.GATE custom Read 1 primer places cycle
1 at the 12-base UMI, followed by the six-base ROI code, fixed `TGAGTTAT`, the eight-base
P domain and transcript cDNA. The authors extracted up to 40 transcript bases from R1;
R2 approaches the transcript from the random tagmentation boundary.

🟢 A GATE*.P5* custom i5 sequencing primer is published. It worked on HiSeq, but the
protocol reports inefficient i5 index sequencing on NovaSeq and therefore recommends a
unique i7 index for every pooled sample. The page displays the complete S502/N701 example;
the protocol also publishes S503 and S505–S508 plus N702–N706.
