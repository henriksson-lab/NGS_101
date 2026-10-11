# Hybrid BAG-seq — joint single-nucleus DNA and RNA sequencing

## Sources read

- 🟢 Li *et al.*, “Joint single-cell DNA and RNA sequencing with hybrid BAG-seq,” *Genome Biology* (2025), DOI `10.1186/s13059-025-03790-5`, PMCID `PMC12476631` (defining paper and methods).
- 🟢 Li *et al.*, “Single-cell DNA sequencing based on hydrogel whole-genome amplification,” *Genome Research* (2020), DOI `10.1101/gr.253047.119`, PMCID `PMC6961581` (original BAG platform and printed capture oligos).
- 🟢 Stepwise hybrid protocol deposited as protocols.io DOI `10.17504/protocols.io.rm7vzqr92vx1/v1`.

## Molecular interpretation

- 🟢 A nucleus is co-encapsulated with acrylamide monomer, a 5′-Acrydite random-TG DNA primer and a 5′-Acrydite poly(dT) RNA primer. Both primers become covalently embedded in the gel.
- 🟢 Hybrid BAG formation uses 85 °C for 5 min rather than the DNA-only protocol's 95 °C for 12 min, avoiding severe RNA destruction.
- 🟢 Reverse transcriptase and DNA polymerase make gel-bound RNA- and genomic-DNA-derived copies. A TSO gives cDNA a double-stranded NlaIII site; Klenow extends captured genomic DNA.
- 🟢 NlaIII cleavage exposes a sticky end. Round 1 ligation-extension adds barcode 1 plus varietal tag 1; round 2 hybridization-extension adds barcode 2 plus varietal tag 2 and a common adapter; split PCR adds barcode 3.
- 🟢 The complete BAG barcode is 17 nt (5 + 6 + 6). The 12-nt varietal tag comprises 4 + 4 random bases plus four genomic bases two bases from the CATG cut site.
- 🟢 Pooled amplicons are tagmented and sequenced 2 × 150 bp on NovaSeq 6000. Read 1 begins with an 80-base tag/adapter region; Read 2 starts from a random tagmentation end.
- 🔴 Complete base-level barcode-family, common-adapter and final sequencing-primer sequences are not all printed in the papers. Their roles and published lengths are preserved without inventing bases.

## Verbatim capture oligos

- 🟢 DNA: `/5ACryd//iSp18/TGTGTTGGGTGTGTTTGGKKKKKKKGKKKKKKKKNN`
- 🟢 RNA: `/5ACryd//iSp18/AAGCAGTGGTATCAACGCAGAGTNNWNNNSTTTTTTTTTTTTTTTTTTTTTTTTTTTT`

These are retained with modification syntax and IUPAC symbols exactly as printed by the original BAG paper.
