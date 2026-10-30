# scG2P: associated-manuscript RNA/DNA co-capture

Defining DOI10.1158/2159-8290.CD-24-0853, *Genotype-to-Phenotype Mapping of Somatic Clonal Mosaicism via Single-Cell Co-Capture of DNA Mutations and mRNA Transcripts*. Current source is NIHMS2135808/PMC12874418.1, an associated author manuscript. Its body is not verified equivalent to the final published article. Registration year is not recipe-version equivalence.

🟢 Methods189–217 defines a physical preparation beyond changing a mutation panel: first Tapestri v2 encapsulation carries lysis/targeted reverse transcription; second encapsulation adds shared cell barcodes and targeted DNA/RNA amplification. RNA targets cross exon junctions; RNA-derived and genomic-DNA amplicons are then enriched into distinct fractions using a biotinylated handle probe.

🟢 RNA forward handle: `GTACTCGCAGTAGTC`; reverse handle: `GTGATACACGACTATGAGCGCTA`. These are partial primer handles, not complete target primers. Whole capture oligo is `/5BioTinTEG/GTGATACACGACTATGAGCGCTA/3C6/`; literal terminal chemistry is retained. A local reverse-complement binding candidate is inferred from primer-handle incorporation, not a complete amplicon proof. SourceResults105 explicitly states noUMIs.

🟢 Original cell/nucleus wash stock already contains1%BSA with DPBS/RNaseinhibitor/DTT; this is co-present input carrier, not a new BSA addition. First-emulsion stocks include10U/ulSSIV,0.5mM EACH dNTP,14.6uM RNA reverse pool, ProteinaseK0.04mg/ml plus stated buffer/detergent. Incubation50°C60min/80°C10min precedes second encapsulation/UVbarcode-primer release and10+12target-PCRcycles. Repeated paragraphs205/207 describe the same second encapsulation, not two repetitions. Emulsion break uses15ul extraction agent; cleanup enzyme37°C60min; two0.72xAMPure cleanups follow.

🟢 Capture oligo0.7uM is annealed99°C5min→ice5min; M270streptavidin capture30minRT separates fractions. Here 'RNA fraction' means RNA-derived DNA amplicons, not an RNA polymer fraction. The source prints10MTris-HCl in binding buffer; this literal is preserved without an asserted correction. DNA fraction gets0.72xcleanup. Each entire recovered fraction enters its own15cycle indexPCR with LibraryMix/kitDNA or customRNAindices. Both are sequenced paired150bp onNovaSeq6000; unreported pooling is not invented.

🟡 Kit molecular components, full target/index primers, barcode identities, whole products/read landings and final-version equivalence remain unassembled. Actual supplemental workbook link is gated/TDM; no claim of sequence absence follows. Companion bulkRNA/Visium/PTA workflows remain separate.

## Selected sources and review scope

- 2159-8290.CD-24-0853_PMC12874418.xml: Associated NIHMS2135808/PMC12874418.1 Methods189–217/Results97/105; author-manuscript identity explicitly retained, final-body equivalence unverified.

This capped pass is not a complete source, linked-information, chemistry or absence audit.
