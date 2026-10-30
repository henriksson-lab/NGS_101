# LM-Seq

[Hou et al., A cost effective RNA sequencing protocol for large scale gene expression studies](https://doi.org/10.1038/srep09570), Scientific Reports 2015.

🟢 Sources read: published PMC4381617 main Results/Discussion/Methods, extracted lines55/73/83–95 and original Figure1, PDFp2; published supplement, Supplemental Methods steps1–8 (PDFp2–4) and Table2 first index orderings. The mandatory bundle scraper was run; relevant main91 and supplement480 pointers were checked against these originals. Investigation closed after two published originals. This is a source-qualified private partial preparation model, not a finished schematic or a complete source review.

## Selected preparation

🟢 The selected branch starts with100ng purified total RNA and uses NEBNext poly(A) isolation with15µl beads per sample, prewashed twice with100µl2×binding buffer and resuspended in50µl2×binding buffer before adding50µl RNA. One wash follows binding. Factory bead manufacture remains opaque; further historical binding/incubation details are delegated to the unread kit protocol. Elution uses8.25µl Elute/Prime/Fragment stock:5µl water,2.5µl5×SmartScribe first-strand buffer and0.75µl20µM HZG883. Heat at85°C7min, cool on ice and transfer an actual6.5µl supernatant off the beads. Beads and untransferred material have separate fates; the selected supernatant carries the pre-added primer into RT.

🟢 Add1µl water,0.25µl100mM DTT,1µl10mM-each dNTPs,0.25µl RNase inhibitor and1µl SmartScribe RT to that6.5µl, giving10µl. Incubate23°C10min then42°C30min. No fresh RT-primer addition is specified here. HZG883 is the reverse complement of canonical small-RNA RA3 plus a random hexamer; its N6 is a priming population, not a UMI. RNA target sequences and actual occupied primer sites are unassigned.

🟢 Add0.5µl RNaseH and0.5µl RNaseA to10µl RT material, making11µl;37°C15min then95°C5min. Add20µl AMPure beads, retain bead-associated material, discard supernatant, wash twice with200µl80% ethanol, dry and elute12µl EB. Transfer10µl; untransferred eluate is unresolved rather than automatically discarded.

## Oligos and a published direction discrepancy

🟢 Main Methods91 gives HZG883 `CCTTGGCACCCGAGAATTCCANNNNNN`; HZG885 `/5phos/AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTddC`; LC056 `AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT`; index primer `CAAGCAGAAGACGGCATACGAGATnnnnnnnnnnGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA`. The model imports matching canonical sequence bodies. The HZG885 terminal ddC is an appended nucleotide, giving{{= chem.ADAPTER_LENGTH }}nt including its terminal C; it is not a modification that replaces the preceding T. Its5′ phosphate and blocked3′ end are retained on the physical donor. The Table2 first ordering ZH001 supplies index `ttggttaggt`; this is an example, not an assigned study sample. Index choice is a10nt sample identifier, and Table2 contains actual alternatives. Index-read orientation is not inferred from the primer ordering alone.

🟢 Published Figure1 caption, original PDFp2, says the modified oligo is ligated to the **5′ of single-stranded cDNA**. Main Methods91 and Results55 instead specify a donor with5′ phosphate and terminal3′ ddC, ligated by T4 RNA ligase1. Both source descriptions are retained here.

🟡 The model's intended junction is cDNA3′OH to donor5′ phosphate, derived from the disclosed donor ends and ligase reaction. It does not silently adopt the caption's5′ cDNA assignment. Source adapter-end geometry is known narrowly; the complete occupied cDNA, its precise termini and whole ligated product remain unassembled. The local intended-substrate constructor marks the cDNA/end interpretation inferred by default. This direction discrepancy remains unresolved by the capped source pass; no further source or correction was opened.

```
{{scene: chem.adapter_scene()}}
```

## Ligation, amplification and reporting

🟢 SI step4 uses10µl cDNA,3µl10×RNA ligase buffer,9.5µl50%PEG8000,3µl10mM ATP,1.5µl20µM HZG885 and3µl T4 RNA ligase1 stock in30µl;22°C12–16h. Add20µl stop solution (10mM Tris-HCl8.0/15mM EDTA), then90µl AMPure for the50µl reaction. Elute22µl/transfer20µl, add36µl AMPure, elute22µl/transfer20µl. Cleanup is a physical workup with uncertain recovery; unreacted primers, adapter, off-products and copied DNA are not assumed perfectly removed.

🟢 SI step4d makes60.1µl PCR:20µl transferred material,8.5µl water,30µl FailSafe PremixE,0.6µl enzyme,0.5µl20µM LC056 and0.5µl20µM sample-index primer. Program95°C1min,18cycles95°C10s/55°C10s/68°C3min, then68°C7min. Original donor ddC chemistry is not copied onto newly synthesized DNA. Full primer landing, whole complementary products and exact run-primer identities remain unassembled; the source's Small RNA index-primer name does not establish an exact nucleotide ordering.

The graph decomposes the one actual30µl PremixE stock into copresent nucleotide and buffer components; these are not independent fresh additions. An N10 index window denotes the unassigned choice of sample-index primer, not a degenerate synthesis pool.

🟢 Add60µl AMPure to the60.1µl PCR (source calls this1:1), elute22µl and transfer20µl. Quantify a separate aliquot with Qubit/PicoGreen. Indexed pooling is an actual physical operation; other prepared library contributions and pooling amounts are unassigned. The selected less-than48-sample branch uses nominal1:1 AMPure, elutes14µl and transfers12µl. The48/96 option instead elutes27µl/transfers25µl and is not merged into this branch. Pooled Qubit/BioAnalyzer QC consumes separate aliquots; remaining stock/loading allocation is explicit. Main Methods reports HiSeq2500 single-end51base sequencing plus10base index. Alignment, deduplication and expression analysis do not physically purify the library.

## Remaining scope

The source-transcription guard has one narrow local formatting exception for HZG885: its parser strips terminal `ddC` entirely and therefore sees only the33nt body. The model preserves the appended physical C. A local check verifies the whole printed `/5phos/...ddC` ordering against main Methods91, and separate checks enforce the34nt donor and blocked acceptor end. Three other whole graph-input oligos remain under strict text matching. This exception does not validate ligation yield, resolve the caption discrepancy or excuse source transcription changes.

🟡 Historical poly(A) kit operations, complete RNA/cDNA molecules and priming occupancy, end heterogeneity, intended ligation junction applicability, adapter/primer off-products, full PCR products, chosen sample indexes/pool membership and exact sequencing-primer landing remain open. The proposed rRNA-depletion and paired-end adaptations, low-input RT/cycle options and short ligation-time comparison stay outside the selected100ng/overnight/single-end branch. Their mention does not create extra method counts. Predecessor and vendor originals remain ordinary unread targets; no absence claim or second absence audit is made.
