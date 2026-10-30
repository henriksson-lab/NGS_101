# TIF-seq — original 2013 paired transcript ends

Pelechano, Wei and Steinmetz, *Extensive transcriptional heterogeneity revealed by isoform profiling*, Nature (2013), DOI [10.1038/nature12121](https://doi.org/10.1038/nature12121).

This is a capped source-positive collection, with a private incomplete local graph. It is not a complete chemistry review or a public schematic. The defining main text read is the associated accepted manuscript PMC3705217, Methods lines 217–290; final main-version applicability remains unreviewed. The published publisher supplement MOESM98 supplies the selected overview and Table 2 (rendered original p38). The associated-manuscript supplement was inspected first for the same selected components; these selected agreements do not establish whole-version equivalence. Mandatory bundle scraper output was checked against source text and rendered original Table S2 p39; the final Table 2 transcription was compared directly.

🟢 The selected main path dephosphorylates non-capped RNA with SAP, purifies, removes caps with TAP, purifies, and ligates a DNA/RNA capping oligo with T4 RNA ligase 1. Column recovery and RNA integrity checking precede a physical split into two aliquots. Each aliquot undergoes oligo-dT SuperScript III reverse transcription, RNase cocktail/RNase H treatment, cleanup, and Phusion PCR with its own matching barcode pair. The published SI places NotI digestion before pooling the sticky-ended aliquots; the selected graph follows that ordering. The associated main instead pools before NotI. This main-version difference is preserved explicitly, not treated as a publication error or proof of whole-version equivalence. Exact pre-pooling digestion conditions/workup are left unassembled across these versions. Dilute T4 DNA ligation, exonuclease I/III treatment, physical workup, sonication, cleanup and biotin capture follow. Historical NEBNext/Illumina forked-adapter preparation, Phusion enrichment and gel size selection create the sequencing preparation. The associated main and published supplement state HiSeq 2000 paired-end 105-base acquisition.

🟢 Published Table 2 selected oligos, verbatim with modifications:

| Oligo | As written |
|---|---|
| 5oligocap | `dCdAdCdTdCdTrGrArGrCrArArUrArCrC` |
| 3cDNANotI_A | `TATAGCGGCCGCCTGGAGTTTTTTTTTTTTTTTTVN` |
| 3cDNANotI_B | `TATAGCGGCCGCCTCCACTTTTTTTTTTTTTTTTVN` |
| 5BioNotI_A | `TATAGCGGCCGCTA[BtndT]CACTCTGAGCAATACC` |
| 5BioNotI_B | `TATAGCGGCCGCTA[BtndT]ATCACTCTGAGCAATACC` |
| 3AmpNotI_A | `CATGTATAGCGGCCGCCTGGAG` |
| 3AmpNotI_B | `ACATGTATAGCGGCCGCCTCCAC` |

The capping oligo physically contains a DNA block followed by an RNA block. `[BtndT]` denotes an internal biotin-dT, not an extra sequence label or terminal biotin. The local constructors keep these distinctions. 🟡 The selected-copy helper computes complement and reverse block order after converting RNA alphabet for a DNA copy; it does not establish a complete RT register or whole PCR backbone. Full primer landing and exact barcode-coordinate feature mapping remain unassembled; the A/B tube controls are distinct from later sample multiplexing and from unique molecule identifiers. No inferred barcode meaning is attached to arbitrary sequence windows.

The RNA-PET comparison is described as conceptual similarity to an independently developed method, not documented recipe inheritance. The 2014 detailed TIF-seq protocol is a possible directly linked follow-up target, not another counted protocol and not imported here. Long-molecule enrichment and non-capped controls remain branches outside this graph.

The published supplement's overview prints `GCCGGCCGC` as the NotI recognition string, while its actual selected oligos contain `GCGGCCGC`. Its figure caption describes intermolecular ligation while the detailed overview and associated main describe dilute intramolecular circularization with intermolecular chimeras as contaminants. The model retains named NotI digestion with unassembled exact cut geometry and explicitly chooses the detailed intramolecular path; it does not silently turn the nine-base printed string into a corrected molecular recipe. No public schematic is promoted by this partial comparison.

No oligo, reagent or readout is claimed confirmed absent. Unassembled historical end preparation, full forked-adapter inventory and PCR/read geometry, main final-version applicability, and the directly cited oligo-capping predecessor remain open work. Physical digestion fragments remain with the reaction until actual cleanup; computational chimera filters operate after acquisition and do not remove physical DNA during preparation.
