# HydEn-seq: sequence a hydrolysis-generated5′DNA end

Clausen et al. (2015), [Tracking replication enzymology in vivo by genome-wide mapping of ribonucleotide incorporation](https://doi.org/10.1038/nsmb.2957).

## Capped sources and qualification

One reader,28October2026; article-associated accepted manuscript PMC4351163, complete HydEn protocol/analysis Methods; publisher Supplementary Information MOESM6, Table1 physically rendered page4 and relevant figure captions. Mandatory full-bundle scraper run. Its concatenated98nt “ARC140” hit is a table-extraction artifact joining ARC49 and ARC76; individual rows/modifications were checked against the rendered original. No sequence is adopted from that concatenated hit.

The complete corrected version-of-record Methods and source-linked originals remain deferred; no complete manuscript/publisher equivalence is assumed. Two chemistry originals inspected; no missing-information acceptance.

## Actual end chemistry and preparation

1. MasterPure extraction from yeast omits RNaseA. Treat1µg DNA with0.3MKOH2h55°C, ethanol precipitate, heat85°C3min. Hydrolysis generates downstream5′OH DNA ends; the hydrolyzed rNMP is immediately adjacent and is not the first base of that downstream fragment. Potential RNA carryover/hydrolysis products are separate species, not covalent DNA extensions.
2. Use **3′-phosphatase-minus T4PNK**,10U37°C30min, then65°C20min and HighPrep purification. This phosphorylates5′OH while preserving the hydrolysis-associated3′phosphate class; do not substitute ordinary PNK with a3′phosphatase activity.
3. Heat85°C3min; T4RNA ligase joins ARC140's3′OH to DNA5′phosphate overnight at room temperature with25%PEG8000 and1mMCoCl3(NH3)6. ARC140 has a **5′amine/C6 linker**, not5′OH. HighPrep purify.
4. Heat85°C3min. Anneal ARC76/ARC77 at room temperature5min; T7DNA polymerase synthesizes the second strand. ARC76's six random bases are priming degeneracy, **not a UMI**; ARC77 has two terminal phosphorothioate bonds. HighPrep purify.
5. KAPA HiFi PCR with ARC49 and a source indexed primer, HighPrep purify and pool. The article reports HiSeq2500 paired-end50base sequencing. Source PCR-cycle counts are not fixed in this Methods paragraph.

## Source oligos and naming discrepancy

| Name | Verbatim5′→3′ sequence/modification |
|---|---|
| ARC140 | /5AmMC6/ACACTCTTTCCCTACACGACGCTCTTCCGATCT |
| ARC76 | GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTNNNN*N*N |
| ARC77 | AGATCGGAAGAGCACACGTCTGAACTCCAGTC*A*C |
| ARC49 | AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT |
| ARC84, represented example | CAAGCAGAAGACGGCATACGAGATACATCGGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT |

The accepted-manuscript Methods text names **ARC79**, whereas the first indexed row in Table1 is **ARC78**. The table also has six- and eight-base index members. The model selects the unambiguous printed ARC84 example and does not silently rename ARC78 as ARC79 or make all indices eight bases long. This source/version difference is shown on the schematic; fullVORMethods remains unread, so it is not classified as a publication error. It does not affect the selected ARC84 sequence.

## Physical product versus analysis

The sequenced end is a captured downstream DNA5′end. Source processing shifts end counts **one base upstream** to assign the hydrolyzed rNMP. Native DNA ends/background can also contribute; the source uses RER/genotype controls and background subtraction, not a chemical guarantee that every read is a newly incorporated ribonucleotide. Nonunique/short and source-oligo-derived read pairs are filtered computationally, not destroyed in the molecular workflow.

Source-local inferred fragment-end assembly and conditional PCR reference can be drawn from the printed primer arms. Exact genomic sequence, random priming positions, selected study index and loaded run primers remain unassigned. RNaseH2-based HydEn adaptations are separate later branches, not silently substituted for this original alkaline preparation.
