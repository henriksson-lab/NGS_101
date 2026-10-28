# CPT-seq: preserve fragment association until physical compartmentalization

Amini et al. (2014), *Haplotype-resolved whole-genome sequencing by contiguity-preserving transposition and combinatorial indexing*, [defining article](https://doi.org/10.1038/ng.3119).

This capped first pass models a representative transposition-well/physical-well route. It remains notes-only: exact historical PCR-primer polymers, source-version applicability and physical index-cycle mapping are not fully established. No complete-protocol marker is emitted.

## Three originals inspected

- Associated accepted manuscript `ng.3119_PMC4409979.xml.txt`: identity/manuscript flag; Results54–66; Methods104–130 covering input DNA, transposome assembly, transposition/EDTA, dilution/SDS, PCR and acquisition. It is not the final published article PDF.
- Publisher Supplemental Material `ng.3119_41588_2014_BFng3119_MOESM1_ESM.pdf`: correction notice; Figures1–4 captions; original Figures3–4 rendered and inspected on PDFpp4–5, including the molecular connector/index layout and reported cycle order.
- Corrected publisher Table4 workbook `ng.3119_41588_2014_BFng3119_MOESM2_ESM.xlsx`: both sheets, including the correction notice and every transferred oligo, run primer and phosphorylated ME complement.

The mandatory complete-bundle scraper was run with `--all`; all24 workbook oligos were checked against their original rows. Its24 matches are source transcriptions, not proof of transposition or PCR. The gated NIHMS Word supplement was not read or declared equivalent to the publisher documents. Investigation stopped at three originals; cited sources and historical kit originals were not followed.

## Correction and apparent connector difference

🟢 Both publisher supplement files carry a20May2015 notice: the first five rows of Supplementary Table4 were originally missing and are restored. The model uses that corrected workbook, not a truncated historical table.

The rendered Figure3 connectorA appears to end `TCCACCC`, whereas the corrected Table4 oligos use `TCCACGC`. The model uses the complete workbook oligos. This apparent graphic-letter difference remains unadjudicated; no claim about the oligos actually manufactured follows. Both locations are identified on the private preview. Do not silently replace the workbook sequence with a diagram reading.

## Source-qualified route

1. 🟢 Anneal each indexed transferred oligo with phosphorylated19nt pMENTS. Load the source EZ-Tn5 stock and combine the8 P5-side and12 P7-side stocks into96 combinations. This does not establish an exact dimer-composition distribution.
2. 🟢 Transpose1ng HMW human DNA per20µL reaction with2.5pmol transposome,55°C10min. Add an equal volume of40mM EDTA,37°C15min. Bound transposase operationally preserves fragment association; EDTA quenching is not modelled as SDS-like protein removal.
3. 🟢 Pool indexed transposition reactions, dilute the25pg/µL pool to1pg/µL in TE, and distribute10pg into each physical well. The source starts with96 transposition combinations and uses96 physical wells. These define possible compartment codes, not an assertion that every code contains one unique molecule or has equal yield.
4. 🟢 Add SDS **after** physical allocation,55°C15min, to dissociate transposase-associated complexes. The source adds i7 PCR primer and BSA before this step; the i7 oligo also acts as a nonspecific surface blocker. Protein-held contiguity is not an unbroken covalent DNA backbone, and no source-specific fragment-size distribution is predicted here.
5. 🟢 Add historical Nextera i5 primer and PCR master mix. Initial extension72°C3min;98°C30s;20cycles98°C10s/63°C30s/72°C3min. The second indexing tier identifies physical compartments. Exact historical FC-121-1012 oligo members and complete PCR geometry remain unassembled; no current kit recipe or same-end amplification prohibition is imported.
6. 🟢 Pool50µL per physical-well reaction, purify on the named Zymo500 column, wash twice and elute2mL; assess library quality. Only actual column flowthrough and washes are discarded.

🟡 `ContiguityState` records association before SDS and dissociation afterward. It does not specify a covalent join between independently drawn fragments, simulate the number of cut sites, calculate collision rates or predict survival. The precursor/well populations are symbolic representatives of the source's repeated operations, not fabricated experimental well identities.

## Corrected Table4 oligos

The model imports canonical shared sequence parts only where the workbook itself matches them. It owns the source-specific connectors, extended mosaic ends and index families locally. A1–A8 correspond to source `P5_i5_1–8_Universal_Connector_A_C15_ME`; B1–B12 correspond to `P7_i7_1–12_Universal_Connector_B_D15_ME`. pMENTS has5′ phosphate; no extra modification is printed on the transferred strands or run primers.

| Source oligo | Sequence5′→3′ |
|---|---|
| A1 | {{= cpt.transposon('A', cpt.T5[0]).top() }} |
| A2 | {{= cpt.transposon('A', cpt.T5[1]).top() }} |
| A3 | {{= cpt.transposon('A', cpt.T5[2]).top() }} |
| A4 | {{= cpt.transposon('A', cpt.T5[3]).top() }} |
| A5 | {{= cpt.transposon('A', cpt.T5[4]).top() }} |
| A6 | {{= cpt.transposon('A', cpt.T5[5]).top() }} |
| A7 | {{= cpt.transposon('A', cpt.T5[6]).top() }} |
| A8 | {{= cpt.transposon('A', cpt.T5[7]).top() }} |
| B1 | {{= cpt.transposon('B', cpt.T7[0]).top() }} |
| B2 | {{= cpt.transposon('B', cpt.T7[1]).top() }} |
| B3 | {{= cpt.transposon('B', cpt.T7[2]).top() }} |
| B4 | {{= cpt.transposon('B', cpt.T7[3]).top() }} |
| B5 | {{= cpt.transposon('B', cpt.T7[4]).top() }} |
| B6 | {{= cpt.transposon('B', cpt.T7[5]).top() }} |
| B7 | {{= cpt.transposon('B', cpt.T7[6]).top() }} |
| B8 | {{= cpt.transposon('B', cpt.T7[7]).top() }} |
| B9 | {{= cpt.transposon('B', cpt.T7[8]).top() }} |
| B10 | {{= cpt.transposon('B', cpt.T7[9]).top() }} |
| B11 | {{= cpt.transposon('B', cpt.T7[10]).top() }} |
| B12 | {{= cpt.transposon('B', cpt.T7[11]).top() }} |
| pMENTS | 5Phos/{{= cpt.nx.ME_RC }} |
| Read1 | {{= cpt.EXTENDED_A }} |
| Read2 | {{= cpt.EXTENDED_B }} |
| Index1 | {{= cpt.revcomp(cpt.EXTENDED_B) }} |

Four separate `linked_read_barcode` components share one group: transpositionT5/T7 and PCRP5/P7. They encode virtual/physical compartment provenance, not UMIs or four independent biological samples. The A1/B1 endpoint is an explicitly inferred example of the source part order; physical-well PCR index members remain unassigned.

## Acquisition and the unresolved physical-cycle map

🟢 Corrected Figure4 reports Read1 cycles1–51; Index1 cycles52–67 (T7 eight, then PCRi7 eight); Index2 cycles68–83 (PCRi5 eight, then T5 eight); Read2 cycles84–134. The associated manuscript specifies HiSeq2000 PE51 and a custom recipe, with custom Read1, Read2 and Index1 primers from Table4.

🟡 Figure3 and the workbook place universal connector sequence between the transposition and PCR index regions. A naive contiguous16-base walk from a single primer must not be labelled as both eight-base codes. How the source custom recipe handles those intervening connectors remains unassembled. No dark-cycle count, index2 primer or four-part cycle map is invented. The sequencing plan therefore records that acquisition limitation instead of claiming computed contiguous index spans. Physical molecular part order and reported acquisition-cycle order are separate facts.

## Lineage and remaining questions

The associated article cites Nextera/transposition predecessors Adey et al.2010 and Caruccio2011. This is limited support for the transposition approach, not evidence of equal source adapters, conditions or products. SCI-seq and later linked-read methods are not imported as parents merely because they also use combinatorial indexing.

Deferred: final published procedural applicability of the associated manuscript; historical FC-121-1012 primer members and exact PCR geometry; precise graphic/table connector comparison; custom recipe and physical-cycle mapping; quantitative cut/association/collision behavior; gated Word supplement identity/contents and linked-original review; independent reader. These are source/model tasks, not confirmed omissions.
