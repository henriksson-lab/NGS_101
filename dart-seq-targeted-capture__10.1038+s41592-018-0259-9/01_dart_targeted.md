# Simultaneous multiplexed amplicon sequencing and transcriptome profiling in single cells

Defining publication: https://doi.org/10.1038/s41592-018-0259-9.

## Identity and scope

🟢 This is Saikia etal.'s targeted-capture DART-seq, DOI10.1038/s41592-018-0259-9, not the distinct APOBEC/m6A DART-seq DOI10.1038/s41592-019-0570-0. Custom targeted primers are added to a subset of Drop-seq bead oligos while ordinary poly(dT) capture remains. B-cell repertoire and whole-transcriptome reads retain the same cell identity/UMI semantics; there is no new antigen barcode or VH:VL molecular fusion.

## Sources inspected and identity checks

- Published article's accepted-manuscript Methods in PMC6378878: bead synthesis, cell preparation, single-cell library preparation; B-cell Results.
- Published DART-seq detailed protocol, publisherMOESM3. Manifest MD5 establishes exact identity with the cached PMC `supplement-2.pdf`; inspected primer design, bead ligation and sections4–5 through sequencing.
- Published SupplementaryTable1, publisherMOESM4. Manifest MD5 establishes identity with PMC `Sup_Table_1.xlsx`. Read exact IgG primer/splint and modification fields.

These are three originals. Publisher supplementary software archive was fetched automatically; its large reference FASTA is not a chemistry original and was not scanned. Its conversion failures are not evidence of a missing protocol document. Superseded preprint routes were not used.

## Base-resolved IgG bead interface

🟢 IgG_CustomPrimer: `/5Phos/GGGGAAGACCGATGGGCCCTTGGTGG`.

🟢 IgG_Splint: `CGGTCTTCCCCAAAAAAAAAAAA`.

🟡 The source IgG primer's first11bases are reverse-complementary to the splint's11-base primer arm. The splint's12A arm pairs with a terminal12T window of the bead oligo. Scene annealing computes both pairings and the adjacency of bead3′OH to custom-primer5′phosphate. The12T window remains a partial view, not the full bead oligo. The published detailed protocol also provides the full source bead baseclasses, which are now declared separately as described below.

🟢 Anneal primer and splint95°C→14°C, then bead ligation with T4DNA ligase/PEG4000 at37°C1h. Heat65°C3min and ice do not discard splint or free oligos; these remain until TE-SDS/TE-TW bead washes. Only primer is covalently appended; splint is not part of the sealed bead primer.

## Downstream local graph

🟢 Single cells and beads are encapsulated for lysis/RNA capture; emulsion is then broken and RNA-bearing beads recovered/washed before bulk MaximaRT plus source template-switch oligo `AAGCAGTGGTATCAACGCAGAGTGAATrGrGrG`. PostRT washes, ExonucleaseI37°C45min and further washes precede amplification. Source specifies one SMARTPCR primer species `AAGCAGTGGTATCAACGCAGAGT`, four initial plus11subsequent cycles with KAPAHiFi, followed by0.6×SPRI cleanup. The graph declares one primer species explicitly while leaving inward-site geometry unavailable; it does not invent a second species.

🟢 NexteraXT FC-131-1096 prepares libraries using the source P5 hybrid primer with terminal phosphorothioates. Both main Methods and detailed protocol specify NextSeq500,20bp Read1 and130bp Read2; detailed protocol additionally supplies customRead1 primer `GCCTGTCCGCGGAAGCAGTGGTATCAACGCAGAGTAC`. Full primer landing and final strands are not yet assembled.

## Deferrals and lineage

Source microfluidic capture/lysis occupancy, occupied bead/RNA junctions and full source-pool product assemblies, exact template-switching geometry, selected Nextera indices, tagmentation end classes and final-strand assembly remain open. All are unassembled, not confirmed absent. Main Methods explicitly says single-cell preparation was carried out as Drop-seq reference2: component ancestry is the Drop-seq capture/RT/amplification backend, with bead-primer modification local to DART. No platform-name recipe import or relation to the m6A homonym is justified.

## Completion and verification

Private incomplete source-local graph. Required source inventory explicitly retains unassembled full-strand chemistry; a valid local reaction contract does not establish source-review completion. Mandatory primer scraper was run against fetched sources and relevant hits checked in originals. No two-reader missing-information pass has been completed.

## Bounded assembly follow-up, 30 October 2026 — model-scope hold

The actual detailed-protocol title/author page and primer page were rechecked. The file's MD5 `287feaffd15ebc1b1b6a48907c7722d2` matches the exact publisher DOI-attached MOESM3 original recorded in the manifest. No printed revision/date identifier appears on its title page. This establishes that supplement-file identity, not complete accepted-main/final-VOR equivalence or applicability of a newer protocol revision.

The model now declares source bead **T7 + PCR handle + J12 + N8 + T30** classes; J12 is the source split-pool cell identifier and N8 the source random UMI, with members unassigned. The selected IgG component appends its known 5′phosphorylated primer. This does not assemble the full biological B-cell primer pool or every occupied RNA/product. Source TSO DNA-prefix plus terminal **rGrGrG**, SMART primer, custom Read1 and custom P5-hybrid **terminal T*A*C** are declared and independently checked against the original protocol. No Drop-seq or modern Nextera constructor supplies those values.

Actual source 6–10 ligation reactions are represented as one representative reaction plus the remaining **5–9** reaction array, with distinct bead/probe allocations and actual product pooling; no arbitrary exact total count is chosen. Other actual B-cell toehold members remain copresent and unassembled, rather than inventing an IgG-only biological experiment. Droplet lysis/capture is followed by actual emulsion break, physical bead recovery/washes and only then bulk RT. One source 200 µl prepared RT mixture supplies copresent primer/TSO, catalyst, nucleotide and buffer components. Selected 2,000-bead PCR contributions, actual 25 µl KAPA ReadyMix components, source cDNA/library Qubit/BioAnalyzer aliquots, selected 600 pg Nextera input and loading/remainders are explicit. The original loaded tagmentation components share a source inventory; free adapters are not invented as an independent fresh addition.

The required `unassembled_source` item remains unchanged and uncovered. Full source-pool, intermediate and final-strand assemblies are still unimplemented, including occupied tagmentation/index interfaces and exact read landing. This is a **blocked model scope**, not a source-absence or fabricated missing-document claim. No complete marker is emitted.
