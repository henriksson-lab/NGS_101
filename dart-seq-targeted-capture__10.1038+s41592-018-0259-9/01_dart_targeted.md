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

🟡 The source IgG primer's first11bases are reverse-complementary to the splint's11-base primer arm. The splint's12A arm pairs with a terminal12T window of the bead oligo. Scene annealing computes both pairings and the adjacency of bead3′OH to custom-primer5′phosphate. The12T window is a partial view, not a claim that the full bead oligo has exactly12T or that upstream barcode/handle geometry has been assembled.

🟢 Anneal primer and splint95°C→14°C, then bead ligation with T4DNA ligase/PEG4000 at37°C1h. Heat65°C3min and ice do not discard splint or free oligos; these remain until TE-SDS/TE-TW bead washes. Only primer is covalently appended; splint is not part of the sealed bead primer.

## Downstream local graph

🟢 Single cells and beads are encapsulated, cells lysed, and cDNA synthesized with MaximaRT plus source template-switch oligo `AAGCAGTGGTATCAACGCAGAGTGAATrGrGrG`. PostRT washes, ExonucleaseI37°C45min and further washes precede amplification. Source specifies one SMARTPCR primer species `AAGCAGTGGTATCAACGCAGAGT`, four initial plus11subsequent cycles with KAPAHiFi, followed by0.6×SPRI cleanup. The graph declares one primer species explicitly while leaving inward-site geometry unavailable; it does not invent a second species.

🟢 NexteraXT FC-131-1096 prepares libraries using the source P5 hybrid primer with terminal phosphorothioates. Both main Methods and detailed protocol specify NextSeq500,20bp Read1 and130bp Read2; detailed protocol additionally supplies customRead1 primer `GCCTGTCCGCGGAAGCAGTGGTATCAACGCAGAGTAC`. Full primer landing and final strands are not yet assembled.

## Deferrals and lineage

Source microfluidic capture/lysis occupancy, full bead barcode/UMI architecture and modification representation, exact template-switching geometry, selected Nextera indices, tagmentation end classes and final-strand assembly remain open. All are unassembled, not confirmed absent. Main Methods explicitly says single-cell preparation was carried out as Drop-seq reference2: component ancestry is the Drop-seq capture/RT/amplification backend, with bead-primer modification local to DART. No platform-name recipe import or relation to the m6A homonym is justified.

## Completion and verification

Private incomplete source-local graph. Required source inventory explicitly retains unassembled full-strand chemistry; a valid local reaction contract does not establish source-review completion. Mandatory primer scraper was run against fetched sources and relevant hits checked in originals. No two-reader missing-information pass has been completed.
