# CoLa-seq: accepted-manuscript preparation and RNA lariat topology

Zeng et al., **Profiling lariat intermediates reveals genetic determinants of early and late co-transcriptional splicing**, [Molecular Cell, DOI10.1016/j.molcel.2022.11.004](https://doi.org/10.1016/j.molcel.2022.11.004).

## Source identity and closed investigation

🟢 The retrieved original is the **accepted author manuscript** NIHMS1854138/PMC10448999, not the definitive published body. Header, introduction line72, workflow lines90–94 and STAR Methods216–230 were inspected in its original JATS text. The definitive DOI route returned HTTP403. No complete final/accepted-version equivalence is established. The automatically retrieved2021preprint and its oligo workbook are superseded and excluded from chemistry authority.

The mandatory bundle scraper was run; its preprint-only oligo hits were not adopted. The accepted Methods explicitly cite Table S4 for depletion probes and Table S5 for library oligos. Those actual documents remain unresolved: the PMC manifest lists unavailable supplemental files without a verified Table-to-filename mapping. Desired primer sequences are not invented unavailable disclosure documents. No confirmed absence or complete source audit is claimed. Source investigation closed before implementation, within the bounded pass.

## Evidence-supported preparation

🟢 K562 chromatin-associated RNA is isolated under cold conditions with EDTA and without added MgCl2 to suppress continued splicing. The fractionation reference uses50million cells, whereas the source library input is RNA from200million cells. Array multiplicity/scaling between these scopes is unassigned; the model must not invent four physical fractionation reactions or one unscaled reagent dose. The accepted JATS prints200mL NRB for the final chromatin-pellet resuspension; intended units remain unassigned without a publication-error claim.

The source retains roughly5% subcellular fractions for western/RNA QC, separates cytoplasmic and nucleoplasmic supernatants from chromatin pellets, extracts RNA using Trizol followed by phenol/chloroform/isopropanol, treats with Turbo DNase and repeats extraction. These physical fractions and actual production/QC stocks remain distinct.

For lariat enrichment, SPRI selection removes short RNA; Ribo-Zero depletes rRNA; biotinylated antisense probes from Table S4 bind abundant ncRNAs before C1-bead removal; RNA is kept from the **supernatant**. Following cleanup and Turbo DNase removal of excess probes, RppH decapping and **Terminator5′→3′ exonuclease** treatment remove susceptible linear RNA. This is not the RNase-R treatment of the separate lariat-seq validation branch. A second SPRI selection removes short material; enrichment is not perfect lariat specificity.

Approximately500ng enriched RNA enters the library route. A pre-adenylated adapter with a source-declared **six-base5′ UMI** is ligated with truncatedK227Q T4RNA ligase2 for3h at22°C. Source conditions are17.5%PEG8000,1x ligase buffer,1µM adapter and10U/µl ligase. Whole adapter order and3′ blocking chemistry are unassembled pending Table S5. Excess adapter is removed by SPRI selection. The source RT primer and SuperscriptIII copy RNA for1h at50°C; ExoSAP-IT treatment is followed by EDTA/NaOH, cDNA cleanup, CircLigase circularization and15–18 cycles of Phusion PCR with Table S5 primers. Reverse PCR primers introduce sequencing indices. A6%TBE gel selects175–750bp DNA; the source reports3nM loading and HiSeq4000 paired-end100nt. Exact index assignment/acquisition, whole oligos and complete cDNA/circle/PCR products remain unassembled.

## A lariat is a branched RNA

🟢 The introduction explicitly specifies a **2′–5′ phosphodiester bond** between the branch nucleotide's2′ oxygen and the intron's5′ splice-site phosphate. The ordinary3′–5′ backbone continues through the branch nucleotide into a free3′ tail. A nascent lariat tail may retain downstream exonic/intronic RNA; an excised lariat tail ends at an intron3′ splice site. These remain distinct biological classes.

The inspected source says RT stops at the2′–5′ linkage, and observed CoLa-seq read starts peak one nucleotide downstream of branch points. The source-specific copied interval therefore excludes the branch nucleotide. This readout rule is local; a generic topology helper does not predict RT behavior, enzyme resistance, splicing or folding.

[Illustrative lariat geometry](lariat_example.svg) is generated from a typed model. **Its short RNA sequence and branch coordinate are invented**, solely to demonstrate bond connectivity and the branch/tail distinction. They are not an actual K562 substrate, not a library oligo, and not a substitute for Table S5. The real preparation keeps unknown genomic lariat sequences and occupied products explicitly unassembled.

## Validation scope

The lariat constructor checks supplied RNA and branch coordinates, conserves each nucleotide once, retains distinct2′ and3′ backbone ports and rejects invalid references. Its geometry must not display a free5′ terminus at the covalently joined splice site. A narrow material-consistency claim checks declared backbone connectivity; it does not establish that splicing occurred or validate the library reaction. The protocol graph separately retains source inventories, biological fractions, aliquots and unassembled molecular transitions. A matching topology is not whole chemical validation.

## Source-qualified preparation page

The public page now uses a validated protocol graph with explicit reviewed source inventories and physical-fate typing. Its full workflow and unresolved chemistry remain visible in an expandable section. Whole experimental molecules and source-unassigned read geometry are not supplied by defaults. The preparation remains accepted-manuscript-qualified. The invented lariat reference is not drawn or used as a real genomic input on this page.
