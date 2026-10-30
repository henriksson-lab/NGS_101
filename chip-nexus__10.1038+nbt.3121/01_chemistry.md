# ChIP-nexus: published detailed protocol scope

He, Johnston and Zeitlinger2015, *ChIP-nexus enables improved detection of in vivo transcription factor binding footprints*, DOI10.1038/nbt.3121.

## Versions and originals

🟢 PublishedMOESM434 SupplementaryProtocol1, datedDecember2014, PartsC–E and oligos-table originalp17 inspected; associatedAMNIHMS647911/PMC4390430 Methods qualified comparison only. The publishedXLSMOESM435 header identifies sample/sequencing results, not an adapter source; only header scope inspected. These are the three originals within cap. Mandatory scraper20candidates/4convertedfiles; source table rendered when checking modifications. No current2019/2021 updated protocol substituted.

PublishedSI controls this preparation: dATP37°C tailing,95°C5min initial denaturation before CircLigase,95°C1min cut-oligo anneal,15–18cyclePCR. AssociatedAM prints ATP for the tailing reaction,95°C5min cut-oligo anneal,80% postcutethanolwash and18cyclePCR; the publishedSI uses cold70% wash. it is not merged into the operative publishedSI recipe or asserted to be a published erratum. Final-main applicability remains unreviewed; its HiSeqSE50 acquisition is qualified explicitly.

## Source components

🟢 AdapterU `/5Phos/GATCGGAAGAGCACACGTCTGGATCCACGACGCTCTTCC`; adapterB `/5Phos/TCAGNNNNNAGATCGGAAGAGCGTCGTGGATCCAGACGTGTGCTCTTCCGATCT`; cut oligo `GAAGAGCGTCGTGGATCCAGACGTG`. Source5randomUMI bases and fixedTCAG are separate; the Scene derives their paired core and actual3′T ligating end. The fixed5′extension is filled before lambda digestion. Read orientation of the selected copied class yields reverse-complementUMI thenCTGA; the entire9nt prefix is not a9ntUMI.

🟢 UniversalPCR primer `AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATC*T`; B01 `CAAGCAGAAGACGGCATACGAGATCGTGATGTGACTGGAGTTCAGACGTGTGCTCTTCCGATC*T`. The terminalT remains a nucleotide and the terminal phosphorothioate is an internucleotide bond, not a3′extension blocker. B01 is an explicit published component, not an assigned study sample; physical primer indexCGTGAT is reverse-complement of table barcodeATCACG. No sample-sheet orientation or full index-read geometry is assumed.

## Operations and physical fates

🟢 Enter crosslinked chromatin already immunoprecipitated onDynabeads; sourcePartA/B full cell/antibody branches are outside selected scope. WashA–D/Tris, endrepairT4polymerase/Klenow/PNK12°C30min, wash, Klenowexo−/dATP37°C30min, wash, adapterligation25°C60min, wash, Klenowfill37°C30min, wash, T4polymeraseendtrim12°C5min, wash, lambda37°C60min, wash, RecJf37°C60min and threeRIPAwashes. Digestion fragments persist until these actual washes.

🟢 Elute65°C20min and collectsupernatant; reversecrosslinks65°C6h; RNaseA37°C2h; PK55°C2h; phenol/phenolchloroform extraction and ethanol precipitation. Split12µl into1µl analyticalqPCR and11µl library. Denature95°C5min/ice, CircLigase60°C1h, annealcutoligo95°C1min→25°C30min, BamHI37°C30min and ethanol precipitate. Circle/relinearization are selected physical classes, not perfect yields; originalDNA, uncutcircles, cutoligo fragments and enzymes are retained until source workup. Phusion15–18cycles,2%gel/MinElute (gel dissolved atRT), Bioanalyzer analyticalallocation and qualifiedHiSeqSE50 readout follow. Computational9nt trimming/deduplication does not remove physicalDNA.

## Boundaries and lineage

PublishedSI PartC explicitly invokes Rhee/Pugh2012 detailedChIP-exo (reference1), while the main also cites original2011ChIP-exo; exact detailed-version scope remains a component ancestry target, not proven complete2011recipe equality. PartE explicitly derives initial library preparation from iCLIP2010 (reference2): ssDNA self-circularization/cut-oligo strategy only, **no RNA or RT substrate inherited**. Full genomic adapter/circle/cut/PCR/read geometry and directly cited historical originals remain unreviewed. Private incomplete graph, no full-source or accepted missing-information audit.
