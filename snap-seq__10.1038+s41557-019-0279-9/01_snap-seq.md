# snAP-seq

Defining DOI: https://doi.org/10.1038/s41557-019-0279-9

Liu et al. (2019), Sequencing abasic sites in DNA at single-nucleotide resolution.

## Source scope

🟢 Associated accepted author manuscript PMC6589398, chemistry/assay rationale and Figures1–2, Methods105–115. Publisher Supplementary Information MOESM1, custom-adapter generation196–200, genomic extraction240–245, oligo Table1 around635–665 and caption713, synthetic N10 control/deduplication282–284/1241–1251. Publisher SI is used rather than the separately cached associated supplement; no version equivalence is presumed. Full VOR comparison remains pending.

## Selected native HeLa AP route

1. 🟢 Extract source HeLa genomic DNA and sonicate around450bp (publisher SI). Treat with hydrazino-iso-Pictet–Spengler probe1,10mM,40mMNaPi pH7.4, room temperature2h; spin-column purify. Perform CuAAC with CuBr/THPTA/biotin-PEG3-azide37°C2h; filter/wash/elute. Source AP and competing reacted aldehyde classes remain distinct; no perfect conversion/selectivity is assumed.
2. 🟢 Use NEBNextUltraII preparation with the source customP7 adapter replacing the kit adapter. Treat with SAP37°C30min, AMPure-purify. Historical kit end-processing details and full duplex assembly are not substituted from another protocol.
3. 🟢 Prewash streptavidin beads, add poly(dI:dC) blocker to ligated DNA and capture15min; wash6times. Room-temperature100mMNaOH10min removes complementary strands; a further alkali wash and3buffer washes follow. DNA blocker and complementary-strand fractions are not covalently joined to the captured strand or treated as enzymatically destroyed.
4. 🟢 HOT100mMNaOH70°C15min cleaves at tagged AP sites and releases truncated DNA with5′ phosphate; quench immediately with Tris. This target fragment is not the complementary-strand eluate from the previous step, and its missing AP/tag region is not drawn as a retained adduct in the sequenced molecule.
5. 🟢 Fresh PREWASHED streptavidin beads plus blocker capture still-biotinylated off-target material, including reacted5-fU. KEEP the supernatant containing AP-cleavage products and purify with the modified Zymo ssDNA protocol. Second capture enriches the unbound product, rather than retaining the biotin fraction again. No complete background removal is claimed.
6. 🟢 Prime with sourceP7, denature95°C1min, anneal65°C30s, hold37°C30min then add Klenow(exo−). Clean dsDNA, ligate customP5 at20°C30min, AMPure-clean, Q5PCR, quantify, MiSeq orNextSeq. Exact run/cycle layout and complete final product remain unassembled.

## Oligos and branch limits

P7 top/extension primer as written: Me-GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT. P7 bottom: p-GATCGGAAGAGCACACGTCTGAACTCCAGTCAC-NNNNNN-ATCTCGTATGCCGTCTTCTGCTTG-spacerC3; the source explicitly calls N6 the TruSeqLT index. Me is preserved as printed; detailed terminal chemistry is not guessed. P5 top as printed: GAATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT; P5 bottom: p-GATCGGAAGAGCG. The leadingG in the P5top is retained, not silently changed to the amplification primer. ForwardPCR: AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGA. ReversePCR: CAAGCAGAAGACGGCATACGAGAT.

The source's synthetic5-hmU controls contain randomizedN10 blocks used to deduplicate CONTROL molecules; no genomic UMI is imported from them. SMUG1/UNG pretreatment routes and TruSeq-indexed input libraries are separate branches, not part of this selected native-AP enrichment graph.

The private model keeps chemical classes, both capture directions, tagged off-target material and physical blocker/reagent carryover explicit. Full atom-level probe/adduct structures, exact adapter/strand assembly and linked historical kit documents remain deferred; no confirmed absence or finished marker is asserted.
