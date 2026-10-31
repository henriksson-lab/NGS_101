# SISSOR

Chu et al., *Ultraaccurate genome sequencing and haplotyping of single human cells*, PNAS (2017), [defining published article](https://doi.org/10.1073/pnas.1707609114).

## Sources and selected scope

🟢 Published PMC5703283 main identity/Results64–74/Figure1, MaterialsAndMethods164, mapping168; published SI AppendixI.1–I.5 and I.7–I.8 overview. Original main PDFp6 was rendered to verify the printed MDA primer-end wording. The mandatory whole-bundle scraper returned eight hits: four from published SITableS9 QC primers and four from the automatically fetched superseded preprint supplement. Only the published source contexts were used; the QC primers are not substituted for library primers. The preprint and its supplement are excluded. Other result figures/statistical tables remain outside this capped chemistry pass.

All selected originals are available in the external source store with assigned URLs/checksums. Detailed predecessor methods and historical kit/primer identities remain unreviewed, not confirmed missing. See [retrieval manifest](ref/MANIFEST.md). This is one capped pass, not a complete two-reader original/linked-material audit.

## Molecular and compartment principle

🟢 SISSOR captures one imaged PGP1 fibroblast cell, uses alkaline lysis/strand separation, mixes and randomly distributes long single-stranded chromosomal fragments among {{= sissor.CHAMBERS }} chambers, then separately neutralizes/amplifies/retrieves each chamber. Chamber libraries receive independent sequencing indexes before same-cell pooling. It is distinct from sequencing one unsplit diploid MDA product or separating daughter cells by biological division.

The local finite-duplex helper derives the complementary physical strand from an explicitly caller-supplied DNA sequence and verifies pairing using Scene. Its caller-supplied sequences are inferred by default; no arbitrary human genomic sequence is drawn as established. Explicit partition assignments retain every supplied strand, including a separate unrecovered class. The helper permits complementary strands in different chambers, co-occurrence in one chamber and multiple strand/homolog classes per chamber. It does not predict occupancy, fragment lengths, recovery, phasing or amplification accuracy.

Source main Results74 and discussion explicitly acknowledge possible losses and mixed-strand fragments. The graph exposes both strand populations in the same lysate, then partitioning; it does not invent a physical pooling operation between denaturation and partitioning. One representative chamber is expanded, the other23 matched source preparations are summarized, and possible unrecovered material stays separate rather than entering a library. Known same-cell library pooling occurs only after individual chamber indexing.

## Source on-chip preparation

🟢 Empty MDA chambers are precoated with BSA/random hexamers/dNTPs, purged/dried and UV-sterilized **before cell loading**. The model does not expose the cell/genome to that sterilization step. Fabrication and cell culture/harvest are outside the graph; source preparation is an explicit input condition.

The source ALS is400mMKOH/10mMDTT/1%Tween20. The lysate is mixed10min at20°C using the rotary pump at5Hz and partitioned. On-chip NS is400mMHCl/600mMTrisHCl/1%Tween20. The source named MDA mixture includes Phi29, dNTPs and dUTP; incubation30°C15h is followed by separate recovery into TE and65°C15min heat inactivation. Heat inactivation is not purification. Original templates, uracil-containing copies and reaction species are not assumed absent from the recovered preparation.

🟢 MainMethods164 and originalPDFp6 print **“3′ phosphorylated random hexamers”** at84µM. The model retains that exact source wording and leaves initiation/end chemistry unresolved. It does not replace it by ordinary3′OH, phosphorothioate protection or an invented deblocking reaction. A source-reported MDA outcome is not computational validation of extension from a guessed primer terminus. The unmodified random hexamers in the separate SI transposon-library second-strand step are a different source reagent, not a correction to this statement.

## CoRE source library route

🟢 SII.3 reports CoRE libraries for all Cell3 libraries and half of Cell1. The expanded CoRE view selects the homogeneous Cell3 route. Its operations are:

1. Source UDG/EndonucleaseIV treatment37°C2h,65°C15min on dU-containing MDA products. Base excision, backbone incision and later physical disposal are separate concepts; source operational language describes gap/nick preparation, while exact AP/end/flap products remain unassembled.
2. Exo-minusKlenow/dNTP treatment37°C1h,75°C15min generates source double-stranded fragmentation products from those positions. Displaced/released DNA and nucleobase species do not disappear merely because extension occurs.
3. Historical NEBNext loop-adapter/T4DNAligase reaction25°C10min; USER37°C15min processes the adapter's uracil-containing loop. Full genomic/adapter ends and other susceptible residual uracil substrates remain unassigned. No modern NEBNextUltraII substrate/end recipe is imported.
4. Actual AMPure cleanup precedes indexed KAPA PCR,12cycles; PAGE selects300–700bp products. The chamber index member, width, complete primer polymers and acquired index read geometry are unassembled.

## Transposon source library route

🟢 SII.4 reports this route for all Cell2 libraries and the other half of Cell1. The independent transposon view selects homogeneous Cell2; it is not physically joined to Cell3's CoRE route.

Source MDA products receive ALS at room temperature3min, ice cooling and the separately printed SI NS (400mMHCl/600mMTrisHCl; no unprinted on-chip Tween addition). DNApolymeraseI plus **unmodified** random hexamers/dNTPs generate transposase-accessible duplex DNA. Ampligase nick sealing precedes ethanol recovery.

Historical Nextera transposase/HMWbuffer55°C5min then source protease50°C10min/70°C20min precede exo-minusKlenow end fill37°C15min/75°C20min. Proteolysis does not itself remove protein fragments. A first7cycle PCR adds adapter/barcode sequences and a second7cycle PCR uses Illumina primers1/2; no unprinted intervening cleanup is inserted. Source PAGEQC/AMPure purification and300–700bp gel selection follow. Actual historical adapter/primer/end geometry remains unassembled.

Cell1's mixed allocation of library routes is recorded as a source option; the model does not invent a same-chamber split into both routes or count the routes as separate protocols.

## Sequencing and computational fates

🟢 Main reports pooled same-cell indexed libraries with100bp paired-end Illumina sequencing. SI separately describes initial single-end36/50bp pilot runs and later single-end100 or paired-end100 acquisitions. These are acquisition options, not new MDA operations. The graphs show the reported deep PE100 readout with explicit incomplete physical primer/index geometry.

The chamber library barcode is a structured sample-index class: its sample is the recovered chamber preparation. It is not a newly assigned per-molecule UMI, a shared unique cell barcode or a biochemical repair marker. No fixed index width or actual index sequence is inferred from a current vendor kit.

Mapping, computationally treating paired reads as separate alignments, HMM fragment segmentation, haplotype assembly and complementary-strand consensus are analytical operations. They do not physically join fragments, repair DNA, remove all MDA errors or enforce one haplotype per chamber. Source accuracy estimates are not constructor guarantees.

## Remaining work

Exact MDA primer/end chemistry and historical formulation; full genomic/fragment/adapter/PCR/index strands; CoRE incision/displacement termini; actual historical Peters/Adey/kit versions and final run-primer assignments; individual parallel-chamber expansion and independent source/linked-material review remain open. No unavailable disclosure document is invented for an unassigned experimental index. The known preparation is modeled operationally, but both source routes remain notes-only/private and chemically incomplete.

## Bounded publication review, 30 October 2026

Published main Fig1/Results64–74/Methods164 and published SI I.3–I.5 were re-read in the closed two-original pass; mandatory bundle scraping was repeated. The main's explicit same-cell library pooling statement was verified at74. Capture and per-chamber retrieval now preserve molecular identity instead of declaring generic covalent changes. The known final same-cell library pool is in a curator-written literal source inventory; source coverage is no longer accumulated from construction calls.

The actual main164 MDA MasterMix supplies are represented as copresent component views: source primer, Phi29, dNTP/dUTP and buffer come from one physical stock. They are not extra independent doses or a physically separated new factory. Source stock volume and per-chamber quantities remain unassigned. The source's literal 3′-phosphorylated random-hexamer wording and unresolved initiation remain intact; the source-reported MDA outcome is still an explicit unvalidated handoff, not successful ordinary extension from an invented 3′OH primer.

Both CoRE and transposon routes now retain exactly the original required `historical_geometry` failure. Exact end/initiation, historical NEBNext/Nextera/primer polymers, complete products and acquired read landing remain unassembled. The scalar representative/other23 chamber treatment is not an individually reconstructed array; downstream source nucleotide/kit/stock/QC details remain local follow-ups. No full-source requirement was renamed, excluded or covered by a pooled-library label. No source absence or final-VOR error is claimed.
