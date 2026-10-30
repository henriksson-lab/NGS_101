# easyCLIP: published workflow and source boundaries

Porter, Miao, Yang et al. (2021), **easyCLIP analysis of RNA-protein interactions incorporating absolute quantification**, Nature Communications12,1569. [Defining article](https://doi.org/10.1038/s41467-021-21623-4).

The positive preparation boundary is on-bead ligation of two fluorescent adapter classes, followed by membrane isolation and capture through the adapter's poly(A) tract before RT/direct PCR. This is a physical library-preparation change, not only crosslink-rate quantification. The poly(A) capture tract belongs to the adapter; this does not restrict targets to naturally polyadenylated RNA.

Three originals inspected: published main PMC7946914 Results78–104 and Methods229–231; published supplementary workbook MOESM3/SupplementaryData1 sheetsLysisToRT,RTtolibraries,L5linkers and selectedOtheroligos rows1–68/109–139/168–192; supplementaryMethods/MOESM1 selected method/read-processing source413–434. Primer scraper run during acquisition:214candidate hits in the then-available bundle; final rerun returned231candidate hits/5files with two large files skipped; those large files and their relevance remain unreviewed, not declared irrelevant or sequence-free. Other workbooks, cloning/standard sheets and figures remain unreviewed. Workbook shared-string XML also checked for selected complete precursor/primer values; no source value was manufactured in a text twin. Source reading stops at these three originals, with no linked-protocol followup.

🟢 Main supports two ordering branches: kinase/5′ligation/L5pooling/phosphatase/3′ligation, or phosphatase/3′ligation/L3pooling/kinase/5′ligation. The selected private graph follows the workbook's detailed3′-first path. Main says preparations are typically combined before a single gel lane. Pool membership and any earlier barcode-specific pooling are unassigned; no exact specimen identity is invented.

🟢 Workbook: UV0.3J/cm² on ice; lysis/sonication/clarified supernatant; antibody-bead purification/washes; RNaseONE30°C3–10min; coldstringency/NT2wash; ATP-freePNK3′reaction37°C30min, magnetically remove reaction without an intervening wash;3′ligation16°C2h-to-overnight; wash;5′PNK/ATP37°C15min;wash;5′ligation16°C2h;wash. Main has distinct sonication and RNase/time settings, kept separate rather than silently assembled into this workbook case. Actual RNase daughters remain until physical separation.

🟢 Selected postligation pool: LDS75°C15min; NuPAGE180V30–50min;nitrocellulose400mA27min;image and cut the minimal region and higher;375µlPKbuffer+25µlProteinaseK50–55°C45min. Transfer extract to20µloligo(dT)beads/600µlBIB4°C20min, thenBIB/NT2/PBSwashes. Main room-temperature capture orKCl/SDS-removal alternatives are unselected. Fluorescence imaging is a measurement, not consumption of the whole library. Optional10percentWestern/dot-blot aliquots are not selected in this graph.

🟢 Workbook heat-elution adds12.4µlwater and3µl1µMRTprimer,95°C3min, then transfers13.3µlsupernatant and adds6.7µlRTmixture. Beads/untransferred liquid remain a separate unresolved fraction. RT feeds PCR directly. The workbook prepares8reactions of master mix but distributes7tubes, each47µlmix plus≤3.3µlRTproduct. These are upper bounds, not seven maximal transfers from20µlRT. RT remainder and unused reagent stock stay unresolved. PCR usesPrimeSTARMax/SYBR,12–18cycles with exponential-phase stopping; optional9+5–8cycling andPhusion substitution are not selected. NucleoSpinGel/PCRcleanup uses85percentethanol instead of kitwash. Whether the seven reactions are pooled before or after cleanup is not assigned: the graph uses an explicit unassembled aggregation/workup boundary, rather than seven invented separate purifications. Source IlluminaHiSeq mention does not assign a particular read length/run.

Selected stock inventories,5′→3′, remain distinct from complete dye-conjugated operational substrates:

- L5Hybrid3 precursor: `/5AzideN/CTTGTGTGTGTGAAGCGAAGGGTATTAGGCrNrNrNrNrNrNrN`.
- L3_12 precursor: `CTGAGATCGGAAGAGCGGTTCAGAAAAAAAAAAAA/iAzideN/AAAAAAAAAAAA/3Bio/-3′`.
- RT/P3.tall.V4: `GGCATTCCTGCTGAACCGCTCTTCCGATCT`.
- P3solexa: `CAAGCAGAAGACGGCATACGAGATCGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT`.
- P6solexa: `AATGATACGGCGACCACCGAGATCTACACTCTTTCCCCTTGTGTGTGAAGCGAAGGGTA`.

The model imports canonical P5/PEPCRprimer2/read2 substrings. Selected L5 sample tag and RNA random window have structured meanings in the precursor inventory only. The final dye attachments/end states, full L3 architecture and complete RNA/cDNA/primer landing geometry remain unassembled; no raw azide precursor is silently treated as an already coupled fluorophore reagent. PCR copies do not inherit nucleotide biotin or a fluorophore by default.

Published alternatives and unresolved discrepancies:

- WorkbookL5linkers literally have seven`rN`units, andOtheroligos explicitly says “random rN L5 (7 nt)”. SupplementaryMethods read-processing text425 describes a “random hexamer UMI”. The model retains the seven-position precursor inventory, with no claim that it establishes the six-base computational handling or exact final read window.
- WorkbookparagraphRT specifies53°C40min/55°C10min/80°C10min, but its adjacentSSIVRTProgram table specifies53°C50min/55°C10min/80°C5min. The selected graph follows the paragraph; main supports40min53°C with or without10min55°C and does not resolve the workbook heat-inactivation difference.
- The selected L5 precursor handle is `CTTGTGTGTGTGAAGCGAAGGGTA`, while P6/RT-standard inventory uses `CTTGTGTGTGAAGCGAAGGGTA`. The original workbook strings preserve thisTG difference. Which source stock revision was paired is unassigned; the model does not silently delete bases or force perfect PCR annealing. This is an unresolved inventory relationship, not a proved universal publication error.

No public schematic is promoted. A future public rendering must retain all these source locations, values, choices and uncertainties. Remaining work includes exact dye-conjugate/ends, operational PCR and read geometry, run and tube assignments, and directly cited ProtocolExchange revision1 DOI10.21203/rs.3.pex-1333/v1. irCLIP's prior adapter is explicitly listed in the workbook, but this does not authorize importing its full recipe. No source information is confirmed absent, and no complete source/linked-material audit has been performed.
