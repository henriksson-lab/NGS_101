# Original SCOPE-Seq (2018)

[Published defining article](https://doi.org/10.1186/s13059-018-1607-x).

🟢 Three originals: published PMC6305572 main39–61; AdditionalFile1 DOCX SupplementaryMethods bead-production/bead-free libraries/cellular libraries/optical decoding; AdditionalFile4 XLSX TableS2 whole modified oligos. Mandatory scraper25candidates/one extracted representation. No additional original may be opened in this cap. Original SCOPE2018 is kept separate from SCOPE2: it requires distinct optical/sequence barcodes and a separately sequenced bead-free lookup library, which the main explicitly discusses changing in future work. No complete source/absence audit.

## Source-positive preparation boundary

Approximately700,000 commercial ChemGenes“Drop-Seq”beads are distributed across64wells. Before each of two optical-ligation rounds, sonicate5min and incubate4°C12h; collect the detached bead-free DNA supernatant to a separate plate and wash beads. Ligate a selected combination from six5′-phosphorylated optical barcode oligos in each round using Rnl1/ATP/PEG24h. SDS/EDTA stop, pool/wash and repeat with the other six optical oligos. Bead-associated capture sites and detached DNA are physical stocks; lookup data are a later computational association, not a new biological analyte.

For each of128 bead-free samples, select5µl with10µl1µM NexterapolyA primer/NEBuffer2,95°C2min/cool; add2µl10mMdNTP andKlenow, room temperature15min/37°C1h/75°C20min. Select4µl into indexedP7/universalP5/KAPA PCR14cycles; select5µl per reaction for pooling. Preserve untransferred DNA/reaction stocks. SPRI2:1,10%PAGE~180bp, BioanalyzerQC;40%PhiX is a physical DNA-library contribution before NextSeq500R1:26/I1:8 custom-read1 sequencing. It is not merely a loading-buffer condition.

Single-cell branch: live image trapped U87/3T3 cells in microwells before loading prepared dual-barcoded beads. Oil-isolatedTCL/mercaptoethanol lysis50°C20min/25°C90min RNA capture, then wash/unseal. MaximaH-minus RT/dNTP/TSO25°C30min/42°C90min. Add twelve unlabeled OBO-complement oligos, then ExoI37°C45min and TE/NaOH washes. Twelve reversibleCy5-probe hybridization/NaOH-stripping cycles decode optical identities on-chip. This does not physically discard computationally unlinked cells.

Cut device into10regions and extract beads by ethanol/vortex/sonication/centrifugation, keeping regions separate. Wash beads, split each region into two50µlPCRs with SMRTpcr/KAPA; source12cycle programme. Pool each region's products/SPRI0.6:1; select0.6ng per NexteraFC-131-1024 tagmentation, use region-specifici7 and source universalP5 instead ofi5 primer. Two SPRI rounds0.6:1 then1:1; physical qualified-library pooling plus20%PhiX. NextSeq500R1:21/R2:66/I1:8 with customRead1. Do not import a modern Nextera/Drop-seq/SCOPE2 recipe, hidden barcode placement or TSO body from names alone.

## Modified oligo inventory and limits

All12OBOs have source prefix `/5Phos/TTTTTTUTTTT/iSp18/`, an18nt optical code, then `/3ddC/`. Retain literalU, the internal non-nucleotideSp18 spacer, and the additional terminal modified C; do not concatenate the spacer as DNA letters or drop ddC to satisfy a parser. Cy5 decoding probes are optical reagents, distinct from unlabeled protection oligos and sequencer primers.

The actual TableS2 also supplies128 lookupP7 primers, NexterapolyA, universalP5 with terminal `T*A*C`, customRead1, andSMRTpcr. Whole source inventories can be checked without pretending all protected/ligated/captured products are assembled. Exact historical bead-bound capture sequence, TSO chemistry, region-index choices, full junctions and read primer landing remain unassembled. The article gives a12nt bead-identifying sequence barcode; other identifier lengths must not be inferred solely from read counts. Label these as model limits, not confirmed publication absence.

## Private implementation scope

The private graph represents the64-well factory as a batch of separate compartments, not a single ligation mixture. One lookup well and one RNA region are detailed; the other127lookup and9region preparations remain explicitly source-qualified compound branches derived from their actual physical inputs. Detached DNA is pooled only after indexedPCR. Unlabeled protection oligo wholebodies remain unassigned; the exact Cy5 probe18mers are not silently substituted.

The model keeps Sp18 as a non-nucleotide part between DNA components and retains the additional terminal ddC. Whole OBO ordering1: {{= chem.OpticalOligo(chem.CODES[0]).ordering }}; nucleotide count excludingSp18: {{= chem.OpticalOligo(chem.CODES[0]).nucleotide_count }}. This inventory does not claim the complete backbone/ligated-product structure.

Selected source primers: universalP5 {{= chem.UNIVERSAL_RAW }}; NexterapolyA {{= chem.NEXTERAPOLYA }}; lookupP7primer1 {{= chem.LOOKUP_P7 }}; SMRTpcr {{= chem.SMRT }}; customRead1 {{= chem.READ1 }}. Exact read landing remains unassembled.
