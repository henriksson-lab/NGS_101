# time-ChIP — pulsed histone-cohort capture

Deaton et al. *Enhancer regions show high histone H3.3 turnover that changes during differentiation*. eLife 5:e15316 (2016). [Published article](https://doi.org/10.7554/eLife.15316).

## Sources and boundary

Read published v2 main Methods343–355/383–385, main reference identity for Bowman2013, and original Supplementary file5, Sheet1 plus empty Sheet2/3. The main PDF p18 was rendered to verify apparent unit problems; PDF/XML are representations of the same main original. Two originals used. Initial main indexing also displayed review correspondence embedded in that JATS original; no correspondence was used as chemistry authority. The mandatory oligo scraper was run and hits reviewed. Other supplement originals and directly cited predecessor recipes remain unread. This is a private partial preparation graph and notes, not a full source or absence audit.

A pre-existing cohort of SNAP-tagged histones is covalently labeled, then its associated chromatin DNA is physically captured. It differs from metabolic labeling of newly synthesized proteins. H3.3/H3.1, ESC/NSC and time points remain branches of one method. The modeled representative arm is H3.3 ESC at0hr; the three later harvests remain a physical unresolved cohort, not invented extra protocols. Genetic integration, cell-line establishment and alternative ChIP/RNA-seq/MNase-accessibility assays are outside this selected library arm.

## Source preparation

🟢 Main351: SNAP-histone expression is induced using500ng/ml doxycycline for four cell doubling times (ESC48hr; NSC72hr). Cells are labeled with10µM CP-biotin,40min37°C, then washed twice with10ml medium, incubated20min37°C in10ml complete medium, and washed again. One-fifth is harvested as0hr; the remainder is plated for3/6/12hr harvest. The literal source cell/volume statement is0.5–2×10^8cells in50µl; it is preserved without asserting a corrected volume.

🟢 Main353: nuclei are extracted/counted and equal cell numbers used; source500µl suspension contains20mM Hepes pH7.7,20mM KCl,0.5mM EDTA,300mM NaCl,1mM DTT and protease inhibitors. Digest with120–140U MNase and3mM CaCl2,20min25°C; EDTA/EGTA quench,0.05% NP40 solubilization,12,000g5min spin. Take400µl actual supernatant into15µl streptavidin magnetic beads, preblocked with literal50mg BSA and200µg/ml yeast tRNA. Incubate1hr4°C. The untransferred supernatant is explicit; no perfect500µl recovery or inferred100µl remainder is asserted.

🟢 Wash beads3×5min, then elute DNA in250µl10mM Tris pH8,100mM NaCl,10mM EDTA,0.5% SDS,10min30°C with DNase-free RNaseA. Phenol/chloroform, then chloroform extraction and Qiagen PCR-purification column precede PicoGreen quantification. No proteinaseK, formaldehyde crosslinking, anti-HA antibody or TruSeq RNA first/second-strand synthesis is imported from the separate comparator assays.

## Apparent published unit problem

Main353 and the actual published PDFp18 both print wash-buffer **500 M NaCl**, beside20mM Hepes/20mM KCl/3mM CaCl2/5mM EGTA/0.5% NP40/1mM DTT/protease inhibitors. This is an apparent unit error; its intended value is unresolved. The model retains the literal source with a warning and does not silently replace it with500mM. The original PDF also prints50mg BSA for bead preblocking; this literal value is retained without a proposed correction. Ordinary PDF text extraction changes someµ symbols to m; the rendered original, not that extraction artifact, establishes200µg/ml yeast tRNA and250µl elution volume. These distinctions prevent a repository extraction problem from being mislabeled as a publication error.

## DNA controls and library boundary

🟢 Main385 and actual S5: synthetic150bp DNA controls with source GC labels30–70 are mixed equimolarly. Fixed absolute control mass is added to each time-point library, equal to one ten-thousandth of the0hr recovered DNA mass; do not make the amount proportional to each later time point. The five literal specified strands are local model constants and independently checked against original XLSX cells. Their partner strands and physical end chemistry remain unassigned; a listed strand is not a claimed single-stranded or double-stranded complete molecule. An actual control-stock aliquot and unused remainder are explicit, with `SPIKE_IN` molecular provenance into the library boundary.

🟢 Main385 delegates DNA-library preparation to Bowman et al.2013, *Multiplexed Illumina sequencing libraries from picogram quantities of DNA*, [DOI10.1186/1471-2164-14-466](https://doi.org/10.1186/1471-2164-14-466). That directly cited original was not read within this pass. The graph keeps this boundary opaque, including the stage at which controls join; it does not reconstruct adapters, enzymes, ends, amplification or sequencing primers from the name “Illumina”. Source PE50 on HiSeq2500 and four barcoded samples/lane are retained. Index cycles and exact geometry remain unassembled.

## Deferred work

Read the directly cited Bowman2013 main and chemistry supplement in a later depth-one pass; verify adapter geometry, index identities, spike-addition stage and reaction ordering before completing this library. Actual native nuclei-extraction details and other temporal arms also remain unassembled. S6 is titled primer sequences and its scraper hits resemble expression qPCR targets; its original has not been read, so its full scope is deferred. S1–4 are source-identified result tables, not inspected as chemistry originals. No actual failed needed URL arose in this pass, and no oligo/reagent/read-layout absence is asserted. Source-specific shared wet-lab ancestry was not established beyond the delegated library citation; the parent is not currently a separate registered protocol.
