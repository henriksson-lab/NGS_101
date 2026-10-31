# DNase Hi-C

Defining DOI: https://doi.org/10.1038/nmeth.3205

Ma et al. (online 2014; journal issue 2015), Fine-scale chromatin interaction maps reveal the cis-regulatory landscape of human lincRNA genes.

## Source scope

🟢 Associated accepted author manuscript PMC4281301, Online Methods, Generation of DNase Hi-C and targeted DNase Hi-C libraries (text lines114–158). Publisher Supplementary Information MOESM45, Supplementary Note2 (7983–8026), explains random fragmentation and in-gel ligation. Supplementary Table23 XLSX (the associated NIHMS export is byte-identical to publisher MOESM57, as recorded by the retrieval manifest), Adaptors-DNase Hi-C-protocol and sequencing adaptors and primers sheets, supplies bridge/blocker oligos and library primers. The manuscript identity is retained; complete publisher-article comparison is deferred. Previous preliminary searches in the associated NIHMS supplement are historical, not an additional chemical authority merged into the current pass.

## Selected whole-genome route

1. 🟢 Select K562: 1% formaldehyde, 10min at room temperature; 0.125M glycine quench10min; PBS wash, detergent lysis. H1 uses2.5% formaldehyde and is a separate input condition.
2. 🟢 SDS treatment, nuclei collection/wash, RNase-containing DNase buffer, DNaseI digestion, then SDS/EDTA stop. Collect chromatin with AMPure, wash twice with80% ethanol. Fragmentation and RNase products remain until physical recovery. Publisher Note2 describes DNase products as5′ phosphate/3′ hydroxyl; genomic cut positions and complete strand geometries are unassigned.
3. 🟢 T4DNA polymerase/Klenow end repair, source SDS/PEG/NaCl bead workup; Klenow(exo−)/dATP tailing and recovery. Ligate T-tailed internally biotinylated bridge adapter together with a nonbiotin blunt blocker; the blocker is deliberately included because not every fragment acquires a dA end. It is not a sequencing index or assumed to vanish after ligation.
4. 🟢 PNK phosphorylation, then ligation in0.4% low-melting agarose,16°C4h then25°C1h. Crosslinked fragments and bridge-bearing/unligated/blocker classes remain a physical population; the model does not predict contact yields or assign every ligation a biological contact.
5. 🟢 Melt gel70°C10min; agarase at42°C overnight, concentration, ProteinaseK/SDS and crosslink reversal65°C overnight then55°C2h. Protein/agarose products persist until filtration/precipitation/QIAquick recovery; protease does not remove covalent DNA blockers.
6. 🟢 Source T4DNA polymerase with dATP+dGTP acts on dangling ends, EDTA stop/concentration, then Covaris fragmentation100–300bp and AMPure recovery. This is not immediate physical disposal of every unligated molecule. A detailed biotin-stripping mechanism and exact strand-resolved products are not reconstructed.
7. 🟢 Fermentas end repair, Klenow(exo−)/dATP, source Table23 sequencing-adapter ligation, each with stated bead workup. Capture biotin-containing DNA on prewashed streptavidin beads and wash. No perfect removal of all background or other DNA classes is assumed.
8. 🟢 Whole-genome branch: KAPA HiFi9–11PCR cycles with Table23 indexed primers, AMPure fractionation, HiSeq2000 paired-end150bp. The targeted branch instead uses pre-capture primers, NimbleGen SeqCap and further PCR; it is outside this graph.

## Oligo and model boundary

🟢 Biotin bridge oligos: /5Phos/GCTGAGGGA/iBiodT/C and CCTCAGCT. Blunt blocker: GCTGAGGGAC and CCTCAGC. The local model computes the antiparallel core from these source-written strands, retaining internal biotin and phosphate. It does not assert a fully sealed bridge/genomic junction from adapter pairing alone. Sequencing Adaptor-Hi includes a terminal phosphorothioate; Adaptor-Lo a5′ phosphate. Table23 PCR/barcode rows were inspected but full library/primer assembly and physical index-cycle orientation remain deferred. No barcode member is invented or copied from another Hi-C method.

🟡 Partial states preserve genomic populations, blocker-associated DNA and reagent debris until actual source separations. Unknown historical end-processing details and exact post-ligation duplexes remain explicit. The private preview has no finished-protocol marker and is blocked by an uncovered source item.

## Lineage and deferred work

The defining article explicitly replaces restriction-enzyme fragmentation of Hi-C with DNaseI; this is strategy-level ancestry, not inheritance of a restriction junction or adapter. Publisher Note2 contrasts TCC with the new in-gel approach, rather than claiming TCC reagent inheritance. Targeted capture is an optional branch and should not silently become the whole-genome route. Full VOR comparison, remaining original tables/figures, source-local bridge/genome ligation chemistry and complete run assembly remain deferred; no missing-information claim is accepted.


## Bounded publication adjudication — 2026-10-30, batch 11

🟢 Source whole-genome K562 route contains actual fixation/DNase/end-repair/dA-tail/biotin bridge+blunt blocker, PNK/in-gel ligation, protein/agarose digestion/workup, dangling treatment/shearing, separate sequencing-adapter ligation, prepared streptavidin beads,9–11cycle library PCR and PE150. Known whole sequencing adapters and primers are available in inspected S23; they must not remain generic placeholders in a future whole-order implementation.

🟡 Original required unreviewed_vor_and_assembly remains uncovered: complete defining VOR comparison and post-ligation strand assembly. Actual local defining body is associated accepted manuscript PMC4281301, not final-version-equivalent. Publisher SI and workbook support known components, but do not complete that required main-version/strand assembly comparison. This is unreviewed scope, not evidence that final publication omits detail.

Source scopes inspected: Associated accepted manuscript PMC4281301 XML/PDF: Online Methods114–158 K562 whole-genome fixation through library; three source end-prep tubes/in-gel handling; source approximately10–100ul PCR wording, Nanodrop and PE150; targeted capture separately excluded from selected whole-genome boundary; Published MOESM45 PDF: Supplementary Note2 original text7983–8026 / printedp36–37 DNase5P/3OH and in-gel proximity rationale; Associated S23 XLSX: Actual original sharedStrings XML and table text: all four bridge/blocker whole orders; sequencing adapter pair with top PS and bottom5P; indexed/universal PCR and separate pre-capture/blocking roles.

The associated accepted-manuscript first chromatin dA-tailing paragraph (Online Methods128) prints **0.15 mM dATP**, and also **3 µL of 20 mM dATP in 300 µL**. The latter quantities calculate to **0.2 mM**. These are different source descriptions within that inspected manuscript; the model does not select an intended concentration or assert a final-publication error. Its later library dA-tailing paragraph separately specifies 2 mM dATP and must not be substituted for the earlier reaction. The required VOR comparison remains incomplete.

Bounded remaining physical/model obligations: Retain existing required final-main comparison/post-ligation assembly hold; no source acquisition chase. Restore actual three-tube branch allocation through400ul end repair/300ul dA-tail/200ul bridge reactions and explicit in-gel collection/pool boundaries. Current scalar graph does not resolve actual array cardinality/fates. Expose physical dNTP/dATP/dGTP supplies separately from catalysts, BSA in actual prepared-bead wash, and GlycoBlue in actual DNA precipitation. Do not invent an independently added DNA-ligase dose from RapidLigationBuffer naming alone. Source dA-tail text prints0.15mMdATP and3ul20mM into300ul; the latter calculates0.2mM. Preserve both values as a source arithmetic/concentration discrepancy if this stage is published, without silently choosing. S23 whole Adapter-Hi ACACTCTTTCCCTACACGACGCTCTTCCGATC*T retains PS before terminalT; Adapter-Lo PO4-GATCGGAAGAGCGGTTCAGCAGGAATGCCGAG retains5P. Known input pairing and bridge internal biotin do not establish complete genomic contact product or sample-index selection. Whole-genome PCR source literal approximately10–100ul reactions is ambiguous; do not infer10actual100ul tubes or borrow16tube targeted-route counts. Actual per-library PCR pool/QC and loading remainders need explicit source-qualified declarations.

The existing private preview and normal publication gate remain unchanged. No required source item is excluded or renamed, no unknown occupied molecule is invented, and this is not a confirmed absence claim. No source acquisition or architecture expansion was undertaken.
