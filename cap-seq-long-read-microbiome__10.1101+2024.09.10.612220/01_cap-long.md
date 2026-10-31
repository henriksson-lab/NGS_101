# CAP-seq (long-read microbiome)

Li et al., *Large-scale single-cell long-read genomics enables high-resolution microbiome profiling*, [bioRxiv10.1101/2024.09.10.612220v3](https://www.biorxiv.org/content/10.1101/2024.09.10.612220v3), **posted13October2025**. This pass uses that explicitly identified preprint revision. A published successor was not identified/qualified in the current API/title checks; that is not a confirmed absence of later publication. The ResearchSquare manuscript and earlier bioRxiv revisions are not merged into this model.

This is the microbial DNA capsule/barcoding/ONT method, separate from CXXC affinity CAP-seq and the Science RNA CapSeq method. A shared word or capsule concept does not establish shared primers, polymers or library chemistry.

## Inspected originals and limits

🟢 Defining v3 PDF header/title and Methods19–23 (extracted text711–887), computational barcode search893–901, Figure1b/SupplementaryFigure1 captions, directly cited source identities and Data/Code availability. Original PDFp22 was rendered to verify the linear-amplification recipe and capsule dissolution; final original pages52–54 were inspected to distinguish result figures from a primer table. The mandatory complete-bundle scraper returned one sequence hit, an analysis motif rather than a full primer. This is a capped original-source pass, not a complete independent source review.

The exact v3 supplemental-material index is available and its title/HTML links were inspected. No referenced SupplementaryTable1 original was located during this pass. The generic fetch returned no attachments; neither an empty attachment list nor the one scraper hit establishes absent disclosure. **SupplementaryTable1 is an actual named source**, not an invented missing-disclosure document. Its identity/location remains unresolved. The main PDF contains source figures/captions but these do not supply a verified full primer inventory. The code repository is identified by the main article as bioinformatics code and was not used as a replacement primer source. See [retrieval manifest](ref/MANIFEST.md).

Historical ABclonal tagmentation chemistry, the barcode-bead fabrication original, material/permeability predecessors, the actual primer table and ONT protocol revisions remain unreviewed. A generic vendor product page or a current kit recipe is not silently substituted.

## Compartment and molecular scope

🟢 The source encloses microbial cells in dextran/PEGDA aqueous two-phase droplets, then UV-gels the outer PEGDA-rich shell. Capsules retain large DNA while allowing reagents to exchange. Figure1b/SupplementaryFigure1 report a dsDNA retention transition around400bp, with100bp/793bp assay examples. This is an empirical permeability description, not a deterministic cutoff or a universal recovery guarantee.

🟢 Single-cell occupancy is probabilistic. Empty capsules, multiplets, cell-dependent lysis, leakage, incomplete amplification and barcode collisions are not ruled out by the model. The graph starts at cell suspension allocated to encapsulation; source culture/stool isolation upstream remains outside it. A selected capsule preparation is not a claim that every genomic molecule is present once or that every capsule contains one cell.

The local preparation type distinguishes intact CAPs, paired droplet-plus-CAP preparations and bulk DNA. DTT dissolves the **barcode bead**, not the PEGDA capsule. Droplet breaking recovers still-intact CAPs. Subsequent NaOH treatment dissolves those capsules and releases DNA into bulk; neutralization itself does not physically discard dissolved polymer or reaction species. Possible carryover remains explicit even after unquantified purification.

## Source operations

1. 🟢 Flow-focus the source two-phase cell suspension in fluorinated oil/surfactant, UV-gel the shell (source200W4min), extract with20%PFO/HFE7500 and wash three times with DPBS/0.1%Tween20. Polymer construction and DNA damage/occupancy outcomes are not predicted.
2. 🟢 Lysozyme/lysostaphin treatment37°C2h is followed by three PBS/EDTA washes. SDS/ProteinaseK55°C30min then precedes three source HEPES/Tween/EDTA/NaCl washes and resuspension. No RNA-depletion step or perfect lysis/purification is invented.
3. 🟢 Source Phi29/pyrophosphatase whole-genome amplification30°C2h,65°C10min, then three washes. The source mixture names these enzymes/dNTPs/buffer/surfactants; exact priming/initiation geometry is not assembled. No default random hexamer or circular starting substrate is imported.
4. 🟢 Source ABclonalRK20239 tagmentation55°C10min, then three washes. Adapter polymers/loading and full insertion/end geometry remain unassembled; a known source operation does not identify a canonical modern Tn5 recipe.
5. 🟢 Pair CAPs with DTT-dissolvable barcode beads and source amplification mix in droplets. The bead primers carry a capsule identity barcode and a genome-fragment-binding adapter, but their exact order, width and sequence are unassembled. Replace emulsification oil with QX200 oil; perform the named source linear-amplification cycling72°C5min/98°C3min, then15cycles98/60/72. Barcoded copies coexist with original/untagged/tagmented DNA and reagents. No new random UMI is inferred from the capsule barcode.
6. 🟢 Break the barcoding droplets with PFO/HFE7500, recover/wash CAPs, then mix20µLCAPs with40µL0.8MNaOH at55°C10min. Add40µL0.8Maceticacid. The stock concentration is not presented as the final reaction concentration. AMPure recovery is the actual physical separation after dissolution/neutralization.
7. 🟢 Bulk KAPA HiFi PCR uses source forward/reverse primers for16cycles, followed by magnetic-bead purification and Bioanalyzer QC. Actual primer species, genomic strands and barcode junctions await Table1. The graph does not synthesize example primers from a short analysis motif.
8. 🟢 The source names ONT LigationSequencingKitV14 SQK-LSK114 and NativeBarcodingKit24V14 SQK-NBD114.24, followed by ONT long-read sequencing. Its generic manufacturer URL does not establish exact historical adapter/end/ligation order. That handoff is explicitly unresolved rather than filled from the site's existing ONT protocol.

## Actual source discrepancies

🟢 Original v3p22 prints a total amplification-mix volume of {{= cap_long.PRINTED_MIX_TOTAL_UL }}µL. Its listed component volumes sum to **{{= cap_long.COMPONENT_TOTAL_UL }}µL**, computed from the rendered source list. The model retains both values, selects no normalized mixture and derives no final enzyme, buffer or DTT concentration from this inconsistent total. The caption/preamble carries this limitation on the private preview too.

🟢 The source's sequencing paragraph combines the labels “GridION MK1”, “P2 solo”, “PromethION” and “FLO-MIN114”. Those source labels are retained as unresolved instrument/flow-cell identity; the model claims only the reported ONT long-read route, not a reconciled hardware configuration or run plan. No current vendor naming rule is imposed as a historical correction.

## Identifiers and computational processing

The capsule-group identity is a structured cell-barcode feature with unassigned width/sequence; it does not establish one cell per barcode. The separate ONT native-barcode kit supplies a sample-index class, with member and molecular geometry unassigned. No Illumina i5/i7 index or random per-molecule UMI is invented.

🟢 Computational barcode detection searches the motif {{= cap_long.ANALYSIS_MOTIF }}, a canonical NexteraS7 sequence imported from the shared sequence library. That match does not establish an entire bead primer, tagmentation adapter, PCR primer or orientation. Alignment, barcode grouping, trimming and genome/plasmid assembly are computational fates; they do not physically remove DNA during preparation.

## Remaining work

Actual SupplementaryTable1 and linked material/protocol identities; full MDA priming and tagmentation polymers; capsule/barcode/primer strand geometry; linear-mix reconciliation; historical ONT handoff/instrument identity; source version/publication reconciliation and independent complete source/linked-material review remain open. The graph follows source-known operations through the sequencing handoff but remains notes-only, incomplete and private. No absent-information claim or complete chemical validation is accepted.


## Bounded publication adjudication — 2026-10-30, batch 15

🟢 Microbial cells occupy UV-gelled PEGDA capsules, are lysed, amplified and tagmented before CAPs are paired with DTT-dissolvable barcode beads in new droplets. Bead dissolution releases primers while capsules remain; emulsion break precedes alkaline capsule dissolution and bulk PCR. CAP barcode/cell-origin grouping and ONT native sample barcode are distinct identities, with actual occupancy/whole primer widths unassigned. Current source readout is long-read nanopore DNA, not RNA CapSeq or CXXC capture.

🟡 Original required supplementary_table1, historical_ont and source_recipe remain uncovered. Real version-specific DNA capsule/barcoding/ONT operations are present, but actual cited primer Table1 identity/location is unresolved; historical ONT kit/reagent/end geometry and instrument identity plus mix/MDA/compartment recipe obligations remain unreviewed. Six generic construct-changing steps/non-library endpoint/missing terminal coverage also fail normal gates. Neither ONT kit names nor a sequencing label fulfill these original obligations.

The selected original remains bioRxiv v3, posted13October2025. No later published-version equivalence or present-day publication-status conclusion is asserted. Actual originalp22 prints a348µL barcoding mix, while its listed120+16+16+12+64+64+64+16µL components sum372µL. Both descriptions remain visible with no intended total chosen. The same source combines GridIONMK1 withP2solo, PromethIONFlowCell andFLO-MIN114 labels; instrument/flow-cell identity remains unresolved, without an independent vendor-error claim. Actual SupplementaryTable1 is cited, but its original file identity/location is unresolved rather than proven inaccessible or absent.

Closed inspected scopes: bioRxiv2024.09.10.612220v3 PDF,posted13October2025: Actual version header; Methods711–754 ATPS/capsule formation and782–887 lysis/MDA/tagmentation/coencapsulation/linearbarcoding/alkaline dissolution/PCR/ONT; originalPDFp22 rendered barcoding mixture, oil exchange and DNA recovery; Exactv3 supplementary-material HTML index: Cached title/anchor inventory only, exactrevision links to main PDF; actual Table1 file identity not resolved in this pass.

Remaining real obligations: Documentation-only source/model hold recommended; preserve all three actual required items without semantic narrowing or exclusion. Source originalp22 prints348ul linear-amplification mix; listed120+16+16+12+64+64+64+16ul components sum372ul. Preserve both, do not normalize quantities or assert intended correction. Source says GridIONMK1 withP2solo, PromethIONFlowCell andFLO-MIN114 together. Record exact source labels as unresolved instrument/flow-cell identity, not as an independently established vendor error. Expose actual molecular dNTPs/PEGDA/dextran/LAP/DTT/stabilizer and source enzyme components instead of catalyst/condition-only labels. Loaded adapter+catalyst must come from the same actual tag supply; do not invent free extra adapter or borrowed barcode bead polymers. Actual300ul CAP MDA/200ul postMDA tag/20ul approximately9000cell CAP dissolution/15ul recovered template and25ul KAPA+5ul each primer allocations are source-specific stocks. Unknown operation cardinalities and postPCRQC/loading remainders require faithful compartment/stock modeling. Actual Table1 is explicitly cited but file origin unresolved; no guessed failed URL. Historical SQK-LSK114/SQK-NBD114.24 manual revisions unreviewed, generic manufacturerURL is not authority for a current/default protocol. A source semi-permeability characterization cutoff around400bp is not a universal molecular sieve, perfect capture or onecellonebead guarantee. The lone analysisS7 motif is not a whole orderable sequencing/barcoding primer.

This is a documentation-only adjudication. The source-supported operational/component work is retained, but no module, source-inventory mapping, terminal or builder is changed. All original required source items stay uncovered. No missing preparation is supplied by an opaque endpoint or by kit/platform resemblance, and no information is claimed confirmed absent. No new originals were acquired and no architecture expansion was undertaken.
