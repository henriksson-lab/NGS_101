# BOLORAMIS — direct RNA-splinted padlocks and nonamer sequencing

[Barcoded oligonucleotides ligated on RNA amplified for multiplexed and parallel in situ analyses](https://doi.org/10.1093/nar/gkab120), published NAR2021. Inspected three published originals: mainPMC8191787, SupplementaryTablesFinal workbook and SupplementaryMaterialsPDF. Source investigation is closed within three originals. The graph models the actual96mRNA+77miRNA multiplex preparation and an exact ACTN1A1 successful-component example, not a synthetic singleplex or an RT-based FISSEQ constructor.

## Source molecules and conditional interfaces

🟢 Original workbookS3C12/A1 gives ACTN1_1932-1957_68_6off_0.9ss:
`/5Phos/AGGTCTGGACAATCTTGGcttcATTAAATCGGTGGCGCGGGCTCTgcaacttcACGTGGT`.

🟢 Its blocks are18nt long arm,4nt copied code,23nt anchor,8nt target code and7nt short arm. The8nt ACTN1 code is `gcaacttc`; the4nt preceding the anchor is its last four bases `cttc`. These are target/gene identifiers, neither molecule-deduplication UMIs nor cell barcodes. The original lower-case distinction is retained in the whole-ordering fixture. Corresponding raw DNA assembly and RCA copy preserve both feature meanings structurally.

🟢 OriginalS3 gives mRNA RCA `AGAGCCCGCGCCACCGATTTAAT`; plus sequencing anchor `ATTAAATCGGTGGCGCGGGCTCT`; minus anchor `/5Phos/ATTAAATCGGTGGCGCGGGCTCT`. The RCA primer's asterisk is part of its source name, not a phosphorothioate token in the actual ordering. The copresent miRNA reagents are also declared exactly: RCA `TCAAGGCGCGGGCAG`; plus anchor `ATAGCGATCTGCCCGCGCCTTGA`; minus anchor `/5Phos/ATAGCGATCTGCCCGCGCCTTGA`. OriginalS3C3/C4/C6–C9/C12 were checked in ZIP/XML, with seven whole ordering fixtures.

🟡 The perfect-matching25nt RNA is derived as the reverse complement of the joined short+long DNA end arms and is **inferred by default**. Its genomic/transcript sequence was not independently verified. A Scene pairs each source DNA end arm to this conditional RNA and marks adjacent5P/3OH ends; the two displayed end tracts belong to one padlock with an off-panel backbone, not two supplied DNA molecules. The shared explicit-end circle primitive preserves the entire source60nt DNA ordering. Its successfully closed subset and shared primer-derived RCA example establish only strand order/primer landing/repeat-copy relations, not ligation specificity, RNA abundance or RCA yield. Three displayed RCA repeats are illustrative, not an experimental copy-number assertion. Aminoallyl-dU substitution positions/occupancy remain unassembled.

## Applied physical preparation

🟢 Main181: cells4%PFA15min/PBS/.25%Triton15min→multiplex10µM padlock pool with RCA primers1:1 molar ratio in10%formamide/6×SSC,1h37°C→three5min hybridization washes/1×SplintR buffer. SplintR210nM(1:50),1hRT closes probes directly on RNA, avoiding reverse transcription. Three hybridization-buffer rinses remove unligated products/enzyme. The retained original RCA primer is not added again. Other956mRNA+77miRNA probes remain physically copresent; their available workbook orderings are unexpanded outside the chosen component, not claimed absent. miRNA sequences excluded from downstream analysis are not retrospectively removed from the physical preparation.

🟢 Actual RCA mix contains1UµlNxGenPhi29,250µMdNTPs,200µgmlBSA,40µM aminoallyldUTP in1×Phi29buffer,2h37°C, then three5min2×SSC washes. Same-tube components are decomposed for inventory bookkeeping, not independent extra additions. Source cells' surrounding protein component remains available for subsequent anchoring; it is not a freshly supplied protein stock.

🟢 Main237: rinse10mMMES6.5; initial172µlwater+4µl500mMMES+4µlBS(PEG)9, begin fixation. Add20µl1Mborate halfway at15min and continue to30minRT. BS(PEG)9 crosslinks incorporated aminoallyl-dUTP to surrounding proteins; full adduct/occupancy remains unassembled. Quench1MTris8 and rinsePBS3×. The optional low-plex fluorescent-FISH readout is not inserted into this selected nonamer-sequencing branch.

🟢 Main239/S8actualPDFp12: each of eight rounds adds sequencing-anchor mix1:100 from100µM stocks in6×SSC/10%formamide,15minRT; wash6×SSC3×1min/T4buffer once. Add200µl mix168µlwater+10µlT4ligase+2µl100µMnonamerstock+20µl10×T4buffer,45–60minRT; wash and image. Between rounds strip with800µlformamide+10µl.01%Triton+190µlwater, preheated80°C,3×3min; rinse6×SSC3× before fresh anchors. Optical imaging does not consume a physical QC aliquot or manufacture an Illumina library.

🟢 ActualS8 shows plus-direction nonamers5P/3′dye and unphosphorylated anchor acceptor, versus minus-direction anchor5P and nonamer3′OH/5′dye. Four source dyes areCy5/TexasRed/Cy3/FAM. Known local anchor binding and terminal polarities are drawn with typed splinted-ligation geometry. Full9nt patterns, cycle offsets, dye/template/reported-base mapping and complete sealed junctions remain unassembled. Source21 polony sequencing10.1126/science.1117389 is an actual cited predecessor, ordinary unread under the cap; its recipe is not inherited. Exact output barcode-read orientation is not inferred solely from raw probe code.

## Limits

Unreacted/off-target classes and residual fractions remain unresolved; no ideal molecular selection/yield is asserted. All concentrations/timing are source-qualified; full occupied cellular molecules, complete nonamer patterns and biochemical fixation products need further directly linked evidence. This private partial graph is not a complete source/chemistry/absence audit and has no public-promotion claim.
