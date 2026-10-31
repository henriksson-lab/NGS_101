# SSiNGLe (Illumina)

Defining DOI: https://doi.org/10.1038/s41467-019-13602-7

Cao et al. (2019), Novel approach reveals genomic landscapes of single-strand DNA breaks with nucleotide resolution in human cells.

## Sources inspected

🟢 Published PMC6925131 main, principle77, Illumina rationale93, common preparation203–215, ILM linear amplification/polyC/library preparation219–223, Illumina read interpretation229. Supplementary Information MOESM1, Note1 beginning (AsiSI validation context, not a new genomic recipe). Other supplementary datasets/figures remain uninspected. Whole-bundle primer hits were checked against main Methods; validation oligos are not substituted for library primers.

## Selected cell-line Illumina route

1. 🟢 Crosslink source K562 cells1%formaldehyde10min; the source prints1.375Mglycine5min for quench (not silently normalized to another concentration), collect/wash/lyse nuclei. Cell-line ILM samples use0.3%SDS37°C1h then1.8%Triton5min and wash; PBMC samples omit that step and remain a separate input branch. Optional Nt.BbvCI/SAP/AsiSI validation treatments are not imposed on native-break samples.
2. 🟢 MNase/RNAIf digestion on ice30min produces a150–500bp preparation. Crucially, MNase fragment ends are3′-phosphate and cannot seed the intended TdT labeling. Accessible native3′OH ends are the eligible class; DSB/repair intermediates are not magically excluded by an SSB label.
3. 🟢 EDTA/SDS stop, ProteinaseK55°C45min then65°C45min, KCl/isopropanol precipitation, ethanol wash and VAHTS cleanup. Reaction products remain until physical workup.
4. 🟢 Denature95°C5min, ice, TdT/dATP37°C30min, then ddNTP37°C30min and70°C10min. PolyA lengths are unassigned. Blocking caps existing accessible3′ ends after tailing; the model does not assume quantitative labeling/blocking. MNase3′P ends remain chemically distinct.
5. 🟢 Clean up; linearly amplify ten cycles with Taq and a chimeric dT50-rT3 primer. The final three residues are ribothymidines, retained as source RNA chemistry, not silently replaced by DNA thymidine or RNA uridine. Free primer's3′rT end is source-described as preventing TdT tailing; an extended product has a newly synthesized DNA3′OH and can receive the next tail.
6. 🟢 Clean up, denature, TdT/dCTP37°C30min, then ddCTP37°C30min and70°C10min. Original ddN-blocked genomic strands and residual free RNA-ended primers remain distinct physical species, not all converted into polyC products.
7. 🟢 PCR with source P5G10/P7T12, the published initial two-cycle annealing program and30further cycles, followed by VAHTS cleanup; paired-end150bp Illumina. Exact complete product/run-primer topology remains unassembled.

Source P5G10: AATGATACGGCGACCACCGAGAtctACACTCTTTCCCTACACGACGCTCTTCCGATCTGGGGGGGGGGHN. Source P7T12: CAAGCAGAAGACGGCATACGAGATcgtgatGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTTTTTTTTTTTTTVN. Source case is retained in this note. HN/VN are degenerate primer-anchor positions; no UMI function is assigned. Index meaning/boundaries are not inferred from lowercase typography alone.

## Readout and boundaries

🟢 The ILM filter expects Read1 startingG10 and Read2 startingT12; the SSB call is the first base after Read2'sT12. Internal-priming, mapping and untailed-control filtering occur after sequencing. The Helicos-specific fill-and-lock one-base correction belongs to SMS, not this Illumina branch. Free3′OH capture is not an exclusive biochemical classifier of all SSB origins.

🟡 Private partial states preserve MNase3′P ends, native3′OH/polyA/ddN history, chimeric linear products and polyC history. Exact primer/template assemblies, RNA/DNA copying details, tail-length distributions and run/index geometry remain deferred; source-positive whole primer sequences are not proof of every molecular transition. Source-associated Helicos preparation/readout is not imported. No finished marker or completed missing-information review is claimed.

## Bounded material corrections; mixed-sugar assembly hold retained

The closed-source implementation adds the actual fragmented-DNA analytical
allocation, 100 ng tailing input and untransferred stock. Source 16.2 µl and
15.2 µl eluates are retained as entire downstream transfers, and the final
20 µl library eluate has distinct Qubit analytical, loading and stock fates.
Fresh dATP/dCTP/ddNTP additions and linear/PCR nucleotide/buffer components
remain physical inputs. The A-tail TdT is carried into the following ddNTP
cap, and the C-tail TdT into the ddCTP cap; neither cap adds an invented fresh
enzyme dose.

The free primer now explicitly preserves 53 thymine nucleotides with separate
sugar tracts: fifty 2′-deoxyriboses followed by three riboses. The terminal
residues are ribothymidines, not RNA U or ordinary DNA T. Its RNA-ended 3′OH
is not chemically erased merely because the source describes free-primer
resistance to TdT. A narrow quantified-source-format exception covers the
prose `oligo-d(T)50-r(T)3`; local checks verify the physical backbone and sugar
tracts rather than flattening the primer to DNA to satisfy the source matcher.

Ten focused checks, private build, source guard and pairing pass. The original
`remaining_source_and_assembly` item remains unchanged and unmet. In particular,
the actual million-nucleus aliquot array and parallel no-MNase HMW preparation,
remaining SI scope and full mixed-sugar copying/tail/anchored-primer/read
assembly are not completed by these contained corrections. Unknown genomic
sequence alone is not the blocker. An opaque final product cannot cover that
actual molecular/source obligation; no completion marker or promotion is added.
No new originals or accepted absence claim follow from this bounded attempt.
