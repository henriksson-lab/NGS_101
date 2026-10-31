# DEPACE-seq

Defining DOI: https://doi.org/10.1093/nar/gkag784.

## Sources inspected

🟢 Published main Methods94–120; SI Adapter preparationp7, genomic control validationpp13–14, rendered TableS1pp17–18. Original cache: `/husky/henriksson/NGS101_papers/sources/depace-seq__10.1093+nar+gkag784/nar_gkag784_PMC13443447.xml.txt`. Mandatory whole-bundle oligo scraper run; source hits checked against original text/table, not treated as independent evidence.

## Molecular scope

Selected published main genomic route. PAB conjugates AP sites; cleavage does not remove fragments from the tube. SI control experiments use NEB/NextSeq conditions and do not replace main VAHTS/HiSeq-or-DNBSEQ recipe. Exact left/right release molecular identity and adduct-bearing end geometry require further assembly; no ideal molecular partition assumed.

1. 🟢 Sonicate genomicDNA200–500bp. Molecular representation: fragmented genomicDNA with AP sites at unassigned positions.

2. 🟢 PAB2mM4xPBS overnight room temperature. Molecular representation: PAB-AP covalent conjugates on DNA plus unlabelled/background DNA.

3. 🟢 Source Zymo-Spin cleanup after PAB. Molecular representation: recovered PAB-conjugated and background DNA.

4. 🟢 Main VAHTS ND607 end repair/dA-tail and customP5 ligation; bead cleanup. Molecular representation: P5-associated PAB-labelled DNA; exact kit/end junction assembly unreviewed.

5. 🟢 EndoIV37C30min at source PAB-labelled sites. Molecular representation: source free3′OH cleavage class and corresponding PAB-associated end fragments, complementary DNA and possible uncleaved DNA; all retained.

6. 🟢 Source Zymo-Spin cleanup after EndoIV. Molecular representation: recovered cleaved/uncleaved DNA population; fragment recovery unquantified.

7. 🟢 Prewashed streptavidin beads, poly dI:dC, binding and six source washes. Molecular representation: captured biotin-associated DNA complexes and possible blocker/background carryover.

## Readout and remaining work

Main HiSeq2000 or DNBSEQ-T7 PE150. Left/right libraries are physically separate; their mapped coordinate intersection is computational. One-library putative AP calls are a separate analysis branch.

Exact modified P5/P7/rhP5 assembly, exact molecular partition through sequential L/R release and full read geometry remain model work; vendor ND607 manual directly named but uninspected.

This capped pass is incomplete and makes no accepted missing-information claim. Local preview has no finished-publication marker.

## Published modified adapter inventory

🟢 SI TableS1 rendered pp17–18; plain-text layout was checked against the original PDF. All written5′→3′:

| Oligo | Sequence and source modifications |
|---|---|
| P5 adapter top | 5′NH2C6-GACTGGAGTTCAGACGTGTGCTCTTCCGATCT |
| P5 adapter bottom | 5′P-GATCGGAAGAGCACACGTCTGAACTCCAGTC-3′C3 Spacer |
| P7 adapter top | 5′P-AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCG-3′C3 Spacer |
| P7 adapter bottom | 5′C3 Spacer-CTACACGACGCTCTTCCGATCTNNNNN*N-3′NH2C6 |
| rhP5 | 5′NH2C6-GACTGGAGTTCAGACGTGTGCTCTTCC/rG/ATCT-3′C3 Spacer |

The six random P7-bottom bases are not labelled UMI by this model. Asterisk is retained as a source backbone modification; C3 spacers add no nucleotide. rhP5 includes the actual riboguanosine, not a DNA G substitution. These inventories are unassembled until exact end geometry is validated. P5-bot1 is separately printed in the table and is not silently assigned as the main genomic adapter bottom.

## Dual sequential elution and library boundary

🟢 Main Methods102–108: L-library uses100mMNaOH at room temperature5min; supernatant collected, neutralized withHCl, column purified;95°C5min denaturation, customP7/Instant Sticky-End ligase25°Covernight, bead cleanup and VAHTS/i5/i7 PCR15cycles. The **same residual bead pellet after L elution**, washed once, yields R-library using95%formamide/10mMEDTA90°C3min and column purification. R-template is extended using5′amine/3′C3 blocked rhP5, thermostableRNaseH2 and Q5 before P7ligation/PCR. RNaseH2 activation of the chimeric blocked primer must be modelled explicitly before asserting a sealed R-library; the graph includes the sequential releases, RNaseH2 activation/extension operation and both ligation/PCR/cleanup paths, while keeping the exact cut/end products unresolved.

Main analysis114–120 assigns complementary strand-specific AP coordinates from Read1 and computationally intersects L/R candidates; this intersection does not physically combine the libraries. Main also defines a separate single-L-library putative-calls path. No exact cleavage-side adduct geometry or perfect strand release is assumed by this capped implementation.

SI pp13–14 control validations use NEB kits, different capture volumes/time and NextSeq; those are named control experiments rather than replacement authority for the main genomic branch. The bounded pass inspected these distinctions, adapter preparationp7 and the modified TableS1; other SI figures/analyses remain unreviewed.

## Bounded publication follow-up (2026-10-30)

The required `pending_assembly` item remains uncovered. Its exact modified P5/P7/rhP5 assembly, cleavage-side adduct/strand identity, L/R molecular partition and full read geometry are available-source modeling work. The directly named VAHTSND607v3 manufacturer instructions remain a relevant unread dependency, with exact document revision/URL unassigned, not a fabricated unavailable disclosure. No endpoint is retyped to bypass these obligations.

The limited private fixes preserve source-known molecules and modifications: rhP5 now contains the actual G at the `/rG/` position, with one RNA-sugar block between the DNA regions and structured5′C6-amino/3′C3 annotations. Its sequence derives from the matching canonical read-primer body; a whole modified-ordering fixture checks the nucleotide string and ribo-G coordinate against published Main106. P7 bottom retains its phosphorothioate linkage between the final two degenerate N residues as structured chemistry, separately from its terminal C3/amino modifications. Exact occupied adapter/primer junctions remain unassembled; N6 is not assigned UMI meaning.

PAB is a molecular conjugation reagent distinct from its incubation conditions, poly dI:dC is an actual DNA blocker contribution, and R-extension dNTPs are a molecular substrate distinct from RNaseH2/Q5. Known adapters and PCR primer classes are typed as adapters/oligos rather than hidden as general molecular material. L-PCR functional polymerase/nucleotide components derive from its actual source25µl VAHTSHiFi stock. Main106 explicitly reuses the L cycling program for R; it does not independently specify the R formulation here, so that formulation/dose is left unassigned rather than assuming the full L mixture.

A literal reviewed subset records the selected genomic route and these verified inputs; it does not mark the broad full-library assembly item fulfilled. Main106's residual-pellet R release follows the L elution. Source-control NEB/NextSeq settings are not substituted for main VAHTS/HiSeq-or-DNBSEQ PE150. This pass uses the two already inspected published originals only, with no acquisition or source expansion. Ten focused checks and the existing private builder pass; notes status remains a model-scope hold.
