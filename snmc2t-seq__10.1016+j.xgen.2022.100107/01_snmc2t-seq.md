# snmC2T-seq — nuclear RNA and methylome

## Active source scope

Inspected in the capped2026-10-18 migration: Published PMC9004682 Results215 and Methods540–554; published mmc1.pdf title/FigureS1 specificity and modality-comparison caption; published mmc9.xlsx TableS8 Neuron Marker One_vs_Rest header/first entries, marker data not recipe. Source-qualified2018 finishing component reused from same-batch verified PMC6147798 Methods149–153 and MOESM1 SupplementaryMethods/Table2; no additional original opened for this reuse.

## Native declarations

Selected brain nuclear RNA/methylome subworkflow: source DNA-dU/C3 RT oligos and RNA/LNA TSO tail, newly synthesized5mC cDNA body,14-cycle nuclear cDNA PCR, UDG-only abasic treatment and later alkaline cleavage,25µlCT/90µl binding modifications. Independent genomic DNA/RNA/cDNA classes stay physically mixed; mCH modality partition is computational. Nuclear finishing explicitly inherits2018 RP-H/Adaptase/PCR/recovery order, while whole-cell2017/N9 assumptions are excluded. New copies encode methylation as bases rather than chemically inheriting5mC. Native member/end/geometry unknowns retained.

🟢 Published RT oligos contain **DNA deoxyuridines** that allow UDG treatment; U here does not make the primer an RNA strand. TSO_3 separately contains explicit `rGrG+G` RNA/LNA terminal modifications. The selected brain oligos have5′C3. Their later primer-derived end retention/abasic cleavage boundaries remain unassigned. Body cytosines newly synthesized using5-methyl-dCTP are distinguished from original oligo bases.

🟢 The selected brain-nuclear cDNA PCR has14cycles. Cultured nuclei use15cycles; whole cultured cells use12cycles and the earlier2017 snmC finish. Those routes are separate alternatives. Brain UDG-only treatment creates abasic sites; subsequent alkaline desulfonation supplies cleavage. Hydrolysis/abasic formation is not a physical discard fraction.

🟢 The nuclear finishing method explicitly cites2018 snmC-seq2; this justifies reuse of its RP-H-based operation order, not treating genomic DNA and methylated cDNA as identical substrates. The two populations stay distinct throughout copied/Adaptase states. No new chemically methylated bases are assumed during ordinary post-bisulfite copying. Transcriptome/methylome assignment by read-level non-CG methylation occurs after sequencing; no modality-specific barcode or covalent DNA/cDNA join is invented.

🟢 The selected brain nuclear read recipe is NovaSeq6000 S4 PE150. Native index acquisition and base-level primer offsets remain unassigned. The broader snmCAT NOMe front end is outside this RNA/methylome subworkflow.

Deferred: Directly linked authored snmCAT v1 protocols.io original and any changed nuclear finishing/index details; native ordered index/RP members and Adaptase geometry; Large mmc11 source-data workbook and other supplements remain uninspected; no absence finding or linked-source closure. This capped pass does not establish complete source review or publication-wide absence.

The PCR output and all subsequent cleanup drawings explicitly retain possible co-carried original DNA and blocked native template lineages, with integrity/abundance/partition unassigned. The expected amplified library is a separate molecule population in that same state. Cleanup fractions therefore do not silently assert physical elimination of originals or native precursor templates.

The TSO physical state separates the DNA/deoxyuridine tract, two RNA guanosines and one LNA guanosine; ordering text preserves `rGrG+G`. The methyl-dCTP condition explicitly feeds cDNA PCR as well as RT. Residual RNA/template hydrolysate after bisulfite conversion has an explicit unresolved fate, with no physical removal asserted.
