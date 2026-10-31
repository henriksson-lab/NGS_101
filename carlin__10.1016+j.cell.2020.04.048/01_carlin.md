# CARLIN — separate genomic, bulk-RNA and single-cell recovery

Bowling et al. (2020), *An Engineered CRISPR-Cas9 Mouse Line for Simultaneous Readout of Lineage Histories and Gene Expression Profiles in Single Cells*, [10.1016/j.cell.2020.04.048](https://doi.org/10.1016/j.cell.2020.04.048).

## Actual source/version scope

One original was read: **accepted author manuscript PMC7529102**, array design/assembly (193–203), **CARLIN Amplification and Bulk Sequencing Protocols** (232–236), and **Single-Cell RNA Sequencing Protocols** (238–244). The manifest marks it an author manuscript/TDM source. A standard attempt at final publisher HTML returned HTTP 403. The actual supplementary originals, including primer/construct Table S1, are behind PMC download gates. The standard bundle scraper finds many primers in a fetched **superseded preprint supplement**; those are excluded and do not fill the missing final-era inventory. Source versions remain distinct; final article/supplement applicability is unresolved, not a publication-error claim.

This capped pass models conditional accepted-manuscript operational routes, not a complete final published chemistry. Recorder construction/editing is described as context before the modeled prepared-DNA/RNA/cDNA inputs. Five displayed variants remain distinct: single-cell v2, single-cell v3, bulk genomic DNA, bulk RNA in vitro and bulk RNA in vivo. They are one method's source variants, not extra protocol-count rows or pooled experiments.

## Recorder and identifier context

🟢 The source CARLIN recorder has ten sense-oriented CRISPR/Cas9 target sites, with source 20-base targets, 3-base PAMs and 4-base linkers. Six guide designs were preserved from the previously tested GESTALT v6/v7 set and four redesigned for mouse compatibility. The recorder is placed in a GFP 3′ UTR before bGH polyA under CAG control. Exact full array and edited alleles are not reconstructed from these dimensions; no fixed edited length or particular indel is predicted.

🟢 The source distinguishes constitutive and inducible guide multiplexes and Dox-inducible Cas9 recording. Most experiments use the constitutive guide system. Induction and repair create inherited lineage alleles before preparation; library PCR does not create those historical edits. A structured lineage barcode is distinct from the Chromium cell barcode/UMI, bulk recovery UMI and sample/library index. Alleles can be shared by descendants and collide between independent lineages.

## Bulk genomic DNA

🟢 Up to 250 ng gDNA enters source UMI tagging with **NGS_UIM_D_F** (source spelling), containing the source CARLIN_fwd2/adaptor region and a 10-base random UMI. It is a **single extension**, ramping annealing/extension five times without another denaturation step, not five independent PCR UMI assignments. After AMPure purification the entire preparation undergoes 35-cycle NGS2_F/NGS1_R PCR, then 10-cycle P5/NGS2R_I# indexing. Exact extension/polymerase/primer/end and final workup details remain unassembled. Cleanup ratios described for the RNA branch are not silently copied to gDNA.

## Bulk RNA variants

🟢 A gene-specific RT_Bulk_CARLIN primer introduces a 10-base UMI during SuperScript III RT from up to 1 µg total RNA.

🟢 In vitro: 1× AMPure cDNA cleanup, 35-cycle NGS2_F/NGS1_Bulk_R PCR with Platinum Taq, 0.8× AMPure cleanup, 10-cycle indexing, then final 0.8× AMPure purification.

🟢 In vivo: purified cDNA enters 15-cycle NGS1_F/NGS1_Bulk_R Q5 PCR, followed by 0.8× AMPure cleanup. **Half** enters 15-cycle NGS2_F/NGS1_Bulk_R Q5 PCR, followed by 0.8× cleanup, 9-cycle indexing and final 0.8× cleanup. The first cDNA-cleanup ratio is not assigned from the specifically in-vitro clause. The untransferred half is unresolved, not presumed discarded.

🟢 The source bulk section reports MiSeq 250/6/250 cycles with 10% PhiX. Exact final molecular read/index/primer binding geometry remains unassembled; the six-base index read is not a UMI.

## Single-cell routes

🟢 Source Chromium 3′ v2 (Figure 5) and v3 (Figure 6) generate amplified cDNA. The graph starts at that source Step 2.3 preparation; actual linked historical capture/RT/library manuals remain unread. Independent aliquots feed whole-transcriptome library continuation and targeted CARLIN enrichment. This is allocation of whole-cDNA mixtures, not selective purification before PCR. Transcriptome acquisition is NextSeq 500 28/8/91 cycles.

🟢 Targeted recovery takes 3–5 ng amplified cDNA into 15-cycle **10X-CARLIN_1-bio/P5-PR1 PCR**, followed by 0.8× SPRIselect. The physical biotin-bearing primer is represented as such, but its attachment position is unassigned pending actual Table S1. Incorporation into its source-derived priming strand does not imply copied biotin in every complementary/new strand.

🟢 Dynabeads kilobaseBINDER captures biotin-tagged products for 3 h at room temperature, followed by actual magnet/wash separation. **Half** the bound beads enter 15-cycle 10X-CARLIN_2/P5-PR1 PCR. Remove Dynabeads by magnet, retain library-containing supernatant and perform 0.8× SPRIselect cleanup. **One tenth** of that product enters eight-cycle SI-primer/Chromium-i7 indexing. Remaining bound beads and unused PCR products have explicit unresolved fates. Source template/primer/capture heterogeneity and carryover are not replaced with perfect purity assertions.

🟢 The source uses polymerase supplied with the **10X v3 kit for all 10X PCRs**, including v2-derived capture material. Targeted acquisition is MiSeq 28/8/350 cycles, with 10% PhiX. Exact cell/UMI subdivisions, sample-index members and source-specific priming geometry are unassembled; current kit specifications do not fill them.

## Component ancestry and open sources

The accepted-manuscript design Methods explicitly starts from GESTALT guide designs and preserves six; the gDNA UMI-extension paragraph cites McKenna et al. (2016). These justify narrowly qualified component-ancestry proposals—guide-design subset and single-extension tagging strategy—not identical arrays, oligos, enzymes, adapters or complete recipes. No source-selected constructor is imported from GESTALT, scGESTALT or generic Chromium code.

Deferred: final published main, actual Table S1/construct/modified oligos, other blocked final-era supplements, source-version applicability, full physical array/strand/modification/end/capture/read geometry, directly cited historical Chromium v2/v3 guide identities and other linked chemistry sources, and independent review. Availability failures do not establish absent disclosure. All previews remain private.

## Bounded publication follow-up (2026-10-30)

The associated accepted manuscript's Methods233–243 supports distinct genomic-DNA, bulk-RNA and single-cell CARLIN recovery procedures, but explicitly refers all primer sequences to TableS1. That actual required original is NIHMS1594358-supplement-8.xlsx, not supplement-7.docx. The required `actual_source_table` item remains uncovered; retyping a generic molecular endpoint as a library cannot substitute for reviewing that original or assembling its actual primer/product scope. The final published main also remains unreviewed. No whole accepted/final equivalence is asserted and no superseded preprint table is adopted.

The source-specific distinctions remain important: genomic DNA uses five anneal/extend ramps without repeated denaturation; bulk in-vitro and in-vivo RNA PCR routes differ; single-cell recovery uses half the bead preparation for nested PCR and one tenth of the resulting product for indexing (Methods243). Biotin belongs to the actual source-modified primer-derived population, not automatically every copied strand. Historical ChromiumV2/V3 guide export identities remain unverified from the cited shortlinks. Source10%PhiX is a physical DNA contribution if a later complete loading trace is modeled, not merely a computational run label.

This pass reread one existing associated-manuscript original, Methods233–243, and ran the mandatory cached-bundle scraper. Final-era TableS1 and the published main were not newly acquired or opened; historical failures are not current access claims. No module/build changes. Full input/primer/read geometry and exact source obligations remain explicit; notes status is a source-scope hold, not an absence finding.
