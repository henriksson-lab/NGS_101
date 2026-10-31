# HydraPsiSeq: protection from uridine-directed chemical cleavage

Marchand et al., *Nucleic Acids Research* (2020), [doi:10.1093/nar/gkaa769](https://doi.org/10.1093/nar/gkaa769).

This is a capped first source pass. The local graph models hydrazine/aniline treatment, physical recovery and end repair, and stops before the historical commercial library kit. It does not assert a finished NGS library, complete source review or absence of adapter disclosures.

## Source identity and actual reading

- `nar_gkaa769_PMC7641733.xml.txt` and corresponding PDF are the **published article**, not an accepted manuscript: source manifest records author-manuscript false and CC BY-NC. Inspected JATS Abstract 57; Introduction 65–69; Methods 87–107; Results 127–135; Discussion 197–211. Key molecular passage is 129.
- `nar_gkaa769_supp_gkaa769_supplemental_file.docx`: score description and Table S1, extracted text 1–45; Fig. S1 caption 43–45. The original embedded Fig. S1 (`word/media/image4.png`) was visually inspected, including the circled 3′ and 5′ phosphate products. Other supplementary figures were not systematically reviewed; the accidentally opened image11 calibration plot does not establish complete panel review.
- The full source-bundle primer scraper was run before reading and returned zero hits. This does not establish that no oligo table exists or that an extraction contains every visible base. No supplementary sequence table was declared inspected from an empty search.

## Molecular distinction and scope

🟢 HydraPsiSeq reads **protection from cleavage**, rather than a CMC-dependent reverse-transcription stop. Hydrazine attacks susceptible U residues and generates RNA abasic sites; aniline promotes backbone scission. Pseudouridine is comparatively resistant. The published schematic qualifies this as no **or very low** reactivity, so the local selected-U event is not a deterministic prediction for all RNA sites or a guarantee that every Ψ remains intact.

🟢 Consecutive β/δ elimination removes the selected residue: the upstream RNA ends at N−1 with 3′ phosphate and the downstream RNA starts at N+1 with 5′ phosphate. The original Fig. S1 explicitly shows ATP-free T4 PNK removing the former without removing the latter. The local `CleavageEvent` holds those distinct ends; it refuses to reinterpret the selected U-scission event as a Ψ-cleavage prediction. Variable flanks are placeholders, not a fabricated tRNA sequence or fixed cleavage probability.

Molecular fragments and released residue-derived products remain molecular material in the reaction mixture. They do not become physical waste merely because the U-containing RNA has been cleaved. The graph marks waste only at the source's ethanol/phase separations, and leaves the exact partition and quantitative recovery of small products unassigned.

## Published preparation

1. 🟢 Treat 50–300 ng RNA with 50% final hydrazine on ice for 30–60 min. Stop by NaOAc/Glycoblue ethanol precipitation, incubate 30 min at −80 °C, pellet and wash twice with 80% ethanol.
2. 🟢 Resuspend in 1 M aniline, pH 4.5; incubate in the dark for 15 min at 60 °C; ethanol-precipitate again.
3. 🟢 Dephosphorylate RNA 3′ ends with 10 U T4 PNK for 6 h at 37 °C in 100 mM Tris-HCl pH 6.5, 100 mM magnesium acetate and 5 mM β-mercaptoethanol. ATP-free operation and preservation of 5′ phosphate are explicit in Fig. S1. Heat-inactivate 20 min at 65 °C, extract with phenol/chloroform and ethanol-precipitate.
4. 🟢 The article names NEBNext Small RNA Library kit and manufacturer instructions, followed by BioAnalyzer/Qubit QC. It reports multiplexed HiSeq 1000 sequencing with 50-base single reads. The applicable historical kit identity/revision and oligos have not been established in this pass; no current kit is substituted.

## Readout and limitations

🟢 Reads beginning at N+1 report cleavage next to the original U position. Protected sites have a reduced end count, interpreted against neighboring U cleavages and calibration controls. This is a computational relation between mapped read ends and chemical sites, not a claim that the missing U is a deletion in an intact RNA molecule. The graph does not fabricate an acquired run-primer sequence or a sample-sheet index.

The source distinguishes Uscore/PsiScore and documents normalization and calibration (Methods 105–107; supplemental score description). It also identifies false protection from residual structure and 5-modified U residues, plus unexpectedly cleavable non-U modified bases (Discussion 211). The local model does not implement a quantitative scoring algorithm, call every protected position Ψ or estimate occupancy without calibration.

Proposed metadata: RNA; bulk RNA; hydrazine/aniline cleavage; 3′-end repair; ligation-based commercial small-RNA library; multiplexing; HiSeq 1000 single-read 50 cycles; published instance. Retain `notes` and `pending_source_review` pending the downstream kit and read-plan curation.

## Ancestry and deferrals

The paper explicitly compares its protection score with RiboMethSeq (Introduction 65; Discussion 197; supplemental score description), with a different cleavage chemistry and a U-only scoring profile. Candidate `ribomethseq-illumina__10.1093+nar+gkw547` → HydraPsiSeq should therefore be scoped to the **inferred analytical protection-score strategy**, not inherited adapters, substrate chemistry or enzymes. The precise registered predecessor directory and its directly cited versions still need lineage review. Similar endpoint sequencing does not establish molecular derivation.

Deferred: identify the applicable historical NEBNext original and exact oligos; model adapter ligation/RT/PCR and acquired primers; finish remaining supplement and linked-original review; represent quantitative cleavage probabilities and scoring only in a separately justified scope. No missing-information assertion follows from this pass.


## Source-qualified public preparation, 2026-10-30

Two originals: published PMC7641733JATS RNAfragmentation/end-repair/library/sequencingMethods87–107; published gkaa769supplementalDOCX originalembeddedFigureS1image4andactualS1caption. Originalfigure/captionexplicitlysupportsN−1RNA3P/N+1RNA5P and ATP-freePNK3dephosphorylation. No historicalvendor/predecessor originalopened or modernkitrecipe inherited.

The graph retains hydrazine/aniline scission and ATP-freePNK3′dephosphorylation, while exposing actualhydrazine/aniline/Glycoblue/β-mercaptoethanol contributions separately from enzymes/conditions. Released residue-derived products remain molecular before actual ethanol/extraction separations; exact carrier/debris partition is unquantified. PublishedMethods explicitly converts recoveredRNA to a NEBNextSmallRNAkit DNA-library product; that actualfinishing interface is represented with kitcomponentpackaging/wholeoligos/RT/PCR/workups unassembled. The kitcollection is not an inventedpremix and no newer kitrecipe is inherited. ActualBioAnalyzerHighSensitivityDNA/Qubit material is selected as QC, followed by selected librarycontributions to the sourcemultiplexpool, otherlibraryinputwithunassignedidentities/counts, final loadingaliquot and untransferred library/poolstocks. PublishedHiSeq1000SE50/10–12pM is preserved. The local selectedU-event model does not predict Ψcleavage, conversion fractions or modification stoichiometry; calibration/scoring remains outside molecular validation.

Source30–60min hydrazinewindow and literalPNKbuffer concentrations are preserved without silently selecting an optimization condition or correcting100mMMgOAc. FigureS1describes relativeΨresistance and explicitly ATP-free3′phosphatase activity, not an absolute guarantee of noncleavage/protection. No new published discrepancy is asserted.

The normal graph publication gate passes for this source-qualified preparation,
with unknown kit geometry/states computed as unresolved. This does not establish
whole-library sequence, assay performance, purity or complete molecular validation.
No required source item was excluded or covered without evidence. Earlier prefix-only
notes describe history. Mandatory bundle scrapers ran during this capped pass; no
negative sequence-disclosure or complete linked-source audit is claimed.
