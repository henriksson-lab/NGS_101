# PAT-seq — capped source-qualified notes

Harrison et al. (2015), PAT-seq: a method to study the integration of 3′-UTR dynamics with gene expression in the eukaryotic transcriptome. DOI10.1261/rna.048355.114.

🟢 This is an NGS adaptation of the ePAT end-extension strategy, not merely a tail-length analysis applied to an ordinary RNA-seq library. The source directly credits Janicke2012 for the Klenow-mediated tagging component. That predecessor original remains unreviewed, so no constructor or complete recipe is inherited.

🟢 Read sources: published main article PMC4509939, methodology/Methods; published `SuppMethods.docx`, steps1–7/reagent list. Mandatory whole-bundle primer scraper completed; adopted oligos were checked in those originals. Figures supplement and expression workbook remain unreviewed. The main article calls the laboratory-ready protocol “Supplemental file3”; the actual downloaded `SuppFile3_DGE_wt_vs_dccr4.xlsx` is expression data, while `SuppMethods.docx` is the workflow. Filename alone did not establish identity.

Published oligos, 5′→3′ unless explicitly indicated:

- Bio-EE: `[Bio]CAGACGTGTGCTCTTCCGATCTTTTTTTTTTTTTTTTTT`.
- SplintA: `CCCTACACGACGCTCTTCCG(rA)(rT)(rC)(rT)`. Both originals explicitly call this terminal block RNA and write rT. The model preserves that notation, without silently substituting rU; the intended ribonucleoside identity is an open question.
- SplintB: `NNNNAGATCGGAAGAGCGTCGTGTAGGG`; the annealing illustration also prints it antiparallel as `3′-GGGATGTGCTGCGAGAAGGCTAGANNNN-5′`. NNNN is a degenerate splint window, without a source-established identifier function. The local Scene constructor verifies that the full known B arm is the reverse complement of A and computes their alignment; this checks base pairing, without resolving the source rT chemical identity.
- Universal: `AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCG`.
- ScriptSeq indexed reverse primer: named product only in the inspected scope; exact historical selection/landing remains unassembled, not confirmed undisclosed.

🟢 Selected detailed workflow: 1µg total RNA → Klenow exo-minus end extension using the separate biotinylated DNA as a template → limited RNaseT1 → immediate phenol/chloroform quench/extraction → streptavidin capture of DNA:RNA duplexes → washes → T4PNK phosphorylation/wash → preannealed mixed-sugar splintA/B ligation by T4RNL2 → wash → SSIII RT primed by the same bead-associated Bio-EE DNA → wash → formamide release/pulse-spin of beads → 6% urea-PAGE, yeast120–300base window → crush/soak/filter/precipitate →30µl →10µl into16cycle AmpliTaqGold360/index PCR. Remaining20µl is explicitly unresolved.

🟢 Ten percent of amplified product is used for an agarose intensity assay; the retained production material is pooled approximately stoichiometrically with other libraries. Facility bead cleanup and quantitation follow that pool. These analytical and production fractions remain separate in the graph. The biotin begins on a separate DNA template and subsequently on the RT primer/cDNA; it is not asserted to become a covalent tag on RNA during end extension. No complete RNA/DNA junction geometry, exact ligation efficiency or perfect capture is predicted.

Published differences retained: main Methods specify RNaseT1 for1min, while detailed workflow step2 uses2min; the selected graph follows the detailed workflow. Main Methods say ethanol precipitation after PAGE, while the detailed recovery specifies ammonium acetate/isopropanol with an ethanol wash; graph uses that detailed branch. Workflow overview says150bp sequencing, while the main actual experiment states HiSeq1500 rapid100base single-end; actual experiment is the readout authority. Supplement step4 labels an8µl assembly, while the displayed2+2+3µl list totals7µl before beads; full volume composition is unassembled. These are source differences, not silently reconciled values.

🟡 Model: private, partial protocol-owned graph. Chemical mixtures, nonproductive species and actual waste are retained to the next physical separation. No exact index/read landing, mixed-rT resolution or complete library sequence is claimed. Software poly(A) masking/alignment is downstream of acquisition, not physical depletion. No public promotion or accepted absence claim.

Open source scope:

- Janicke2012 ePAT predecessor original, cited by main Methods; exact DOI/version to verify before use.
- ProtocolExchange1889 linked by supplementary PAGE step6; same-scope original not followed within cap.
- Historical ScriptSeq indexed primer product/revision and selected indices; source selection not yet assigned.
