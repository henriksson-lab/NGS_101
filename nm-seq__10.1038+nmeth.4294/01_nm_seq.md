# Nm-seq

Defining original: https://doi.org/10.1038/nmeth.4294

Dai et al. (2017), *Nm-seq maps 2′-O-methylation sites in human mRNA with base precision*.

🟢 Publisher Supplementary Protocol1 original `nmeth.4294_41592_2017_BFnmeth4294_MOESM3_ESM.pdf`, steps1–19, printed DOI verified. Available `nmeth.4294_PMC5712428.xml.txt` is explicitly an accepted author manuscript (NIHMS922173), not the VOR: chemistry overview71–77, Methods117–133 and Fig1 caption341 were inspected. Published article procedure remains to be reconciled before promotion. Full-bundle scraper run; model-oligo/fluorescent adapter hits checked in manuscript Methods117–121, not transferred into the library recipe.

The discovery ledger previously pointed to nmeth0318-226b as a purported Nm-seq correction. The actual publisher original identifies an unrelated protein-crosslink-structure paper, DOI10.1038/nmeth.3838. It is excluded, and establishes no Nm-seq correction. This is a discovery identity error, not a published chemistry discrepancy.

## Terminal chemistry

🟢 Mechanistic description: periodate oxidizes a terminal2′,3′vicinal diol to dialdehyde; beta elimination removes the terminal nucleoside, leaving the shortened RNA3′monophosphate. Dephosphorylation restores3′OH for another cycle. Terminal2′OMe lacks the required vicinal diol and resists the selected pathway. These are conditional molecular relations; no cycle efficiency, quantitative Nm fraction or absence of other products is inferred.

🟢 Publisher protocol steps7–13 performs eight OED cycles, with rSAP directly after quench and a column purification each cycle. Steps14–15 uses final OE without phosphatase. Thus susceptible RNA becomes3′phosphate and incompatible with3′adapter ligation, whereas terminalNm retains3′OH. Released products stay in molecular material until actual column cleanup. The local constructor enforces reaction order, modified-site positions, terminal chemistry and residue conservation before separation.

🟢 Steps16–19: T4PNK3′phosphatase-minus plus ATP,37°C30min, column cleanup, NEBNextSmallRNAE7330S kit with overnight3′ligation16°C. The manuscript Methods125 says60min for PNK. This is a source-version comparison pending final publication reconciliation; no within-publication contradiction has been established. The module uses the publisher supplementary protocol's30min condition. Within that original, step16 is printed as a “5’ dephosphorylation reaction” beneath the5′PHOSPHORYLATION heading, with ATP and T4PNK3′phosphatase-minus. The model follows the heading/reagents as5′phosphorylation; the contradictory step label is recorded here rather than silently adopted as a different reaction.

## Oligos, readout and boundary

🟢 Manuscript model assay: `ACGCGCACCCCCU` / terminal `Um`; fluorescent ligation assay adapter `5′-rAppAGATCGGAAGAGCACACGTCT-Alexa488-3′`. These are calibration/model assay reagents, not established as the preparative library adapters. They are deliberately not substituted for the E7330S library.

The source reports paired-end HiSeq2500 and an asymmetric read pileup whose RNA3′end corresponds to Nm (manuscript Methods125–129/Fig1). The model does not invent a fully acquired read plan from trimming strings. It covers eight OED cycles, finalOE, physical cleanups and5′phosphorylation, with an explicit downstream boundary at post-PNK cleanup/historical kit/adapters/PCR/read geometry. Notes status. Publisher SI figures, workbook, cited ProtocolExchange DOI10.1038/protex.2017.088 and applicable historical kit guide remain deferred; no information is claimed absent.

No explicit catalogue-parent lineage established. Manuscript Introduction cites periodate/nucleoside chemistry predecessors and compares2OMe-seq/RiboMeth-seq, but conceptual comparison is not derivation of their library components. Do not infer family inheritance from shared small-RNA adapters.


## Source-qualified public preparation, 2026-10-30

Two originals: associated NIHMS922173manuscript nmeth.4294_PMC5712428.xml.txt chemistry71–77/Methods117–133; published SupplementaryProtocol1 MOESM3 fullSteps1–19, actualprintedPDFpp2–3 tables/Step16wording inspected. Publishedsupplementcontrols selectedtreatedpreparation; associatedmanuscriptprovidesqualifiedmechanism/inputcontrol/HiSeqcontext, notcompleteAM/VOR procedural equivalence. OtherSI/ProtocolExchange/correction pointers notopened/adopted thispass.

The graph now starts with source10µgDNase-treated totalRNA, fragmentation and Antarctic3′repair/columns, then eight OED rounds and a ninth finalOEwithout dephosphorylation. Actual periodate/lysine/ethylene-glycol supplies are molecular inputs, rSAP and5′PNK are catalysts, ATP is separate from PNK. Released residue-derived products stay copresent until actual column separation, without a quantitative debris/recovery claim. SourceStep18post-PNKcleanup elutes8µL before the actualStep193′adapterligation16°Covernight and sourcekitfinishing. The source-positive kit continuation constructs DNA-library material; sequencing is not attached to the RNAintermediate. Historical adapters, downstream dual-ligation/RT/PCR/workups, exact kitstockcount and wholeDNA remain unassembled. Synthetic fluorescent control-adapter sequence is not borrowed into the sequencinglibrary. RNA3′OH/3′Pclasses are preserved through finalPNK/columnmetadata. A separately reported associated-manuscript untreated-input library branch is retained as an opaque sourcekit/controlloading interface; its earlier end handling is unknown and does not borrow treated OED/PNKconditions. Selected library-to-run contributions and untransferred stocks are explicit; actualtreated/controlrun grouping and rawlength/indexmask remain unknown.

PublishedProtocolStep16literally says5′dephosphorylation, while its sectionheading says5′PHOSPHORYLATION and the table suppliesATP plusT4PNK3′phosphataseminus. The model uses sourceheading/ATP-supported5′kinase direction and preserves the printed wording difference. Selectedpublishedprotocoluses20mMperiodate and30minPNK. The associatedmanuscript gives10mMperiodate in the synthetic-model context and60minPNK in Nm-seq; those values are not merged into this selectedpublishedsupplementary route or called a confirmed final-paper error. CompleteAM/VORreconciliation remains open.

The normal graph publication gate passes for this source-qualified preparation,
with unknown kit geometry/states computed as unresolved. This does not establish
whole-library sequence, assay performance, purity or complete molecular validation.
No required source item was excluded or covered without evidence. Earlier prefix-only
notes describe history. Mandatory bundle scrapers ran during this capped pass; no
negative sequence-disclosure or complete linked-source audit is claimed.
