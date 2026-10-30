# SPACE-seq: associated-manuscript selected spatial lineage preparation

Defining DOI: 10.1016/j.stem.2026.05.017. The catalogue identifies the final publication as **SPACE-seq integrates spatial transcriptomics and lineage tracing in native tissues**. The actual original inspected is PMC13401910.1/NIHMS2190436, an associated author manuscript headed **SPACE-seq Enables High-Resolution Lineage Tracing in Native Tissue Contexts**. This is an explicit version limit: final published body and supplements have not been verified, and matching selected components would not establish complete equivalence.

## Source scope

🟢 Associated manuscript main Methods246–254 and analysis264 were read, alongside the document identity. One original, source investigation closed. The mandatory source-bundle oligo scraper ran, and the printed bead ordering was checked in the main original. The exact TableS1 and SupplementaryFigures files are identified by the manuscript/manifest, but the PMC download routes were gated and those files were not available in its TDM Cloud bucket. This is positive associated-manuscript preparative evidence, not a complete published-source or missing-information audit.

## Preparation boundary and ancestry limits

🟢 The source modifies predecoded Slide-seq bead pucks by ligating DARLIN-specific capture oligos, captures tissue RNA, performs bead-bound RT and second-strand synthesis, and physically divides amplified cDNA between transcriptome libraries and separately amplified spatial-lineage libraries. The selected private graph follows a Tigre spatial-DARLIN library, with other actual allocations retained as unresolved branches. Engineered DARLIN/Cas9-Dntt mouse stocks are source inputs; complete recording cassettes are not invented.

The existing targeted-capture DART-seq modifies Drop-seq beads before droplet processing. SPACE-seq's declared physical boundary is the spatial puck capture and library split, not a new tissue, target panel or analysis software alone. Similar ligation vocabulary does not establish identical strands, ends or recipes. Methods246/250/252 explicitly delegates puck preparation, second-strand details and transcriptome finishing to reference23; the local graph imports neither a DART nor a whole Slide-seq constructor.

## Literal bead ordering and within-manuscript discrepancies

🟢 Methods246 prints:
`5′-TTT_PC_GCCGGTAATACGACTCACTATAGGGCTACACGACGCTCTTCCGATCTJJJJJJJJTCTTCAGCGTTCCCGAGAJJJJJJNNNNNNNNVT30-3′`.

The PC marker is retained as a non-base linker specification with detailed chemistry unassembled. J represents a bead/spatial barcode split by the fixed spacer; both blocks carry the same structured feature identity. The eight printed N positions carry the UMI feature. The V is printed before T30 and is not silently moved or incorporated into the UMI.

🟢 **Within the associated manuscript**, Methods246's literal bead has eight Ns, while its prose calls the UMI seven bases. Methods246 prints eight-plus-six J positions (fourteen total); analysis264 instead prints eight-plus-seven J positions in its filter string. The local inventory preserves Methods246's literal eight-N/eight-plus-six-J ordering. It does not silently substitute a canonical Slide-seq ordering. Final publication/SI version identity remains unverified: these are not labelled errors in the final published article, and complete read spans are not asserted.

🟢 Methods248 names 5′-phosphorylated CS1, CS2 and CS2L plus SplintOligo and T4DNA ligase37°C1h. Their whole sequences are not assigned locally because TableS1 is gated. The ligation-ready geometry is explicitly unavailable, rather than borrowed from DART's IgG splint. No source absence is claimed.

## Source order and physical fates

🟢 Methods250: ten-micrometre section → hybridize on puck in6×SSC/RNase inhibitor, at least15minRT → dip into RT buffer → MaximaHMinus/TSO reaction30minRT then52°C90min → add equal volume tissue-clearing buffer and PK37°C30min → release individual beads by pipetting, wash/centrifuge → Klenow/dN-SMRT second-strand synthesis → centrifuge → Terra/primer mixture → divide into50µl amplification reactions →1.2×SPRIselect → physically split cDNA into library branches.

The capture primers remain bead-associated through RT; they are represented as produced same-tube retained components rather than new free primers. The TSO is a new template-switch acceptor in a combined RT/TS boundary. The clearing buffer concentrations are source stock concentrations before equal-volume addition, not silently assigned final concentrations. ProteinaseK is typed as deproteinization. Actual supernatants/washes are physical fates; unresolved puck residues and parallel reactions are not declared empty or discarded. New PCR/RT bulk mixes contribute same-stock enzyme, primer and buffer components without duplicate additions.

🟢 Methods254, selected Tigre branch: ten ng amplified cDNA →10-cycle P5-CARLIN-F1/TA-RA-Spatial-PCR1-F →1.2×SPRI →P5-PR1/Tigre_Bulk_PCR2_R (cycle count unassigned) →1.2×SPRI →half product transferred →10-cycle P5/Spatial_Lib_Prep_i7_index →0.5×/1.2×double-sided SPRI →analytical Qubit/BioAnalyzer allocation →physical10%PhiX spike →sequence. Whole PCR oligos/products and exact capture/RT/second-strand geometry remain unassembled; names alone do not authorize importing canonical full primer orderings.

## Readout and open source obligations

🟢 Associated manuscript254 specifies spatial-DARLIN R1=50, i7=8, R2=350, MiSeq500-cycle v2 or NextSeq1000/2000 P1 600-cycle kits and10%PhiX. Methods252's transcriptome50/8/50 and Methods240's bulk250/8/250 with5%PhiX are separate branches. They are not merged or repaired by total-cycle arithmetic. All full primer landing, final strands, actual sample-index values and the disputed barcode/UMI layout remain unassigned.

Outstanding originals are the final published body and actual identified TableS1/SupplementaryFigures. Other relevant linked kit/precursor/version details remain unreviewed within the cap; none is claimed absent. This private incomplete graph has no public schematic or finished-status promotion.
