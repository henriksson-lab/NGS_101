# Epi-Decoder: HO/KanMX barcode input and TAP-ChIP libraries

Published defining paper: [PLOS Biology 2018](https://doi.org/10.1371/journal.pbio.2005542).

## Inspected sources and scope

🟢 Published main PMC6059479 Results55/Fig1caption62, Methods170/174/178; actual original Fig1 PDF page3. Published S7 Table s012 original XLSX cells (six library/custom-read orders). Published S1 Figure s001 actual TIFF (the standard source manifest stores this valid TIFF as file.html). Mandatory whole-bundle scraper:40 hits in2files; pointers verified against originals, not accepted as facts themselves.

Selected preparation is pooled yeast HO/KanMX reporter strains, input and IgG-Sepharose TAP-ChIP, separate UP/DN PCRs and their indexed sequencing pool. ADE2 capture, construct engineering, qPCR and alternative sonication times remain outside this graph. This capped three-original pass is not a complete source or absence audit.

## Molecular identities

🟢 Each clone has paired genomic reporter barcodes and an endogenous TAP-tagged protein. BC_UP and BC_DN report clone-linked protein binding at the HO locus; they are neither single-cell identifiers nor deduplication UMIs. The selected HO computational extraction uses20barcode bases and a six-base sample index. Selected actual barcode and index values, genomic alleles and complete occupied primer geometry are unassembled; S1 Table barcode lookup was not read within the cap.

🟢 Whole source orders, computed directly from the module:

{{scene: chem.oligo('2B_UP_Fw')}}
{{scene: chem.oligo('2B_UP_Rv_all')}}
{{scene: chem.oligo('2B_DN_Fw')}}
{{scene: chem.oligo('2B_DN_Rv_all')}}
{{scene: chem.oligo('2B_UP_Seq')}}
{{scene: chem.oligo('2B_DN_Seq')}}

S1 Figure explicitly places P5–seq–U1/D1 forward primers and P7–index–U2/D2 reverse primers on the two genomic reporter regions. Custom read primers are suffixes of the corresponding whole forward order. This does not resolve a selected complete genomic amplicon or authorize importing another yeast-library construct.

## Physical preparation

🟢 Log culture OD660~0.4 is fixed20min with one-tenth culture volume freshly prepared FixSolution, printed as1%formaldehyde/50mMHEPES/100mMNaCl/1mMEDTA. The graph preserves this stock and ratio wording; the source final-concentration convention is unassigned. Glycine125mM final quenches5min. Washed cells are frozen; ~1.5e9 cells are bead-beaten in200µl breaking buffer, washed into FA and sonicated10min at30s intervals, then clarified.

🟢 Allocate100µl chromatin as INPUT and1ml separately for IgG-Sepharose IP, retaining untransferred stock. Beads are washed with PBS containing5mg/mlBSA; IP incubates6h4°C, followed by FA/high-salt/RIPA/TE washes and100µl SDS eluate65°C10min. Input and IP remain separate through RNaseA/ProteinaseK, reverse crosslink and QIAquick purification.

🟢 Each purified fraction supplies independent UP and DN PCR allocations: four assay categories, not invented technical replicates. Exact primer bodies are declared; polymerase identity, PCR formulation and occupied allele geometry remain unassembled. Equimolar contributions are selected from all four PCR stocks, leaving remainders, pooled and100–150bp agarose-selected. Indexed sample-pool contributions enter a sequencing masterpool; a final loading aliquot is separate from the untransferred masterpool.

🟢 Single-end sequencing >50bp on HiSeq2500/MiSeq uses custom UP/DN primers, one or a mixture. No exact integer run length or index-cycle allocation is invented. Computational constant-region search, barcode counting and ChIP/input normalization are not physical DNA processing.

## Open boundaries

The full genomic reporter and selected index/barcode lookup require S1 Table and applicable strain/construct sources. Historical PCR formulation and exact run landing are unassembled, not confirmed absent. No predecessor/vendor recipe is inherited. The public source-qualified page retains these limits; its typed publication gate accounts for reviewed operations and material fates without asserting complete molecular validation.
