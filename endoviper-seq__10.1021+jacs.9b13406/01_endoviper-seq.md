# EndoVIPER-seq

Knutson et al. (2020), *Selective enrichment of A-to-I edited transcripts from cellular RNA using Endonuclease V*. [Defining article](https://doi.org/10.1021/jacs.9b13406).

This capped pass uses the associated author manuscript NIHMS1587417 / PMC7286354. The final published body and SI applicability remain unverified; this is a private partial preparation, not a complete chemistry or missing-information review.

## Inspected source and boundary

🟢 Main Methods96 describes a sequencing preparation distinct from the synthetic oligo binding/cleavage assays76 and gene-specific Sanger/qPCR validation92. EndoV–MBP recognizes inosine-containing RNA in calcium-supported conditions rather than deliberately generating EndoV cleavage ends. Glyoxal modification suppresses RNA structure; magnetic anti-MBP capture enriches RNA fragments before glyoxal removal and RNA-seq. Biological samples and supporting protocol versions are not additional methods.

## Selected preparation

🟢 Begin with 2µg human brain mRNA. Fragment with the NEBNext Magnesium RNA Fragmentation Module for one minute at94°C, then ethanol precipitate. React fragmented RNA in100µl containing50% DMSO and6% glyoxal at50°C for one hour and precipitate again. Quantify the resuspended pellet by NanoDrop. The source allocates500ng to each of two frozen RNA-seq comparison tubes and500ng to each of two EndoVIPER tubes. These are actual parallel aliquots; recovery and any remaining stock are unassigned.

🟢 For each enriched tube, incubate500ng in250µl with100nM eEndoV and120U RNasin Plus in1×EV binding buffer for30minutes at room temperature. Main76 defines the optimized binding buffer as19mM Tris,100mM NaCl,1mM CaCl2,15µM EDTA,150µM DTT,0.025% TritonX-100,30µg/ml BSA,7% glycerol,pH7.4. The MBP-tagged protein is a recognition reagent here, not an inferred nuclease reaction.

🟢 Wash300µl anti-MBP magnetic bead slurry separately with EV wash buffer; combine with the RNA/protein preparation and rotate for two hours at room temperature. Discard the actual unbound supernatant. Wash the beads three times with500µl EV wash buffer. Elute bound RNA in200µl EV elution buffer at95°C for10minutes. Main76 defines the optimized wash buffer as19mM Tris,100mM NaCl,1mM CaCl2,7% glycerol,pH7.4 and elution buffer as0.5M TEAA,pH8.6,47.5% formamide,0.01% SDS.

🟢 Collect the elution supernatant and remove residual beads with the source-described0.22µM microfuge spin filter (unit retained as printed, not repaired). Purify using Monarch RNA Cleanup. To remove remaining glyoxal adducts, incubate in100µl50% DMSO/137mM NaCl/2.7mM KCl/8mM Na2HPO4/2mM KH2PO4,pH7.4 for two hours at65°C, ethanol precipitate and resuspend. Exact adduct occupancy, complete removal and recovery are not established by the graph.

🟢 Assess starting, fragmented comparator and enriched RNAs with Bioanalyzer/6000 RNA Pico. Use8ng of each RNA-seq and EndoVIPER replicate for SMARTer Stranded Total RNA-Seq Kitv2–Pico Input library preparation with standard8-base i5/i7 indexes. Sequence on NextSeq550 with paired150-base reads. Whole kit oligos, RT/PCR intermediates, actual index values and complete read geometry remain unassembled; no predecessor constructor or present-day kit recipe is imported.

🟢 Main100 describes computational HEADCROP3 on read2, attributed to the PicoV2 template-switching oligo, followed by mapping/deduplication. This is read processing, not physical removal of three RNA bases or evidence of a molecular UMI. Main92 gene-specific primers are analytical validation oligos and are not substituted for kit sequencing primers.

## Open scope

The final published article and actual named SI files remain subject to source verification. Synthetic oligo/hairpin assays, glyoxal chemical atom maps, complete RNA/adapter/library molecules, per-step off-products, exact historical kit internals and read landing are outside this capped model. No information is declared absent.

The associated manuscript Results52 describes approximately4ng for library preparation, while detailed Methods96 specifies8ng from each replicate. The model uses the detailed Methods allocation; the actual value and final published wording remain unverified. This is a discrepancy within the inspected associated manuscript, not a claimed final-publication error.
