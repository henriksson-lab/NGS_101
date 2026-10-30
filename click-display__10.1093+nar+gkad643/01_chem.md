# Click display: corrected linker and mock-library sequencing

Zeng et al. (2023), **Click display: a rapid and efficient in vitro protein display method for directed evolution**. DOI: [10.1093/nar/gkad643](https://doi.org/10.1093/nar/gkad643). Mandatory published correction (2024): [10.1093/nar/gkae1172](https://doi.org/10.1093/nar/gkae1172).

This private incomplete model follows the published mock-library experiment and known downstream NGS route. It does not equate every DARPin/ClfA selection application with that exact run. Source display creates stable covalent protein–cDNA association through a modified two-arm linker; neither reporter chemicals nor amplified DNA copies are mistaken for reported codes.

## Originals and published discrepancies

🟢 Published main PMC10484664: Table1 and Methods134–174. 🟢 Published correction PMC11662653: corrected table/descriptions, independently rendered originalp1. 🟢 Original publishedSI TableS1/FigS3 selected. Three-original cap reached; corrected SI is available but uninspected. Mandatory bundle scraper11 candidates; correction-table full modified strings verified in original and tested with whitespace-only normalization.

The original2023Fig1A/Fig2B/TableS1/Introduction/FigS1 exchange oligo names or azide/DBCO placement. The2024correction expressly states the originalVOR remains unrewritten. The model uses these corrected source strings:

- Oligo1: `/5AzideN/AAAA/i5-TAMK/AAAAAAAAAAAAAAAAAAAAACC/3Puro/`
- Oligo2: `/5DBCON/AAATTCCA/i3Cyan/GCCGCCCCCCGTCCT`

TAMK is the source internal TAMRA label and CnvK is the photoreactive noncanonical nucleotide token, not ordinary A/C/G/T. The modeled source-token arms retain both original DNA bodies. Consumed5′azide/DBCO groups become a nonphosphodiester **5′-to-5′ triazole bridge**; they are not a conventional linear3′-5′ligase junction. Exact complete modification chemistry is not inferred from a base projection, and these scenes remain whole-polymer-unassembled for sourceguard purposes. The independent original-transcription tests preserve all modified strings.

A further visible published discrepancy remains: correctedFig1C labels the template-cutting step “BstI”, while mainTable1 and Methods explicitly use BsaI-HFv2/BsaI-HF. The model follows those explicit recipe instructions (BsaI), retains this caveat on the private preview and leaves the figure-label conflict unresolved. This is a source discrepancy, not a repository transcription correction.

## Operational graph and retained species

1. Corrected Oligos1/2 are mixed at250µM each,10min room temperature in the dark. The source assumes100%for later concentration calculations; the model makes no quantitative conversion claim and retains unreacted arm/reagent classes.10folddilution/storage precedes separate per-round ML aliquots; the same physical aliquot is not consumed three times.
2. Mock-library template contains naiveN3CDARPinDNA plus0.01%D3G86. R0 reference and production allocations remain separate; exact upstream template amplification/construct chemistry is unassembled. The graph follows three source panning rounds and separately prepares R0–R3NGS libraries.
3. Table1:37°C20min transcription;30sUV crosslinks linkerCnvK to associated RNA; BsaI treatment37°C10min follows. DNA digestion fragments remain in the same reaction.37°C30min shaking translation and20min room-temperature peptide transfer follow. The conditional successfully transferred class retains both linker arms, RNA crosslink and puromycin–peptide linkage; unconverted classes are not deleted by the model.
4. EDTA inactivates ribosomes without physically removing them, followed by MgCl2/buffer/dNTP/DTT and SuperScript RT42or50°C30min. Newly synthesized cDNA extends arm2; original RNA and other companions remain until physical workup. The source RT variant is unresolved for the selected mock experiment.
5. IMAC preparation has two printed alternatives (ultrafiltration/PBSexchange or dilution intoPBS+10mMimidazole). Allocation to the actual mock preparation is opaque. Source one-step Ni-NTA purification precedes200nMbiotinylatedeGFP binding, magnetic pull-down and3–8washes; exact per-round beads/wash counts are not invented.
6. Beads are resuspended in StandardTaqbuffer and PCR copied with1938/2643. No invented cleared-elution step intervenes. PCR daughter DNA does not inherit the template's puromycin, TAMRA, CnvK or triazole chemistry. Gel purification physically removes fractions before library aliquots and next-round templates are separated. Third-round and library remainders retain explicit unresolved fates.
7. Q5PCR with source round-specific NGS primers adds Illumina sequence and2nt aa/ac/at/ag round barcodes. Those identify samples/rounds, not molecule UMIs. Each library receives ZymoClean gel recovery;150ng per R0–R3library is physically allocated and combined. Remaining library material is separate.

🟢 Source acquisition is Amplicon-EZ/Azenta MiSeq paired-end sequencing. Read cycles and complete primer landing remain unassembled. FLASH merging, translation and counting are computational analyses and do not change physical DNA. No public finished-library claim is made.

## Ancestry and deferred reading

Main Introduction contrasts mRNA/ribosome/TRAP approaches and explains its new readily obtainable click linker. Such comparison alone does not establish exact inherited recipes or a parent edge to PLATO/CFMS. Direct predecessor components, complete template assembly and corrected supplementary figures remain pending. No confirmed-absence audit has been performed.
