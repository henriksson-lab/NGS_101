# Nano-tRNAseq — source-qualified partial preparation

[Defining paper](https://doi.org/10.1038/s41587-023-01743-6), *Quantitative analysis of tRNA abundance and modifications by nanopore RNA sequencing*.

Two published originals were inspected: mainPMC10791586 (publishedCCBY, not AM), actualFig1PDFp3 and Methods202/206; actualTableS25 inMOESM3. One source slot remains unused. No modernRNA004 constructor or failed strategy is inherited.

🟢 Both mature-tRNA ends receive preannealed RNA padding. The selected3′RNA:DNA adapter has{{= chem.sugar_count('rna3', 'RNA') }}RNA residues followed by **{{= chem.sugar_count('rna3', 'DNA') }}DNA A residues**; tokenize the rawrA-prefixed tail rather than counting printed terminalAletters in isolation. The5′RNAadapter is{{= len(chem.bases('rna5')) }}nt and endsrN; it is a complementary matureCCA/discriminator binding class, not a UMI. All-RNA3′padding and polyadenylation/DNA-only strategies are separate failed/optimization branches. SourceFig1b distinguishes nativeRNA transport from protective cDNA.

🟢 Deacylate biologicalyeasttRNA in10µlwater+95µl100mMTris9/37°C30min; sourceZymoR1016 retention≥17nt with1.3×ethanol, TapeStationRNAHSQC. Anneal5′/3′padding1:1 in10mMTris7.5/50mMNaCl/1µlRNasin,50ngµl,75°C15s then.1°Cs to25°C. Prepare sourceRTAoligoA/B similarly. Add selected preannealed padding to tRNA at1.2:1 (average90nt assumption source-specific),50µl20%PEG/Rnl2buffer/4µl6mgml in-houseRnl2/1µlRNasin;RT2h. RNAclean2×/15min/70%ethanol/elute10min;RNAHSQubit.

🟢 Select200ng paddedRNA plusRTA1:2 (source~4.3pmolRNA/8.6pmolRTA),15µlQuickLigationbuffer/1.5µlT4DNAligase/.5µlRNasin,RT30min. Add25µl actualRTmaster (13µlwater/2µl10mMdNTP/8µl5×Maximabuffer/2µlMaximaHMinus) directly,60°C1h/85°C5min/4°C. RNAclean2×. RMXligation40µlQuickLigationbuffer/3µlT4DNAligase/6µlRMX,RT30min;RNAclean2×10min/two150µlWSBwash/elute20µlELB10min. Add17.5µlwater+37.5µlRRB;keepice. FLO-MIN106QC/priming/loading follow sourceSQK-RNA002. Actual stock selections/remainders and opaque RMX/motor/flowcell boundaries remain explicit; no sequence body or active enzyme dose borrowed from a newer kit.

🟢 TableS25 provides full5′RNA,3′RNA:DNA,RTA-A andRTA-B orders; the all-RNA3′adapter and3′RTtestadapter are comparators. All mixed-sugar residues and5P states remain literal. MaximaRT protects/linearizesRNA; reads remain nativeRNA, not cDNA sequencing.

Main250 describes custom1secondadapter/2secondstrand thresholds; ExtendedData4caption300 describes2secondadapter/1secondstrand thresholds. Both are retained as published computational wording, not a chemistry change or proven actual setting. Exact configuration used for a given run is unassembled; no physical70/140nt size-selection step is inferred. Full occupiedRNA/RTA/RMX/cDNA products, in-house enzyme production and exact historical kit internals remain unassembled, not confirmed absent.

The private graph selects the source commercial total yeast tRNA branch, beginning with its actual deacylation input. Cell culture/extraction, synthetic tRNA manufacture, PAGE for LC–MS/MS and RT-validation RNase digestion are outside this branch. Full source oligos have typed sugar identities and complete ordering fixtures. Computed adapter pairing establishes only the local core: the RTA divergent arms remain unpaired, and the exact poly(A)/T10 register is unassembled. Known adapter 5P donor geometry is retained without claiming a complete sealed native tRNA product; the RNA5 5′ chemistry remains unassigned. QC and unused stocks have separate fates. One actual25µl master supplies the RT enzyme, nucleotides and buffer; no extra primer or fresh duplicate components are invented.

Source TableS25 orders (5′→3′), retaining modifications and sugar notation:

| Component | Verbatim ordering |
|---|---|
| 5′ RNA padding | `{{= chem.RAW['rna5'] }}` |
| 3′ RNA:DNA padding | `{{= chem.RAW['rna3'] }}` |
| RTA A | `{{= chem.RAW['rta_a'] }}` |
| RTA B | `{{= chem.RAW['rta_b'] }}` |

Local padding geometry; the target tRNA and sealed products are unassembled:

```
{{scene: chem.padding_scene()}}
```

Cropped RTA paired-stem fragments only; divergent arms and poly(T) tail are omitted from this view. Complete source oligos are listed above. The RNA poly(A) binding register is unassembled:

```
{{scene: chem.rta_core_scene()}}
```

## Source-qualified public preparation page

The page now publishes the selected commercial total yeast tRNA route, with a curator-written literal inventory linking the source-supported stages to the model. The two existing SourceItems describe the selected branch and computational-setting difference; neither was a broad missing-library obligation. Both remain declared and covered. Whole occupied tRNAs, ligation junctions, RMX/motor and modified-base read behavior remain unassembled, without an absence claim. The complete workflow and unresolved chemistry remain visible in a collapsed section.

Analytical TapeStation/Qubit allocations are typed as QC. Used and remaining padding/RTA stocks retain the same known source oligos and modifications. The three source RNasin contributions (1 µl padding anneal, 1 µl padding ligation, 0.5 µl RTA ligation) are explicit molecular inputs. The residue parser remains unchanged: terminal `rN` is one RNA residue, making RNA5 entirely RNA; it is not a molecular identifier. Source main250 versus Extended Data Fig4caption300 read-classification settings are disclosed on the page. No inferred occupied tRNA or example genome sequence is drawn.
