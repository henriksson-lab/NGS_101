# LABS — selected extracted-cfDNA preparation

Defining paper: [LABS: linear amplification-based bisulfite sequencing for ultrasensitive cancer detection from cell-free DNA](https://doi.org/10.1186/s13059-024-03262-2), 2024.

Private incomplete source-qualified model. The bounded pass inspected published main129/211–237, SI FigS1caption and TableS1 sequencing-information workbook including cfDNA-comparison rows. Three originals read; no complete linked-source/missing-information review. Historical kit internals remain unassembled.

🟢 The chosen input is already extracted cfDNA. Main215 describes plasma workup and QIAamp circulating-nucleic-acid extraction from0.3–2mlplasma, but these are outside the selected graph start. The150–400bp gDNA preparation from MethylC-seq comparison is not imposed on native cfDNA.

🟢 Main219 gives the complete long-strand ordering and blocked helper:

```
/5Phos/iMe-dC/iMe-dC/iMe-dC/TATAGTGAGT/iMe-dC/GTATTAATTT/iMe-dC/G/iMe-dC/GGGG/iMe-dC/T
CGACTCACTATAGGGT/3Phos/
```

All literal long-strand cytosines are5mC; the helper's3′phosphate is a distinct block. Do not methylate the unmodified helper by analogy. Equimolar50µM oligos are annealed95°C5min, cooled−0.25°C/min to4°C, diluted15µM and stored−20°C. Preannealed adapter is a produced reagent consumed at ligation, not an extra later stock addition. Exact whole cfDNA/adapter junction geometry is unassembled.

🟢 KAPAHyperPrepKK8504 endrepair/A-tailing acts20°C30min→65°C30min. Source ligation alternatives are20°C4h or4°Covernight; the selected route is20°C4h. Postligation1.4×AMPure cleanup precedes MethylCodeMECOV50 bisulfite conversion, producing9–20µl eluate. Unmodified DNA C is converted to chemical U while5mC is protected; damaged/nicked/incompletely converted classes remain possible. This is distinct from the eventual DNA T readout.

🟢 EpiMarkHotStartTaqM0490S extends with the complete T7 primer:

```
AGCCCCGCGAAATTAATACGACTCACTATAGGG
```

No extension thermal programme is assigned by this selected source pass. Qiagenprotease50°C2h→75°C30min follows extension. No invented intermediate physical cleanup removes the protease reaction before HiScribeT7HighYieldE2040S/SUPERaseIn transcription37°C12–16h.

🟢 DNaseI37°C10min cleaves DNA templates; RNA CleanandConcentratorR1013 is the subsequent physical cleanup. QubitRNAHS material is separate from production. At most100ngRNA is allocated to KAPARNAHyperPrepKK8541; any remainder has an unresolved fate. Source overview129 explicitly includes reverse transcription, second-strand synthesis, adapter ligation and amplification. They remain one source-specific compound boundary because exact historical enzyme additions, oligos and substrate/geometry are not assembled; another KAPA/Illumina recipe is not substituted.

🟢 Published TableS1 marks the cfDNA comparison librariesSE75. This is kept separate from the main227/231 NextSeq500 SR80 descriptions of the MethylC/EpiGnome comparators. Instrument, index geometry and complete primer/read landing are not filled from those comparators. Computational methylation calls and tissue classification do not physically select or remove molecules.

No documented parent is asserted merely from shared T7/linear-amplification concepts. Cited prior methods and historical vendor guides remain unreviewed targets, not confirmed absent information.
