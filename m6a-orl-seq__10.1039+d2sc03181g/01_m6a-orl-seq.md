# m6A-ORL-Seq

[Xie et al., Transcriptome-wide profiling of N6-methyladenosine via a selective chemical labeling method](https://doi.org/10.1039/d2sc03181g), Chemical Science2022.

🟢 Read the published defining main, PMC9600483, Figure1 and Methods73–121; published ESI part10, S20–S22, with actualPDFp20/21 tables and TableS1 oligos. The third original, supplementary XLSX-s001, was inspected only for its first `MeRIP_1` sheet header/coordinate rows: it is analytical data, not a chemistry source for this model. Mandatory whole-bundle scraper pointers were verified against the selected originals. The source investigation is closed at three originals. This is a private, source-qualified partial model.

## Selected chemical preparation

🟢 The selected branch begins with5µg HEK293T poly(A)+RNA, fragmented70°C15min to the source60–200nt range. RNA extraction/DNase/poly(A) capture is source-described upstream; complete historical kit operations and actual RNA sequences remain outside this graph. Model oligo, total-RNA, MeRIP and SELECT controls are not merged into this branch.

🟢 Main Figure1 and Methods describe sequential N-nitrosation of m6A to Nm6A, TDO reduction to Am6A, then BA-probe labelling to biotin-bearing Bm6A, with RCC purification after each chemical treatment. The RNA N-nitrosation route uses NO exposure and a4°C20h incubation in3M acetate4.5/2mM EDTA; the selected total/poly(A) RNA TDO route uses100mM TDO/300mM carbonate-bicarbonate10/2mM EDTA,50°C10min in30µl. BA labelling uses2mM BA/100mM MES/10%DMF5.5,37°C1h in50µl. Prepared NO, TDO and BA are physical molecular reactants, not catalyst defaults. Reagent manufacture is outside the graph.

🟡 The local typed residue markers follow that source order without asserting a complete atom map, conversion probability or occupied transcript. Unreacted RNA and side changes A→I, G→X and C→U remain possible/source-described classes; the nucleoside conversion yield is not assigned to transcript RNA. The source also discusses m6Am, so the preparation does not establish exclusive m6A specificity.

🟢 Add50µg yeast tRNA carrier **after** the three chemical treatments. Source25µl C1 streptavidin beads are prewashed twice with200µl100mMNaOH/50mMNaCl and twice with RNase-free water; resuspend in100µl binding solution (10µl source high-salt buffer plus90µl water). Bind1h, wash three times with1ml high-salt buffer. Elute first with50µl100mM DTT37°C15min, collect supernatant, then apply a second50µl DTT50°C5min to the residual beads. Pool the two actual eluates, RCC-clean and quantify with Qubit. Select50ng for library construction; stock and analytical aliquots remain separate.

🟢 Figure1B labels the released species Tm6A (thiolated m6A). DTT cleaves the **BA disulfide linker**, leaving a thiolated adduct attached to the RNA. It does not restore unmodified m6A or cleave the noncovalent biotin–streptavidin bond. Biotin-side material can remain with beads; whole fragment/adduct and carryover populations remain unassembled.

## Source library construction

🟢 ESI10a treats5µl RNA with1µl10×PNK buffer,1µlPNK,1µl inhibitor and2µl water,10µl total;37°C1h/65°C20min. The table lists no ATP addition. The model therefore does not claim this creates RNA5′ phosphate; the source-compatible intended3′OH ligation class remains unassembled rather than being assigned to every molecule.

🟢 ESI10b RNA3′ linker is `5′rAPP-AGATCGGAAGAGCGTCGTG-3SpC3`:{{= chem.RNA_LINKER_LENGTH }}DNA nucleotides, a5′ adenylate activation group and a non-nucleotide3′C3 blocker. The adenylate is an activation/leaving group, not an extra nucleotide in the linker backbone. Add1µl10µM linker to10µl RNA,65°C5min/ice; add14.5µl ligation mix (2.5µl buffer,2.5µl100mM DTT,7.5µl50%PEG8000,1µl Rnl2truncatedKQ,1µl inhibitor), making25.5µl. Incubate25°C2h/16°C overnight. Apply5′ deadenylase30°C1h **then** RecJf37°C1h. OCC-elute12µl and transfer10µl, retaining the untransferred remainder separately. Workup does not guarantee complete adapter/off-product removal.

🟢 RT primer is `ACACGACGCTCTTCCGATCT`,{{= len(chem.RT_PRIMER) }}nt. Its3′{{= chem.RNA_LINKER_LENGTH }}nt pair to the linker; the first5′A remains unpaired. Add1µl2µM primer to10µl RNA,65°C5min/ice, then9µl SSIII reaction mix (4buffer/2dNTP/1DTT/1enzyme/1inhibitor),20µl total. Source program is25°C3min/42°C7min/52°C30min/4°C hold. The immediately following instruction says not to raise samples above37°C to avoid denaturing conditions. Both are retained: its handling scope is ambiguous, and the model does not rewrite the stated RT program or impose a global37°C ceiling.

The following local duplex uses a partial postligation template: upstream RNA is inferred/unassembled, and the known19nt linker retains its3′C3 blocker. The donor's activating AMP has already left during ligation. Released AMP stays copresent with RNA/linker/off-products through the deadenylase/RecJf reaction mixture until OCC workup; its partition and recovery remain unassigned.

```
{{scene: chem.rt_binding_scene()}}
```

🟢 Add1µl ExoI37°C30min; then10µl1MNaOH/10µl0.5MEDTA/10µlwater65°C15min to hydrolyse RNA. OCC-elute7µl/transfer6µl. Add0.3µl80µM DNA5′ adaptor `5Phos-NNNNNNNNNNAGATCGGAAGAGCACACGTCTG-3SpC3` plus1µl DMSO,65°C5min/ice. This is{{= chem.DNA_ADAPTER_LENGTH }}backbone nucleotides; SpC3 is a non-nucleotide blocker, not an appended C. Add12.7µl mix (2buffer/.2ATP/9PEG/1.5Rnl1),20µl total;25°C overnight, OCC-elute30µl. Known donor ends are represented narrowly; whole cDNA acceptor geometry and product remain inferred/unassembled.

## Identifiers, QC and sequencing

🟢 Main119 calls the source N10 a random decamer barcode ligated to fragments, but also says preprocessing removes a **random hexamer**. ESI10d prints ten Ns. The physical model retains N10 and assigns only the stated random inline-barcode role; it does not infer a UMI, deduplication rule, cell identity or sample identity. The N10/N6 processing discrepancy remains unresolved and visible in the graph evidence notice; no silent normalization is applied.

🟢 ESI10e takes3µl of the30µl cDNA for10µl qPCR with5µl SYBR master mix and1µl each10µM forward/reverse primer. TableS1 gives FP `CAGACGTGTGCTCTTCCGATCT` ({{= len(chem.QPCR_FORWARD) }}nt) and RP equal to the RT primer. This QC aliquot is separate from production. ESI10f transfers23µl from the remaining27µl to50µl PCR with25µl Q5 UltraII master stock plus1µl each named NEB universal/index primer; nominal4µl remains untransferred. The one master stock supplies nucleotide/enzyme/buffer components, not separate fresh additions. Whole NEB primer bodies/index choices are unassembled rather than borrowed by product name. Cycles refer to qPCR; source98°C30s, cycling98°C10s/65°C75s,65°C5min final/4°C hold.

The one5µl prepared SYBR qPCR stock supplies copresent enzyme, nucleotide, fluorescent-dye and buffer components. Their separate doses and full identities remain unassigned.

🟢 Add45µl beads to50µl PCR (0.9×), elute23µl/transfer21µl and store at−20°C before sequencing. The actual loading allocation is separated from retained stock. Main reports HiSeqXTen PE150. Exact run primers, index identities/lengths and binding/read geometry remain unassembled; readout/peak-calling/trimming are computational, not physical cleanup. Physical RNA adducts, RT-stop signatures and copied DNA are distinct classes.

## Limits

🟡 Whole atom-resolved chemical products, RNA occupancy/end heterogeneity, reaction yields, all side products, RNA/cDNA ligation junctions, complete PCR/library molecules, kit operations and sequencing-primer landing remain unresolved. Source-labelled comparisons, chemical-reagent synthesis and predecessor recipes were not imported. The3′OH acceptor views are inferred by construction; narrow source-known donor ends do not establish all ligation products. No missing-information claim, full chemistry audit or second absence pass is made.
