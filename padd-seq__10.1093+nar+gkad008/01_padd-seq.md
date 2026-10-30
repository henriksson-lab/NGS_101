# PADD-seq — selected Pol II/CPD preparation

[Genome-wide mapping of protein–DNA damage interaction by PADD-seq](https://doi.org/10.1093/nar/gkad008),2023.

This source-positive private partial model uses the published main91–115/129 and published supplement TableS1 p2 plus Supplementary Protocol160–413 (printedpp13–17). Original oligo table p2 and PCR table p16 were rendered and inspected. Mandatory scraper pointers were checked there. Two originals were counted before opening; the source investigation is closed. Historical kit manuals, Damage-seq predecessor originals and full molecular products remain unassembled. No confirmed missing-information or full reaction-validation claim is made.

🟢 The source first selects protein-associated chromatin, reverses crosslinks and isolates DNA, then selects damage-containing DNA and maps stalled primer extension. This serial physical selection distinguishes PADD-seq from independent ChIP and damage profiles. The selected detailed supplement route is PolII/CPD with NVP-2, not the separate H3K9ac/cisplatin preparation or photolyase-repaired Input/ChIP-seq comparators.

🟢 Three15cm dishes of approximately3×10^7 fibroblasts receive250nMNVP-2 for5min before20J/m²254nmUV-C; recovery time is experiment-specific and remains unassigned. Crosslink1%formaldehyde10minRT; quench150mMglycine5min; wash/pellet. Sequential lysis buffers1/2 retain pellets, then630µl1.5%SDS-RIPA produces mixed lysate. Source Q800 sonication targets300–600bp. The mixture remains physically copresent until20000g clarification transfers supernatant; detergent treatment is not an invented soluble/insoluble separation.

🟢 Select167µg chromatin; dilute with14volumes SDS-free RIPA. Anti-RPB1 binding1–2h4°C precedes addition of separately washed/blocked60µlProteinAagarose; capture overnight4°C. Source washes three times each RIPA/RIPA500/LiCl and once TE. Elute240µl65°C15min; spin and transfer eluate; RNaseA37°C30min, ProteinaseK55°C2h, reverse crosslink65°C~10h, phenol/chloroform and ethanol. Main97 specifies25µg anti-RPB1, while detailed SI217 specifies25µl of the same antibody. Stock concentration connecting the mass and volume is unassigned; the model retains the detailed-volume recipe and both source statements.

🟢 Select50–100ng ChIP DNA (SI minimum25ng) for NEBNext UltraII end repair/dA tailing. Whole source orderings and modifications are:

```
Ad1T /5Phos/GATCGGAAGAGCACACGTCTGAACTCCAGTCA/3SpC3/
Ad1B NNNNNGACTGGTTCCAATTGAAAGTGCTCTTCCGATC*T
O3P GACTGGAGTTCAGACGTGTGCTCTTCCGATCT
Ad2T /5Phos/AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT/3SpC3/
Ad2B ACACTCTTTCCCTACACGACGCTCTTCCGATCTNNNNN/3SpC3/
```

Source100µM each oligo in50µl hybridization buffer (0.5mMTris8/100mMNaCl/0.05mMEDTA) is boiled2min and cooled naturally, with stock stored−20°C. Scene computes the local Ad1 terminal stem/3′T and full Ad2 donor/splint pairing. N5 pools are random splint/adapter windows, with no UMI meaning assigned. C3 blocks and the Ad1B terminal phosphorothioate are retained. These local interfaces do not establish complete sealed genomic product geometry or yield.

🟢 Source Ad1 ligation uses50pmol (1.25µl40µM stock)4°Covernight. Ad2 later uses1µl40µM stock. Each model consumes a selected aliquot from its independent prepared stock and leaves an explicit untransferred remainder. The100µM annealing-to40µM working-stock dilution bridge is unassembled. Entire50µl stocks are not silently consumed by these ligations.

🟢 After Ad1,0.8×FragSelect cleanup gives12.2µl; add5µl8Murea,98°C2min/ice1min. Add2.5µl8×IP buffer and5.5µl separately prepared anti-CPD beads overnight4°C. Preparation combines1.25µlProteinG and1.25µlanti-rabbit beads, rabbit anti-mouse IgG, anti-CPD and salmon-DNA blocking, with actual bead washes before addition. Sequential U/I/II/III/IV/TE washes retain captured material. Elute100µl1%SDS/TE65°C10min; isolate eluate, phenol/chloroform/ethanol, dissolve6µl. Binding is not assumed perfectly specific or quantitatively complete.

🟢 Add7.5µlQ5 master mix and30pmolO3P for extension98°C50s/65°C5min. Lesion-stalled, full-length and unextended classes remain alongside original templates; stalled product is not a perfect yield assertion. Add1µlExoI37°C15min; digestion does not physically remove molecules. Add25µlTE,1.1×FragSelect45.6µl, elute100µl; ethanol concentrate6.5µl. Denature98°C3min/ice1min; Ad2/InstantStickyEndsLigase4°Covernight;1.1×FragSelect55µl gives15.5µl.

🟢 Allocate15µl to production,0.3µl to analytical22-cycleQ5PCR/5%PAGE, and leave the nominal0.2µl untransferred remainder explicit. Production PCR uses15cycles, not the QC cycle count:98°C30s;15×98°C10s/65°C75s;65°C5min. SI originalp16 literally specifies10mM i5/i7 stocks for both QC and production; this suspicious unit is retained as an unresolved source value, not silently changed toµM. Exact historical NEB multiplex primer orderings and selected index identities remain unassembled.

🟢 Final size selection adds50µlTE then40µlbeads, transfers the supernatant, then adds50µlfresh beads, discards that second supernatant, washes twice70%ethanol and elutes16µlTris8, nominal300–1000bp. The first bead-bound fraction's later fate is unassigned. QubitHS and NovaSeq PE150 follow. The two-nucleotide-upstream damage-site inference and read filters are computational, not physical DNA loss.

Main77/99 explicitly describes a modified Damage-seq component; ancestry is scoped to that component, without importing a predecessor constructor or unreviewed historical recipe. Independent critical review checked the stated main/SI scope and literal original tables; it is not a complete linked-source or absence audit.
