# Five- and six-letter epigenetic sequencing

Defining publication: [Simultaneous sequencing of genetic and epigenetic bases in DNA](https://doi.org/10.1038/s41587-022-01652-0), 2023.

One protocol entry contains two source preparation branches. The private partial model selects80ng genomic DNA. The published main article121/123/185–195/203 and published SI (original Fig7 p8 control, captions and Tables1–3 identities) were inspected. Two originals were counted before opening; the source investigation is closed. Mandatory primer-scraper pointers were checked against the actual main orderings. Superseded preprints, comparator protocols and other originals are not authority for this model. Full linked-source/missing-information review remains unfinished.

🟢 Sonicated genetic DNA receives hairpins, then a linked unmodified copy before methylation protection and APOBEC conversion. Comparison of original and copied strands provides genetic/epigenetic information. Six-letter sequencing adds early5hmC glucosylation and DNMT5 methylation copying; it is a branch of this preparation, not another standalone entry. This is distinct from MASD-seq linking original Watson/Crick strands and overwriting sample CpG methylation. No MASD constructor is imported.

🟢 Source controls include unmethylated pUC19, methylated lambda and a synthetic5hmC duplex (SI Fig7 original inspected; green5hmC and gray CpG annotations are meaningful). The selected model starts with the source ground-truth spike stock; composition and component ratios for an individual library are unassigned. It does not assert that every library receives all three controls. Study DNA is not subjected to the lambda control's M.SssI methylation recipe.

🟢 The selected80ng gDNA is sonicated to approximately250bp; cfDNA preserves its native profile. Add1µl spike at0.5ng/µl for>=80ng input. The lower-input spike stock is a separate branch. End preparation uses3.5µl EndPrep buffer and1.5µl enzyme mix (NEB E7647),20°C30min then65°C30min. Detailed historical enzyme formulation and resulting substrate ends remain unassembled.

🟢 Adapter1 ordering (all C are5mC; U denotes deoxyuridine):

```
ATGACGATGCGTTCGAGCATCGUCAUT
```

No5′phosphate is invented. The candidate local hairpin fold is computed from the whole ordering and marked inferred. The shared Scene checks complementary fixed stem arms and shows the terminal3′T interface. Both U residues lie in the paired arm, not the loop. Accordingly the model does not reuse a generic USER-loop-opening recipe. Exact post-USER strand paths, released fragments and genomic junction/new-end geometry remain unassembled; chemical cleavage does not itself physically remove fragments.

🟢 Add3.75µl adapter1,0.5µl ligation enhancer and15µl ligation master mix E7647,20°C15min. SPRI, elute23.75µl. Add3µl10×rCutSmart and3.25µl USER,37°C30min. SPRI, elute12µl. Add2µl each10×NEB buffer4,10mM-each dNTPs,Klenow exo-minus and T4PNK,37°C30min. The source3′end-copying relation is retained as unassembled; no external primer stock or fabricated intramolecular primer port is added.

🟢 Six-letter branch adds2.5µl UDP-glucose and1µl T4BGT in that same PNK/Klenow reaction, before methylation copying. Both branches then heat95°C2min and cool−0.1°C/s to reanneal.

🟢 Adapter2 is a ready Y-adapter stock with whole orderings:

```
ACACTCTTTCCCTACACGACGCTCTTCCGATC*T
GATCGGAAGAGCACACGTCTGAACTCCAGTCA
```

All C are5mC; the top terminal phosphorothioate is retained, and source5′phosphorylation remains unassigned. The computed local terminal stem and3′T interface are checked with Scene. This establishes only adapter local compatibility; exact substrate ends and complete sealed products remain unassembled. Add2.5µl TOTAL ready adapter stock,0.5µl enhancer,15µl ligation mix and12µl water,20°C15min. Do not add2.5µl of each free strand or a second adapter stock. SPRI follows; five-letter eluate is30µl.

🟢 Six-letter branch then uses15µl DNA,3.95µl water,0.75µl100mMATP,1.5µl100mMMgCl2,0.4µlSAM,6µl5×DNMT buffer (250mMTris8/10mMDTT) and2.4µlDNMT5,37°C3h. DNMT5 occurs after adapter2 ligation, not before. Post-DNMT5 transfer/cleanup bookkeeping into the30µl protection input is unassembled. The model inserts no invented SPRI or perfect enzyme clearance.

🟢 Both branches protect30µl DNA with10µl TET2 fivefold supplement stock (10mMαKG/0.25MTris8/10mMATP),1µl UDP-glucose,1µl T4BGT,1µl100mMDTT and2µlTET2. Add5µl fresh1:1250 dilution of500mMFeSO4 and incubate37°C60min. The literal iron stock/dilution is preserved. This protection-stage BGT/UDP-glucose is a fresh source addition, distinct from the six-letter early reaction.

🟢 SPRI/elute31µl. Add12.75µlwater,17.5µl4×APOBEC buffer (200mMBisTris6.1/0.4%Tween stock),1.75µl100mMATP,3.5µl100mMMgCl2,2µlAPOBEC-A3A and2.5µlUvrD,37°C1.5h. SPRI/elute20µl. Protected cytosine derivatives and residual unprocessed substrate coexist with C-to-U products; no conversion yield is assumed. Chemical U is distinct from PCR-reported T.

🟢 PCR uses20µlDNA,5µlAbclonal UDI primer pair RK21624SetA and25µl2×KAPA HiFi Uracil+ KK2802. The two primer components come from one supplied5µl pair stock; their exact historical orderings and selected indexes are unassembled. Programme98°C30s, then6cycles98°C10s/62°C30s/65°C60s, then65°C5min for80nggDNA. The10ng and2ng cfDNA routes use7 and8cycles respectively and are not imposed on the selected genomic route. SPRI/elute15µl10mMTris8; TapeStationD5000 aliquot remains separate from production.

🟢 Dilute1.2–1.8nM; add8%PhiX to converted libraries. NovaSeqS4 run alternatives are111/8/8/111 or151/8/8/151; the selected run remains unassigned. Source R1 is P7-primed original, R2 P5-primed copy. Full molecular primer/cycle geometry remains unassembled. Computational N masking or filtering does not physically remove DNA.

Remaining work is source-qualified assembly of USER/copy/junctions, historical Abclonal primer orderings and six-letter DNMT5 transfer bookkeeping. SI Tables1–3 are results/fragment statistics in the inspected scope, not substitute primer recipes. Unreviewed historical kits and linked predecessor methods are deferred reading, not evidence that information is absent. MASD's source-documented adaptation can be curated as a scoped ancestry edge separately; matching adapter handles never authorizes shared constructor inheritance.
