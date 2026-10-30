# SMM-seq

Single-molecule, quantitative detection of low-abundance somatic mutations by high-throughput sequencing. Science Advances (2022), [doi:10.1126/sciadv.abm3259](https://doi.org/10.1126/sciadv.abm3259).

🟢 Sources: published PMC8993124 main76/82/134–140; original PDFp4; publishedSI tablesS1–S3/FigS1–S4captions. TWO originals counted before opening; source investigation closed. SI tables are donor/results tables, not an accepted absence finding for linked oligo disclosures. Mandatory scraper reviewed6hits/2representations of main. Historical NEB/BIORON and cited predecessor originals remain unreviewed.

🟢 SMM-seq ligates hairpins at both ends to make dumbbell-like substrates, uses thermostable pulse-RCA to copy representations of both strands, and PCR-amplifies individual copied units. This differs from single-strand circle sequencing or vector-based CypherSeq. Exact copy number, independence, yield and error elimination are not guaranteed by this private incomplete graph; strand-family filtering remains computational.

🟢 Genomic DNA is Quick-gDNA-extracted/Qubit-quantified, digested AluI/MluCI37°Covernight and cleaned1.5×AMPure. PublishedFig1A/main82 explicitly identifies end-repaired/A-tailed DNA; this supports the broad NEBNextUltraII endprep operation, while its formulation/conditions/geometry remain unassembled. Ligation uses customP5_HP6N instead of kitadapter. Literal5′→3′ ordering:

```
TCTTC TACAGT NNNNNN AGATCG GAAGAG CACACG TCTGAA CTCCAG TC /ideoxyU/ ACACTC TTTCCC TACACG ACGCTC TTCCGA TCT
P5-RCA: GTAGGGAAAGAGTGTAGACTGGAGTTC
```

The6nt UID is structured as a UMI; `/ideoxyU/` remains DNAU, not T or N. Adapter preparation:86µl10µMoligo in0.1×TE,95°C5min/37°C5min;add10µlCutSmart/2µl10mMdNTP/Klenowexo-minus10U37°C30min; QIAquickNucleotideRemoval; HpyCH4III10U37°C1h; secondNucleotideRemoval/elute100µlEB. No undisclosed phosphorylation, phosphorothioate or precise cut/end geometry is invented.

🟢 After genomic adapterligation, double-sidedAMPure size selection precedes separateQubit/BioanalyzerHSQC. Source dilutionD=M×N_A×2/10^25 targets~8amol(~5million molecules)/µl underPE150/30Gb assumptions. Only1µl diluted preparation enters20µlpulse-RCA; unselected material remains outside this path. Mix:1µlP5-RCA,BIORONstrand-displacementpolymeraseHS25U(0.5µl),2µlbuffer,1µl10mMdNTP,0.6µl100mMMgCl2,13.9µlwater.

🟢 **Printed pulse program differs from the narrative:** main76 describes repeated denaturation/annealing/extension, but main138 and originalPDFp4 print92°C2min(step1),92°C30s(step2),60°C30s(step3),65°C150s(step4),“go to(3)ninetimes”,4°Chold. The model preserves that literal jump to step3; it does not silently repeat step2 or claim a corrected experimental program. The intended loop remains an open source discrepancy.

🟢1.5×AMPure/elute23µl; wholeRCAproduct goes to50µlPCR with25µlNEBNextUltraIIQ5 and listed1µlP5/P7dualindexoligos. Preserve source wording rather than assuming undocumented index identities or concentration. PCR98°C30s,98°C10s/65°C75s,go-to-step2eighttimes,65°C5min.0.7×AMPure/elute30µlTE; separateQubit, pool and NovaSeqPE150. The conventional control uses undiluted ligation material and USER before PCR; USER is deliberately not inserted into the selected SMM path, whose dumbbell U remains intact during pulse-RCA.

Whole hairpin text matching is exempted narrowly because the shared parser splits `TC /ideoxyU/` across hits. A local test normalizes the complete published modified ordering and checks every base; the explicit raw fixture and original page locator are in the module. This is not an exemption from source verification. Open: historical kit components, exact fold/Hpy cleavage/sealed junctions/RCA/PCR/read geometry and intended pulse loop. No confirmed missing disclosure, full validation or public schematic is claimed.
