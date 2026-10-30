# LISH-seq

[Credle et al., published defining article](https://doi.org/10.1093/nar/gkx471), *Multiplexed analysis of fixed tissue RNA using Ligation in situ Hybridization* (2017).

🟢 Sources read: published PMC5737328 Methods76–112; published supplement Additional Figure4 (originalp6), TableS1 hRPS19_2 (originalp7), TableS3 whole primers (originalp11). Mandatory whole-bundle oligo scraper run; original rendered table columns checked, rather than trusting interleaved extracted text. Two originals in this capped pass, no absence audit. Fixed-tissue immune-panel Plus-slide branch selected; LCM thermal release, analytical qPCR and fungal ONT adaptation remain alternatives of this method.

## Preparation

🟢 Source10µm formalin-fixed paraffin sections onPlus slides:60°C1h bake, Trilogy95°C2×30min,4water rinses,3×10minPBS/.1%Triton,2×5minprehyb2×SSC/20%formamide/.1%Tween. Hybridize250µlbuffer with5pMeach of an equimolar probe panel,60°C30min then45°C2h. Buffer includes2×SSC/20%formamide/.2mgmlBSA/2mMvanadyl/1mgmlE.colitRNA/.1gml dextransulfate. The tRNA is an actual molecular carrier contribution, not an inert omitted label.

🟢 Wash45°C2×SSC/2mMvanadyl/1mgmltRNA then37°CRnl2buffer. Ligate200µl1.5UµlRnl2,37°C2h,50mMTris/10mMMg/5mMDTT/1mMATP. Aspirate reaction mix and wash2×SSC. Release in150µlRNaseHmix6U,50mMTrispH8.3/75mMKCl/3mMMg/fresh10mMDTT,37°C30min. Aspirate released solution,20,000g15min4°C, collect cleared supernatant; tissue residue/pellet have distinct physical fates. OligoClean&Concentrator purifies released material. Exact recovery, RNA-cleavage spectrum and reporter preservation are unassigned.

## Whole source oligos and chemical distinction

🟢 TableS1 hRPS19_2 left3′-Diribo: `GGAGCTGTCGTTCACTCTAGAACCAGTTCTCATCGrUrA`; right5′-Phospho: `/5Phos/GGGAGCAAGCTCTTTGTGCTAGATCGGAAGAGCACAC`. Both probes contain37 nucleotides. The intended sealed reporter has two internal RNA residues, with a terminal RNA-A3′OH joined to the right DNA-G5′P. RNA-template cleavage does not itself turn these reporter RNA residues into DNA. The PCR-copied DNA sequence and original chimeric probe class must be distinguished. Full target occupancy, other panel members, off-target ligations and residual unligated probes remain populations.

🟢 TableS3 OF `AGATCGGAAGAGCACAC`; OR `GGAGCTGTCGTTCACTC`. P5/i5 ordering `AATGATACGGCGACCACCGAGATCTACACBBBBBBBBACACTCTTTCCGATCTGGAGCTGTCGTTCACTC`; P7/i7 ordering `CAAGCAGAAGACGGCATACGAGATBBBBBBBBGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT`. B8 is explicitly “index nts”, not an IUPACB degenerate pool. Actual sample index values remain unassigned. The model retains the literal orderings and represents the index window as a typed sample-index placeholder.

🟢 CustomRead1 `ACACTCTTTCCGATCTGGAGCTGTCGTTCACTC`; customNextSeqi5 `GAGTGAACGACAGCTCCAGATCGGAAAGAGTGT`; standardi7 primer `GATCGGAAGAGCACACGTCTGAACTCCAGTCAC`.

## Source orientation discrepancy and sequencing

🟡 TableS1 prints the right probe3′tail identical to TableS3 OF. These literal strands do not provide opposed outside primers on the expected sealed hRPS19_2 reporter: an opposed primer at that end requires the reverse complement of the printed OF. TableS3 indexed P7 primer does end in that reverse complement. AdditionalFig4 is a schematic, and does not resolve the full effective preamplification ordering. We retain the printed sequences and leave preamp geometry unresolved; no source primer is silently repaired, and no claim is made that the experimental method failed. A future public schematic must display this source comparison. Rechecking the two originals/predecessor design and author clarification remain open, not an accepted absence claim.

🟢 SourceHerculaseII preamplification20µl:95°C2min,20cycles95°C20s/59°C30s/72°C30s,72°C3min (30cycles for very small/LCM specimens, not merged into selected routine branch). Select2µl into20µlindexPCR10cycles. Pool actual indexed products, then twoQIAquickpurifications. SourceNextSeq500orHiSeq2500SE50, customRead1 and customi5 forNextSeq500. Source describes40nttarget-derived products; trim3′10readbases and align to possible probe-pair reporters. Mismatched ligation (~5%) is computationally filtered, not physically removed by an invented purification.

## Limits and linked sources

🟡 Source reference11 delegates probe design; directly linked predecessor and historical manufacturer column/polymerase instructions remain ordinary unread originals outside this two-original pass. Main's ONTNSK-007 hairpin/leader fungal adaptation, LCM alternative and qPCR validation were inspected as boundaries but not modeled or split into extra protocols. Full reporter/template Scene, effective preamp primer, index landing, complete panel/other library allocations and manufacturer internals remain unassembled. No missing-information claim has passed an absence audit.
