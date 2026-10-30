# o2n-seq

Ultrasensitive and high-efficiency screen of de novo low-frequency mutations by o2n-seq. Nature Communications (2017), [doi:10.1038/ncomms15335](https://doi.org/10.1038/ncomms15335).

Source scope: 🟢 published main PMC5458117 Methods121–151, Results65–67, analysis165–172; published supplement PDF Fig1 caption and remaining figure captions read. Two originals read; the Fig1 original graphical layout remains deferred because text extraction omits its strand drawing. Mandatory scraper reviewed (8 candidates/3 files). No complete linked-source or absence audit. The private graph selects the uncaptured o2n route, not comparator Cir-seq or AmpliSeq.

🟢 o2n-seq uses nick-controlled copying of a single-strand circle. It does not assert independent agreement of the two original duplex strands. The source contrasts controlled nick-directed strand displacement with unbounded RCA; importing Phi29/GenomiPhi from its separate Cir-seq comparator would change the method. Reported consensus is from the two tandem copies within PE reads, not cross-cluster barcode families.

🟢1–4µgDNA in130µlAE is shearedCovarisS220,10%duty/intensity5/100cyclesburst/600s below4°C, purifiedOligoCleanConcentrator,4%agarose80V70min,60–120bp bands extractedMinElute6×QG. The alternative900s/no-gel route remains separate. E7442S endprep:55.5µlDNA+6.5µlbuffer+3µlenzyme,20°C30min/65°C30min. AddE7445Sligasemix15µl+adapter2µl+enhancer1µl,20°C30min/65°C5min/ice. MinElute cleanup,65°C5min/ice,1.8×AMPure; QubitdsDNAHS is independent QC.

🟢 First-adapter orderings, both5′ phosphorylated:

```
/5phos/5′-GATCAGTCGTACGTGCTTACTCTCAATAGCAGCTT-3′
/5phos/5′-GTGGGCAGTCGGTGAACGACTGAUCT-3′
```

The internalU is DNA deoxyuridine, retained as a nicking substrate. Equimolar strands are annealed to50µM in10mMTrispH7.5/1mMEDTA/0.1mMNaCl,94°C5min then0.1°C/s cooling toRT. The published Y structure/full sealed ligation geometry is not reconstructed from the caption alone.

🟢30µlDNA<300ng is denatured95°C3min/ice3min;3µl10×CircLigaseIIbuffer,1.5µl50mMMnCl2,1.5µlCircLigaseII(stated100U) are added.60°C14h/25°C2h/95°C2min/ice3min. ExoI1µl+ExoIII1µl+Fpg1µl37°C1h; MinElute3×ERC, then independentQubitssDNAQC. No perfect linear/damaged-DNA digestion or recovery is imposed.

🟢 Second-strand annealing:15.8µlcircleDNA+1µlprimer+2µlNEBuffer4,95°C3min/45°C5min/ice3min. Add0.5µl2.5mMdNTP/0.2µl100×BSA/0.5µlT4DNApolymerase,25°C30min/75°C20min. USER1µl37°C30min/50°C4min/ice,1.5×AMPure. Bst2WarmStart strand displacement:16µlDNA+2µlisothermalbuffer+1µl2.5mMdNTP+1µlenzyme,60°C1h;1×AMPure. Exact primer identities, nick sites and products remain unassembled rather than replacing them with generic random priming or unbounded concatemer predictions.

🟢 Standard NGS finishing is another endprep/ligation:18.5µlDNA+2.17µlbuffer+1µlenzyme20°C30min/65°C30min, then5µlligasemix+0.33µlenhancer+0.17µlNEXTflex514102barcodeadapter20°C30min;0.8×AMPure.50µlPCR contains24µlDNA,1µlNEXTflexPrimerMix,25µlKAPAHiFiHotStartReadyMix.98°C45s;5–11cycles98°C15s/65°C30s/72°C1min;72°C4min/4°Chold. Twice0.8×AMPure or2%gel size selection precedes optional capture/sequencing. Source o2n acquisition is HiSeq2500PE125; comparatorCir-seqPE250 and validationAmpliSeqPE150 are separate. Computational first5base/internaladapter trimming and PE-copy consensus do not remove physical DNA.

🟢 Optional HCC0.42MbRocheNimbleGen capture adds an o2n-blocker1000µM to the manufacturer's capture conditions; it is not imposed on every library. Exact second-strand primer/blocker and historical NEXTflex adapter/index/primer geometry remain open. The main intro explicitly contrasts Cir-seq and Droplet-CirSeq, but this does not authorize chemistry inheritance; source-specific curated parent relations need further scoped evidence.
