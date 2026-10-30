# TF-chRDP: shared p53 chromatin capture, separate DNA and RNA libraries

Defining published original: Ning et al., [Frontiers2025 DOI10.3389/fcell.2025.1561540](https://doi.org/10.3389/fcell.2025.1561540).

## Inspected source and boundary

🟢 Published main PMC11925928, actual header/reagents/Methods265–479/Results497–503 and originalFigure1 PDFp9. This is one original, with XML/PDF representations, already selected and closed by the coordinator in `/tmp/ngs43-root-originals.tsv`; implementation rereads the same original only. Mandatory bundle scraper returned0hits, which is not an oligo-absence claim. Historical Tn5/index-kit/mix manuals and any linked supplements remain ordinary unread. No whole canonical Nextera/Illumina oligos are borrowed.

The selected HCT116/5-FU/p53 route starts with same-treatment cells, applies FA then EGS dual-crosslinking, and captures intact protein/DNA/RNA complexes. Three physical post-IP allocations supply DNA, RNA and LC-MS/MS; these are separate fractions, not a single molecule with joined DNA/RNA/protein sequence tags. Protein output is retained as an analytical fate, not a fabricated NGS library.

## Shared preparation and real fractions

🟢 Source culture receives375µM5-FU for9h. Per-dish references are20ml1%FA20min, glycine1.74ml2.5M10min, washes, then20ml2mMEGS45min and the same glycine quench/washes. Approximately20million cells per dish are pooled by treatment, then30million-cell samples are selected. Dish multiplicity and aggregate reagent amounts are unassigned rather than inventing a one-dish30million-cell preparation.

🟢 Prepare150µlProteinGbeads per antibody and50µlpreclear beads per source sample; prewash each three times. Anti-p5310µg binds washedIPbeads in500µl0.1%Triton/DPBS8h4°C. Nuclear lysis/permeabilization precedes three actual1ml sonication tubes per sample; source suspensions are recombined before clarification. Recovered chromatin1ml-tube distribution and subsequent aggregate IP allocation convention remain unassigned.

🟢 After preclearing, allocate **two10µl controls**, one DNA and one RNA, before IP. Actual remaining chromatin enters washed antibody-bound beads overnight4°C; excluded supernatant and low-salt/high-salt/TE wash fractions stay explicit. Only after washing are actual bead-associated complexes allocated into the three DNA/RNA/protein tubes, with unassigned fraction sizes.

## Separate downstream routes

🟢 DNA and DNAQC use65°Covernight ProteinaseK decrosslinking, magnetic removal of production beads, Zymo DNA purification,15µlelution and14µlcollection. Production and QC DNA remain separate. DNAQC supplies source Tn5-testing reactions and analytical fragment/quantity checks.

🟢 RNA and RNAQC instead use50°C45min ProteinaseK,95°C10min, chill and500µlTRIzol. Chloroform extraction retains400µlaqueous fraction, followed by600µlEtOH and miRNeasy micro purification. Actual80µlDNaseincubationmix is the same physical stock containing DNase andRDDbuffer, not an added fresh enzyme plus duplicated whole mix. Source RWT flow-through is reapplied once, then discarded; RPE/EtOH washes and14µlelution precede65°C15min DNase inactivation.

🟢 Select13µlRNA, add1µl50ng/µl random hexamers and1µl10mMdNTPs. RandomN6 is a primer family, **not a deduplication UMI**. First-strand SSIII reaction totals30µl; second-strand synthesis uses actual E.coliDNApolymeraseI/RNaseH/DNAligase at16°C2h, then3.3µlT4DNApolymerase12°C15min and10µl0.5MEDTA. Source294µlAMPure1.8×,14µlelution/12µlcollection and Qubit/Qsep precede RNA-control testing or RNA-library preparation.

🟢 Section3.18 explicitly reuses local3.13 tagmentation/PCR/size selection for RNA-derived double-stranded cDNA. The DNA-testing recipe specifies0.33µlTn5 and8µl4ngDNA. Optimized production dose/transfer conventions and RNA testing amounts remain unassigned; those reference numbers are not asserted as measured optimized amounts. Catalyst and functional DNA-adapter inventory come from the same actual selected transposome stock; unknown whole adapters are not separate fresh reagents. Source25µlHigh-Fidelity2×PCRMasterMixM0541 supplies enzyme/nucleotide/buffer components from that same physical allocation. Indexone(i7)/Indextwo(i5)1µl each have unassembled whole orders/selected values. Source DNA-testing reference programme:55°C10min tagmentation;3µl0.2%SDS37°C5min/RT5min/ice>5min; PCR72°C3min/98°C30s/eightcycles98°C15s60°C30s72°C3min/72°C5min.

🟢 Testing uses input-control DNA or input-derived cDNA and1×AMPure cleanup, ending as QC. Production uses the corresponding captured DNA or RNA-derived cDNA. After adjusting PCR product100µl, add61µlAMPure and RETAINsupernatant; add19µl more to that fraction, giving0.8× relative to the original100µl, and RETAINbeads. Elute8µl, collect6µl, quantify and select a loading aliquot, retaining unused stocks. Both independent libraries have NovaSeq6000PE150 acquisition; they are not physically combined with each other or their controls in this graph.

## Molecular limits

Whole loaded Tn5 adapters, indexed primer bodies, historical kit formulations, exact substrate sequences/crosslinks, transcript priming sites, genomic target junctions and complete read landing remain unassembled. Protein fraction is a non-NGS assay continuation outside this model, not a QC control. Source programme reuse is local and explicit; it does not establish equal DNA/cDNA substrates or products. Testing thresholds are source empirical calibration criteria, not imposed perfect yields. No full-molecule/absence audit or source-missing claim.
