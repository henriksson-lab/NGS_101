# MASD-seq — selected Appendix cfDNA preparation

Defining publication: [Methyltransferase-assisted single duplex sequencing for detecting circulating tumor DNA](https://doi.org/10.1038/s44321-026-00509-2), 2026.

This is a source-qualified private partial model. Published main90/96–100/298–310, Appendix PartI126–502(printed14–20), and the original TableEV1 workbook were inspected. Three-original combined-reader cap is closed; the peer-review file and data-only workbooks were not opened. Mandatory primer scraping returned source pointers, followed against these originals. Full cfDNA product paths, conversion yields and read landing are not validated or publicly promoted.

🟢 The selected input is20ng extracted plasma cfDNA. Tissue MspI-RRS, WGS and the separate bisulfite comparison remain distinct preparations. MspI-RRS's CIP/ddGTP/end-repair/VAHTS-UMI recipe is not inherited into MASD-seq.

🟢 Appendix PartI uses one50µl reaction containing42µl cfDNA,5µl10×rCutSmart,1µl8mMSAM,1µlMspI100U/µl and1µlM.SssI20U/µl at37°C1h. Restriction and CpG methylation occur together, not in an invented sequential intermediate. The intended CpG methylation protects CpGs during conversion, and overwrites original CpG methylation information. Incomplete digestion/methylation remains possible.150µlAMPure cleanup and two200µl70%ethanol washes give41µl eluate;40.5µl enters ligation, with the remainder explicit.

🟢 Selected Appendix hairpin ordering:

```
5'-p-CGGCCGCAGTCDDDD/Int Biotin dT/DDDDGACTGCGGC-3'
```

TableEV1 lists this same ordering and a second `CGGCCGCTCAG...CTGAGCGGC` hairpin variant; the selected Appendix explicitly uses the former. D denotesA/G/T, not unrestrictedN. Eight random positions flank the internal biotin-dT and have source UMI meaning. Main310 records a4-base computational UMI tag when available; its complete correspondence to the eight source random positions remains unassigned. No random bases are complemented by eye.

🟡 Local fold is computed from that whole source ordering: the two terminal arms contain an exactly complementary9bp stem, a9nt D4–biotin-dT–D4 loop and a2nt5′CG cohesive interface. The paired-end Scene is a local view of two arms of the same covalent oligo; the complete ordering is retained separately. The fold can pair; no thermodynamic occupancy or homogeneous folding is predicted.

🟢 Whole Y-strand orderings:

```
5'-ACACTCTTTCCCTACACGACGCTCTTCCGATCT*G-3'
5'-p-CGCAGATCGGAAGAGCACACGTCTGAACTCCAGTC-3'
```

The top strand has the stated terminal phosphorothioate; the bottom strand has5′phosphate. Appendix says all adapter Cs are methylation-modified. TableEV1's original rich-text runs mark body Cs in red and the initial bottomCGC in blue as part of the reconstituted NotI site. Blue is a recognition-site annotation and does not establish that those Cs are unmodified. The selected model follows Appendix's all-C instruction and retains this annotation scope.

🟡 Exact source Y orderings compute a14bp terminal stem with unpaired fork arms and a bottom-strand5′CG interface. Both adapter interfaces are compatible with the MspI5′CG end derived from its source C|CGG cut. This narrow known cohesive geometry is modeled; entire sealed fragment sequences, ligation alternatives and their proportions remain unassembled.

**Published main/Appendix differences are not silently repaired:**

| Component | Main source | Appendix selected here |
| --- | --- | --- |
| Methylation order | Fig1caption96 puts methylation after hairpin ligation | Methods298/Results98 and Appendix PartI perform concurrent methylation/digestion before ligation |
| Ligase amount | Main298:2000U |2µl at2000U/µl, i.e.4000U |
| Ligation duration | Main298:1h at20°C |30min at20°C, then65°C10min |
| Formamide elution | Main300:10µl |20µl,65°C5min |
| Deamination MgCl2 | Main300:5µM |3.5µl100mM stock in70µl, i.e.5mM |

The Appendix protocol supplies the selected operational values because it is the explicit step-by-step preparation. The source discrepancies remain unresolved. Appendix's DEA buffer is literally labeled1× while7µl is used in70µl; no corrected stock strength is invented. Post-NotI150µl beads are also retained literally instead of recalculating a3× volume from51µl.

🟢 In the selected50µl ligation, hairpin8ng and Y-adapter4ng are source-ready stocks. Their preparation products feed ligation once, with no new free strands added again. NotI20U37°C15min cuts adapter-dimer classes; physical cleanup follows, rather than declaring all dimers destroyed. Elute40µl after150µlAMPure.

🟢 Capture uses MyOne StreptavidinC1(typically10µl), washed three times and resuspended40µl2×wash buffer before40µlDNA addition. Bind15min at room temperature; the noncaptured supernatant is explicitly discarded, then retained beads are washed three times with200µl1×buffer.95%formamide/10mMEDTA(pH8.2),65°C5min releases DNA into transferred supernatant. This is not a covalent biotin-cleavage reaction. Residual bead-associated DNA remains an unresolved physical fraction.60µlAMPure cleanup gives28µl eluate, of which27µl enters protection.

🟢 EM-seqE7125 oxidation/protection uses the declared45µl mixture plus5µl freshly diluted1:1000Fe active reagent,37°C1h, then1µlstop reagent37°C30min. The source describes TET2 oxidation and beta-glucosyltransferase protection; exact kit component identities are not imported from another protocol.150µlAMPure yields54µl;53.55µl enters the70µl APOBEC/UvrD reaction. Appendix uses1.75µl100mMATP,3.5µl100mMMgCl2,1µlAPOBEC mix,2.5µlUvrD(20µg/ml,50ng) and0.7µlBSA,37°C90min. Its separate initial denaturation step is omitted. Helical strand separation and C→chemicalU conversion occur together while successfully ligated original strands remain covalently linked.

🟢210µlAMPure yields23.5µl;22.5µl enters50µlPCR with25µlKAPAHiFiHotStartUracilPlusKK2802 and2.5µlUDI primer mix. Source forward/reverse complete fixed arms and two8nt sample-index windows are modeled, with actual index values unassigned. PCR98°C1min;14cycles98°C20s/65°C30s/72°C30s;72°C1min. Final60µlAMPure yields25µl. Qubit/fragment-analyzer QC is separate from production. NovaSeqXPlus uses150bp paired-end reads. Chemical U versus the reported DNA T base is explicitly distinguished; the graph does not predict a duplex error rate or correction of all transition classes.

Main90 explicitly describes modification from five-letter sequencing(10.1038/s41587-022-01652-0), and main98 scopes inheritance to its UvrD-assisted EM-seq component while changing original-strand linkage and adapters. That method is currently only a queued candidate, so no imported constructor or registered-parent link is manufactured. CODEC is a comparison, not an asserted ancestor. Historical kits and remaining geometry are unreviewed, not confirmed absent.
