# HiCuT — bounded published-source preparation

[Defining publication](https://doi.org/10.1371/journal.pgen.1010121). Two closed published originals: PMC8979432 main135–159 (chemistry143/147/151/153/157) and originalS4Table`s010.pdf` actualp1. Whole-bundle scraper `/tmp/ngs39-hicut-scraper.log` run before handoff; actualS4render independently inspected. No vendor/predecessor/new original followed.

🟢 HiCuT combines in-situ contact ligation with antibody-directed pAG-Tn5 targeting/tagmentation and indexed library amplification. Plain **dATP**, not biotin-dATP, is used in contact fill-in; no biotin-junction pulldown is inserted. ConA attaches nuclei for washing; it is not a physical antibody IP. Selected CTCF branch represents the preparation; PolII/H3K27ac and cell types are branches rather than extra methods. Complete occupied contact molecules, native fragment ends and loaded transposome oligo bodies remain unassembled.

🟢 Source100,000cells fixed.5ml1%formaldehyde10minRT/.2Mglycine5minRT+15minice; wash.05%BSA/PBS and2000g5min. Then3mMDSG500µlPBS40minRTrotation/.4Mglycine5min, repeatBSA/PBSwash. Hi-Clysis10mMTris8/10mMNaCl/.2%Igepal/proteaseinhibitor30minice;2500g5min/washNEBuffer3.1. BSA is an actual molecular carrier/additive, not an unexplained condition label.

🟢 Nuclei161µlNEBuffer3.1+19µl1%SDS65°C10min/ice;21.5µl10%Triton37°C15min900rpm. Add20µl10UµlDdeI/4µl50UµlDpnII/2µl1×NEBuffer3.1,37°C3h orovernight with30son/4minoff900rpm;65°C20min. Source calls35µl fill master:18.75µl.4mMdATP;.75µl dTTP/dGTP/dCTP at10mMeach;5µlDNApolI,37°C2–3hrotation. Listed nucleotide-volume wording/remaining balance is retained without inventing water or normalizing total. Actual332.5µl ligationmaster60µl10×T4buffer/50µl10%Triton/6µl10mgmlBSA/2.5µl400UµlT4ligase/214µlwater,RT2–3h. Pellet/wash200µlexchangebuffer, resuspend100µl. Select25µl; remaining75µl is not silently consumed or discarded.

🟢 Activate10µlConA beads with2×100µlactivationbuffer, resuspend100µl and add to25µlnuclei;RT10min, magnet/remove supernatant. Add50µlantibodybuffer/1µgCTCFprimary,2hRT orovernight4°C; magnet/remove supernatant. Add50µllowsaltdigitonin/.25µgsecondary,RT30min; wash2×200µl. Add50µlhighsaltdigitonin/2.5µl20×CUTANApAG-Tn5,RT1h/wash2×200µl. This loaded source complex contains enzyme and oligos but their whole sequences/ends are not inferred from commercial brand or PCR tails.

🟢 Tagmentation50µlbuffer10mMMgCl2/37°C1h, **remove this supernatant**. Retained material is resuspended in50µlReleaseBuffer, which literally lists10mMTAPS8.5, **0.2mMEDTA and22.5mMEDTA**, .5%SDS and1µl20mgmlproteinaseK. Both EDTA descriptors remain explicit, with intended concentration unassigned.58°C1h/68°C2h, then **collect the new release supernatant** for ZymoD5210/elute21µl. The removed tagmentation supernatant and collected release supernatant are different fractions.

🟢 PCR source adds2µluniversali5, barcodedi7 at10µM with **volume unspecified**, and25µlCUTANAHighFidelity2×master. The one actualmaster supplies polymerase/nucleotide/buffer components; no duplicate fresh reagents are invented. Initial58°C5min/72°C5min/98°C45s;18cycles98°C15s/60°C10s/72°C1min. AMPure size selection has no assigned ratio in this pass; elute15µlQiagenbuffer. SeparateactualQubit/qPCR analytical aliquots, retained production/loadingstock and untransferred stock. NextSeqPE75; complete readlanding/index-cycle/mature-library structures remain unassembled.

🟢 S4 universal i5 sequence wraps within its **same Sequence cell**: first55nt plus`GTGTAT` gives61nt. The final6nt are neither a separate barcode-column value nor an asserted publication error. Whole source order:

```
AATGATACGGCGACCACCGAGATCTACACTAGATCGCTCGTCGGCAGCGTCAGATGTGTAT
```

🟡 Its computed parts are canonicalP5 + fixedi5`TAGATCGC` + NexteraS5 + first10ntME. Each i7 isP7 + selected8ntsampleindex + NexteraS7 + first7ntME. These are local sequence decompositions, not proof of complete loaded-Tn5 oligos or an inherited vendor recipe. S4 contains12unique i7 orderings; printedi7_11 andi7_12 rows repeat identically. The model permits source-listed indexes without assigning a study-specific member; N8 means unassigned sampleindex, not degenerate synthesis or UMI. A narrow universal-i5 text-guard exception covers extraction of the wrapped cell, backed by the complete printed source fixture and actualp1 image; shared parser unchanged.

Full reagent factories, exact loaded pAG-Tn5 structures, occupied contact/junction/library products, final primer landing and cleanup yields remain unassembled. This capped private partial model is not a complete source/chemical review or accepted missing-information audit.
