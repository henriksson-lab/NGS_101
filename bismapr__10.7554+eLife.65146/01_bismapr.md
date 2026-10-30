# BisMapR

Defining DOI: https://doi.org/10.7554/eLife.65146.

## Actually inspected sources

🟢 Published PMC7901872 Fig1/Results52–76 and Methods142–156; source-cited2019 MapR frontend inspected separately. Whole-bundle mandatory primer scraper run; hits are reading aids and not independent evidence. Superseded preprints, when automatically fetched, are not read or used.

## Scope and boundaries

Selected2021 BisMapR route omitsRNaseA from source-cited MapR stop buffer. Do not use original MapR RNaseA digestion or infer protocol-wide equality between2019 and2020 MapR versions. C-to-U on exposed DNA and dUTP second-strand marking create distinct physical species; UDG does not physically remove DNA from the tube. Exact enzymes/end cleavage/adapter geometry and complete uracil depletion remain unassembled.

## Readout

Source NextSeq500PE38; Read1 mapping strand reports hybrid-template DNA. Computational strand separation is not physical molecular purification.

The private graph records source-known operations, while unresolved exact molecular geometry is an explicit model boundary. This capped pass is not a complete source review or a missing-information audit.

## Full source-known route

🟢 Methods144 explicitly uses MapR (Yan2020 andYan2019) on5million cells **withoutRNaseA in the stop buffer**. The local frontend is independently declared from source-cited2019operations with version-specific quantities left qualified; it does not import theMapR recipe. Native RNA/DNA hybrid RNA persists through release/extraction for subsequent protection. TheMapR-only/MNase control andRNaseA/NBSS controls are separate controls, not silently substituted into the selected chemistry.

🟢 Tenµlsample+10µlwater+130µlZymoD5005CTreagent at25°C3h performs non-denaturing bisulfite conversion. Exposed displaced-DNA C-to-U is a selected chemical class; hybrid-protected DNA remains distinct. Kit desulfonation/column purification yields20µlM-Elution. NEBE7760buffer8µl/enzyme4µl/water48µl,16°C1h, performsRNaseH/dUTP second-strand synthesis, followed by1.8×AMPure. RNA digestion products and the new dU-containing DNA are independent physical species.

🟢 Methods152: EnzymaticsEnd-RepairMix→Klenowexo-minus A-tail→MinElute→NEBE7600 adapters/T4DNAligase→AMPure greater150bp→Q5dual-indexPCR withUDG **in that PCRmix**→MinElute→NextSeq500PE38. The graph representsUDG before copying as an ordered view of the coupled PCRmix, not a source-described separately purified reaction. Uracil excision/depletion is not physical disappearance: AP-DNA/possible secondary fragmentation, excised uracil and original templates remain until actual column cleanup. Adapter hairpin processing is source-reported but the exact historical molecules and mechanism are unassembled. No exact tag/UMI/index architecture or perfectUDG specificity/depletion is inferred.

## Scoped ancestry

Explicit `mapr__10.1016+j.celrep.2019.09.052`→`bismapr__10.7554+eLife.65146`: Methods144 directly citesYan2019 andYan2020 (10.1002/cpmb.113) for theMapR frontend, with the explicitRNaseA omission. Results52/Fig1 then adds non-denaturing bisulfite anddUTP/UDG selection. Scope is nuclease-guided nativeR-loop release, not identical strand handling or entire library construction. Both named precursor versions remain distinguishable;2020fulloriginal deferred.
