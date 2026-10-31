# hiCLIP

Defining original: https://doi.org/10.1038/nature14280

Sugimoto et al.(2015), *hiCLIP reveals the in vivo atlas of mRNA secondary structures recognized by Staufen1*.

## Source identity and capped scope

🟢 `nature14280_PMC4376666.xml.txt` explicitly identifies an accepted author manuscript (provenance44), not the VOR. Read Introduction46–48 and detailed hiCLIP Methods104–174, sequencing206 and barcode-processing240/290, including source oligos and complete downstream outline. Associated original `Guide.docx` opening source-file inventory and `Note.docx` opening discussion/identity were inspected. The latter contains embedded reference-code noise; it was not a complete linked-source review. Original publisher supplement mirrors arrived later and are deferred under the three-original cap. Final published procedure applicability remains unverified; manuscript detail is used explicitly as manuscript-qualified evidence, not silently promoted to VOR.

Preliminary2026-10-28 fetched files were unread before pause. Resumed2026-10-29 pass inspected actual originals and ran the mandatory full-bundle scraper, checking modified adapter and RT-primer hits in Methods104/164. No absence claim or independent second pass.

## Ordered nucleotide/adapter chemistry

🟢 AdaptorA as written: `/5Phos/rArGrArUrCrGrGrArArGrArGrCrGrGrUrUrCrArG/3ddC/`. Its body is20RNA nucleotides; terminal ddC is an additional blocked cytidine nucleotide, not a non-nucleotide spacer. AdaptorB: `/5Phos/rCrUrGrUrArGrGrCrArCrArUrArCrArArUrG/3Phos/`,19RNA nucleotides with a3′phosphate block, no extra terminal nucleotide. Source Mth RNA ligase/ATP preadenylates both5′ends, followed by PAGE purification.

🟢 After capture and3′dephosphorylation of RNA, source T4RNA ligase2 truncatedK227Q attaches preadenylated A/B donors to RNA3′OH. The source's nominal50% hetero-adapter fraction explicitly assumes perfect/random ligation and is not treated as measured recovery. The local graph follows the selected A/B-bearing class with residual classes unquantified.

T4PNK plus ATP then phosphorylates RNA5′ends and removes B's3′phosphate; it does not restore an OH to A's ddC. After physical washes, source T4RNA ligase1 joins B3′OH to the other RNA5′phosphate. Conditional order is RNAarm1–B–RNAarm2–A. The constructor requires B already attached at5′and deblocked at3′; it preserves both RNA arms, the bridging linker and terminal A ddC through joining. Exact3D stem geometry and fractions remain unassembled. Free adapters/material remain in the reaction until actual source bead washes.

## RT, identifiers and remaining boundary

🟢 High-RNase RT primer `/5Phos/NNTATTNNNAGATCGGAAGAGCGTCGTGGATCCTGAACCGC`; other local sample barcodes AACC, CGCC, GTGG. Methods240 identifies random read positions1–2 and7–9 versus experimental positions3–6; Methods290 uses matching random barcodes to remove PCR duplicates. Accordingly the two random windows form a composite identifier with distinct stable partIDs and a shared group; the intervening4bases are a separate inline sample barcode. A source-qualified computed alignment shows the RT3′nine-base seed pairing inside A's RNA body before the terminal ddC; the5′primer tail is unpaired at the adjacent ddC. This is a local alignment relation, not a general RT-binding/efficiency predictor; terminal ddC is retained on the original RNA template.

🟢 Source Methods160–174: protease release from selected RNP membrane region; SSIII; frequent crosslink-associated truncation; NaOH RNA hydrolysis; ethanol and100–300nt cDNA selection; CircLigaseII60°C1h; Cut_oligo `GTTCAGGATCCACGACGCTCTTCAAAA`/BamHI reopening; source P5/P3 PCR primers and24cycles. These downstream steps are recorded but not assembled into the prefix graph. No iCLIP module is imported merely because the circle/BamHI design resembles it. The ribosome-profiling/mRNA-seq branches in this paper do not supply missing hiCLIP chemistry.

🟢 Methods206 reports MiSeq ReagentKitv2(300cycle), acquired196cycles for hiCLIP. The50cycle HiSeq statement applies to separate ribosome/mRNA-seq branches; it is not substituted here.

Graph boundary: upstream capture/workup, stringent purification, proteolysis/RT/hydrolysis, cDNA circle/cut/PCR and final sequencing geometry remain unassembled, alongside VOR procedural reconciliation. No finished schematic.

## Component lineage

Introduction46 explicitly compares direct-arm CLASH ligation with hiCLIP's inserted bridge strategy: candidate CLASH→hiCLIP is a strategy refinement, not complete recipe inheritance. Main Methods calls denaturing purification adapted from Kiel etal.ref31 and iCLIP protocolsref19/32; exact editions/components require independent predecessor review before an import/registry claim. Shared RT motifs or matching primers alone do not establish ancestry.

## Bounded unlinked-page adjudication: source/model scope retained

The 30 October 2026 read-only pass revisited only the existing associated
accepted manuscript's detailed Methods 104–174 and sequencing paragraph 206.
No additional original, publisher retry, vendor source or reference was opened.
The required `unreviewed_final_vor` source item remains unresolved; this pass
did not establish a new version-of-record access failure. Its applicability
is unperformed source/version scope, not a claim of missing disclosure.

The present six-step adapter-gating graph still omits the actual adapter factory,
cell/capture preparation, urea denaturation and second immunoprecipitation,
20% radiolabel allocation/recombination with the other 80%, membrane selection,
proteinase extraction, RT/hydrolysis, cDNA size selection, circularization,
BamHI and PCR producer path. Those actual preparation and physical-fate
obligations cannot be excluded merely to qualify a coarse RNA endpoint.

The second-ligation paragraph states a 20 µl added mixture but lists component
volumes totaling 100 µl. A possible bulk-stock versus aliquot interpretation
remains unassigned here; neither normalization nor a publication-error claim is
justified by that referent alone. The later source 11 µl cDNA stock contributes
5 µl to PCR, with 6 µl nominal untransferred stock. The 24-cycle hiCLIP PCR and
MiSeq 196-cycle report must remain distinct from the separate ribosome/mRNA
protocol's modifications and cycles. The model is unchanged and remains a
private source-qualified prefix, with no accepted absence claim.
