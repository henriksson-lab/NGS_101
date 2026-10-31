# fastGRO — standard HPDP-biotin nascent RNA preparation

[Rapid and Scalable Profiling of Nascent RNA with fastGRO](https://doi.org/10.1016/j.celrep.2020.108373), Cell Reports33,108373(10November2020). Three originals were inspected: associated author manuscriptNIHMS1645786; the exact cited standard protocol export; and the published publisher-layout article plus SupplementalInformation, cached as **NIHMS1645786-supplement-2.pdf**. The latter's printed title/DOI/date/copyright and actual MethodsPDFp16/e3 were verified; its published Methods have priority. The external manifest now identifies this actual document rather than implying it is only a supplement. The blocked ScienceDirect duplicate is not a remaining document request.

## Distinct physical preparation

🟢 Standard fastGRO performs a nuclear run-on with4-thio-UTP, extracts RNA, fragments it mechanically, chemically biotinylates thiol-containing RNA with HPDP-biotin, captures it on M280 streptavidin beads, and releases enriched RNA with DTT before directional RNA-library preparation. This differs physically from biotin-NTP PRO-seq and from4sU base-conversion assays; speed/low input alone is not the admission boundary. The selected model is the **standard HPDP** route with a later prepared RNA spike. Low-input fastGRO-LI uses MTS-biotin and a different workup; it is not silently merged or counted separately.

🟢 Published MethodsPDFp16/e3:5–20million cells→two coldPBS washes→swelling5minice in10mMTris7.5/2mMMg/3mMCa/2UmlSUPERasein→swellingbuffer+10%glycerol wash→1%Igepal lysis with the same salts/glycerol/protection→two lysis-buffer washes→freezing buffer40%glycerol/5mMMg/printed“0.1mM0.5MEDTA”/50mMTris8.3,2×10^7nuclei per100µl→dryice/−80°C→thawice. The EDTA wording is retained without inferring a new dose.

🟢 Add an **equal volume** of prewarmed NRO buffer:10mMTris8/5mMMg/300mMKCl/1mMDTT/500µMATP/500µMGTP/500µM4-thio-UTP/2µMCTP/200UmlSUPERasein/1%Sarkosyl. Incubate7min30°C. These listed concentrations describe the added buffer, not automatically the final equal-volume mixture. The catalyst is resident engaged nuclear polymerase, not freshly supplied purified enzyme. ATP/GTP/4sUTP/CTP are contributing molecular substrates from one actual buffer stock.

🟢 TRIzolLS extraction/ethanol precipitation/water→QubitRNAHS→select up to150µg RNA into a new tube. Add5–10%prepared labelled Drosophila RNA if nuclei were not spiked. The selected route has no nuclei-spike addition. The source spike factory names **Drosophila S9**, **50mM4sU5minRT**, Direct-zol and2µg frozen aliquots; these literals are not changed to S2 or50µM. Actual spike amount/percentage denominator/number of stock aliquots used remain unassigned. Untransferred source RNA remains an unresolved stock.

🟢 BioruptorUCD200 high1–5cycles30sON/30sOFF; fragmented and unfragmented RNA are compared by TapeStationRNAHS. Actual analytical aliquots and retained preparation are separate in the graph. Biotinylate in **20mMTris7.5/2mMEDTA8/40%DMF/200µgmlHPDP-biotin for2h dark25°C800rpm**. Ethanol precipitate/resuspend. Prepare100µl M280 beads/sample with two washes of two volumes100mMTris7.5/10mMEDTA8/1MNaCl/.1%Tween20 and resuspend one volume. Published Results500–503 explicitly identifies a reversible disulfide bridge between biotin and the thiolated uracil base. DTT reductively cleaves this biotin-adduct linkage and releases thiol-labelled RNA; it is not cleavage of a biotin–streptavidin bond or conversion to unmodified uridine. Full atom mapping/occupancy remain unassembled. Capture15min4°Crotation; remove unbound material; wash three times65°C then three timesRT. Release RNA with100mMDTT. Residual beads/remaining RNA stay unresolved rather than automatically discarded.

🟢 ZymoRNA Clean and Purification with **in-column DNase** follows release; QubitRNAHS precedes NEBNextUltraIIDirectional barcoded RNA-library preparation and **NextSeq500**. Exact kit oligo bodies, sample-index values, internal RT/PCR recipe/order, run lengths and read landing are unassembled. The kit's actual oligo contributor inventory is represented, not hidden as buffer. No analytical RT-qPCR primer is substituted for a library primer: the published appended SupplementaryTable2 explicitly identifies its twelve orders as iPSC/NPC RT-qPCR primers.

## Detailed protocol identity and differences

🟢 The directly cited protocol DOI is **10.17504/protocols.io.bbmgik3w**. The fetched export URL is fastgro-yxmvmxe4ol3p, printed ProtocolIntegerID32136, created22January2020, lastmodified12October2020, top date8October2020/footer7October2020. Matching cited DOI and selected components establishes identity/applicability at that scope; it does not prove complete historical revision equivalence.

The following source values are retained separately. The selected model uses the published article, without silently importing a conflicting protocol value:

| Location | Published articlePDFp16/e3 | Retrieved cited protocol export |
|---|---|---|
| Freezing buffer / nuclei | pH8.3;2×10^7nuclei per100µl | pH8; step8 says2×10^7 per10µl, while step13 later uses100µl nuclei+100µl2×NRO |
| Biotinylation |20mMTris/2mMEDTA/25°C | steps28–29 add100µl buffer containing100mMTris/10mMEDTA into1ml total, implying10mM/1mM;24°C |
| RNA workup | ethanol precipitation summarized | explicit isopropanol/75%ethanol workup; no assumed identical detailed schedule |
| DTT release |100mMDTT, no fixed volume/count stated | steps56–58:100µl addition, then self-reference “Repeatsteps57–58”, intended200µlcombined eluate; exact second addition/order unresolved |

The export also prints **10006×g** for a wash spin; no automatic10000×g repair is made. Its originalPDFp15 establishes on-columnDNase order despite floated text-extraction columns:5µlDNaseI+35µlbuffer,15minRT, thenRNAPrepBuffer/washes. It gives6µl finalwater elution. These are recorded protocol-specific facts, not silently substituted into the selected published recipe. Optional glycogen and storage branches are not mandatory additions to this graph.

## Remaining limits

Whole RNA/engaged-complex/adduct structures, exact4sU occupancy, specificity/recovery, nonlabelled/off-target products and complete kit/read geometry remain unassembled. Superseded preprint material was automatically fetched but excluded from authority. The low-input export and separate supplement1 are ordinary unread originals under the cap. No new vendor or predecessor recipe is inherited solely by kit name. Source-review history records a late first manuscript-identity header counter entry; total unique originals remained three. No complete missing-information audit or accepted absence claim is made.


## Public page scope

The public page follows the selected published standard preparation and displays the source/protocol differences above. Four analytical allocations are typed QC materials. A literal curator-written source inventory covers each external input, operation and sequencing endpoint. Passing that gate establishes reviewed workflow coverage and terminal fates, while the computed molecular relations remain unresolved. Historical kit/index/readout facets are labelled unassigned in the model rather than filled from a modern kit recipe.
