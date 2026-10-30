# PADSE-seq

[Published defining article](https://doi.org/10.1073/pnas.2425000122), *Proximity-activated DNA scanning encoded sequencing for massive access to membrane proteins nanoscale organization* (2025).

🟢 Two originals: publishedPMC12012555 mechanism121–135/165/177–179 and actualpublishedSI Methods133–187/NGSPCR124–131/TableS1originalp28/p29 rendered. Mandatory scraper run; extractedtable columns interleave names/bodies, so originalimages verified. Sourcefixed-cellroute is selected; syntheticDNAorigami, in-vitro fluorescent orSanger validations are not merged into it. Proteinpanel/probedistance/cellline alternatives are branches, not extra methods.

## Physical preparation

🟢 Sourceantibody10µl1mg/ml reacts with.5mMDBCO-sulfo-NHS1hRT andAmicon30kcleanup; add15µl15µMazideDNA4°Coneweek/rotation andAmicon30kcleanup/store4°C. DNA contributes physical barcode/scanning structure; DBCO andazide linkage are retained as modifications, not DNA bases. N₃ in TableS1 is an azide group, not three randomN nucleotides. HER2SPcel1 andHER1/EpCAM/PDL1BPcel1–3 define the selected~30nm sourcebranch; exactposition/ligandoccupancy and absolute reach are not established by a barcode map.

🟢10,000cells percollagenAcoatedcoverslip overnight;4%formaldehyde15minRT. Perwell20µl reaction, blockQuickBlock1hRT; HER2Ada highsalt dilutionstock(literal5MNaCl descriptor, finalconcentrationunassigned) overnight4°C. SourceAssori500nM,1µlspermDNA.5mg/ml,4µlformamide and2µl20xSSC mixture5h37°C blocksSPcel; otherAdas recognizeproximityproteins under sourceconditions. Threewashes withnextstepbuffer for each reaction are retained as physicalwashfractions.

🟢 Activate5UNt.BbvCI1xCutSmart3h37°C; TPcel1µM/spermDNA/formamide/2xSSC5h37°C; thenSDA10UNt.BbvCI/5UKlenow3′→5′exo−/500nMdNTP1xCutSmart3h37°C. Source scanningprobe is reused cyclically; activation/releasedBPprimer populations andfailed/free/backgroundclasses remain mixed, not perfect proximity-specific conversion. Collectactualsupernatant, leaving cell/surfacefraction, then80°C20min inactivation.

🟢 NGS uses same applicablePCR asorigamiNGS:1UHotStartTaq,100nMeachRPcel/FP142s,500nMdNTP;95°C10min,30×95°C15s/60°C30s/72°C20s,72°C5min. Select5ngPCRDNA forSangon sequencinglibraryprep/NGS, retaining untransferredPCRstock. This is a known outsourcedlibraryboundary, not permission to invent Illumina ends/indexes/runmask. Sanger-specifictwo-roundextensionPCR and fluorophore-qPCR recipes are excluded from this NGSroute.

## Inventory and limits

🟢 TableS1 giveswholemodifiedSPcel/BPcel probes,Assori,TPcel,RPcel/FP142s. TPcel contains a randommolecularID region; its intendedidentifierboundary is sourcequalified and not automatically an ordinaryUMI orsequencerindex. Knownlocalcomplementarysubstrands may be computed; whole nicked/SDA/bidirectionallyextendedproducts remain unassembled unless locallyproven. This protein-proximitymethod has DNA reporter libraries; it does not assay native genomicDNA simply because its barcode substrates areDNA. No confirmedmissing or fullchemistryclaim.

Computedwholemodifiedfixed-cellinputorderings: SPcel1 {{= chem.SOURCE_ORDERINGS["sp"] }}; BPcel1 {{= chem.SOURCE_ORDERINGS["her1"] }}; BPcel2 {{= chem.SOURCE_ORDERINGS["epcam"] }}; BPcel3 {{= chem.SOURCE_ORDERINGS["pdl1"] }}; TPcel {{= chem.TP }}; Assori {{= chem.ASSORI }}; FP142s {{= chem.FP }}; RPcel {{= chem.RP }}. The N3prefix is a chemicalazide; TP N30is a source molecularID window, not an inferredsequencerindex orordinaryUMI.
