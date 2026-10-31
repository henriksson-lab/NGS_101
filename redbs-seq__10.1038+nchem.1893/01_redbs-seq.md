# redBS-seq:5fC reduction changes a bisulfite readout

Booth et al. (2014), [Quantitative sequencing of5-formylcytosine in DNA at single-base resolution](https://doi.org/10.1038/nchem.1893).

## Capped sources and study branch

One reader,28October2026; associated accepted manuscript PMC4188980, complete chemical rationale/Methods; publisher Supplementary Information MOESM165, experimental methods pp2–5 and control/primer Table3. Mandatory scraper run; Table3 control-template/primer hits were distinguished from the genomic TruSeq primer mix. The publisher supplement supplies the positive genomic preparation; no later vendor recipe is imported. Complete version-of-record Methods and linked originals remain deferred.

The graph selects the actual **reduced-representation genomic redBS arm and untreated BS comparator**. The study also has an oxBS arm, synthetic-control sequencing and Sanger experiments; their chemistry is not merged into this route.

## Genomic library preparation

1. Extract mouse ES-cell DNA with proteinaseK/RNaseA, phenol/chloroform and ethanol. Digest twelve1µg preparations with MspI overnight37°C, quenchEDTA, GeneJET purify. Add the source synthetic control mixture. MspI reduced representation is an experimental front end, not a compulsory property of every possible redBS application.
2. Perform Illumina TruSeq library preparation with1:10adapters. Following finalAMPure purification, dry30min and elute in water. The study's exact adapter/primer members and kit revision remain unassigned here.
3. Four libraries are reduced, four oxidized and four untreated. For redBS, add5µL fresh1MNaBH4 to15µL DNA:250mM final,1h at room temperature in the dark. Quench slowly with10µL750mMsodium acetate,pH5; wait at least10min/until gas release ends. The quench does not physically purify DNA or imply complete chemical conversion.
4. redBS and untreatedBS arms undergo EpiTect bisulfite:30µLDNA,80µLbisulfitemix,30µLprotectbuffer; twoFFPEthermal cycles and kit workup. Source oxidation arm is not modeled here.
5. KAPA HiFiUracil+ and TruSeq primer mix amplify12cycles (98°C45s;12×98°C15s/60°C30s/72°C30s;72°C60s). AMPурe purification, Qubit/TapeStation, HiSeq2000 **single-end100bases**.

## Chemical product is not the reported letter

| Starting base | redBS reduction | After redBS bisulfite | Expected reported base |
|---|---|---|---|
| C | C | U | T |
|5mC |5mC |5mC | C |
|5hmC |5hmC |CMS | C |
|5fC |5hmC |CMS | C |

UntreatedBS instead converts5fC toU, reportingT. The source explicitly identifies **cytosine-5-methylenesulfonate(CMS)** as the5hmC-bisulfite product; the graph does not label CMS as chemically unchanged cytosine. Difference inC readout between redBS andBS estimates5fC; it is not a single-read unambiguous base identity. Source controls show incomplete conversion/background, including residual5fC readingC in ordinaryBS; idealized classes do not establish quantitative yield or a perfect per-molecule transformation. Source says5caC readout is unchanged between arms; its exact chemical path is not assigned in this model.

## Scope and ancestry

The source explicitly selects RRBS/MspI enrichment and compares with its published oxBS method. These are source-qualified study design/component relations; similarity does not authorize adopting another protocol's library adapters or precise kit revision. Synthetic-control “Illumina” primers in Table3 do not automatically establish genomic TruSeq primer members. Deferred: historical kit originals, complete VOR, control/oxidation branches, exact final library/run primers, full linked-source checklist and independent second reader. No absence claim accepted.

## Publication preparation review

The generated page remains a private draft. A green graph gate and publication marker do not resolve the source/model mismatch: the genomic study has twelve independently prepared MspI/TruSeq libraries, with four allocated to each chemical arm. The current graph instead prepares one collection and splits it into arms. A follow-up must represent independent preparation inputs/reactions, preserve the excluded oxidation scope and avoid borrowing synthetic-control pooling into the genomic route. Source-known NaBH4/acetate, prepared KAPA/primer-mix components and actual QC/stock allocations also need explicit inventory/fate treatment. No source absence or inaccessible-original conclusion follows from this model issue.
