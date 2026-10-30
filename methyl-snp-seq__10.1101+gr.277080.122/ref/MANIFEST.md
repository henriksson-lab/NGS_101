# Source manifest

Originals live outside the checkout at `/husky/henriksson/NGS101_papers/sources/methyl-snp-seq__10.1101+gr.277080.122`. Re-obtain with `python3 tools/get_sources.py --slugs methyl-snp-seq__10.1101+gr.277080.122`; then mandatory `tools/scrape_primers.py`. No third-party originals or extracted text are checked in.

Selected published originals/sections:

- [gr.277080.122_PMC9808626.xml](https://pmc-oa-opendata.s3.amazonaws.com/PMC9808626.1/PMC9808626.1.xml): PublishedVORPMC9808626 overview46–50/discussion130–138/Methods148–164/readout200–210; linkedcopy andnative5mC readouts separate.
- [gr.277080.122_PMC9808626.pdf](https://pmc-oa-opendata.s3.amazonaws.com/PMC9808626.1/PMC9808626.1.pdf): Samepublishedarticle originalp2088 rendered toverifyprintedEndoVI(M0304) andhairpinsequence; selectedoriginal not a second article.
- [gr.277080.122_supp_supp_gr.277080.122_Supplemental_Materials.pdf](https://pmc-oa-opendata.s3.amazonaws.com/PMC9808626.1/supp_gr.277080.122_Supplemental_Materials.pdf): PublishedSupplementalProtocolsteps1–7, text352–491; SupplementalText233–235; originalp8 rendered verifiesEndoIV(M0304), directadditionPolI/Taq, SIadapter/Abcam identifiers. Commands/figures not exhaustivelyreviewed.
- [gr.277080.122_supp_supp_gr.277080.122_Supplemental_Tables.xlsx](https://pmc-oa-opendata.s3.amazonaws.com/PMC9808626.1/supp_gr.277080.122_Supplemental_Tables.xlsx): PublishedSupplementalTables selectedsheetheader/inventory and samplemetric scope; no claim of all numericalcells reviewed.

Deferred directly named sources:

- Main156/160andSIprotocol NEBE7546/M0367/M0280/M0304/E7536/E7165/E7140/M0597; historicalmanufacturer originals: Actual vendor originals unreviewed; exact kit/version identity andfullhairpin/Yprimerends pending.
- Main156AbcamFastBisulfiteab117127 versusSIp8ab11727: Preserveprintedcatalogueidentifiers; no manufactureridentity resolutionwithin cap.
- Main160TwistHumanCoreExome102025 guide andNEBM0544/E7603: Optionaltargetedsourceknownoperation, internalcapture/primergeometry unreviewed.
- SupplementalProtocol commands linked Read_Processing/README.pdf GitHuboriginal: Computationalhairpin/readcycleassembly unreviewed; sourceoperations known without inventingcyclegeometry.

The source-review JSON retains actual failed-manifest URLs and additional uninspected inventory. A download or scraper match alone does not establish a chemical fact. No accepted missing-information claims.
