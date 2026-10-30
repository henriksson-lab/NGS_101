# ProtSeq: Toward high-throughput, single-molecule protein sequencing via amino acid conversion into DNA barcodes

Defining DOI: 10.1016/j.isci.2021.103586. Bounded source pass;no accepted absence or finished public schematic.

Sources and exact selected scopes:

- https://pmc-oa-opendata.s3.amazonaws.com/PMC8717419.1/PMC8717419.1.xml — `j.isci.2021.103586_PMC8717419.xml.txt`: Published VOR Results262–278,STARMethods588–609/643–657/669–673;general DNA1/6cycle route,bridge modification,NGS finish
- https://pmc-oa-opendata.s3.amazonaws.com/PMC8717419.1/mmc1.pdf — `j.isci.2021.103586_supp_mmc1.pdf.txt`: Published SI TableS2 foundation/bridge and TableS4 selected binder identifiers;FigS3 sixDNA-DNAcycle caption;other tables deferred

🟢 Published identity: PMC8717419 manifest author-manuscriptFalse (published iScience original),Epub11December2021/issue21January2022. This is a different method from bulk ProteinSeq2011.

🟢 Scope is **experimentally demonstrated DNA-DNA barcode-cycle sequencing**:sourceResults reports six sequential barcode incorporation cycles on chip with imperfect incorporation/misassignment. Peptide–nanobody and thrombin capture are separately described demonstrations;complete amino-acid/protein sequencing is an intended future platform. The graph contains no peptide degradation step and no deterministic inference of six successful ligations per molecule.

🟢 Source1/6cycle operational order:hybridize FC/RC/barcoded foundation/barcoded target in staged95→25°C →source foundationmix/T4BluntTA ligation toMiSeqchip28°C15min →100%formamidewash40°C90s →surface/DNAblocking37°C15min →eachround binder/bridge incubation,unbound-binder wash,T4ligation,EcoRI+complementarybridge cleavage,andphysical release wash →adapter bridge-ligation →two500µLwater washes →direct modifiedMiSeq acquisition. Full endpoint/junction geometry andsource assignments remain unassembled;no perfect hybridization/specificity/yield predicted.

🟢 EcoRI cuts the source DNA binder body from retained barcode/consensus DNA. Released fragments remain in the reaction until actual washes. It does not delete nucleotides magically or remove an amino acid. Foundation/original binder DNA persists through ligation;NGS reads report a chain of recorded binding events and cycles,not a directly measured protein sequence.

Verbatim selected oligos/modifications:

- TableS2 Fd7: `/5Phos/CGACTGCGAGCTGATGGCCTTGATGATAACG`
- Universal bridge: `CTGCGCCTATAGGAATTCGTTATC/i5NitInd//i5NitInd//i5NitInd//i5NitInd//i5NitInd//i5NitInd//i5NitInd//i5NitInd//i5NitInd//i5NitInd//i5NitInd//i5NitInd/GGACACGGCCGTTATC`
- Universal NGS adapter: `/5Phos/AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT`
- NGS adapter bridge9/5: `TTCCGATCTCGTTA`
- Restriction bridge: `CTGCGCCTATACGAATTCGTTATC`

🟢 Each5NitInd is an actual5-nitroindole universal-base analogue. The local typed bridge retains an ordered12residue inventory and fixed flanks;it is not N12 or12inosines. Its complete drawn molecular/pairing state remains unassembled and is honestly reported as unrepresented by the generic DNA-strand sourceguard. A local whole-modified-notation transcription test independently checks it against the published original. Other exact DNA inputs are sourceguard-checked. Source-foundation identifier is a structured feature barcode,not a cell barcode or UMI.

🟢 Source adapter completion uses 2:1 universalNGSadapter:bridge,BluntTAligase,two30µLloads40°C2min45s;waterwash90s apart. PublishedMiSeqV2 uses a45–600cycle read depending on experiment;the graph leaves the specific choice and primer sites unassigned. PhiX loaded to cartridge is a separate control,not ligated into every chain. The distinct TargetSwitchSELEX177bpPippinPCR library is not merged into this on-chip BCS route.

🟡 CanonicalP5 reverse complement imported only as a sequence component,whole adapter checked verbatim. Chain observations may skip cycles;code records observed binding events rather than filling missing incorporation. No importable predecessor recipe inheritance established in these selected sections.

Deferred targets:

- Whole FC/RC,source target–foundation/binder pair assignments and full per-cycle barcode/ligation/EcoRI geometry;TablesS2–S8 full original scope
- Peptide–nanobody and thrombin physical conjugation branches,POC/modifier manufacturing,independent TargetSwitchSELEX branch and its original oligos
- Specific run45–600cycle choice,full primer landing/read geometry and actual measured six-cycle sequences
- Other source-linked chemistry originals/buffer/manufacturing details deferred;no complete version/source audit

Mandatory standard source retrieval and primer scraper completed. Original pointers/modifications checked only in scopes above;uninspected originals and source hits are not evidence. No two-pass absence audit performed.
