# FICC-Seq — endogenous TRMT2A selected preparation

[FICC-Seq: a method for enzyme-specified profiling of methyl-5-uridine in cellular RNA](https://doi.org/10.1093/nar/gkz658), published NAR2019. Inspected published PMC6821191 main identity, introduction/Fig1 mechanism, Methods FICC-Seq and computation, plus the original SupplementaryMethods DOCX. Two originals; complete current-source/absence review is not claimed.

## Physical principle and scope

🟢 HEK293 or HAP1 cultures receive100µM5-fluorouracil for24h. Its FUTP metabolite competes with uridine incorporation into nascent RNA. Endogenous TRMT2A methylation proceeds through a covalent enzyme–base intermediate;5FU prevents normal release. This is a chemical catalytic trap, not UV crosslinking, cytidine5aza incorporation or an imported mutant-miCLIP construct. Select HEK293/endogenousTRMT2A. The source compares HAP1/iCLIP/mutant-miCLIP and LCMS controls separately; those are not additional counted methods or silently merged preparations. Covalent protein/RNA adduct classes and unreacted/off-target molecules remain unresolved.

🟢 The published supplement gives the full selected library recipe. Lysis0.8ml/4°C30min→2µlTurboDNase+2µlRNaseI1:200/37°C3min/ice5min→clarify13000rpm15min4°C. Separately prepare100µlProteinG beads with5µgantiTRMT2A OTI1C3/0.5mlbuffer30min4°C, wash/resuspend in lysate/IP2h4°C. Wash5×lysis/2×PNK. ATP-free listed20µlPNKpH6.5 setup healsRNA3ends37°C20min; wash3×PNK. AddpreAppL3 withT4RNA ligase/PEG40016°Covernight. The listed ligation components total12µl, without a further water addition supplied here; no normalizing volume is invented.

🟢 Afterwash, addfreshPNK/gamma32P-ATP37°C5min, remove supernatant, heat20µlNuPAGEbuffer70°C5min. Run4–12%BisTris/MOPS, transfer to nitrocellulose, exposefilm and excise radioactiveRNP membrane regions. PK200µg/200µl37°C20min precedes addition200µl7Mureamix/further20min (3.5M finalurea). Phenol/chloroform400µl, retainaqueousphase;1µlGlycoBlue/40µl3MNaOAc5.5/1mlEtOH/−20°Covernight→pellet/wash/resuspend5µl. The RNA may retain a crosslinked peptide remnant; it is not asserted to be completely unmodified naked RNA.

🟢 Add1µlRTprimer and1µldNTP to the actual RNA before65°C5min/ice5min. AddSSII/buffer/DTT/RNaseinhibitor/water, then25°C5min/42°C20min/50°C40min/80°C5min. The original external RT primer and dNTPs are retained from this annealing mixture rather than freshly added twice. Ethanolprecipitation with1µlGlycoBlue; resuspend6µl, add2µlTBEureabuffer/70°C3min, run6%TBEurea40min180V and select70–150ntcDNA. Recovergel-elutedflowthrough, add1µlglycogen/NaOAc/EtOH andresuspend6.5µl.

🟢 CircLigaseII0.3µl/buffer0.8µl/Mn0.4µl yields8µlreaction,60°C1h. Add26µlwater/3µlFastDigestbuffer/1µl10µMBamHIoligo;95°C2min then70cycles95°Cdown1°C/min. Add2µlFastDigestBamHI37°C30min. GlycoBlue/NaOAc/EtOH recovery/resuspend20µl. Select10µlforPCR and retain the other stock;1µleach10µMP5/P3/20µlAccuprimeSupermix/8µlwater,94°C2min then25cycles94°C15s/65°C30s/68°C30s,68°C3min. Run6%TBEgel, SYBRGreenI staining, select150–200nt amplified library, elute/filter,1µlglycogen/NaOAc/EtOH, resuspend20µl. Final sourceHiSeq4000singleend50cycles; loading allocation and full read landing remain unassembled.

## Exact source oligos and limits

🟢 L3-App: `rAppAGATCGGAAGAGCGGTTCAG/ddC/`. Preadenylation is an activated5prime end; terminalddC is an additional physical C nucleotide, not an omitted nonnucleotide label.

🟢 RT: `5Phos/NNXXXXNNNAGATCGGAAGAGCGTCGTGGATCCTGAACCGC`. N positions are the five random UMI bases, split2+3; X is a fixed4nt sample-demultiplexing code with actual sequence unassigned. The fixed adapters and BamHIsite are local source components; complete occupied cDNA/circle products remain unassembled.

🟢 BamHIoligo: `GTTCAGGATCCACGACGCTCTTCAAAA`.

🟢 P5Solexa: `AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT`.

🟢 P3Solexa: `CAAGCAGAAGACGGCATACGAGATCGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT`.

The original DOCX XML was checked for all five whole orderings and modifications. No oligos from separate iCLIP/miCLIP protocols were inherited. Shared canonical sequence bodies are used only where exact local whole-ordering fixtures match.

🟢 The supplement twice prints **400mlTE** forgelcrushing/extraction, followed by collection into a1.5mltube; other workup volumes are inµl. The literal400ml wording is retained as an apparent unit inconsistency in this original, with executed extraction volume unassigned. No silent400µl repair or complete source-error audit is claimed. The listed12µlligation mixture is likewise not silently normalized.

The results workbook and supplement figures are ordinary unread originals under the capped pass; linked iCLIP/miCLIP/vendor originals were not followed. Unknown whole RNA/cDNA molecules, end distributions, trapped-adduct geometry, conversion/yield and exact read landing are unassembled, not confirmed absent. The partial graph is private and does not pass the missing-information audit gate.
