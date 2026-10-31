# PCL-seq

**High-resolution spatial transcriptomics in fixed tissue using a cost-effective PCL-seq workflow** (2025), [Genome Research](https://doi.org/10.1101/gr.279906.124).

🟢 Published sourcePMC12401048 main69/Fig2 and Methods159–201, plus published Supplemental Figures and Tables originalTable2p10. Two originals; cap closed. The SI firstheader was opened immediately before its counter update; this timing oversight is recorded, without exceeding the cap. Mandatory scraper15candidates inoneSIfile, selected literal rows verified in rendered original. SupplementalDataworkbook is downloaded but uninspected; empty searches would not establish absence.

🟢 PCL-seq performs in-situRT with photocapped primers, blocks synthesis-related uncaged background primers, then exposes5′phosphates by region-selective illumination. A biotinylated partially double-stranded adapter is ligated; stringentformamide washes remove noncovalently associated material before tissue release and streptavidin capture. A second RT/template-switch reaction precedes amplification and single-adapterTn5 library finishing. Main69 explicitly calls the TSO adapted fromSMART-seq2; that supports a component link, not inheritance of its entire reaction recipe or TSO modification pattern.

🟢 Selectedfresh-frozen route:10µmcryosection,4%formaldehyde10minRT,1.25Mglycine5min, .5%Triton3min, .1NHCl5min,1MTris8.0neutralization10minRT. FFPEdeparaffinization/collagenase/heat/pepsin and cell-coverslip conditions are separate alternatives. The model does not mix those pretreatments into the fresh-frozen branch.

🟢 RTprimer hybridization uses2×SSC/10%formamide/**printed50%dextransulfate**/1%RNaseinhibitor30minRT, followed by3×2×SSC/1×RTbuffer washes. MaximaHMinus25U/µl/100µMdNTP/1%RNaseinhibitor/1×RTbuffer,42°C90min. The same retained primer is tracked as a produced, same-tube component after washing and consumed as an external primer; no second stock addition or intramolecular self-primer is fabricated. SourceTable2 orderings5′→3′:

```
RT-PC: /PC linker/AGGCCAGAGCATTCGNNNNNNNVTTTTTTTTTTTTTTTTTTTTTTTTVN
Linker: CGAAUGCUCUGGCCUCTCAAGCACGTGGAT
Block: ATCCACGTGCTTGAG
Adapter1: /Biotin/AAGCAGTGGTATCAACGCAGAGTGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTAAACCGGATCCACGTGCTTGAG
Adapter2: /Biotin/AAGCAGTGGTATCAACGCAGAGTGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTCAAGCCTATCCACGTGCTTGAG
TSO: AAGCAGTGGTATCAACGCAGAGTGAATrGrGrG
ISPCR: AAGCAGTGGTATCAACGCAGAGT
```

RT-PC has7NsourceUMI, adjacentV,24T andterminalVN. The adjacentV is not silently declared an extraUMIbase. The mixedU/TLinker ordering is retained; sugar chemistry of all other positions is unassigned rather than converting everyUtoT or assuming uniformRNA/DNA. Table2 colors adapter1AAACCGG andadapter2CAAGCCT as7ntsample/replicatecodes. Those inline codes differ from Illumina index-read sequences. The selected graph usesadapter1 and retainsadapter2 as a verified alternative inventory. ExactPCchemicalordering/cleavage-junction product is unassembled; mainFig2 establishes exposed5′phosphate. The local adapter/splint interface is source-known: the first15Linkerbases retain literalU/T and pair with the primer-derived15ntcDNAhandle, and the last15 pair with the adapterBLOCK tail. A checkedScene places the full5′Biotinadapter and cDNA5′phosphate at an adjacentnick, with a typed splinted3′adapter end and5′substrate end. It does not assign complete cDNA, detailed linker sugar partition, sealed product or ligation yield. Background blocking geometry remains unassembled.

🟢 Backgroundblocking:10µMlinker+10µMblock hybridize15minRT,2×SSC/1×T4buffer wash, T4DNAligase20U/µl30min37°C.65%formamidewash removes noncovalent partners. Region illumination uses100%UVpower10sec with face-upsection on Nikon/AndorMosaic DMD; actualPBS washing removes cleavedcap. Then1µMlinker+1µMbiotinadapter hybridize15minRT,2×SSC/T4buffer wash, ligase20U/µl30min37°C; **65%formamide wash precedes tissue release**. Unilluminated/background/incomplete species are not assumed perfectly absent.

🟢 Initialrelease:50mMTris/1mMEDTA/.2MNaCl/ProteinaseK1mgml55°C2min. Collect initiallysate and recover3additionalslidewashsolutions from the same remaining section, using fresh additions of the same ProteinaseK-containing buffer. The graph derives those washes from the retained source section, rather than introducing them as a new external DNA stock; residual washed-slide material remains unresolved. Pool the collections, addSDSfinal1%,55°Covernight. ProteinaseK is already carried into this overnight mixture; it is not added anew at that incubation. .6×AMPure, thenMyOneStreptavidinC1 capture:10mMTris7.5/1mMEDTA/1MNaCl,1600rpm30minRT. Physically remove unretained fraction; washonce1×SSC/.1%SDS15minRT and3×.1×SSC/.1%SDS10min65°C. Washes do not establish perfect binding specificity.

🟢 Continuedbead-boundcDNAextension/template-switching:1×RTbuffer/1mMdNTP/1µMTSO/1%RNaseinhibitor/MaximaHMinus2U/µl,42°C90min. The source calls this a secondRT reaction, but the graph consumes retained nascentcDNA/RNA context with a newTSO acceptor and RT catalyst; the TSO is not reinterpreted as a new RNA-priming oligo. ThenKAPAHIFIHotStartKB2500/.4µMISPCR/**printed.3MdNTP**/.2mMMgCl2/1U/µlHiFi.98°C3min,23cycles98°C20s/67°C10s/72°C6min,72°C5min. .6×AMPure/elute10µl10mMTris8/Qubit. No normalization of the literaldNTPunits is made.

🟢 DilutecDNA1ng/µl;Tn5A55°C10min, quench.2%SDS,1×AMPure. The source-prepared loadedTn5A is represented as **contributing components of one stock**, not addition of free oligos to tissue or a newly invented transposome assembly. SourceTable2DNA inventory is computed as a pairedMEduplex, with the short strand5′phosphate; complete activeprotein/adaptergeometry and stoichiometry remain unassembled. Literal source orderings:

```
Tn5A long: TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG
Tn5 short: /PHO/CTGTCTCTTATACACATCT
P7xx: CAAGCAGAAGACGGCATACGAGAT[6nt]GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT
S5xx: AATGATACGGCGACCACCGAGATCTACAC[8nt]TCGTCGGCAGCGTC
```

FinalKAPAHIFI KK2602/.2µMP7/.2µMS5/**printed.3MdNTP**/.2mMMgCl2/1UµlHiFi,72°C5min then15cycles95°C20s/55°C30s/72°C30s,72°C5min. AMPureratio/finalelution unassigned. NovaSeqatNovogene; actual run/index masks and primer landings unresolved.

One narrow source-text exemption coversadapter1: Table2 extraction interleaves the row name between the handle and wrapped redbarcode/linker continuation. Wholebiotinadapter fixture was checked againstoriginalp10. It does not waive source transcription. Open: actual sample-index assignments, fullPC/blocked/ligated/RT/TS/Tn5/PCRproduct geometry and relevantlinked originals. The schematic is source-qualified and molecularly partial; it does not establish full validation or accepted missing information.

## Public preparation scope

The public schematic covers the selected published preparation through its physical sequencing contribution. Its source inventory is a literal curator-written set of reviewed inputs and operations; adding graph targets cannot automatically extend review. Existing scoped source findings remain required and mapped locally, without exclusions or deletion. Full preparation, stock/analytical fractions and chemistry warnings remain visible in a collapsed workflow.

Methods195/199literal0.3MdNTPsand16750percentdextransulfateretainedwithoutmM/currentkitrepair. FullPCchemistry/occupiedjunctions/linkersugarpartitionunknown. Full occupied molecular products, primer landing and unspecified read masks remain unassembled rather than silently supplied from another method or current kit.

Source RT and template-switching supplies expose Maxima, dNTPs, RNase inhibitor and buffer functions. Source KAPA amplification supplies expose polymerase, dNTP and buffer functions without extra additions or a new formulation claim. Both source 0.3 M dNTP statements remain literal; whole stock packaging and concentration convention stay unresolved.
