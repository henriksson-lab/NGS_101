"""10x Chromium Next GEM Single Cell Multiome ATAC + Gene Expression."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
from batch_ngs import seg
from chemdraw import Construct, Row, feature, revcomp
import illumina as il
import nextera as nx
import seqprimers as sp
from multimodal_spatial import finish_pcr_branch, modality_split
from pcr import PCRGeometry, PCRPrimer
from adapter_ligation import AdapterLigationGeometry
from protocol_graph import (Disposition, InputRole, MaterialKind, ProtocolGraph,
                            SourceClaim, SourceDiscrepancy, StepKind, Terminal)
from protocol_model import molecule, role

TITLE = "10x Chromium Single Cell Multiome ATAC + Gene Expression"
NOTES = "01_multiome.html"
SOURCE = 'Commercial protocol: 10x Genomics <a href="https://www.10xgenomics.com/support/instruments/chromium-x-series/chromium-next-gem-single-cell-multiome-atac-plus-gene-expression-reagent-kits-user-guide">CG000338 Rev G</a>.'
SUMMARY = "Nuclei are tagmented before GEM formation. Inside each GEM, one gel-bead barcode is transferred to both accessible-DNA fragments and poly(A)-derived cDNA, yielding linked ATAC and gene-expression libraries from the same nucleus."
CAVEAT = ("CG000338 Rev G prints both gel-bead structures and both complete final-library "
          "layouts. Variable barcode, UMI and sample-index windows remain placeholders, "
          "but every published constant is represented exactly. Rev G does not print the "
          "transposon, Template Switch Oligo, Pre-Amp Primers, cDNA Primers, Adaptor Oligos or "
          "sample-index primers, nor the mechanism joining the ATAC gel-bead oligo to s5; "
          "those steps are marked INFERRED or unresolved. SI-PCR Primer B is assigned from the "
          "same part number printed in CG000168 Rev A.")
CELL_BARCODE = feature("cell_barcode", "cell_barcode", "whitelist",
                       whitelist="10x-multiome-gem-barcode")
UMI_FEATURE = feature("umi", "umi", "random")
I7_FEATURE = feature("sample_index_i7", "sample_index", "unknown")

GEX = Construct([
    seg("P5",il.P5,"p5"),
    seg("i5 sample index","J"*10,"cbc",placeholder=True,
        feature=feature("sample_index_i5","sample_index","unknown")),
    seg("TruSeq Read 1",il.TRUSEQ_READ1,"r1"),
    seg("cell barcode","B"*16,"cbc",placeholder=True,feature=CELL_BARCODE),
    seg("UMI","U"*12,"umi",placeholder=True,feature=UMI_FEATURE),
    seg("poly(dT)VN","T"*30+"VN",placeholder=True),
    seg("cDNA insert","X"*36,placeholder=True),
    seg("Read 2 arm reverse complement",revcomp(il.TRUSEQ_READ2),"r2"),
    seg("i7 sample index reverse complement","I"*10,"cbc",placeholder=True,
        feature=I7_FEATURE),
    seg("P7 reverse complement",il.P7_RC,"p7")],
    name="Multiome gene-expression library")
GEX_PRIMERS=tuple(sp.TRUSEQ[k] for k in ("R1","I1","I2","R2"))
_gex_problems=sp.verify(GEX,GEX_PRIMERS,required_roles=tuple(p.role for p in GEX_PRIMERS))
if _gex_problems:
    raise ValueError("Multiome gene-expression library: "+"; ".join(_gex_problems))
ATAC = Construct([
    seg("P5", il.P5, "p5"),
    seg("10x cell barcode", "B" * 16, "cbc", placeholder=True,
        feature=CELL_BARCODE),
    seg("10x spacer","CGCGTCTG"), seg("S5", nx.S5, "s5"),
    seg("left mosaic end", nx.ME, "me"),
    seg("accessible genomic DNA", "X" * 34, placeholder=True),
    seg("right mosaic end reverse complement", nx.ME_RC, "me"),
    seg("S7 reverse complement", nx.S7_RC, "s7"),
    seg("i7 sample index", "I" * 8, "cbc", placeholder=True, feature=I7_FEATURE),
    seg("P7 reverse complement", il.P7_RC, "p7")], name="Multiome ATAC library")
ATAC_PRIMERS = tuple(sp.NEXTERA[k] for k in ("R1", "I1", "I2", "R2"))
_atac_problems = sp.verify(ATAC, ATAC_PRIMERS,
                           required_roles=tuple(p.role for p in ATAC_PRIMERS))
if _atac_problems:
    raise ValueError("Multiome ATAC library: " + "; ".join(_atac_problems))

FINAL_LIBRARIES = (
    ("Final gene-expression library", GEX, GEX_PRIMERS,
     "Source-exact CG000338 Rev G layout. Read 1 reports the GEM cell barcode and UMI.",
     "The declared Illumina primers are located on the final gene-expression duplex."),
    ("Final ATAC library", ATAC, ATAC_PRIMERS,
     "Source-exact CG000338 Rev G layout. Paired genomic reads flank the accessible fragment; the 24-cycle Index 2 read covers the 8-base spacer and the 16-base GEM barcode, and Index 1 reads the sample index.",
     "The declared Nextera sequencing primers are located on the final ATAC duplex."),
)

def oligos():
    return [
        "GEX gel-bead oligo: 5′-CTACACGACGCTCTTCCGATCT-[16-bp cell barcode]-[12-bp UMI]-T30VN-3′",
        "ATAC gel-bead oligo: 5′-AATGATACGGCGACCACCGAGATCTACAC-[16-bp cell barcode]-CGCGTCTG-3′",
        "Final GEX shell: P5-[10-bp i5]-TruSeq Read 1-[cell barcode]-[UMI]-T30VN-cDNA-Index 1/Read 2 arm-[10-bp i7]-P7′",
        "Final ATAC shell: P5-[cell barcode]-CGCGTCTG-s5-ME-insert-ME′-s7′-[8-bp i7]-P7′",
    ]

def sections():
    return [
        ("Tagment accessible chromatin before partitioning",
         [Row(chunks=[("nucleus chromatin → Tn5 cleavage + adapter transfer → tagged open-DNA fragments", "me", False)])],
         "Transposition (ATAC Buffer B + ATAC Enzyme B, 37 °C 60 min) occurs in bulk nuclei; intact nuclei retain both tagged chromatin and RNA."),
        ("Transfer one GEM barcode to RNA and ATAC molecules",
         [Row(chunks=[("poly(A) RNA + barcoded oligo-dT + Template Switch Oligo → cell barcode + UMI + cDNA", "cbc", False)]),
          Row(chunks=[("tagged genomic DNA + P5/barcode/spacer gel-bead oligo → cell barcode + ATAC fragment", "cbc", False)])],
         "Both reactions occur in the same GEM incubation (37 °C 45 min, 25 °C 30 min, then Quenching Agent) and therefore use the same published 16-base barcode; CG000338 prints each gel-bead structure but not the joining mechanism for the ATAC oligo."),
        ("Pre-amplify both modalities together, then split", modality_split("gene-expression cDNA library", "ATAC DNA library"),
         "After Dynabeads and SPRIselect cleanup, one 7-cycle Pre-Amplification PCR fills gaps and amplifies both products; aliquots then go to ATAC Sample Index PCR and to cDNA amplification."),
    ]

GEX_BEAD = Construct([
    seg("partial TruSeq Read 1", il.TRUSEQ_READ1[-22:], "r1"),
    seg("cell barcode", "B" * 16, "cbc", placeholder=True, feature=CELL_BARCODE),
    seg("UMI", "U" * 12, "umi", placeholder=True, feature=UMI_FEATURE),
    seg("poly(dT)VN", "T" * 30 + "VN", placeholder=True)], name="Single Cell Multiome Gel Beads A GEX oligo")
ATAC_BEAD = Construct([
    seg("P5", il.P5, "p5"),
    seg("10x cell barcode", "B" * 16, "cbc", placeholder=True, feature=CELL_BARCODE),
    seg("10x spacer", "CGCGTCTG")], name="Single Cell Multiome Gel Beads A ATAC oligo")
# CG000338 Rev G Appendix p. 80 prints these two gel-bead oligos verbatim.
assert GEX_BEAD.get("partial TruSeq Read 1").top == "CTACACGACGCTCTTCCGATCT"
# Rev G Sequencing tables: ATAC 50/8/24/49 (i5 24 = 8-nt spacer + 16-nt barcode;
# custom 8-dark-cycle recipe on NextSeq 500/550); GEX 28/10/10/90.
ATAC_READ_LENGTHS = {"Read 1": 50, "Index 1 (i7)": 8, "Index 2 (i5)": 24, "Read 2": 49}
GEX_READ_LENGTHS = {"Read 1": 28, "Index 1 (i7)": 10, "Index 2 (i5)": 10, "Read 2": 90}
# SI-PCR Primer B PN-2000128 is listed in Multiome Amp Kit A (Rev G p. 6) without bases;
# the same part number is printed as AATGATACGGCGACCACCGAGA in CG000168 Rev A Appendix
# (10x scATAC v1). Catalog-linked assignment, not a Rev G reprint.
SI_PCR_PRIMER_B = il.P5[:22]

def protocol_graph():
    g=ProtocolGraph(TITLE,authority=SOURCE)
    nucleus=g.input("nucleus","nucleus containing accessible chromatin and poly(A) RNA",MaterialKind.BIOLOGICAL,state=molecule("intact nucleus",(role("accessible chromatin"),role("poly(A) RNA"))),primary=True)
    transposome=g.input("transposome","ATAC Enzyme B PN-2000265 (loaded transposome)",MaterialKind.ADAPTER,inferred=True,source=SOURCE+" (transposon oligos not printed; Read 1N/Read 2N arms appear only in the final ATAC layout)")
    tn5=g.input("tn5","Tn5 transposase",MaterialKind.CATALYST)
    tagged=g.product("tagged_nucleus","nucleus retaining tagmented chromatin and RNA",MaterialKind.MOLECULAR,state=molecule("tagged nucleus",(role("adapter-tagged accessible DNA"),role("poly(A) RNA"))))
    waste=g.product("tn5_waste","candidate same-adaptor tagmentation products",MaterialKind.WASTE,inferred=True,source=SOURCE+" (candidate end classes; not described in CG000338)")
    g.react("tagment","Tagment accessible chromatin before partitioning",StepKind.TAGMENTATION,((nucleus,InputRole.TEMPLATE),(transposome,InputRole.STRUCTURAL),(tn5,InputRole.CATALYST)),(tagged,waste))
    g.mark_terminal(Terminal(waste.id,Disposition.UNRESOLVED,"Candidate insertion products stay in the nuclei through partitioning; no physical separation of them is described and their recovery is not assigned."))
    droplets=g.input("droplets","Chip J GEM partitioning condition",MaterialKind.CONDITION)
    occupied=g.product("occupied","occupied GEM products",MaterialKind.MOLECULAR,state=tagged.state)
    failed=g.product("failed","empty or unresolved GEM products",MaterialKind.WASTE)
    g.react("partition","Partition tagged nuclei with barcoded gel beads",StepKind.PARTITION,((tagged,InputRole.TEMPLATE),(droplets,InputRole.CONDITION)),(occupied,failed),changes_molecule=False)
    g.mark_terminal(Terminal(failed.id,Disposition.UNRESOLVED,"Empty or multi-nucleus GEMs are pooled and recovered with the rest; exclusion is computational, not a physical discard."))
    rtprimer=g.input("rt_primer","GEX gel-bead oligo (Gel Beads A PN-2000261)",MaterialKind.OLIGO,state=GEX_BEAD,source=SOURCE+" (Rev G Appendix p. 80)")
    rt=g.input("rt","Barcoding Enzyme Mix PN-2000266 (RT activity)",MaterialKind.CATALYST)
    first=g.product("first_strand","GEM mixture with first-strand cell-barcoded cDNA",MaterialKind.MOLECULAR,state=molecule("GEM products",(role("cell barcode + UMI + first-strand cDNA"),role("adapter-tagged accessible DNA"))))
    g.react("reverse_transcribe","Capture poly(A) RNA and reverse-transcribe in the GEM",StepKind.REVERSE_TRANSCRIPTION,((occupied,InputRole.TEMPLATE),(rtprimer,InputRole.STRUCTURAL),(rt,InputRole.CATALYST)),(first,))
    tso=g.input("tso","Template Switch Oligo PN-3000228",MaterialKind.OLIGO,inferred=True,source=SOURCE+" (in the GEM Master Mix; bases not printed in Rev G)")
    switched=g.product("switched","GEM mixture with template-switched cDNA",MaterialKind.MOLECULAR,state=molecule("GEM products",(role("cell barcode + UMI + cDNA + TSO complement"),role("adapter-tagged accessible DNA"))),inferred=True,source=SOURCE+" (template-switch end not printed)")
    g.react("template_switch","INFERRED — copy the Template Switch Oligo at the cDNA 3′ end",StepKind.TEMPLATE_SWITCHING,((first,InputRole.TEMPLATE),(tso,InputRole.STRUCTURAL),(rt,InputRole.CATALYST)),(switched,))
    atac_bead=g.input("atac_bead","ATAC gel-bead oligo (Gel Beads A PN-2000261)",MaterialKind.OLIGO,state=ATAC_BEAD,source=SOURCE+" (Rev G Appendix p. 80)")
    barcoded=g.product("barcoded_gem","GEM mixture with barcoded cDNA and barcoded ATAC fragments",MaterialKind.MOLECULAR,state=molecule("GEM products",(role("cell-barcoded cDNA"),role("P5 + cell barcode + spacer joined to s5 fragments"))))
    g.react("barcode_atac","INFERRED — join the P5/barcode/spacer oligo to the s5 end of tagged fragments",StepKind.MAGIC,((switched,InputRole.TEMPLATE),(atac_bead,InputRole.STRUCTURAL)),(barcoded,),reason="CG000338 Rev G states only that GEM incubation produces 10x Barcoded DNA from transposed DNA; the final library places the CGCGTCTG spacer directly 5′ of s5, but the joining reagent (any splint) and its geometry are not disclosed.")
    quench=g.input("quench","Quenching Agent PN-2000269; Dynabeads MyOne SILANE and SPRIselect cleanup",MaterialKind.CONDITION)
    recovered=g.product("recovered","recovered barcoded GEM products",MaterialKind.MOLECULAR,state=barcoded.state)
    gem_waste=g.product("gem_waste","partitioning oil, leftover reagents and primers removed by bead cleanup",MaterialKind.WASTE)
    g.react("recover","Quench, break GEMs and purify (Dynabeads, SPRIselect)",StepKind.CLEANUP,((barcoded,InputRole.TEMPLATE),(quench,InputRole.CONDITION)),(recovered,gem_waste),changes_molecule=False)
    g.discard(gem_waste,"Rev G Step 3: silane beads purify products from leftover biochemical reagents and primers; supernatants are removed.")
    pre_f=g.input("preamp_primers","Pre-Amp Primers PN-2000271",MaterialKind.OLIGO,inferred=True,source=SOURCE+" (bases not printed)")
    pre_r=g.input("preamp_primers_2","Pre-Amp Primers PN-2000271 (second end)",MaterialKind.OLIGO,inferred=True,source=SOURCE+" (bases not printed; one mix, primer count unknown)")
    pre_pol=g.input("preamp_mix","Pre-Amp Mix PN-2000270",MaterialKind.CATALYST)
    preamp=g.product("preamp","pre-amplified ATAC and cDNA products",MaterialKind.MOLECULAR,state=recovered.state)
    g.react("preamp_pcr","Pre-Amplification PCR: 72 °C 5 min gap fill, then 7 cycles",StepKind.PCR,((recovered,InputRole.TEMPLATE),(pre_f,InputRole.STRUCTURAL),(pre_r,InputRole.STRUCTURAL),(pre_pol,InputRole.CATALYST)),(preamp,),geometry=PCRGeometry.unavailable(ATAC,reason="Rev G names Pre-Amp Primers PN-2000271 but does not print them; the guide states that this PCR fills gaps and amplifies both modalities."))
    atac=g.product("atac_branch","40 µl pre-amplified aliquot for ATAC",MaterialKind.MOLECULAR,state=preamp.state)
    rna=g.product("rna_branch","35 µl pre-amplified aliquot for cDNA amplification",MaterialKind.MOLECULAR,state=preamp.state)
    g.react("split","Split the pre-amplified product into ATAC and cDNA aliquots",StepKind.SPLIT,((preamp,InputRole.TEMPLATE),),(atac,rna),changes_molecule=False)
    cf=g.input("cdna_forward","cDNA Primers PN-2000089 (forward)",MaterialKind.OLIGO,inferred=True,source=SOURCE+" (bases not printed in Rev G)")
    cr=g.input("cdna_reverse","cDNA Primers PN-2000089 (reverse)",MaterialKind.OLIGO,inferred=True,source=SOURCE+" (bases not printed in Rev G)")
    camp=g.input("cdna_amp","Amp Mix PN-2000047",MaterialKind.CATALYST)
    cdna=g.product("cdna","amplified cell-barcoded cDNA",MaterialKind.MOLECULAR,state=molecule("Multiome cDNA",(role("cell barcode"),role("UMI"),role("cDNA"))))
    g.react("cdna_amplification","Amplify barcoded cDNA",StepKind.PCR,((rna,InputRole.TEMPLATE),(cf,InputRole.STRUCTURAL),(cr,InputRole.STRUCTURAL),(camp,InputRole.CATALYST)),(cdna,),geometry=PCRGeometry.unavailable(GEX,reason="cDNA Primers PN-2000089 are named but not printed in CG000338 Rev G."))
    frag=g.input("fragmentase","Fragmentation Enzyme/Buffer PN-2000090/2000091",MaterialKind.CATALYST)
    tailed=g.product("tailed","fragmented, end-repaired, A-tailed cDNA",MaterialKind.MOLECULAR)
    g.react("frag_er_at","Fragmentation, end repair and A-tailing (one reaction)",StepKind.FRAGMENTATION,((cdna,InputRole.TEMPLATE),(frag,InputRole.CATALYST)),(tailed,))
    adaptor=g.input("adaptor","Adaptor Oligos PN-2000094",MaterialKind.ADAPTER)
    ligase=g.input("ligase","DNA Ligase PN-220110",MaterialKind.CATALYST)
    adapted=g.product("adapted","adaptor-ligated cDNA fragments",MaterialKind.MOLECULAR)
    g.react("ligate","Adaptor ligation",StepKind.LIGATION,((tailed,InputRole.TEMPLATE),(adaptor,InputRole.STRUCTURAL),(ligase,InputRole.CATALYST)),(adapted,),geometry=AdapterLigationGeometry.unavailable("CG000338 Rev G names Adaptor Oligos PN-2000094 but prints neither its strands nor its ligating ends; the final library alone does not establish reagent geometry."))
    reviewed=["nucleus","transposome","tn5","tagment","droplets","partition","rt_primer","rt","reverse_transcribe","tso","template_switch","atac_bead","barcode_atac","quench","recover","gem_waste","preamp_primers","preamp_primers_2","preamp_mix","preamp_pcr","split","cdna_forward","cdna_reverse","cdna_amp","cdna_amplification","fragmentase","frag_er_at","adaptor","ligase","ligate"]
    gex_geometry=PCRGeometry.exact(
        GEX,
        PCRPrimer("SI-TT i5",il.P5+"A"*10+il.TRUSEQ_READ1,anneal_nt=len(il.TRUSEQ_READ1)),
        PCRPrimer("SI-TT i7",il.P7+"A"*10+il.TRUSEQ_READ2,anneal_nt=len(il.TRUSEQ_READ2)),
        statement="Dual Index Kit TT Set A supplies index pairs, not printed full primers; the Rev G final-library layout fixes the 3′ TruSeq Read 1 and Read 2 landing regions, and representative index windows do not affect annealing.")
    atac_geometry=PCRGeometry.exact(
        ATAC,PCRPrimer("SI-PCR Primer B PN-2000128",SI_PCR_PRIMER_B),
        PCRPrimer("Sample Index N, Set A",nx.n7xx_primer("A"*8),anneal_nt=len(nx.S7)),
        statement="Rev G names SI-PCR Primer B and Single Index Kit N Set A without printing primer bases; Primer B is catalog-linked to its CG000168 Rev A Appendix sequence, and the i7 primer's s7 3′ end is fixed by the Rev G ATAC layout.")
    # Rev G Steps 5.1 and 7.5 name the index reagents but print no primer bases, so
    # only source names are passed; there is no printed state to check.
    finish_pcr_branch(g,"gex",adapted,GEX,GEX_PRIMERS,reviewed,label="gene-expression",geometry=gex_geometry,read_lengths=GEX_READ_LENGTHS,
                      step_label="Sample Index PCR with one Dual Index TT Set A i5/i7 pair (Step 7.5)",
                      forward_label="Dual Index TT Set A i5 sample-index primer (bases not printed)",
                      reverse_label="Dual Index TT Set A i7 sample-index primer (bases not printed)",
                      polymerase_label="Amp Mix PN-2000047/2000131")
    finish_pcr_branch(g,"atac",atac,ATAC,ATAC_PRIMERS,reviewed,label="ATAC",geometry=atac_geometry,read_lengths=ATAC_READ_LENGTHS,
                      step_label="Sample Index PCR with SI-PCR Primer B and a Single Index Kit N Set A i7 index (Step 5.1)",
                      forward_label="SI-PCR Primer B PN-2000128",
                      reverse_label="Single Index Kit N Set A i7 sample-index primer mix (bases not printed)",
                      polymerase_label="Amp Mix PN-2000047/2000103")
    g.annotate(SourceDiscrepancy(
        "gex_library", "CG000338 Rev B omits four P5 bases in its gene-expression layout",
        (SourceClaim("CG000338 Rev B, Appendix p. 80, Gene Expression Library",
                     "The 5′ P5 arm is printed as AATGATACGGCGACCACCGAGATCT before the N10 i5 index."),
         SourceClaim("CG000338 Rev B, Appendix p. 80, ATAC gel bead and library; "
                     "Rev F/G, Appendix pp. 80–81, Gene Expression Library",
                     "The P5 arm is printed as AATGATACGGCGACCACCGAGATCTACAC before the index or barcode.")),
        "The gene-expression construct uses the full P5 arm ending ACAC, as printed in Rev F/G.",
        "Rev B's own ATAC layout and later guide revisions agree on the full P5 arm; "
        "the canonical Illumina P5 sequence agrees as well.",
        "Whether Rev B's shortened gene-expression line was a diagram error or represented "
        "a reagent variant is not stated by 10x."))
    g.declare_reviewed_scope(reviewed,SOURCE)
    return g.require_publishable()
