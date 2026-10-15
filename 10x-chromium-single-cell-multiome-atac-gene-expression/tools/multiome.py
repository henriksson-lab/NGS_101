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
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SourceClaim,
                            SourceDiscrepancy, StepKind)
from protocol_model import molecule, role

TITLE = "10x Chromium Single Cell Multiome ATAC + Gene Expression"
NOTES = "01_multiome.html"
SOURCE = 'Commercial protocol: 10x Genomics <a href="https://www.10xgenomics.com/support/instruments/chromium-x-series/chromium-next-gem-single-cell-multiome-atac-plus-gene-expression-reagent-kits-user-guide">CG000338 Rev G</a>.'
SUMMARY = "Nuclei are tagmented before GEM formation. Inside each GEM, one gel-bead barcode is transferred to both accessible-DNA fragments and poly(A)-derived cDNA, yielding linked ATAC and gene-expression libraries from the same nucleus."
CAVEAT = ("CG000338 Rev G prints both gel-bead structures and both complete final-library "
          "layouts. Variable barcode, UMI and sample-index windows remain placeholders, "
          "but every published constant is represented exactly. The internal bridge reagent "
          "and some reaction-mixture oligos remain proprietary.")
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
     "Source-exact CG000338 Rev G layout. Paired genomic reads flank the accessible fragment; Index 2 reads the 16-base GEM barcode and Index 1 reads the sample index.",
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
         "Transposition occurs in bulk nuclei; intact nuclei retain both tagged chromatin and RNA."),
        ("Transfer one GEM barcode to RNA and ATAC molecules",
         [Row(chunks=[("poly(A) RNA + barcoded oligo-dT → cell barcode + UMI + cDNA", "cbc", False)]),
          Row(chunks=[("tagged genomic DNA + barcoded bridge → cell barcode + ATAC fragment", "cbc", False)])],
         "Both reactions occur inside the same GEM and therefore use the same published 16-base barcode; CG000338 prints each gel-bead structure."),
        ("Separate and amplify the two modalities", modality_split("gene-expression cDNA library", "ATAC DNA library"),
         "After GEM cleanup, modality-specific amplification preserves the common cell identity while producing independently sequenceable libraries."),
    ]

def protocol_graph():
    g=ProtocolGraph(TITLE,authority=SOURCE)
    nucleus=g.input("nucleus","nucleus containing accessible chromatin and poly(A) RNA",MaterialKind.BIOLOGICAL,state=molecule("intact nucleus",(role("accessible chromatin"),role("poly(A) RNA"))),primary=True)
    transposome=g.input("transposome","Multiome loaded transposome",MaterialKind.ADAPTER,inferred=True,source=SOURCE+" (complete oligo sequence not disclosed)")
    tn5=g.input("tn5","Tn5 transposase",MaterialKind.CATALYST)
    tagged=g.product("tagged_nucleus","nucleus retaining tagmented chromatin and RNA",MaterialKind.MOLECULAR,state=molecule("tagged nucleus",(role("adapter-tagged accessible DNA"),role("poly(A) RNA"))))
    waste=g.product("tn5_waste","non-library tagmentation products",MaterialKind.WASTE)
    g.react("tagment","Tagment accessible chromatin before partitioning",StepKind.TAGMENTATION,((nucleus,InputRole.TEMPLATE),(transposome,InputRole.STRUCTURAL),(tn5,InputRole.CATALYST)),(tagged,waste))
    g.discard(waste,"does not retain an amplifiable accessible-DNA fragment")
    droplets=g.input("droplets","Multiome GEM partitioning condition",MaterialKind.CONDITION)
    occupied=g.product("occupied","occupied GEM products",MaterialKind.MOLECULAR,state=tagged.state)
    failed=g.product("failed","empty or unresolved GEM products",MaterialKind.WASTE)
    g.react("partition","Partition tagged nuclei with barcoded gel beads",StepKind.PARTITION,((tagged,InputRole.TEMPLATE),(droplets,InputRole.CONDITION)),(occupied,failed),changes_molecule=False)
    g.discard(failed,"does not yield a resolvable nucleus barcode")
    rna=g.product("rna_branch","GEM RNA branch",MaterialKind.MOLECULAR,state=occupied.state)
    atac=g.product("atac_branch","GEM ATAC branch",MaterialKind.MOLECULAR,state=occupied.state)
    g.react("split","Separate gene-expression and ATAC products",StepKind.SPLIT,((occupied,InputRole.TEMPLATE),),(rna,atac),changes_molecule=False)
    rtprimer=g.input("rt_primer","cell-barcode/UMI oligo(dT) gel-bead primer",MaterialKind.OLIGO,source=SOURCE+" (Rev G Appendix p. 80 prints its constant structure and variable windows)")
    rt=g.input("rt","reverse transcriptase",MaterialKind.CATALYST)
    cdna=g.product("cdna","cell-barcoded, UMI-tagged cDNA",MaterialKind.MOLECULAR,state=molecule("Multiome cDNA",(role("cell barcode"),role("UMI"),role("cDNA"))))
    g.react("reverse_transcribe","Transfer the GEM barcode to poly(A)-derived cDNA",StepKind.REVERSE_TRANSCRIPTION,((rna,InputRole.TEMPLATE),(rtprimer,InputRole.STRUCTURAL),(rt,InputRole.CATALYST)),(cdna,))
    bridge=g.input("atac_bridge","cell-barcode ATAC bridge oligo",MaterialKind.OLIGO,inferred=True,source=SOURCE+" (complete bridge sequence not disclosed)")
    pol=g.input("atac_extension","ATAC barcode-extension polymerase",MaterialKind.CATALYST)
    barcoded_atac=g.product("barcoded_atac","cell-barcoded accessible-DNA fragments",MaterialKind.MOLECULAR,state=molecule("Multiome ATAC precursor",(role("cell barcode"),role("accessible DNA"))))
    g.react("barcode_atac","Transfer the GEM barcode to tagged genomic DNA",StepKind.EXTENSION,((atac,InputRole.TEMPLATE),(bridge,InputRole.STRUCTURAL),(pol,InputRole.CATALYST)),(barcoded_atac,))
    reviewed=["nucleus","transposome","tn5","tagment","droplets","partition","split","rt_primer","rt","reverse_transcribe","atac_bridge","atac_extension","barcode_atac"]
    gex_geometry=PCRGeometry.exact(
        GEX,
        PCRPrimer("SI-TT i5",il.P5+"A"*10+il.TRUSEQ_READ1,anneal_nt=len(il.TRUSEQ_READ1)),
        PCRPrimer("SI-TT i7",il.P7+"A"*10+il.TRUSEQ_READ2,anneal_nt=len(il.TRUSEQ_READ2)),
        statement="The final-library layout in CG000338 fixes the 3′ TruSeq Read 1 and Read 2 landing regions; representative index windows do not affect annealing.")
    atac_geometry=PCRGeometry.exact(
        ATAC,PCRPrimer("ATAC forward",il.P5[:22]),
        PCRPrimer("N7xx",nx.n7xx_primer("A"*8),anneal_nt=len(nx.S7)),
        statement="CG000338 prints the P5-side ATAC primer and the P7-index/s7 endpoint; their 3′ regions face inward on the completed ATAC molecule.")
    finish_pcr_branch(g,"gex",cdna,GEX,GEX_PRIMERS,reviewed,label="gene-expression",geometry=gex_geometry)
    finish_pcr_branch(g,"atac",barcoded_atac,ATAC,ATAC_PRIMERS,reviewed,label="ATAC",geometry=atac_geometry)
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
