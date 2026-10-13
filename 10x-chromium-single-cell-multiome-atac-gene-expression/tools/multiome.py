"""10x Chromium Next GEM Single Cell Multiome ATAC + Gene Expression."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
from batch_ngs import seg, truseq_library
from chemdraw import Construct, Row, feature
import illumina as il
import nextera as nx
import seqprimers as sp
from multimodal_spatial import finish_pcr_branch, modality_split
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, StepKind
from protocol_model import molecule, role

TITLE = "10x Chromium Single Cell Multiome ATAC + Gene Expression"
NOTES = "01_multiome.html"
SOURCE = 'Commercial protocol: 10x Genomics <a href="https://www.10xgenomics.com/support/instruments/chromium-x-series/chromium-next-gem-single-cell-multiome-atac-plus-gene-expression-reagent-kits-user-guide">CG000338 Rev G</a>.'
SUMMARY = "Nuclei are tagmented before GEM formation. Inside each GEM, one gel-bead barcode is transferred to both accessible-DNA fragments and poly(A)-derived cDNA, yielding linked ATAC and gene-expression libraries from the same nucleus."
CAVEAT = "10x does not disclose the complete gel-bead oligos. Published barcode, UMI and adapter roles are shown as length-preserving placeholders; no proprietary bases are inferred."
CELL_BARCODE = feature("cell_barcode", "cell_barcode", "whitelist",
                       whitelist="10x-multiome-gem-barcode")
UMI_FEATURE = feature("umi", "umi", "random")
I7_FEATURE = feature("sample_index_i7", "sample_index", "unknown")

GEX, GEX_PRIMERS = truseq_library([
    seg("cell barcode", "B" * 16, "cbc", placeholder=True, inferred=True,
        feature=CELL_BARCODE),
    seg("UMI", "U" * 12, "umi", placeholder=True, inferred=True,
        feature=UMI_FEATURE),
    seg("poly(dT) junction", "T" * 12, placeholder=True, inferred=True),
    seg("cDNA", "X" * 36, placeholder=True)], "Multiome gene-expression library",
    inferred_adapters=True)
ATAC = Construct([
    seg("P5", il.P5, "p5"),
    seg("10x cell barcode", "B" * 16, "cbc", placeholder=True, inferred=True,
        feature=CELL_BARCODE),
    seg("S5", nx.S5, "s5"), seg("left mosaic end", nx.ME, "me"),
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
     "INFERRED — Read 1 reports the GEM cell barcode and UMI; undisclosed adapter and bead-oligo bases remain visibly inferred.",
     "The declared Illumina primers are located on the final gene-expression duplex."),
    ("Final ATAC library", ATAC, ATAC_PRIMERS,
     "INFERRED — Paired genomic reads flank the accessible fragment; Index 2 reads the undisclosed 16-base GEM barcode and Index 1 reads the sample index.",
     "The declared Nextera sequencing primers are located on the final ATAC duplex."),
)

def sections():
    return [
        ("Tagment accessible chromatin before partitioning",
         [Row(chunks=[("nucleus chromatin → Tn5 cleavage + adapter transfer → tagged open-DNA fragments", "me", False)])],
         "Transposition occurs in bulk nuclei; intact nuclei retain both tagged chromatin and RNA."),
        ("Transfer one GEM barcode to RNA and ATAC molecules",
         [Row(chunks=[("poly(A) RNA + barcoded oligo-dT → cell barcode + UMI + cDNA", "cbc", False)]),
          Row(chunks=[("tagged genomic DNA + barcoded bridge → cell barcode + ATAC fragment", "cbc", False)])],
         "INFERRED — Both reactions occur inside the same GEM and therefore use the same published barcode role; undisclosed bases are placeholders."),
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
    rtprimer=g.input("rt_primer","cell-barcode/UMI oligo(dT) gel-bead primer",MaterialKind.OLIGO,inferred=True,source=SOURCE+" (complete gel-bead sequence not disclosed)")
    rt=g.input("rt","reverse transcriptase",MaterialKind.CATALYST)
    cdna=g.product("cdna","cell-barcoded, UMI-tagged cDNA",MaterialKind.MOLECULAR,state=molecule("Multiome cDNA",(role("cell barcode"),role("UMI"),role("cDNA"))))
    g.react("reverse_transcribe","Transfer the GEM barcode to poly(A)-derived cDNA",StepKind.REVERSE_TRANSCRIPTION,((rna,InputRole.TEMPLATE),(rtprimer,InputRole.STRUCTURAL),(rt,InputRole.CATALYST)),(cdna,))
    bridge=g.input("atac_bridge","cell-barcode ATAC bridge oligo",MaterialKind.OLIGO,inferred=True,source=SOURCE+" (complete bridge sequence not disclosed)")
    pol=g.input("atac_extension","ATAC barcode-extension polymerase",MaterialKind.CATALYST)
    barcoded_atac=g.product("barcoded_atac","cell-barcoded accessible-DNA fragments",MaterialKind.MOLECULAR,state=molecule("Multiome ATAC precursor",(role("cell barcode"),role("accessible DNA"))))
    g.react("barcode_atac","Transfer the GEM barcode to tagged genomic DNA",StepKind.EXTENSION,((atac,InputRole.TEMPLATE),(bridge,InputRole.STRUCTURAL),(pol,InputRole.CATALYST)),(barcoded_atac,))
    reviewed=["nucleus","transposome","tn5","tagment","droplets","partition","split","rt_primer","rt","reverse_transcribe","atac_bridge","atac_extension","barcode_atac"]
    finish_pcr_branch(g,"gex",cdna,GEX,GEX_PRIMERS,reviewed,label="gene-expression",inferred_primers=True)
    finish_pcr_branch(g,"atac",barcoded_atac,ATAC,ATAC_PRIMERS,reviewed,label="ATAC",inferred_primers=True)
    g.declare_reviewed_scope(reviewed,SOURCE)
    return g.require_publishable()
