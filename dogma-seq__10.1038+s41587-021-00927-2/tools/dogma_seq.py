"""DOGMA-seq molecular model."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
from antibody_tags import feature_library
from batch_ngs import nextera_library, seg, truseq_library
from chemdraw import Row, feature
from multimodal_spatial import finish_pcr_branch, modality_split
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, StepKind
from protocol_model import molecule, role

TITLE = "DOGMA-seq — RNA, ATAC and protein from one cell"
NOTES = "01_dogma-seq.html"
SOURCE = 'Defining source: <a href="https://doi.org/10.1038/s41587-021-00927-2">Mimitou et al., <i>Nature Biotechnology</i> (2021)</a>.'
SUMMARY = "A 10x Multiome-compatible workflow jointly captures poly(A) RNA, accessible chromatin and antibody-derived tags from each permeabilized cell."
CAVEAT = ('The three outputs share cell identity but are distinct sequencing libraries. '
          'The cited <a href="https://assets.ctfassets.net/an68im79xiti/5tnYerYC4FdDwYK0hlcfWj/cc46120ef02327249b5c42f24f44d699/CG000338_ChromiumNextGEM_Multiome_ATAC_GEX_User_Guide_RevA.pdf">10x CG000338 Rev A guide</a> '
          'prints both Multiome gel-bead primer architectures; diagrams use role tokens for the variable barcode and UMI members and antibody-tag conjugates.')
CELL = feature("dogma_cell", "cell_barcode", "whitelist", whitelist="10x Chromium Next GEM bead barcode")
UMI = feature("dogma_rna_umi", "umi", "random")
RNA, RNA_P = truseq_library([seg("cell barcode", "C"*16,"cbc",placeholder=True,feature=CELL),seg("RNA UMI","U"*12,"umi",placeholder=True,feature=UMI),seg("3-prime cDNA","X"*40,placeholder=True)], "DOGMA-seq RNA library")
ATAC, ATAC_P = nextera_library([seg("cell barcode association", "C"*16,"cbc",placeholder=True,feature=CELL),seg("accessible genomic DNA","X"*40,placeholder=True)], "DOGMA-seq ATAC library")
ADT, ADT_P = feature_library("DOGMA-seq antibody-tag library", barcode_id="dogma_antibody",cell_id="dogma_cell",
                            barcode_len=15, umi_len=10, capture_len=18,
                            cell_barcode_len=16,
                            bead_whitelist="single-cell bead barcode set",
                            panel_whitelist="published antibody-tag panel")
FINAL_LIBRARIES = (
 ("RNA library",RNA,RNA_P,"Poly(dT) capture, GEM barcoding and cDNA library construction recover gene expression.","Read 1 carries cell barcode/UMI; the insert read reports cDNA."),
 ("ATAC library",ATAC,ATAC_P,"Tagmented chromatin is barcoded and PCR-completed as the accessibility library.","Paired genomic reads use Nextera sites."),
 ("Antibody-tag library",ADT,ADT_P,"Antibody-derived oligos are selectively amplified from the same emulsion product.","Feature barcode and UMI are read separately from genomic and cDNA inserts."),
)
READ_LENGTHS={"RNA library":{"Read 1":28,"Read 2":90},"ATAC library":{"Read 1":50,"Read 2":50},"Antibody-tag library":{"Read 1":28,"Read 2":50}}

def oligos():
    return [
        "<b>10x Multiome RNA gel-bead primer (CG000338 Rev A, PN-2000261)</b> — 5′-CTACACGACGCTCTTCCGATCT-N16-N12-TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTVN-3′ (N16: cell barcode; N12: UMI)",
        "<b>10x Multiome ATAC gel-bead primer (same guide)</b> — 5′-AATGATACGGCGACCACCGAGATCTACAC-N16-CGCGTCTG-3′",
        "<b>ADT additive for Multiome preamplification</b> — 5′-CCTTGGCACCCGAGAATT*C*C-3′ (*: phosphorothioate linkage)",
        "<b>Antibody-tag SI-PCR forward</b> — 5′-AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTC-3′",
        "<b>Antibody-tag RPI-x P7 index family</b> — 5′-CAAGCAGAAGACGGCATACGAGATxxxxxxxxGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA-3′ (x: sample index)",
    ]

def sections(): return [
 ("Stain cells with DNA-barcoded antibodies",[Row(chunks=[("antibody — feature barcode — UMI — capture handle","cbc",False)])],"Antibody tags encode the detected epitope."),
 ("Tagment accessible chromatin",[Row(chunks=[("Tn5 inserts adapters at accessible DNA","me",False)])],"Tagmentation occurs while RNA and antibody tags remain associated with the cell."),
 ("Partition cells with barcoded gel beads",[Row(chunks=[("one GEM barcode → RNA + ATAC + antibody tag","cbc",False)])],"The shared bead identity links the three measurements computationally."),
 ("Recover three modality-specific libraries",modality_split("RNA library","ATAC library","antibody-tag library"),"Separate enrichments preserve a common cell barcode while using modality-specific primers."),
]

def protocol_graph():
 g=ProtocolGraph(TITLE,authority=SOURCE)
 cells=g.input("cells","permeabilized cells containing RNA and accessible chromatin",MaterialKind.BIOLOGICAL,state=molecule("DOGMA-seq cells",(role("poly(A) RNA"),role("accessible chromatin"),role("protein epitopes"))),primary=True)
 antibodies=g.input("antibodies","DNA-barcoded antibodies",MaterialKind.SELECTION)
 stained=g.product("stained","cells retaining antibody-derived tags",MaterialKind.MOLECULAR,state=molecule("stained cells",(role("poly(A) RNA"),role("accessible chromatin"),role("antibody tags"))))
 unbound=g.product("unbound","unbound antibody tags",MaterialKind.WASTE)
 g.react("stain","Stain proteins with DNA-barcoded antibodies",StepKind.BINDING,((cells,InputRole.TEMPLATE),(antibodies,InputRole.SELECTION)),(stained,unbound)); g.discard(unbound,"does not report a cell-associated protein")
 transposome=g.input("transposome","Multiome-compatible loaded transposome",MaterialKind.ADAPTER)
 tn5=g.input("tn5","Tn5 transposase",MaterialKind.CATALYST)
 tagged=g.product("tagged","cells retaining tagmented chromatin, RNA and antibody tags",MaterialKind.MOLECULAR,state=stained.state)
 twaste=g.product("tn5_waste","non-library tagmentation products",MaterialKind.WASTE)
 g.react("tagment","Tagment accessible chromatin",StepKind.TAGMENTATION,((stained,InputRole.TEMPLATE),(transposome,InputRole.STRUCTURAL),(tn5,InputRole.CATALYST)),(tagged,twaste)); g.discard(twaste,"does not retain an amplifiable accessibility fragment")
 droplets=g.input("droplets","10x Multiome GEM partitioning condition",MaterialKind.CONDITION)
 occupied=g.product("occupied","cell-barcoded multimodal products",MaterialKind.MOLECULAR,state=tagged.state)
 failed=g.product("failed","empty or unresolved GEM products",MaterialKind.WASTE)
 g.react("partition","Partition cells with barcoded gel beads",StepKind.PARTITION,((tagged,InputRole.TEMPLATE),(droplets,InputRole.CONDITION)),(occupied,failed),changes_molecule=False); g.discard(failed,"does not retain a resolvable cell barcode")
 rna=g.product("rna_branch","RNA branch",MaterialKind.MOLECULAR,state=occupied.state); atac=g.product("atac_branch","ATAC branch",MaterialKind.MOLECULAR,state=occupied.state); adt=g.product("adt_branch","antibody-tag branch",MaterialKind.MOLECULAR,state=occupied.state)
 g.react("split","Separate RNA, ATAC and antibody-tag products",StepKind.SPLIT,((occupied,InputRole.TEMPLATE),),(rna,atac,adt),changes_molecule=False)
 bead=g.input("bead_primer","cell-barcode/UMI oligo(dT) bead primer",MaterialKind.OLIGO,inferred=True,source=SOURCE+" (complete historical commercial sequence not yet verified in source audit)")
 rt=g.input("rt","reverse transcriptase",MaterialKind.CATALYST)
 cdna=g.product("cdna","cell-barcoded, UMI-tagged cDNA",MaterialKind.MOLECULAR,state=molecule("DOGMA cDNA",(role("cell barcode"),role("RNA UMI"),role("cDNA"))))
 g.react("reverse_transcribe","Capture and reverse-transcribe poly(A) RNA",StepKind.REVERSE_TRANSCRIPTION,((rna,InputRole.TEMPLATE),(bead,InputRole.STRUCTURAL),(rt,InputRole.CATALYST)),(cdna,))
 tag_primer=g.input("tag_primer","antibody-tag capture primer",MaterialKind.OLIGO,inferred=True,source=SOURCE+" (complete historical commercial sequence not yet verified in source audit)")
 tag_pol=g.input("tag_polymerase","antibody-tag barcode-transfer polymerase",MaterialKind.CATALYST)
 copied=g.product("copied_tags","cell-barcoded antibody-derived tags",MaterialKind.MOLECULAR,state=molecule("DOGMA antibody-tag precursor",(role("cell barcode"),role("UMI"),role("antibody barcode"))))
 g.react("copy_tags","Transfer cell identity to antibody-derived tags",StepKind.EXTENSION,((adt,InputRole.TEMPLATE),(tag_primer,InputRole.STRUCTURAL),(tag_pol,InputRole.CATALYST)),(copied,))
 reviewed=["cells","antibodies","stain","transposome","tn5","tagment","droplets","partition","split","bead_primer","rt","reverse_transcribe","tag_primer","tag_polymerase","copy_tags"]
 finish_pcr_branch(g,"rna",cdna,RNA,RNA_P,reviewed,label="DOGMA-seq RNA",read_lengths=READ_LENGTHS["RNA library"],inferred_primers=True)
 finish_pcr_branch(g,"atac",atac,ATAC,ATAC_P,reviewed,label="DOGMA-seq ATAC",read_lengths=READ_LENGTHS["ATAC library"],inferred_primers=True)
 finish_pcr_branch(g,"adt",copied,ADT,ADT_P,reviewed,label="DOGMA-seq antibody tag",read_lengths=READ_LENGTHS["Antibody-tag library"],inferred_primers=True)
 g.declare_reviewed_scope(reviewed,SOURCE)
 return g.require_publishable()
