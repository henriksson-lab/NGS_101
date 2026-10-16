"""SEC-seq molecular model."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"lib"))
from chemdraw import Construct, Row, Scene, Segment, feature, revcomp
import illumina as il
import nextera as nx
import seqprimers as sp
from pcr import PCRGeometry, PCRPrimer
from adapter_ligation import AdapterLigationGeometry
from endprep import repair_and_dA_tail, dA_tailed_scene
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind
from protocol_model import molecule, role

TITLE="SEC-seq — secretion-linked single-cell RNA sequencing"
NOTES="01_sec-seq.html"
SOURCE='<a href="https://doi.org/10.1038/s41467-023-39367-8">Cheng et al., <i>Nature Communications</i> (2023)</a>.'
SUMMARY="A cell is held in an antibody-coated nanovial that captures its secreted protein. An oligo-barcoded detection antibody and the cell transcriptome are then read together with 10x Chromium."
GUIDE = "10x CG000317, Appendix pp71–72 (Rev E): https://www.10xgenomics.com/support/user-guides/chromium-next-gem-single-cell-3-reagent-kits-v31-dual-index-with-feature-barcode-technology-for-cell-surface-protein-user-guide"
CAVEAT = ("This schematic represents the main human ASC experiment: Chromium 3′ v3.1 with "
          "CG000317 and 28/10/10/90 sequencing. The paper’s mixed-species validation uses a "
          "distinct 5′ v2/TotalSeq-C branch. The <a href='https://cdn.10xgenomics.com/image/upload/v1660261285/support-documents/CG000317_ChromiumNextGEMSingleCell3_v3.1_CellSurfaceProtein_RevE.pdf'>public CG000317 Rev E guide</a> supplies the dual-index "
          "GEX and Capture Sequence 1 feature architectures shown here; the exact guide revision, "
          "anti-PE barcode member and selected sample-index pair remain unassigned in this model.")
CAVEAT += (" <strong>Apparent vendor-guide inconsistency:</strong> CG000317 Rev E p.72 prints "
           "the TT forward primer and GEX product with the 29-base P5-side prefix "
           "AATGATACGGCGACCACCGAGATCTACAC before the i5 index; "
           "<a href='https://assets.ctfassets.net/an68im79xiti/5tnYerYC4FdDwYK0hlcfWj/cc46120ef02327249b5c42f24f44d699/CG000338_ChromiumNextGEM_Multiome_ATAC_GEX_User_Guide_RevA.pdf'>CG000338 Rev A p.80</a> "
           "prints a 25-base GEX prefix ending AGATCT. Both identify Dual Index Kit TT Set A "
           "PN-1000215 / plate PN-3000431. This page retains the CG000317 Rev E 29-base "
           "layout shown by its primer and product; the historical physical prefix used in "
           "SEC-seq has not been independently established, so a reagent variant is not asserted.")
CAPTURE_SEQ1 = "TTGCTAGGACCGGCCTTAAAGC"
GEX_BEAD_HANDLE = "CTACACGACGCTCTTCCGATCT"
FEATURE_AMPLIFICATION_HANDLE = "GCAGC" + "GTCAGATGTGTATAAGAGACAG"
TSO = "AAGCAGTGGTATCAACGCAGAGTACAT"
ROLES = ("Read 1", "Index 1 (i7)", "Index 2 (i5)", "Read 2")

def _identity_segments(branch):
    return (
        Segment("cell barcode", "B" * 16, "cbc", placeholder=True,
                feature=feature("sec_cell", "cell_barcode", "whitelist",
                                whitelist="10x Chromium 3′ v3.1")),
        Segment("UMI", "N" * 12, "umi", placeholder=True,
                feature=feature("sec_" + branch + "_umi", "umi", "random")),
    )

def _index(name):
    return Segment(name, "N" * 10, "cbc", placeholder=True,
                   feature=feature("sec_" + name, "sample_index", "unknown"))

def gex_precursor():
    return Construct((Segment("GEX bead handle", GEX_BEAD_HANDLE, "r1"),
                      *_identity_segments("gex"),
                      Segment("poly(dT)", "T" * 30),
                      Segment("VN anchor", "VN", placeholder=True, length_bp=2),
                      Segment("cDNA insert", "X" * 40, placeholder=True),
                      Segment("terminal dA", "A"),
                      Segment("ligated TruSeq Read 2 arm", il.INDEX1_PRIMER, "r2")),
                     name="CG000317 GEX adaptor-ligated PCR precursor")

def secretion_precursor():
    return Construct((Segment("feature amplification handle", FEATURE_AMPLIFICATION_HANDLE, "r1"),
                      *_identity_segments("secretion"),
                      Segment("Capture Sequence 1", CAPTURE_SEQ1),
                      Segment("feature upstream N9", "N" * 9, placeholder=True, length_bp=9),
                      Segment("anti-PE feature barcode", "N" * 15, "cbc", placeholder=True,
                              feature=feature("sec_antibody", "feature_barcode", "unknown")),
                      Segment("feature downstream N10", "N" * 10, placeholder=True, length_bp=10),
                      Segment("TruSeq Read 2 arm", revcomp(il.TRUSEQ_READ2), "r2")),
                     name="CG000317 secretion feature PCR precursor")

def gene_expression_library():
    precursor = gex_precursor()
    return Construct((Segment("P5", il.P5, "p5"), _index("i5"),
                      Segment("TruSeq Read 1", il.TRUSEQ_READ1, "r1"),
                      *precursor.segments[1:], _index("i7"),
                      Segment("P7 reverse complement", il.P7_RC, "p7")),
                     name="SEC-seq human ASC gene-expression library")

def secretion_feature_library():
    precursor = secretion_precursor()
    return Construct((Segment("P5", il.P5, "p5"), _index("i5"),
                      Segment("Nextera Read 1", nx.READ1_PRIMER, "r1"),
                      *precursor.segments[1:], _index("i7"),
                      Segment("P7 reverse complement", il.P7_RC, "p7")),
                     name="SEC-seq human ASC Capture Sequence 1 secretion library")

RNA = gene_expression_library()
SEC = secretion_feature_library()
RNA_SEQ_PRIMERS = sp.standard_primers(RNA, ROLES)
# CG000317 explicitly identifies this mixed Nextera-R1 / TruSeq-R2 shell.
# Nextera R2 also partially matches the R1 mosaic end, so use the vendor-declared
# canonical identities instead of treating that partial match as an acquired R2.
SEC_SEQ_PRIMERS = (sp.NEXTERA["R1"], sp.TRUSEQ["I1"],
                   sp.NEXTERA["I2"], sp.TRUSEQ["R2"])
if problems := sp.verify(SEC, SEC_SEQ_PRIMERS, required_roles=ROLES):
    raise ValueError("; ".join(problems))
FINAL_LIBRARIES = (
    ("Gene-expression library", RNA, RNA_SEQ_PRIMERS,
     "Main human ASC branch: CG000317 dual-index GEX, TruSeq Read 1.",
     "Read 1 records the cell barcode and UMI; Read 2 enters the cDNA insert."),
    ("Secretion-tag library", SEC, SEC_SEQ_PRIMERS,
     "Main human ASC branch: CG000317 dual-index Capture Sequence 1 feature library, Nextera Read 1.",
     "Read 1 records the shared cell barcode and a molecule-specific UMI; Read 2 traverses the feature’s N10 flank before the antibody barcode."),
)
READ_LENGTHS = {heading: dict(zip(ROLES, (28, 10, 10, 90)))
                for heading, *_ in FINAL_LIBRARIES}
READ_USES = {heading: {"Read 1": "extract", "Index 1 (i7)": "demultiplex",
                      "Index 2 (i5)": "demultiplex", "Read 2": "include"}
             for heading, *_ in FINAL_LIBRARIES}
READ_USES["Secretion-tag library"]["Read 2"] = "extract"

def _amplify_capture(graph, prefix, template, final, forward_seq, reverse_seq, reviewed):
    forward = graph.input(prefix + "_capture_forward", "CG000317 capture-product amplification forward primer",
                          MaterialKind.OLIGO, source=GUIDE)
    reverse = graph.input(prefix + "_capture_reverse", "CG000317 capture-product amplification reverse primer",
                          MaterialKind.OLIGO, source=GUIDE)
    polymerase = graph.input(prefix + "_capture_polymerase", "cDNA amplification polymerase", MaterialKind.CATALYST)
    amplified = graph.product(prefix + "_amplified", "amplified " + prefix + " capture product",
                              MaterialKind.MOLECULAR, state=final, source=GUIDE)
    geometry = PCRGeometry.exact(template.state,
        PCRPrimer("capture-product forward primer", forward_seq, anneal_nt=22),
        PCRPrimer("capture-product reverse primer", reverse_seq, anneal_nt=21),
        statement="CG000317 Appendix p71 publishes the separate GEX and feature amplification primers.")
    graph.react(prefix + "_capture_pcr", "Amplify the " + prefix + " capture product", StepKind.PCR,
                ((template, InputRole.TEMPLATE), (forward, InputRole.STRUCTURAL),
                 (reverse, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)),
                (amplified,), geometry=geometry)
    reviewed.extend((forward.id, reverse.id, polymerase.id, prefix + "_capture_pcr"))
    return amplified

def _indexed_branch(graph, prefix, template, precursor, final, primers, reviewed, heading):
    forward_tail = il.TRUSEQ_READ1[:24] if prefix == "rna" else nx.READ1_PRIMER
    forward_sequence = il.P5 + "N" * 10 + forward_tail
    reverse_sequence = il.P7 + "N" * 10 + il.TRUSEQ_READ2[:21]
    forward = graph.input(prefix + "_forward", "CG000317 indexed forward primer family",
                          MaterialKind.OLIGO, source=GUIDE)
    reverse = graph.input(prefix + "_reverse", "CG000317 indexed reverse primer family",
                          MaterialKind.OLIGO, source=GUIDE)
    polymerase = graph.input(prefix + "_polymerase", "library PCR polymerase", MaterialKind.CATALYST)
    library = graph.product(prefix + "_library", heading, MaterialKind.LIBRARY, state=final, source=GUIDE)
    geometry = PCRGeometry.exact(
        precursor, PCRPrimer("indexed forward primer family", forward_sequence,
                             anneal_nt=13 if prefix == "rna" else 26),
        PCRPrimer("indexed reverse primer family", reverse_sequence, anneal_nt=21),
        statement="CG000317 Appendix p72 publishes both primer families and their sites; ten-base index members remain symbolic.")
    graph.react(prefix + "_pcr", "Amplify and dual-index the " + heading, StepKind.PCR,
                ((template, InputRole.TEMPLATE), (forward, InputRole.STRUCTURAL),
                 (reverse, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)),
                (library,), geometry=geometry)
    graph.sequence(library, SequencingPlan(final, primers, ROLES,
                   read_lengths=READ_LENGTHS[heading], read_uses=READ_USES[heading]))
    reviewed.extend((forward.id, reverse.id, polymerase.id, prefix + "_pcr", library.id))

def oligos():
    return (
        "CG000317 GEX bead-primer family — 5′-" + GEX_BEAD_HANDLE + "-N16-N12-" + "T" * 30 + "VN-3′",
        "CG000317 feature bead-primer family — 5′-" + FEATURE_AMPLIFICATION_HANDLE[5:] + "-N16-N12-" + CAPTURE_SEQ1 + "-3′",
        "CG000317 template-switch oligo — 5′-" + TSO + "rGrGrG-3′",
        "The N15 secretion-antibody barcode and selected ten-base i5/i7 members remain symbolic; these are guide families, not complete ordering sequences for the SEC-seq anti-PE conjugate.",
    )

def sections(): return [
 ("Capture one cell in a nanovial",[Row(chunks=[("antibody-coated cavity ← cell",None,False)])],"The cavity retains the cell and concentrates secreted protein on its surface."),
 ("Label captured secretion",[Row(chunks=[("nanovial — secreted IgG — detection antibody — DNA feature barcode","cbc",False)])],"An oligo-barcoded secondary antibody converts secretion into a sequenceable feature."),
 ("Partition cell-loaded nanovials",[Row(chunks=[("one GEM barcode → cellular mRNA + nanovial antibody tags","cbc",False)])],"The nanovial fits through the Chromium microfluidic channel and does not replace the standard bead chemistry."),
 ("Amplify two libraries",[Row(chunks=[("barcoded cDNA ├─ gene expression\n              └─ secretion feature tags",None,False)])],"The paper pools 80% RNA and 20% antibody-tag library for sequencing."),]

def protocol_graph():
 g=ProtocolGraph(TITLE,authority=SOURCE)
 cells=g.input("cells","cells containing poly(A) RNA and secreted protein analytes",MaterialKind.BIOLOGICAL,state=molecule("SEC-seq cells",(role("poly(A) RNA"),role("secreted proteins"))),primary=True)
 nanovials=g.input("nanovials","capture-antibody-coated nanovials",MaterialKind.SELECTION)
 loaded=g.product("loaded","cell-loaded nanovials retaining secreted protein",MaterialKind.MOLECULAR,state=molecule("loaded nanovial",(role("cell"),role("captured secretion"),role("poly(A) RNA"))))
 uncaptured=g.product("uncaptured","uncaptured secretion and empty nanovials",MaterialKind.WASTE)
 g.react("capture","Capture cells and their secretion in antibody-coated nanovials",StepKind.BINDING,((cells,InputRole.TEMPLATE),(nanovials,InputRole.SELECTION)),(loaded,uncaptured)); g.discard(uncaptured,"does not associate secretion with a retained cell-loaded nanovial")
 detection=g.input("detection_antibody","oligo-barcoded secretion-detection antibody",MaterialKind.SELECTION,inferred=True,source=SOURCE+" (selected anti-PE conjugate sequence not declared in this model)")
 labelled=g.product("labelled","nanovials bearing sequenceable secretion tags",MaterialKind.MOLECULAR,state=molecule("labelled nanovial",(role("cell"),role("poly(A) RNA"),role("secretion antibody tag"))))
 unbound=g.product("unbound","unbound detection-antibody tags",MaterialKind.WASTE)
 g.react("label","Label captured secretion with oligo-barcoded antibodies",StepKind.BINDING,((loaded,InputRole.TEMPLATE),(detection,InputRole.SELECTION)),(labelled,unbound)); g.discard(unbound,"does not report nanovial-captured secretion")
 droplets=g.input("droplets","10x Chromium GEM partitioning condition",MaterialKind.CONDITION)
 occupied=g.product("occupied","cell-barcoded transcript and secretion-tag products",MaterialKind.MOLECULAR,state=labelled.state)
 failed=g.product("failed","empty or unresolved GEM products",MaterialKind.WASTE)
 g.react("partition","Partition cell-loaded nanovials with barcoded beads",StepKind.PARTITION,((labelled,InputRole.TEMPLATE),(droplets,InputRole.CONDITION)),(occupied,failed),changes_molecule=False); g.discard(failed,"does not retain a resolvable nanovial-associated cell barcode")
 rna=g.product("rna_branch","gene-expression branch",MaterialKind.MOLECULAR,state=occupied.state); sec=g.product("sec_branch","secretion-tag branch",MaterialKind.MOLECULAR,state=occupied.state)
 g.react("split","Separate transcript and secretion-tag products",StepKind.SPLIT,((occupied,InputRole.TEMPLATE),),(rna,sec),changes_molecule=False)
 bead=g.input("bead_primer","CG000317 TruSeq-handle cell-barcode/UMI oligo(dT) bead primer",MaterialKind.OLIGO,source=GUIDE)
 rt=g.input("rt","reverse transcriptase",MaterialKind.CATALYST)
 cdna=g.product("cdna","cell-barcoded transcript cDNA",MaterialKind.MOLECULAR,state=Construct((*gex_precursor().segments[:-2], Segment("template-switch handle", revcomp(TSO + "GGG"))), name="CG000317 GEX reverse-transcription product"))
 g.react("reverse_transcribe","Capture and reverse-transcribe cellular mRNA",StepKind.REVERSE_TRANSCRIPTION,((rna,InputRole.TEMPLATE),(bead,InputRole.STRUCTURAL),(rt,InputRole.CATALYST)),(cdna,))
 tagprimer=g.input("tag_primer","CG000317 Nextera-handle Capture Sequence 1 bead primer",MaterialKind.OLIGO,source=GUIDE)
 tagpol=g.input("tag_extension","secretion-tag barcode-transfer polymerase",MaterialKind.CATALYST)
 copied=g.product("copied","cell-barcoded secretion-antibody tags",MaterialKind.MOLECULAR,state=Construct((Segment("GEM feature handle", FEATURE_AMPLIFICATION_HANDLE[5:], "r1"), *secretion_precursor().segments[1:]), name="CG000317 secretion GEM product"))
 g.react("copy_tags","Transfer cell identity to secretion-antibody tags",StepKind.EXTENSION,((sec,InputRole.TEMPLATE),(tagprimer,InputRole.STRUCTURAL),(tagpol,InputRole.CATALYST)),(copied,))
 reviewed=["cells","nanovials","capture","detection_antibody","label","droplets","partition","split","bead_primer","rt","reverse_transcribe","tag_primer","tag_extension","copy_tags"]
 amplified_rna = _amplify_capture(g,"gex",cdna,cdna.state,GEX_BEAD_HANDLE,TSO[:22],reviewed)
 amplified_sec = _amplify_capture(g,"feature",copied,secretion_precursor(),FEATURE_AMPLIFICATION_HANDLE,il.TRUSEQ_READ2[:21],reviewed)
 fragmentation = g.input("fragmentation", "CG000317 cDNA fragmentation enzymes", MaterialKind.CATALYST)
 fragment_state = Construct(gex_precursor().segments[:-2], name="selected GEX 3′ cDNA fragment")
 fragmented = g.product("fragmented_gex", "fragmented GEX cDNA", MaterialKind.MOLECULAR, state=fragment_state)
 other_fragments = g.product("gex_other_fragments", "cDNA fragments without the 3′ cell-barcode handle", MaterialKind.WASTE)
 g.react("fragment_gex", "Fragment amplified GEX cDNA", StepKind.FRAGMENTATION,
         ((amplified_rna,InputRole.TEMPLATE),(fragmentation,InputRole.CATALYST)),(fragmented,other_fragments))
 g.discard(other_fragments,"lacks the barcode-side primer site required for the final GEX amplicon")
 endprep = g.input("gex_endprep", "CG000317 end-repair and A-tailing enzymes", MaterialKind.CATALYST)
 prepared_state = repair_and_dA_tail(fragment_state)
 prepared = g.product("prepared_gex", "end-repaired and dA-tailed GEX fragments", MaterialKind.MOLECULAR,
                      state=dA_tailed_scene(prepared_state))
 g.react("prepare_gex_ends", "End-repair and dA-tail GEX fragments", StepKind.END_REPAIR,
         ((fragmented,InputRole.TEMPLATE),(endprep,InputRole.CATALYST)),(prepared,))
 adaptor = g.input("gex_adaptor", "CG000317 PN-2000094 adaptor oligos", MaterialKind.ADAPTER, source=GUIDE)
 ligase = g.input("gex_ligase", "adaptor ligation enzyme", MaterialKind.CATALYST)
 sealed = Scene()
 sealed.strand("top",gex_precursor().segments,label="GEX ligated strand")
 sealed.junction("top","terminal dA","ligated TruSeq Read 2 arm")
 ligated = g.product("gex_ligated", "GEX adaptor-ligated PCR precursor", MaterialKind.MOLECULAR, state=sealed)
 g.react("ligate_gex_adaptor", "Ligate the GEX Read 2 adaptor", StepKind.LIGATION,
         ((prepared,InputRole.TEMPLATE),(adaptor,InputRole.STRUCTURAL),(ligase,InputRole.CATALYST)),(ligated,),
         geometry=AdapterLigationGeometry.unavailable("The PN-2000094 reagent and sealed junction are published in CG000317 p72; its complete ligation-ready end geometry has not been reconstructed in this model."))
 reviewed.extend(("fragmentation","fragment_gex","gex_endprep","prepare_gex_ends","gex_adaptor","gex_ligase","ligate_gex_adaptor"))
 _indexed_branch(g,"rna",ligated,gex_precursor(),RNA,RNA_SEQ_PRIMERS,reviewed,"Gene-expression library")
 _indexed_branch(g,"sec",amplified_sec,secretion_precursor(),SEC,SEC_SEQ_PRIMERS,reviewed,"Secretion-tag library")
 g.declare_reviewed_scope(reviewed,SOURCE)
 return g.require_publishable()
