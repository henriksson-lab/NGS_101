"""Conventional ligation-first MethylC-seq / WGBS model."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
import illumina as il
import seqprimers as sp
from adapter_ligation import AdapterLigationGeometry
from chemdraw import Construct, MolecularState, Row, Scene, Segment, Workflow, feature
from pcr import PCRGeometry
from protocol_graph import (Disposition, Input, InputRole, Material, MaterialKind,
                            ProtocolGraph, SequencingPlan, Step, StepKind, Terminal)

TITLE = "MethylC-seq — conventional whole-genome bisulfite sequencing"
NOTES = "01_methylc-seq.html"
SOURCE = 'Protocol source: <a href="https://doi.org/10.1038/nprot.2014.114">Urich et al., <i>Nature Protocols</i> (2015)</a>.'
SUMMARY = "Fragment genomic DNA, repair and dA-tail it, ligate fully methylated adapters, then denature and bisulfite-convert before limited uracil-tolerant PCR."
CAVEAT = ("The paper specifies NEXTflex Bisulfite-Seq Barcodes–12 (Bioo Scientific, "
          "cat. 511912), not an unspecified adapter. It is a single-index, 6-nt i7, "
          "TruSeq-compatible library: standard Read 1, Index 1 and Read 2 sequencing "
          "primers bind; there is no i5 index.")

# The 12 six-base indices supplied with the product family.  They are kept as data rather
# than collapsed to an 8-base generic barcode so the library constructor enforces the
# actual single-index geometry used by this protocol.
NEXTFLEX_I7 = (
    "CGATGT", "TGACCA", "ACAGTG", "GCCAAT", "CAGATC", "CTTGTA",
    "ATCACG", "TTAGGC", "ACTTGA", "GATCAG", "TAGCTT", "GGCTAC",
)


def nextflex_library(index: str = NEXTFLEX_I7[0], *,
                     insert_label: str = "C→T-converted genomic insert") -> tuple[Construct, tuple]:
    """Exact NEXTflex/TruSeq-compatible single-index endpoint used by MethylC-seq.

    P5 and the Read 1 site share their terminal/initial ACAC.  Modelling that overlap
    explicitly avoids inventing four bases, while the dA junction supplies the leading A
    of the reverse-complemented Read 2 site on the opposite adapter arm.
    """
    if index not in NEXTFLEX_I7:
        raise ValueError("index is not one of the 12 NEXTflex Bisulfite-Seq barcodes")
    lib = Construct([
        Segment("P5 before shared ACAC", il.P5[:-4], "p5"),
        Segment("Read 1 site (shares ACAC with P5)", il.TRUSEQ_READ1, "r1"),
        Segment(insert_label, "X" * 36, placeholder=True),
        Segment("dA junction / Read 2 site start", "A"),
        Segment("Index 1 / Read 2 arm", il.INDEX1_PRIMER, "r2"),
        Segment("6-nt i7", index, "cbc",
                feature=feature("sample_i7", "sample_index", "whitelist", whitelist="NEXTflex Bisulfite-Seq Barcodes-12")),
        Segment("P7 reverse complement", il.P7_RC, "p7"),
    ], name="NEXTflex MethylC-seq library")
    primers = (sp.TRUSEQ["R1"], sp.TRUSEQ["I1"], sp.TRUSEQ["R2"])
    problems = sp.verify(lib, primers, required_roles=tuple(p.role for p in primers))
    if problems:
        raise ValueError("invalid NEXTflex MethylC-seq endpoint: " + "; ".join(problems))
    return lib, primers


FINAL_LIBRARY, SEQ_PRIMERS = nextflex_library()
FINAL_CAPTION = ("PCR-restored NEXTflex library: a single 6-nt i7 and the standard "
                 "TruSeq-compatible sequencing-primer sites. Unmethylated genomic C is "
                 "read as T; protected adapter 5mC remains C through bisulfite treatment.")
SEQUENCING_INTRO = ("Standard TruSeq Read 1, Index 1 and Read 2 primers bind to the "
                    "NEXTflex library. There is no Index 2 read. Paired reads traverse "
                    "the converted insert, whose C/T outcome is compared with the "
                    "unconverted reference sequence.")


def ligated_nextflex_scene() -> Scene:
    ligated, _ = nextflex_library(insert_label="unconverted genomic insert")
    sc = Scene.duplex(list(ligated), label="NEXTflex adapter-ligated molecule")
    sc.junction("top", "Read 1 site (shares ACAC with P5)",
                "unconverted genomic insert", "adapter ligation")
    sc.junction("top", "dA junction / Read 2 site start",
                "Index 1 / Read 2 arm", "adapter ligation")
    sc.labels("top")
    return sc

def sections():
    wf=Workflow(MolecularState("Fragmented genomic DNA",(Row(chunks=[("~200-bp genomic fragments",None,False)]),)))
    wf.react("End-repair and dA-tail",[Row(chunks=[("blunt 5′-phosphorylated ends → 3′ dA",None,False)])],note="Mechanical fragmentation is followed by conventional end preparation.")
    wf.ligate_adapter("Ligate protected NEXTflex adapters before conversion",
        AdapterLigationGeometry.unavailable("The protocol identifies NEXTflex Bisulfite-Seq Barcodes–12 and dA-tailed inserts, but the accessible source does not print the ligation adapter duplex or explicitly report its overhang."),
        ligated_nextflex_scene().rows(),note="The PCR-restored endpoint is TruSeq-compatible; ** marks the final adapter/insert boundaries, not evidence for an unprinted reagent end.")
    wf.react("Convert and amplify",[Row(chunks=[("unmethylated C → U → T in PCR     5mC → C","w1",False)]),Row(chunks=[("adapter-ligated duplex → denatured non-complementary single strands → PCR duplex",None,False)])],note="Bisulfite conversion occurs after adapter ligation and can break already tagged molecules.")
    return [(r.action,r.after.rows,r.note) for r in wf.reactions]


def protocol_graph() -> ProtocolGraph:
    """Conventional pre-bisulfite adapter-ligation provenance graph."""
    source = "Urich et al. 2015, doi:10.1038/nprot.2014.114"
    graph = ProtocolGraph(TITLE, authority=source)

    def add_input(key, label, kind, *, state=None, primary=False):
        value = Material(key, label, kind, state=state, primary=primary, source=source)
        graph.add_input(value)
        return value

    dna = add_input("genomic_dna", "genomic DNA", MaterialKind.BIOLOGICAL,
                    state=(Row(chunks=[("high-molecular-weight genomic DNA", None, False)]),),
                    primary=True)
    shearing = add_input("mechanical_shearing", "mechanical shearing condition",
                         MaterialKind.CONDITION)
    fragmented = Material(
        "fragmented_dna", "mechanically fragmented genomic DNA", MaterialKind.MOLECULAR,
        state=(Row(chunks=[("~200-bp genomic fragments", None, False)]),), source=source)
    graph.add_step(Step(
        "fragment", "Mechanically fragment genomic DNA", StepKind.FRAGMENTATION,
        (Input(dna.id, InputRole.TEMPLATE), Input(shearing.id, InputRole.CONDITION)),
        (fragmented.id,), source=source), (fragmented,))

    repair_enzymes = add_input("end_repair_enzymes", "end-repair enzyme mix",
                               MaterialKind.CATALYST)
    repaired = Material(
        "end_repaired_dna", "blunt 5-prime-phosphorylated genomic fragments",
        MaterialKind.MOLECULAR,
        state=(Row(chunks=[("blunt 5′-phosphorylated genomic fragments", None, False)]),),
        source=source)
    graph.add_step(Step(
        "end_repair", "End-repair fragmented DNA", StepKind.END_REPAIR,
        (Input(fragmented.id, InputRole.TEMPLATE),
         Input(repair_enzymes.id, InputRole.CATALYST)),
        (repaired.id,), source=source), (repaired,))

    tailing_enzyme = add_input("a_tailing_enzyme", "dA-tailing enzyme mix",
                               MaterialKind.CATALYST)
    tailed_scene = Scene.duplex([
        Segment("left 3′ dA", "A"),
        Segment("genomic fragment", "X" * 36, placeholder=True),
        Segment("right 3′ dA", "A"),
    ], label="dA-tailed genomic DNA")
    tailed = Material("da_tailed_dna", "dA-tailed genomic fragments",
                      MaterialKind.MOLECULAR, state=tailed_scene, source=source)
    graph.add_step(Step(
        "da_tail", "Add one 3-prime dA to each end", StepKind.TAILING,
        (Input(repaired.id, InputRole.TEMPLATE),
         Input(tailing_enzyme.id, InputRole.CATALYST)),
        (tailed.id,), source=source), (tailed,))

    adapters = add_input("nextflex_adapters", "protected NEXTflex Bisulfite-Seq adapters",
                         MaterialKind.ADAPTER)
    ligase = add_input("dna_ligase", "DNA ligase", MaterialKind.CATALYST)
    ligated = Material("ligated_dna", "NEXTflex-adapter-ligated genomic DNA",
                       MaterialKind.MOLECULAR, state=ligated_nextflex_scene(), source=source)
    missing = AdapterLigationGeometry.unavailable(
        "The protocol identifies NEXTflex Bisulfite-Seq Barcodes–12 and dA-tailed inserts, but does not print the ligation adapter duplex or its overhang.")
    graph.add_step(Step(
        "adapter_ligation", "Ligate protected NEXTflex adapters before conversion",
        StepKind.LIGATION,
        (Input(tailed.id, InputRole.TEMPLATE),
         Input(adapters.id, InputRole.STRUCTURAL),
         Input(ligase.id, InputRole.CATALYST)),
        (ligated.id,), geometry=missing, source=source), (ligated,))

    bisulfite = add_input("bisulfite_reagent", "bisulfite conversion reagent",
                          MaterialKind.CONDITION)
    converted_template = Construct([
        Segment("protected left adapter", "[adapter]", placeholder=True),
        Segment("C→U-converted genomic strand", "X" * 36, placeholder=True),
        Segment("protected right adapter", "[adapter]", placeholder=True),
    ], name="bisulfite-converted adapter-ligated template")
    converted = Material("converted_template", "bisulfite-converted ligated template",
                         MaterialKind.MOLECULAR, state=converted_template, source=source)
    graph.add_step(Step(
        "bisulfite_conversion", "Denature and bisulfite-convert the ligated DNA",
        StepKind.CHEMICAL_CONVERSION,
        (Input(ligated.id, InputRole.TEMPLATE),
         Input(bisulfite.id, InputRole.CONDITION)),
        (converted.id,), source=source), (converted,))

    forward = add_input("nextflex_forward_primer", "NEXTflex library PCR forward primer",
                        MaterialKind.OLIGO)
    reverse = add_input("nextflex_reverse_primer", "NEXTflex indexed PCR reverse primer",
                        MaterialKind.OLIGO)
    polymerase = add_input("uracil_tolerant_polymerase", "KAPA HiFi uracil+ Readymix",
                           MaterialKind.CATALYST)
    library = Material("methylc_library", "MethylC-seq sequencing library",
                       MaterialKind.LIBRARY, state=FINAL_LIBRARY, source=source)
    graph.add_step(Step(
        "library_pcr", "Four-cycle uracil-tolerant library PCR", StepKind.PCR,
        (Input(converted.id, InputRole.TEMPLATE),
         Input(forward.id, InputRole.STRUCTURAL),
         Input(reverse.id, InputRole.STRUCTURAL),
         Input(polymerase.id, InputRole.CATALYST)),
        (library.id,), geometry=PCRGeometry.unavailable(
            converted_template,
            reason="The protocol names the commercial NEXTflex primer mix but does not print one selected primer pair."),
        source=source), (library,))
    graph.mark_terminal(Terminal(library.id, Disposition.SEQUENCING), SequencingPlan(
        FINAL_LIBRARY, tuple(SEQ_PRIMERS), tuple(p.role for p in SEQ_PRIMERS)))

    graph.declare_reviewed_scope((
        "genomic_dna", "mechanical_shearing", "fragment", "end_repair_enzymes",
        "end_repair", "a_tailing_enzyme", "da_tail", "nextflex_adapters",
        "dna_ligase", "adapter_ligation", "bisulfite_reagent",
        "bisulfite_conversion", "nextflex_forward_primer", "nextflex_reverse_primer",
        "uracil_tolerant_polymerase", "library_pcr", "methylc_library",
    ), source)
    return graph.require_publishable()
