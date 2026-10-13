"""Original 22-nt-UMI CELLO-seq chemistry from the defining paper."""
from adapter_ligation import AdapterLigationGeometry
from chemdraw import (Construct, MolecularState, Scene, Segment, Workflow,
                      duplex_rows, feature, revcomp)
from pcr import PCRGeometry
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind)


TITLE = "CELLO-seq — long-UMI single-cell long-read RNA-seq"
NOTES = "01_cello-seq.html"
SOURCE = ('Defining source: <a href="https://doi.org/10.1038/s41587-021-01093-1">'
          'Berrens et al., <i>Nature Biotechnology</i> (2022)</a>. An '
          '<a href="https://doi.org/10.1038/s41596-025-01203-2">updated 2025 protocol</a> '
          'replaces the original 22-nt UMI with a 50-nt UMI.')
SUMMARY = ("A phosphorylated oligo(dT) and amino-blocked TSO make full-length first-strand "
           "cDNA. An exact two-strand splint then ligates a 24-nt cell barcode and "
           "22-nt RYN-patterned UMI onto its 5′ end before deliberate duplicate generation "
           "and Oxford Nanopore sequencing.")
CAVEAT = ("This page represents the defining 2022 CELLO-seq protocol and its 22-nt UMI. "
          "The 2025 protocol publishes a revised 50-nt-UMI bottom splint and should not be "
          "silently mixed into the original chemistry. ONT motor-adapter bases remain proprietary.")

NEXTERA_R1 = "CGTCGGCAGCGTCAGATGTGTATAAGAGACAG"
NEXTERA_R2 = "GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG"
DT_HANDLE = "TGGCGTAGCGGGTTCGAGCGCACCGCAGGGTATCCGGCTA"
DT_TRACT = "T" * 30 + "VN"
DT_OLIGO = DT_HANDLE + DT_TRACT
TSO_DNA = "TCGTCGGCAGCGTCAGATGTGTATAAGAGACA"
TSO_WRITTEN = TSO_DNA + "rGrGrGAmn"
TOP_SPLINT = "TAGCCGGATACCCTGCGGTGCGCTCGAACCCGCTACGCCATGTAAAACGACGGCCAGT"
BARCODE_1 = "AAGAAAGTTGTCGGTGTCTTTGTG"
SPLINT_CONSTANT = "AGTGGTATC"
UMI_PATTERN = "NRY" * 7 + "N"
BOTTOM_SPLINT = (NEXTERA_R2 + BARCODE_1 + SPLINT_CONSTANT + UMI_PATTERN
                 + revcomp(TOP_SPLINT) + "T" * 14)

CELL = feature("cello_cell_barcode", "cell_barcode", "whitelist",
               whitelist="Berrens et al. Supplementary Table 5, splintlig_spc_idx1–96")
UMI = feature("cello_umi", "umi", "random",
              note="original 22-nt NRY repeat design")

FINAL_LIBRARY = Construct([
    Segment("proprietary motor-loaded ONT adapter", "[motor adapter]", "r2",
            placeholder=True, inferred=True),
    Segment("Nextera R2 priming site", NEXTERA_R2, "r2"),
    Segment("cell barcode 1 representative", BARCODE_1, "cbc", feature=CELL),
    Segment("barcode/UMI boundary", SPLINT_CONSTANT),
    Segment("22-nt RYN UMI", UMI_PATTERN, "umi", placeholder=True, feature=UMI),
    Segment("splint-derived junction", revcomp(TOP_SPLINT)),
    Segment("splint terminal dT14", "T" * 14),
    Segment("oligo(dT) handle", DT_HANDLE),
    Segment("oligo(dT)-derived tract", DT_TRACT, placeholder=True),
    Segment("full-length cDNA", "X" * 42, placeholder=True),
    Segment("template-switch / Nextera R1 site reverse complement",
            revcomp(TSO_DNA), "tso"),
], name="ONT-adapted original CELLO-seq cDNA")
FINAL_CAPTION = ("Original CELLO-seq identifier end, with the exact idx1 splint member and "
                 "22-nt NRY UMI, joined to full-length cDNA. Only the terminal ONT motor "
                 "adapter is inferred because its bases are proprietary.")
SEQ_PRIMERS = ()
SEQUENCING_ENDING = ("There is no synthesis sequencing primer. The ONT motor meters the "
                     "adapted strand through the pore; the exact cell barcode and patterned "
                     "UMI are recovered inline and UMI families form corrected consensuses.")


def oligos():
    return (
        f"<b>Phosphorylated oligo(dT)</b> — 5′-/5Phos/{DT_OLIGO}-3′",
        f"<b>3′-amino-blocked TSO</b> — 5′-{TSO_WRITTEN}-3′",
        f"<b>Top splint</b> — 5′-{TOP_SPLINT}/3InvdT/-3′",
        f"<b>Bottom splint, idx1</b> — 5′-{BOTTOM_SPLINT}-3′",
        f"<b>Nextera R1 oligo</b> — 5′-{NEXTERA_R1}-3′",
        f"<b>Nextera R2 oligo</b> — 5′-{NEXTERA_R2}-3′",
    )


def _rt_scene() -> Scene:
    mrna = [Segment("transcript", "X" * 34, placeholder=True),
            Segment("poly(A)", "A" * 30)]
    primer = [Segment("phosphorylated dT handle", DT_HANDLE),
              Segment("dT30", "T" * 30), Segment("VN anchor", "VN", placeholder=True)]
    sc = Scene(); sc.strand("mRNA", mrna, label="polyadenylated full-length RNA")
    sc.anneal("dT oligo", primer, to="mRNA", pair=("dT30", "poly(A)"),
              unpaired=("phosphorylated dT handle", "VN anchor"),
              label="phosphorylated oligo(dT)", mod5="/5Phos/")
    sc.arrow("dT oligo", "high-temperature reverse transcription")
    sc.labels("dT oligo")
    return sc


def _template_switch_scene() -> Scene:
    cDNA = [Segment("oligo(dT) primer and cDNA", "X" * 38, placeholder=True),
            Segment("MMLV-added CCC", "CCC")]
    tso = [Segment("Nextera R1 / TSO handle", TSO_DNA, "tso"),
           Segment("three RNA G", "GGG")]
    sc = Scene(); sc.strand("first-strand cDNA", cDNA, label="full-length first strand")
    sc.anneal("TSO", tso, to="first-strand cDNA",
              pair=("three RNA G", "MMLV-added CCC"),
              unpaired=("Nextera R1 / TSO handle",), label="3′-amino-blocked TSO",
              mod3="/3Amn/")
    sc.arrow("TSO", "template switching")
    sc.labels("TSO")
    return sc


def splint_scene() -> Scene:
    top = [Segment("cDNA-handle splint arm", TOP_SPLINT[:40]),
           Segment("identifier splint arm", TOP_SPLINT[40:])]
    bottom = [
        Segment("Nextera R2 tail", NEXTERA_R2, "r2"),
        Segment("idx1 cell barcode", BARCODE_1, "cbc", feature=CELL),
        Segment("barcode/UMI boundary", SPLINT_CONSTANT),
        Segment("22-nt NRY UMI", UMI_PATTERN, "umi", placeholder=True, feature=UMI),
        Segment("identifier-arm complement", revcomp(TOP_SPLINT[40:])),
        Segment("cDNA-handle-arm complement", revcomp(TOP_SPLINT[:40])),
        Segment("terminal dT14", "T" * 14),
    ]
    sc = Scene(); sc.strand("top splint", top, label="top splint", mod3="/3InvdT/")
    sc.anneal("bottom splint", bottom, to="top splint",
              pair=("identifier-arm complement", "identifier splint arm"),
              unpaired=("Nextera R2 tail", "idx1 cell barcode", "barcode/UMI boundary",
                        "22-nt NRY UMI", "terminal dT14"), label="bottom splint idx1")
    sc.ligation_end("bottom splint", "terminal dT14", kind="splinted",
                    polarity="3-prime", owner="adapter",
                    text="3′ splinted ligation end presented to phosphorylated cDNA")
    sc.labels("bottom splint")
    return sc


def _ligated_scene() -> Scene:
    parts = [Segment("Nextera R2 tail", NEXTERA_R2, "r2"),
             Segment("idx1 cell barcode", BARCODE_1, "cbc", feature=CELL),
             Segment("barcode/UMI boundary", SPLINT_CONSTANT),
             Segment("22-nt NRY UMI", UMI_PATTERN, "umi", placeholder=True, feature=UMI),
             Segment("splint-derived junction", revcomp(TOP_SPLINT)),
             Segment("terminal dT14", "T" * 14),
             Segment("phosphorylated oligo(dT) and full-length cDNA", "X" * 40,
                     placeholder=True),
             Segment("template-switch end", revcomp(TSO_DNA), "tso")]
    sc = Scene(); sc.strand("ligated cDNA", parts, label="identifier-ligated full-length cDNA")
    sc.junction("ligated cDNA", "terminal dT14",
                "phosphorylated oligo(dT) and full-length cDNA", "HiFi Taq ligation")
    sc.labels("ligated cDNA")
    return sc


def workflow() -> Workflow:
    wf = Workflow(MolecularState("Poly(A)-primed RNA", tuple(_rt_scene().rows())))
    wf.react("Reverse-transcribe at high temperature and template-switch",
             _template_switch_scene().rows(), name="Full-length first-strand cDNA",
             note="The TSO carries the exact Nextera R1-derived handle and a 3′ amino block that suppresses TSO-primed artefacts.")
    geometry = AdapterLigationGeometry.known(
        "splinted", splint_scene(),
        statement="Supplementary Table 5 prints both splint strands, including idx1, the 22-nt NRY UMI, the inverted-dT block and terminal dT14.")
    wf.ligate_adapter("Anneal the two splint strands and ligate identifiers to cDNA",
                      geometry, _ligated_scene().rows(),
                      name="Cell- and molecule-indexed cDNA",
                      note="HiFi Taq ligase joins the splint-presented 3′ end to the phosphorylated oligo(dT)-derived cDNA end.")
    wf.react("PCR-amplify deliberately to create UMI families",
             _ligated_scene().rows(), name="Duplicate family for consensus correction",
             note="Nextera R1 and R2 sites flank the full-length molecule. High duplicate depth is intentional rather than treated as duplication noise.")
    proprietary = AdapterLigationGeometry.unavailable(
        "The defining paper identifies an Oxford Nanopore ligation kit, but its motor-adapter end sequence and overhang geometry are proprietary.")
    wf.ligate_adapter("Attach the ONT motor adapter", proprietary,
                      duplex_rows(FINAL_LIBRARY, label="ONT library",
                                  unpaired=("proprietary motor-loaded ONT adapter'",)),
                      name="ONT-adapted CELLO-seq library",
                      note="Only the proprietary terminal adapter remains inferred; the internal CELLO-seq oligos are source-exact.")
    return wf


def protocol_graph() -> ProtocolGraph:
    source = "Defining source: https://doi.org/10.1038/s41587-021-01093-1"
    graph = ProtocolGraph("CELLO-seq", authority=source)
    rna = graph.input("rna", "polyadenylated single-cell RNA",
                      MaterialKind.BIOLOGICAL, primary=True)
    dt = graph.input("dt_primer", "5'-phosphorylated oligo(dT) primer",
                     MaterialKind.OLIGO)
    rt = graph.input("reverse_transcriptase", "high-temperature reverse transcriptase",
                     MaterialKind.CATALYST)
    first = graph.product("first_strand", "phosphorylated full-length first-strand cDNA",
                          MaterialKind.MOLECULAR)
    graph.react("reverse_transcription", "Reverse-transcribe full-length polyadenylated RNA",
                StepKind.REVERSE_TRANSCRIPTION,
                ((rna, InputRole.TEMPLATE), (dt, InputRole.STRUCTURAL),
                 (rt, InputRole.CATALYST)), (first,))
    tso = graph.input("tso", "3'-amino-blocked template-switch oligo",
                      MaterialKind.OLIGO)
    switched = graph.product("switched_cdna", "template-switched full-length cDNA",
                             MaterialKind.MOLECULAR)
    graph.react("template_switch", "Install the Nextera R1-derived handle by template switching",
                StepKind.TEMPLATE_SWITCHING,
                ((first, InputRole.TEMPLATE), (tso, InputRole.STRUCTURAL),
                 (rt, InputRole.CATALYST)), (switched,))
    splint = graph.input("identifier_splint", "cell-barcode/22-nt-UMI splint duplex",
                         MaterialKind.ADAPTER, state=splint_scene())
    ligase = graph.input("hifi_ligase", "HiFi Taq DNA ligase",
                         MaterialKind.CATALYST)
    ligated = graph.product("identifier_ligated_cdna", "cell- and molecule-indexed full-length cDNA",
                            MaterialKind.MOLECULAR, state=_ligated_scene())
    geometry = AdapterLigationGeometry.known(
        "splinted", splint_scene(),
        statement="Supplementary Table 5 prints both splint strands, the idx1 barcode, 22-nt NRY UMI, inverted-dT block and terminal dT14.")
    graph.react("identifier_ligation", "Splint-ligate the cell barcode and UMI to cDNA",
                StepKind.LIGATION,
                ((switched, InputRole.TEMPLATE), (splint, InputRole.STRUCTURAL),
                 (ligase, InputRole.CATALYST)), (ligated,), geometry=geometry)
    r1 = graph.input("nextera_r1", "Nextera R1 PCR primer", MaterialKind.OLIGO)
    r2 = graph.input("nextera_r2", "Nextera R2 PCR primer", MaterialKind.OLIGO)
    polymerase = graph.input("pcr_polymerase", "PCR polymerase", MaterialKind.CATALYST)
    families = graph.product("umi_families", "deliberately amplified UMI duplicate families",
                             MaterialKind.MOLECULAR, state=_ligated_scene())
    graph.react("family_pcr", "PCR-amplify to create UMI families for consensus correction",
                StepKind.PCR,
                ((ligated, InputRole.TEMPLATE), (r1, InputRole.STRUCTURAL),
                 (r2, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)),
                (families,), geometry=PCRGeometry.unavailable(
                    FINAL_LIBRARY, reason="The ligated cDNA is drawn as a scene; the two source-printed Nextera primer sites flank it."))
    ont = graph.input("ont_adapter", "motor-loaded Oxford Nanopore ligation adapter",
                      MaterialKind.ADAPTER, inferred=True,
                      source=source + " (motor-adapter bases and end geometry are proprietary)")
    library = graph.product("library", FINAL_LIBRARY.name, MaterialKind.LIBRARY,
                            inferred=True,
                            source=source + " (terminal ONT motor adapter is proprietary)")
    graph.react("ont_ligation", "Attach the Oxford Nanopore motor adapter",
                StepKind.LIGATION,
                ((families, InputRole.TEMPLATE), (ont, InputRole.STRUCTURAL),
                 (ligase, InputRole.CATALYST)), (library,),
                geometry=AdapterLigationGeometry.unavailable(
                    "The defining paper identifies the ONT ligation kit, but its motor-adapter bases and overhang geometry are proprietary."))
    graph.sequence(library, SequencingPlan(
        limitation="Nanopore sequencing uses a proprietary motor adapter rather than a synthesis sequencing primer; its exact bases are unavailable.",
        acquired_roles=("Nanopore strand",),
        read_uses={"Nanopore strand": "include"}))
    graph.declare_reviewed_scope((
        "rna", "dt_primer", "reverse_transcriptase", "reverse_transcription",
        "tso", "template_switch", "identifier_splint", "hifi_ligase",
        "identifier_ligation", "nextera_r1", "nextera_r2", "pcr_polymerase",
        "family_pcr", "ont_adapter", "ont_ligation", "library"), source)
    return graph.require_publishable()
