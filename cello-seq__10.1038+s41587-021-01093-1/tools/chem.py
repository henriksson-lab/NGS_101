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
          "silently mixed into the original chemistry. ONT motor-adapter bases remain proprietary. "
          "INFERRED: the Methods state only that the annealed splint is ligated to first-strand "
          "cDNA (Fig. 1A draws it at the 5′-phosphorylated oligo(dT) end); which bases template "
          "the nick between the bottom-splint 3′ T15 and the cDNA 5′ phosphate is not described, "
          "so the drawn junction is a model interpretation. The CELLO-seq PCR primers are named "
          "only as 'forward and reverse PCR primers'; the Nextera R1/R2 pair matching the TSO "
          "and bottom-splint handles is inferred.")

NEXTERA_R1 = "CGTCGGCAGCGTCAGATGTGTATAAGAGACAG"
NEXTERA_R2 = "GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG"
DT_HANDLE = "TGGCGTAGCGGGTTCGAGCGCACCGCAGGGTATCCGGCTA"
DT_TRACT = "T" * 25 + "VN"  # Table S5 oligo_dT: 40-nt handle + T25VN
DT_OLIGO = DT_HANDLE + DT_TRACT
TSO_DNA = "TCGTCGGCAGCGTCAGATGTGTATAAGAGACA"
TSO_WRITTEN = TSO_DNA + "rGrGrGAmn"
TOP_SPLINT = "TAGCCGGATACCCTGCGGTGCGCTCGAACCCGCTACGCCATGTAAAACGACGGCCAGT"
BARCODE_1 = "AAGAAAGTTGTCGGTGTCTTTGTG"
SPLINT_CONSTANT = "AGTGGTATC"
UMI_PATTERN = "NRY" * 7 + "N"
BOTTOM_SPLINT = (NEXTERA_R2 + BARCODE_1 + SPLINT_CONSTANT + UMI_PATTERN
                 + revcomp(TOP_SPLINT) + "T" * 15)  # Table S5 idx1: terminal T15

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
    Segment("splint terminal dT15", "T" * 15),
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
SEQUENCING_MECHANISM = "nanopore_motor"
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
            Segment("poly(A)", "A" * 25)]
    primer = [Segment("phosphorylated dT handle", DT_HANDLE),
              Segment("dT25", "T" * 25), Segment("VN anchor", "VN", placeholder=True)]
    sc = Scene(); sc.strand("mRNA", mrna, label="polyadenylated full-length RNA")
    sc.anneal("dT oligo", primer, to="mRNA", pair=("dT25", "poly(A)"),
              unpaired=("phosphorylated dT handle", "VN anchor"),
              label="phosphorylated oligo(dT)", mod5="/5Phos/")
    sc.arrow("dT oligo", "SuperScript IV RT with betaine and ET SSB, 57 °C 10 min then 42 °C 120 min")
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
        Segment("terminal dT15", "T" * 15),
    ]
    sc = Scene(); sc.strand("top splint", top, label="top splint", mod3="/3InvdT/")
    sc.anneal("bottom splint", bottom, to="top splint",
              pair=("identifier-arm complement", "identifier splint arm"),
              unpaired=("Nextera R2 tail", "idx1 cell barcode", "barcode/UMI boundary",
                        "22-nt NRY UMI", "terminal dT15"), label="bottom splint idx1")
    sc.ligation_end("bottom splint", "terminal dT15", kind="splinted",
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
             Segment("terminal dT15", "T" * 15),
             Segment("phosphorylated oligo(dT) and full-length cDNA", "X" * 40,
                     placeholder=True),
             Segment("template-switch end", revcomp(TSO_DNA), "tso")]
    sc = Scene(); sc.strand("ligated cDNA", parts, label="identifier-ligated full-length cDNA")
    sc.junction("ligated cDNA", "terminal dT15",
                "phosphorylated oligo(dT) and full-length cDNA", "HiFi Taq ligation")
    sc.labels("ligated cDNA")
    return sc


def workflow() -> Workflow:
    wf = Workflow(MolecularState("Poly(A)-primed RNA", tuple(_rt_scene().rows())))
    wf.react("Reverse-transcribe at high temperature and template-switch",
             _template_switch_scene().rows(), name="Full-length first-strand cDNA",
             note="The TSO carries the exact Nextera R1-derived handle and a 3′ amino block that suppresses TSO-primed artefacts. Exonuclease I (37 °C 1 h) then digests residual oligos and 95 °C 3 min denatures enzyme and cDNA.")
    geometry = AdapterLigationGeometry.known(
        "splinted", splint_scene(),
        statement="Supplementary Table 5 prints both splint strands, including idx1, the 22-nt NRY UMI, the inverted-dT block and terminal T15. The nick template joining the T15 3′ end to the cDNA 5′ phosphate is not described in the Methods (model interpretation of Fig. 1A).")
    wf.ligate_adapter("Anneal the two splint strands and ligate identifiers to cDNA",
                      geometry, _ligated_scene().rows(),
                      name="Cell- and molecule-indexed cDNA",
                      note="Splints are pre-annealed (95 °C, 1 °C steps to 4 °C); HiFi Taq ligase with NAD+ ligates at 55 °C 1 h. INFERRED junction: the bottom-splint 3′ T15 is drawn joined to the phosphorylated oligo(dT)-derived cDNA end; the nick template is not described. Proteinase K, pooling and 0.5× AMPure XP clean-up follow.")
    wf.react("PCR-amplify deliberately to create UMI families",
             _ligated_scene().rows(), name="Duplicate family for consensus correction",
             note="KAPA HiFi Uracil+ with betaine and ET SSB; cycle number from a proxy qPCR (25–29 cycles, 10 min extension), then 0.5× AMPure XP. INFERRED primers: the Methods say only forward and reverse PCR primers; the Nextera R1 and R2 sites flank the molecule. High duplicate depth is intentional.")
    proprietary = AdapterLigationGeometry.unavailable(
        "200 fmol prepared with ONT Ligation Sequencing Kit SQK-LSK109 (hiPSC, 2-cell blastomeres) or SQK-LSK110 (mESC) and run on MinION R9.4.1; the kit's motor-adapter end sequence and overhang geometry are proprietary and the kit guide was not read.")
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
    rt = graph.input("reverse_transcriptase", "SuperScript IV reverse transcriptase",
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
        statement="Supplementary Table 5 prints both splint strands, the idx1 barcode, 22-nt NRY UMI, inverted-dT block and terminal T15; the nick template joining the T15 3′ end to the cDNA 5′ phosphate is not described in the Methods (model interpretation of Fig. 1A).")
    graph.react("identifier_ligation", "Splint-ligate the cell barcode and UMI to cDNA",
                StepKind.LIGATION,
                ((switched, InputRole.TEMPLATE), (splint, InputRole.STRUCTURAL),
                 (ligase, InputRole.CATALYST)), (ligated,), geometry=geometry)
    primer_note = source + " (Methods name only 'forward and reverse PCR primers'; Nextera R1/R2 from Table S5 inferred from the TSO and bottom-splint handles)"
    r1 = graph.input("nextera_r1", "Nextera R1 PCR primer (inferred identity)",
                     MaterialKind.OLIGO, inferred=True, source=primer_note)
    r2 = graph.input("nextera_r2", "Nextera R2 PCR primer (inferred identity)",
                     MaterialKind.OLIGO, inferred=True, source=primer_note)
    polymerase = graph.input("pcr_polymerase", "KAPA HiFi Uracil+ hot-start polymerase",
                             MaterialKind.CATALYST)
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
