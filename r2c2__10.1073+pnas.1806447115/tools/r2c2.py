"""Molecular model for the original R2C2 full-length-cDNA protocol."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))

from chemdraw import (Construct, MolecularState, Scene, Segment, Workflow,
                      annotation_rows, feature)
from adapter_ligation import AdapterLigationGeometry
from circular import CircularizationGeometry, ClosedCircle
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind)

TITLE = "R2C2 — rolling circle to concatemeric consensus"
CAVEAT = "The illustrated ONT 1D ligation finish is the SIRV E2 branch. The paper used RAD4 rapid preparation for single B cells; that different adapter installation is not represented here."
SOURCE = ('Defining source: <a href="https://doi.org/10.1073/pnas.1806447115">'
          'Volden et al., <i>PNAS</i> (2018)</a>.')


def seg(name, top, tag=None, **kw):
    return Segment(name, top, tag, **kw)


def duplex_rows(con, label=""):
    unknown = tuple(s.name + "'" for s in con if s.is_role_token())
    return [*Scene.duplex(list(con), label=label, unpaired=unknown).rows(),
            *annotation_rows(con)]


TSO_INDEX = feature("r2c2.tso-index", "sample_index", "combinatorial",
                    note="7-nt index introduced by the Tn5Prime TSO")
PCR_INDEX = feature("r2c2.pcr-index", "sample_index", "combinatorial",
                    note="8-nt Nextera A index introduced during cDNA PCR")

# Supplementary Table S3 prints the orderable RT, TSO, cDNA PCR and splint-PCR
# oligos (transcribed in the research note). The complete lambda amplicon and
# exact strand placement of the indexed cDNA ends have not been reconstructed;
# retain role tokens for those assembled molecules.
CDNA = Construct([
    seg("ISPCR end", "[ISPCR site]", "p5", placeholder=True),
    seg("7-nt TSO index", "NNNNNNN", "cbc", placeholder=True, feature=TSO_INDEX),
    seg("full-length cDNA", "X" * 34, placeholder=True, length_bp=34),
    seg("8-nt Nextera A index", "NNNNNNNN", "cbc", placeholder=True, feature=PCR_INDEX),
    seg("Nextera A end", "[Nextera A site]", "p7", placeholder=True),
], name="Tn5Prime full-length cDNA")

SPLINT = Construct([
    seg("ISPCR overlap", "[ISPCR overlap]", "p5", placeholder=True),
    seg("lambda splint backbone", "[about 200-bp lambda-DNA splint]", placeholder=True,
        length_bp=200),
    seg("Nextera A overlap", "[Nextera A overlap]", "p7", placeholder=True),
], name="Gibson circularisation splint")

# R2C2 uses overlap (Gibson) assembly of dsDNA, not ssDNA ligase circularisation. One
# strand is shown to keep the closed topology readable; every RCA unit is derived from it.
CIRCLE_TEMPLATE = Construct([
    *CDNA.segments,
    seg("splint junction", "[DNA splint]", "r2", placeholder=True, length_bp=200),
], name="R2C2 circular template strand")
CIRCLE = ClosedCircle(CIRCLE_TEMPLATE, "splint junction", "ISPCR end")


def circle_rows_model():
    return CIRCLE.rows("Gibson-assembled covalent circle")


def rca_construct(copies=3):
    """A schematic concatemer whose repeat units come from the circular template."""
    # Exact hexamer landing is stochastic, so show units rather than fabricate one site.
    parts = []
    for i in range(copies):
        for s in CIRCLE_TEMPLATE:
            parts.append(seg(f"copy {i + 1}: {s.name}", s.top, s.tag,
                             placeholder=s.placeholder, inferred=s.inferred,
                             bottom=s.bottom, feature=s.feature, length_bp=s.length_bp))
    return Construct(parts, name=f"{copies}-repeat R2C2 RCA concatemer")


def branched_rca_rows():
    """Model-derived displaced products represent random-primed phi29 branches."""
    sc = Scene()
    sc.strand("primary RCA product", list(rca_construct(3)), label="primary RCA strand")
    sc.strand("displaced branch 1", list(rca_construct(2)), label="displaced branch")
    sc.strand("displaced branch 2", list(rca_construct(1)), label="short displaced branch")
    return sc.rows()


def debranched_rca_rows():
    """A single linear concatemer released from the branched RCA network."""
    product = rca_construct()
    sc = Scene()
    sc.strand("debranched concatemer", list(product), label=product.name)
    sc.labels("debranched concatemer")
    return [*sc.rows(), *annotation_rows(product)]


def sections():
    linear = duplex_rows(CDNA, "full-length cDNA")
    circ = circle_rows_model()
    return [
        ("Gibson assembly with terminally complementary DNA splint", circ,
         "NEBuilder HiFi joins both cDNA ends to the roughly 200-bp lambda-DNA splint."),
        ("Exonuclease I, III and lambda exonuclease", circ,
         "Linear DNA is digested; closed circles are retained."),
        ("Random-hexamer phi29 rolling-circle amplification", branched_rca_rows(),
         "Each circle yields a branched high-molecular-weight tandem-repeat product."),
        ("T7 endonuclease debranching and >2-kb selection", debranched_rca_rows(),
         "Debranched concatemers are the input to the ONT 1D ligation library workflow."),
    ], linear

INITIAL_NAME = "Indexed full-length ds cDNA"
INITIAL_ROWS = tuple(duplex_rows(CDNA, "full-length cDNA"))

def workflow():
    steps, _ = sections()
    wf = Workflow(MolecularState(INITIAL_NAME, INITIAL_ROWS))
    for action, rows, note in steps:
        wf.react(action, rows, note=note)
    geometry = AdapterLigationGeometry.unavailable(
        "The R2C2 paper names the ONT 1D ligation kit, but neither the paper nor its "
        "oligo table discloses the motor and distal adapter terminal bases.")
    wf.ligate_adapter("Attach ONT 1D adapters", geometry, final_rows(),
                      name="motor-loaded R2C2 nanopore library",
                      note="This is a proprietary single-strand nanopore adapter attachment, not an inferred sticky-end sequence.")
    return wf


def final_library():
    core = rca_construct()
    return Construct([
        seg("motor-loaded ONT leader", "[ONT 1D motor adapter]", "r1", placeholder=True),
        *core.segments,
        seg("distal ONT adapter", "[ONT 1D distal adapter]", "r2", placeholder=True),
    ], name="R2C2 nanopore library")


def final_rows():
    lib = final_library()
    unknown = tuple(s.name + "'" for s in lib if s.is_role_token())
    sc = Scene.duplex(list(lib), label="R2C2 concatemer", unpaired=unknown)
    sc.junction("top", "motor-loaded ONT leader", "copy 1: ISPCR end", "adapter ligation")
    sc.junction("top", "copy 3: splint junction", "distal ONT adapter", "adapter ligation")
    return [*sc.rows(), *annotation_rows(lib)]


def protocol_graph():
    graph = ProtocolGraph(TITLE, authority=SOURCE)
    cdna = graph.input("indexed_cdna", "Tn5Prime indexed SIRV E2 full-length ds cDNA",
                       MaterialKind.MOLECULAR, state=CDNA, primary=True)
    splint = graph.input("dna_splint", "terminally complementary lambda-DNA splint",
                         MaterialKind.ADAPTER, state=SPLINT)
    assembly = graph.input("hifi_assembly", "NEBuilder HiFi assembly enzyme mix",
                           MaterialKind.CATALYST)
    circle = graph.product("circle", "Gibson-assembled full-length cDNA circle",
                           MaterialKind.MOLECULAR, state=CIRCLE)
    graph.react("assemble_circle", "Join both cDNA ends to splint (NEBuilder, 55°C, 60 min)",
                StepKind.CIRCULARIZATION,
                ((cdna, InputRole.TEMPLATE), (splint, InputRole.STRUCTURAL),
                 (assembly, InputRole.CATALYST)), (circle,),
                geometry=CircularizationGeometry.known(
                    "overlap assembly",
                    statement="The source uses NEBuilder HiFi and a roughly 200-bp DNA splint with terminal overlaps matching both cDNA ends."))
    exonucleases = graph.input("exonucleases", "exonuclease I, III and lambda exonuclease",
                               MaterialKind.CATALYST)
    selected_circle = graph.product("selected_circle", "exonuclease-resistant R2C2 circle",
                                    MaterialKind.MOLECULAR, state=CIRCLE)
    linear_waste = graph.product("linear_waste", "digested non-circular DNA",
                                 MaterialKind.WASTE)
    graph.react("digest_linear", "Digest non-circular DNA, then 0.8x SPRI-purify circles above 500 bp", StepKind.DIGESTION,
                ((circle, InputRole.TEMPLATE), (exonucleases, InputRole.CATALYST)),
                (selected_circle, linear_waste))
    graph.discard(linear_waste, "is not protected by covalent circle closure")
    random_hexamers = graph.input("random_hexamers", "exonuclease-resistant random hexamers",
                                  MaterialKind.OLIGO)
    phi29 = graph.input("phi29", "phi29 strand-displacing polymerase",
                        MaterialKind.CATALYST)
    branched = graph.product("branched_rca", "branched R2C2 rolling-circle product",
                             MaterialKind.MOLECULAR, state=tuple(branched_rca_rows()))
    graph.react("rolling_circle", "Random-prime rolling-circle amplification (30°C, overnight)",
                StepKind.EXTENSION,
                ((selected_circle, InputRole.TEMPLATE),
                 (random_hexamers, InputRole.STRUCTURAL),
                 (phi29, InputRole.CATALYST)), (branched,))
    pre_size = graph.input("pre_debranch_selection", "0.5x SPRI selection above 2 kb before debranching", MaterialKind.SELECTION)
    retained_branched = graph.product("retained_branched", "size-selected branched RCA", MaterialKind.MOLECULAR, state=branched.state)
    pre_short = graph.product("pre_debranch_short", "short RCA products before debranching", MaterialKind.WASTE)
    graph.react("select_before_debranch", "Select RCA products above 2 kb (0.5x SPRI) before T7 treatment", StepKind.SELECTION,
                ((branched, InputRole.TEMPLATE), (pre_size, InputRole.SELECTION)), (retained_branched, pre_short))
    graph.discard(pre_short, "physical SPRI size-selection fraction")
    t7 = graph.input("t7_endonuclease", "T7 endonuclease I", MaterialKind.CATALYST)
    debranched = graph.product("debranched", "debranched R2C2 concatemers",
                               MaterialKind.MOLECULAR, state=rca_construct())
    graph.react("debranch", "Debranch RCA with T7 endonuclease (37°C, 2 h)",
                StepKind.DIGESTION,
                ((retained_branched, InputRole.TEMPLATE), (t7, InputRole.CATALYST)),
                (debranched,))
    size_matrix = graph.input("size_selection", ">2-kb SPRI size-selection condition",
                              MaterialKind.SELECTION)
    long_concatemers = graph.product("long_concatemers", ">2-kb debranched concatemers",
                                     MaterialKind.MOLECULAR, state=rca_construct())
    short_waste = graph.product("short_products", "short RCA and debranching products",
                                MaterialKind.WASTE)
    graph.react("select_long", "Select debranched products above 2 kb",
                StepKind.SELECTION,
                ((debranched, InputRole.TEMPLATE), (size_matrix, InputRole.SELECTION)),
                (long_concatemers, short_waste))
    graph.discard(short_waste, "falls below the retained concatemer size range")
    ont_adapter = graph.input(
        "ont_adapter", "historical ONT 1D motor and distal adapters",
        MaterialKind.ADAPTER, inferred=True,
        source=SOURCE + " (historical adapter bases unassigned)")
    ligase = graph.input("adapter_ligase", "ONT library-adapter ligase",
                         MaterialKind.CATALYST)
    library_construct = final_library()
    library = graph.product(
        "sequencing_library", library_construct.name, MaterialKind.LIBRARY,
        state=library_construct, inferred=True,
        source=SOURCE + " (historical ONT adapter bases unassigned)")
    graph.react("attach_ont_adapters", "Attach ONT 1D sequencing adapters",
                StepKind.LIGATION,
                ((long_concatemers, InputRole.TEMPLATE),
                 (ont_adapter, InputRole.STRUCTURAL),
                 (ligase, InputRole.CATALYST)), (library,),
                geometry=AdapterLigationGeometry.unavailable(
                    "The source identifies the ONT 1D kit for SIRV; its selected historical adapter-end geometry has not been verified in this model."))
    graph.sequence(library, SequencingPlan(
        acquired_roles=("Sequencing output",),
        limitation="Nanopore entry is motor-mediated and primer-free; selected historical adapter bases remain unassigned."))
    graph.declare_reviewed_scope((
        "indexed_cdna", "dna_splint", "hifi_assembly", "assemble_circle",
        "exonucleases", "digest_linear", "random_hexamers", "phi29",
        "rolling_circle", "pre_debranch_selection", "select_before_debranch", "t7_endonuclease", "debranch", "size_selection",
        "select_long", "ont_adapter", "adapter_ligase", "attach_ont_adapters",
        "sequencing_library"), SOURCE)
    return graph.require_publishable()
