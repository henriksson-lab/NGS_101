"""Protocol-owned, source-scoped molecular graph; unresolved coordinates stay explicit."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
from chemdraw import Construct, Row, Scene, complement_segments
from circular import CircularizationGeometry
from pcr import PCRGeometry
from page import head, info
from protocol_graph import Disposition, InputRole as R, MaterialKind as K, ProtocolGraph, StepKind as S, SequencingPlan
from protocol_model import Drawing, closed_circle, role
from protocol_views import WorkflowSection, PlanLimitationSection, assemble_page


def duplex(label, *names):
    # Role tokens identify molecular regions, not bases or measurable lengths.
    parts = tuple(role(name, inferred=True) for name in names)
    bottom = complement_segments(parts)
    sc = Scene()
    sc.strand("top", parts, label=label)
    sc.anneal("bottom", bottom, to="top", pair=(bottom[-1].name, parts[0].name),
        unpaired=tuple(seg.name for seg in bottom))
    return sc


def circle(label, *names):
    # One reference strand of the duplex circle is opened at the drawing origin.
    # Its role ordering is inferred; exact closure bases are unassigned.
    return closed_circle(label + "; reference-strand schematic", tuple(role(name, inferred=True) for name in names),
        closure_label="circle closure; exact boundary bases unassigned")

PAGE_FILE = "inverse_pcr.html"
TITLE = "Inverse PCR — flanking-sequence recovery"
SOURCE = "Ochman, Gerber and Hartl, Genetics120:621–623 (1988), doi10.1093/genetics/120.3.621. Current graph covers the defining abstract only; full published methods/Figure1 remain deferred."


def protocol_graph():
    g = ProtocolGraph(TITLE, authority=SOURCE, required_output=Disposition.CONSTRUCT)
    dna = g.input("dna", "genomic DNA containing a known region", K.MOLECULAR,
        state=duplex("INFERRED schematic locus", "unknown left flank", "known region", "unknown right flank"), primary=True)
    nuclease = g.input("nuclease", "restriction enzyme; original methods not assigned", K.CATALYST)
    fragment = g.product("fragment", "restriction fragment with known sequence and flanks", K.MOLECULAR,
        state=duplex("INFERRED restriction-fragment order", "left flank", "known region", "right flank"))
    g.react("digest", "INFERRED fragment layout — restrict genomic DNA", S.DIGESTION,
        [(dna,R.TEMPLATE),(nuclease,R.CATALYST)], [fragment])
    ligase = g.input("ligase", "circle-closing ligase; identity deferred", K.CATALYST)
    circular = g.product("circular", "self-ligated restriction-fragment circle", K.MOLECULAR,
        state=circle("INFERRED circular fragment layout", "left flank", "known region", "right flank"))
    g.react("circularize", "INFERRED circle layout — self-ligate the restriction fragment", S.CIRCULARIZATION,
        [(fragment,R.TEMPLATE),(ligase,R.CATALYST)], [circular],
        geometry=CircularizationGeometry.unavailable("The defining abstract establishes self-ligation into a circle. End chemistry and exact junction coordinates await the original full methods; no blunt/cohesive end is invented."))
    forward = g.input("forward", "known-region outward primer A; sequence deferred", K.OLIGO)
    reverse = g.input("reverse", "known-region outward primer B; sequence deferred", K.OLIGO)
    polymerase = g.input("polymerase", "PCR polymerase; original identity deferred", K.CATALYST)
    product = g.product("amplicon", "linear flanking-sequence amplicon", K.MOLECULAR,
        state=duplex("INFERRED circle-spanning amplicon order", "known-region boundary A", "flanking DNA across restriction-fragment closure", "known-region boundary B"))
    g.react("ipcr", "INFERRED amplicon layout — inverse PCR across the circularized flanks", S.PCR,
        [(circular,R.TEMPLATE),(forward,R.STRUCTURAL),(reverse,R.STRUCTURAL),(polymerase,R.CATALYST)], [product],
        geometry=PCRGeometry.unavailable(duplex("INFERRED opened-circle reference", "known-region primer sites; coordinates deferred", "circularized flanking DNA"),
            reason="The original abstract specifies reverse-oriented primers on a self-ligated restriction fragment. Their actual sequences/sites and the retained known-region gap await the full original; no exact primer-placement claim is made."))
    g.constructed(product, "Scoped amplification product, not a sequencing-ready library. Downstream analysis and complete original method remain outside the abstract-scoped draft.")
    g.declare_reviewed_scope(SOURCE_INVENTORY,SOURCE)
    return g.require_publishable()


SOURCE_INVENTORY = ("dna","nuclease","digest","ligase","circularize","forward","reverse","polymerase","ipcr","amplicon")


def render_page():
    return assemble_page(protocol_graph(), (head(TITLE), '<div class="wrap">',
        '<h1>'+TITLE+'</h1><p class="research-notes"><a href="01_inverse_pcr.html">Research notes</a></p>',
        info('<a href="https://doi.org/10.1093/genetics/120.3.621">Ochman et al., Genetics (1988)</a>'),
        info("DRAFT — abstract-scoped architecture. INFERRED region order is schematic; exact bases, primer sites, end chemistry and full original methods remain unassigned. Generic locus placeholders do not select the original organism or insertion sequence."),
        WorkflowSection("Restriction fragment → circle → flanking amplicon"),
        info("The outward-facing primers become convergent around the circle's flanking arc. The graph ends at constructed DNA; sequencing adapters, indexes and a sequencing run are not inferred."), '</div>'))


def selftest():
    from checks import Check, Source, have
    from storage import sources_root
    from chemistry_validation import Status, check_graph
    check=Check()
    check("no contradictory molecular claims", [c for c in check_graph(protocol_graph()) if c.status is Status.CONTRADICTORY], [])
    source=Source(sources_root()/Path(__file__).resolve().parents[1].name/"genetics_120.3.621_PMC1203539.html.txt", "https://pmc.ncbi.nlm.nih.gov/articles/PMC1203539/")
    if have(check,source,label="defining abstract"):
        text=' '.join(source.path.read_text().split())
        check("original abstract establishes self-ligated restriction-fragment template", "restriction fragment that has been ligated upon itself to form a circle" in text, True)
        check("original abstract establishes reverse primer orientation", "primers oriented in the reverse direction" in text, True)
    check.report()
