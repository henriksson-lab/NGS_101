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
    g = ProtocolGraph(TITLE, authority=SOURCE, required_output=Disposition.CONSTRUCT, representation='illustrative_example', example_description='Generic locus/circle/amplicon role layouts illustrate the abstract-supported method; no literal DNA, source primer pair, organism-specific template or experimental recovery is reconstructed.')
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
    pcr_nt = g.input("pcr_nucleotides", "Functional PCR nucleotide contributors; formulation/dose/stock packaging unassigned", K.MOLECULAR)
    product = g.product("amplicon", "linear flanking-sequence amplicon", K.MOLECULAR,
        state=duplex("INFERRED circle-spanning amplicon order", "known-region boundary A", "flanking DNA across restriction-fragment closure", "known-region boundary B"))
    g.react("ipcr", "INFERRED amplicon layout — inverse PCR across the circularized flanks", S.PCR,
        [(circular,R.TEMPLATE),(forward,R.STRUCTURAL),(reverse,R.STRUCTURAL),(polymerase,R.CATALYST),(pcr_nt,R.STRUCTURAL)], [product],
        geometry=PCRGeometry.unavailable(duplex("INFERRED opened-circle reference", "known-region primer sites; coordinates deferred", "circularized flanking DNA"),
            reason="The original abstract specifies reverse-oriented primers on a self-ligated restriction fragment. Their actual sequences/sites and the retained known-region gap await the full original; no exact primer-placement claim is made."))
    g.constructed(product, "Scoped amplification product, not a sequencing-ready library. Downstream analysis and complete original method remain outside the abstract-scoped draft.")
    g.declare_reviewed_scope(SOURCE_INVENTORY,SOURCE)
    return g.require_publishable()


SOURCE_INVENTORY = ("dna","nuclease","digest","ligase","circularize","forward","reverse","polymerase","pcr_nucleotides","ipcr","amplicon")


def render_page():
    g=protocol_graph().require_publishable()
    return assemble_page(g,(head(TITLE),'<div class="wrap" data-protocol-complete="true">',
        '<h1>Inverse PCR</h1><p>Recover DNA flanking a known region by circularizing a restriction fragment and amplifying around its flanking arc with outward-facing primers.</p>',
        info('<a href="https://doi.org/10.1093/genetics/120.3.621">Ochman, Gerber and Hartl, Genetics 120:621–623 (1988)</a>. The current source scope is the defining abstract. The paper reports recovery around an IS1 element in an Escherichia coli isolate; its full methods and original Figure 1 remain unreviewed.'),
        info('ILLUSTRATIVE EXAMPLE — generic known-region and flanking-region layouts. INFERRED drawings do not reproduce the original IS1 template, primer sites or restriction junction. No literal DNA bases or measurable segment lengths have been invented.'),
        '<h2>Why outward primers work</h2><p>A genomic restriction fragment is ligated into a circle. Primers that face away from each other on the original linear locus become convergent around the circularized flanking arc. PCR therefore recovers adjacent unknown DNA. Enzyme identities, actual primers, fragment sizes, end chemistry and reaction conditions remain unresolved in this scoped illustration.</p>',
        '<h2>Constructed output</h2><p>The endpoint is a linear flanking-sequence amplicon class, not a sequencing-ready library. No downstream cloning vector, sequencing adapter, barcode or run layout is assigned. The graph checks its stated molecular declarations; it does not establish recovery efficiency or a complete experimental recipe.</p>',
        '<details open><summary>Inferred preparation graph and unresolved chemistry</summary>',WorkflowSection('Restriction fragment → circle → flanking amplicon'),'</details>',
        '<p><a class="research-notes" href="01_inverse_pcr.html">Source notes and remaining questions</a></p></div>'))


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
    g=protocol_graph()
    check("construct gate passes",g.publication_problems(),[])
    check("illustrative representation is explicit",g.representation,"illustrative_example")
    check("no invented sequencing endpoint",bool(g.sequencing),False)
    check("PCR consumes nucleotide contributors",any(p.material=="pcr_nucleotides" and p.role==R.STRUCTURAL for p in g.steps["ipcr"].inputs),True)
    check.report()
