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

PAGE_FILE = "rca_gip.html"
TITLE = "RCA-GIP — rolling-circle assisted genome walking"
SOURCE = "Tsaftaris, Pasentzis and Argiriou, Biotechnol Lett32:157–161 (2010), online17Sep2009, doi10.1007/s10529-009-0128-9; defining author abstract in PMID19760115. Full methods remain inaccessible in this bounded pass."


def protocol_graph():
    g=ProtocolGraph(TITLE, authority=SOURCE, required_output=Disposition.CONSTRUCT)
    dna=g.input("dna","genomic DNA surrounding a known region",K.MOLECULAR,
        state=duplex("INFERRED generic locus", "unknown left flank", "known region", "unknown right flank"),primary=True)
    circles=g.product("circles","circular genomic DNA fragments",K.MOLECULAR,
        state=circle("INFERRED circular genomic fragment", "known-region class", "unknown flanking sequence"))
    g.react("make_circles","INFERRED fragment layout — generate circular genomic fragments",S.MAGIC,[(dna,R.TEMPLATE)],[circles],
        reason="The defining abstract establishes circular genomic fragments but does not specify the restriction/ligation recipe. Full original methods, enzymes, end chemistry and fragment selection remain unreviewed; no later application recipe is imported.")
    polymerase=g.input("phi29","phi29 DNA polymerase",K.CATALYST)
    primer=g.input("rca_primer","RCA priming reagent; identity not assigned from abstract",K.OLIGO)
    concatemer_state=Drawing((Row(chunks=[("Linear tandem-copy concatemer class; strand composition, copy count and exact repeat boundaries unassigned",None,False)]),))
    concatemer=g.product("concatemer","linear genomic concatemers",K.MOLECULAR,state=concatemer_state)
    g.react("rca","Rolling-circle amplification of the circular genomic templates",S.EXTENSION,
        [(circles,R.TEMPLATE),(primer,R.STRUCTURAL),(polymerase,R.CATALYST)],[concatemer])
    forward=g.input("forward","known-region inverse primer A; bases unassigned",K.OLIGO)
    reverse=g.input("reverse","known-region inverse primer B; bases unassigned",K.OLIGO)
    pcr_enzyme=g.input("pcr_enzyme","PCR polymerase; identity deferred",K.CATALYST)
    amplicon=g.product("amplicon","genome-walking amplicon class",K.MOLECULAR,
        state=duplex("INFERRED intended inverse-PCR amplicon", "known-region boundary A", "unknown flanking DNA", "known-region boundary B"))
    g.react("ipcr","INFERRED amplicon layout — inverse PCR from the RCA concatemer",S.PCR,
        [(concatemer,R.TEMPLATE),(forward,R.STRUCTURAL),(reverse,R.STRUCTURAL),(pcr_enzyme,R.CATALYST)],[amplicon],
        geometry=PCRGeometry.unavailable(Scene(),reason="The defining abstract identifies the linear concatemer as inverse-PCR template, not an Illumina adapter ligation or an exact primer/site declaration. Selected primers and repeat-traversal boundaries await the full original."))
    g.constructed(amplicon,"Scoped genome-walking amplicon. The abstract names sequencing or cloning as downstream options; neither is silently converted into an NGS-ready endpoint.")
    g.declare_reviewed_scope(SOURCE_INVENTORY,SOURCE)
    return g.require_publishable()


SOURCE_INVENTORY=("dna","make_circles","circles","phi29","rca_primer","rca","concatemer","forward","reverse","pcr_enzyme","ipcr","amplicon")


def render_page():
    return assemble_page(protocol_graph(),(head(TITLE),'<div class="wrap">',
        '<h1>'+TITLE+'</h1><p class="research-notes"><a href="01_rca_gip.html">Research notes</a></p>',
        info('<a href="https://doi.org/10.1007/s10529-009-0128-9">Tsaftaris et al., Biotechnology Letters (online2009; issue2010)</a>'),
        info("DRAFT — defining-abstract scope. INFERRED locus/amplicon order remains schematic. Circle generation, priming, exact repeat sequence and inverse-PCR sites are unresolved; no target-specific promoter design is selected. The RCA class does not assign a copy number or claim successful recovery."),
        WorkflowSection("Circular genomic fragments → RCA concatemers → inverse-PCR amplicon"),'</div>'))


def selftest():
    from checks import Check, Source, have
    from storage import sources_root
    from chemistry_validation import Status, check_graph
    check=Check()
    check("no contradictory molecular claims",[c for c in check_graph(protocol_graph()) if c.status is Status.CONTRADICTORY],[])
    source=Source(sources_root()/Path(__file__).resolve().parents[1].name/"19760115.html.txt","https://pubmed.ncbi.nlm.nih.gov/19760115/; full original doi10.1007/s10529-009-0128-9")
    if have(check,source,label="defining author abstract"):
        text=' '.join(source.path.read_text().split())
        check("original abstract independently identifies phi29", "29 DNA polymerase" in text, True)
        check("original abstract independently identifies linear concatemers", "linear concatemers" in text,True)
    check.report()
