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

PAGE_FILE = "tla.html"
TITLE = "TLA — targeted locus amplification"
SOURCE = "de Vree et al., Nat Biotechnol32:1019–1025 (2014), doi10.1038/nbt.2959. Published main Online Methods: TLA sample preparation and Bioinformatics; Supplementary Fig1 and generic tissue-culture Supplementary Protocol through TLA PCR. Generic locus only, not a target-specific demonstration."
MAIN = "DeLaat_Key_2014_DeVree_Targeted-sequencing-by-proximity-ligation-for-comprehensive-variant-deletion-and-local-haplotyping.pdf.txt"
SUPPLEMENT = "nbt.2959_41587_2014_BFnbt2959_MOESM1_ESM.pdf.txt"


def restriction_hierarchy():
    primary, secondary = "CATG", "RCATGY"
    if secondary[1:-1] != primary:
        raise ValueError("declared secondary recognition core must retain the primary motif")
    return primary, secondary


def protocol_graph():
    restriction_hierarchy()
    g=ProtocolGraph(TITLE,authority=SOURCE)
    cells=g.input("cells","cells; generic locus, no selected study target",K.BIOLOGICAL,
        state=Drawing((Row(chunks=[("Cells containing genomic DNA; no covalent joining between genomic regions is implied",None,False)]),)),primary=True)
    formaldehyde=g.input("formaldehyde","formaldehyde",K.CONDITION)
    fixed=g.product("fixed","crosslinked genomic material",K.MOLECULAR,
        state=Drawing((Row(chunks=[("Physically proximal genomic DNA held by formaldehyde crosslinks; not a DNA ligation junction",None,False)]),)))
    g.react("crosslink","Crosslink physically proximal genomic regions",S.CROSS_LINKING,[(cells,R.TEMPLATE),(formaldehyde,R.CONDITION)],[fixed])
    glycine=g.input("glycine","glycine quench",K.CONDITION)
    quenched=g.product("quenched","quenched fixed cells",K.MOLECULAR,state=fixed.state)
    g.react("quench","Quench crosslinking",S.QUENCHING,[(fixed,R.TEMPLATE),(glycine,R.CONDITION)],[quenched])
    lysis=g.input("lysis_buffer","source detergent lysis buffer; SDS/Triton preparation",K.CONDITION)
    accessible=g.product("accessible","restriction-accessible crosslinked material",K.MOLECULAR,state=fixed.state)
    supernatant=g.product("lysis_supernatant","removed lysis/wash supernatant",K.WASTE)
    g.react("lyse","Lyse and expose chromatin to restriction digestion",S.SOLUBILIZATION,[(quenched,R.TEMPLATE),(lysis,R.CONDITION)],[accessible,supernatant])
    g.discard(supernatant,"Supplementary Protocol fixation/lysis steps9–11 retain the pellet and remove the supernatant. Exact DNA losses are not assigned.")
    nla=g.input("nla","NlaIII; primary CATG recognition motif",K.CATALYST)
    fragments=g.product("fragments","crosslinked primary restriction fragments",K.MOLECULAR,
        state=Drawing(tuple(Row(chunks=[(label,None,False)]) for label in (
            "Anchor-containing restriction-fragment class; exact bases, cuts and primer coordinates unassigned",
            "Physically associated captured-fragment classes; separate molecules before proximity ligation"))))
    g.react("digest_primary","Digest crosslinked DNA with NlaIII",S.DIGESTION,[(accessible,R.TEMPLATE),(nla,R.CATALYST)],[fragments])
    ligase=g.input("ligase","Roche ligase; primary and secondary ligation",K.CATALYST)
    contact_state=duplex("INFERRED example chimeric DNA order", "anchor-containing primary fragment", "captured primary fragments; count and order unassigned")
    contact=g.product("contact","proximity-ligated chimeric DNA",K.MOLECULAR,state=contact_state)
    g.react("proximity","INFERRED molecule order — ligate proximal primary fragments",S.PROXIMITY_LIGATION,
        [(fragments,R.TEMPLATE),(ligase,R.CATALYST)],[contact])
    protease=g.input("protease","Proteinase K and heat; reverse crosslinks",K.CATALYST)
    released=g.product("released","de-crosslinked chimeric DNA",K.MOLECULAR,state=contact_state)
    g.react("de_crosslink","Reverse crosslinks and remove protein",S.DEPROTEINIZATION,[(contact,R.TEMPLATE),(protease,R.CATALYST)],[released])
    rnase=g.input("rnase","RNase A",K.CATALYST)
    clean=g.product("clean","purified chimeric DNA",K.MOLECULAR,state=contact_state)
    extraction_waste=g.product("extraction_waste","removed extraction/precipitation fractions",K.WASTE)
    g.react("purify","RNase treatment, extraction and precipitation",S.CLEANUP,[(released,R.TEMPLATE),(rnase,R.CATALYST)],[clean,extraction_waste])
    g.discard(extraction_waste,"Supplementary Protocol steps29–37 transfer the aqueous fraction and retain the DNA pellet. No molecule-specific recovery or depletion is asserted.")
    nsp=g.input("nsp","NspI; secondary RCATGY recognition motif",K.CATALYST)
    trimmed=g.product("trimmed","secondary restriction-fragment classes",K.MOLECULAR,
        state=duplex("INFERRED secondary fragment order", "anchor fragment context", "captured NlaIII fragments; individual cut coordinates unassigned"))
    g.react("digest_secondary","INFERRED fragment arrangement — trim at NspI sites",S.DIGESTION,[(clean,R.TEMPLATE),(nsp,R.CATALYST)],[trimmed])
    diluted=g.input("dilution","dilute ligation conditions promoting intramolecular closure",K.CONDITION)
    circular=g.product("circular","circular chimeric template classes",K.MOLECULAR,
        state=circle("INFERRED chimeric circle order", "anchor fragment context", "captured fragments"))
    g.react("circularize","INFERRED circle layout — re-ligate secondary fragments",S.CIRCULARIZATION,
        [(trimmed,R.TEMPLATE),(ligase,R.CATALYST),(diluted,R.CONDITION)],[circular],
        geometry=CircularizationGeometry.unavailable("The original identifies NspI fragment self-ligation. A generic locus does not specify selected cleavage coordinates or exact junction bases; exact end geometry is not assigned merely from the intended circle."))
    purified=g.product("purified","purified circular TLA templates",K.MOLECULAR,state=circular.state)
    cleanup_waste=g.product("cleanup_waste","removed precipitation/column fractions",K.WASTE)
    g.react("cleanup_circles","Precipitate and column-purify circular templates",S.CLEANUP,[(circular,R.TEMPLATE)],[purified,cleanup_waste])
    g.discard(cleanup_waste,"Supplementary Protocol steps46–53 retain the precipitated/column-purified DNA. Exact off-product depletion is unassigned.")
    forward=g.input("forward","target-specific outward anchor primer A; Table2 member unselected",K.OLIGO)
    reverse=g.input("reverse","target-specific outward anchor primer B; Table2 member unselected",K.OLIGO)
    polymerase=g.input("polymerase","Phire DNA polymerase",K.CATALYST)
    pcr_state=duplex("INFERRED linear chimeric amplicon order", "anchor primer boundary A", "captured fragments across the circle", "anchor primer boundary B")
    amplified=g.product("amplified","source pooled population of TLA inverse-PCR products",K.MOLECULAR,state=pcr_state)
    g.react("ipcr","INFERRED amplicon boundaries — inverse PCR on anchor-containing circles",S.PCR,
        [(purified,R.TEMPLATE),(forward,R.STRUCTURAL),(reverse,R.STRUCTURAL),(polymerase,R.CATALYST)],[amplified],
        geometry=PCRGeometry.unavailable(duplex("INFERRED opened circle reference", "outward anchor primer sites", "captured circular DNA"),
            reason="Supplementary Table2 contains study-specific primer pairs; none is selected for this generic locus view. Their exact sites and junction-spanning amplicons need the matching target map. Source inverse primers lack dangling NGS adapters (Supplementary Fig1)."))
    # Individual physical aliquots are not separate molecular classes here.
    # The source pools parallel PCRs; represent their combined product population
    # without inventing a second independent substrate to satisfy a POOL signature.
    pooled=amplified
    prep=g.input("ngs_prep","manufacturer Illumina library-preparation components; selected kit unassigned",K.MOLECULAR)
    library=g.product("library","Illumina libraries from randomly fragmented TLA amplicons",K.LIBRARY,
        state=duplex("INFERRED library class; adapter identity/order unassigned", "unassigned sequencing-library end", "random fragment of chimeric amplicon", "unassigned sequencing-library end opposite"))
    g.react("prepare_library","INFERRED end arrangement — random fragmentation and Illumina library preparation",S.MAGIC,
        [(pooled,R.TEMPLATE),(prep,R.STRUCTURAL)],[library],
        reason="Supplementary Fig1 specifies random fragmentation followed by manufacturer library preparation; main Methods identifies Illumina platforms but no selected kit or exact adapter/index geometry in the reviewed generic recipe. No Nextera, TruSeq, tagmentation enzyme or ligation recipe is silently substituted.")
    g.sequence(library,SequencingPlan(acquired_roles=("Read1","Read2"),limitation="Published Methods identifies paired-end Illumina sequencing on MiSeq, HiSeq2000 and HiSeq2500. Exact kit, sequencing-primer geometry, indexes and cycle counts are unassigned in this generic graph. Paired ends were mapped independently in the general analysis owing to sequence reshuffling; this is computational handling, not physical read removal."))
    g.declare_reviewed_scope(SOURCE_INVENTORY,SOURCE)
    return g.require_publishable()


SOURCE_INVENTORY=("cells","formaldehyde","crosslink","glycine","quench","lysis_buffer","lyse","nla","digest_primary","ligase","proximity","protease","de_crosslink","rnase","purify","nsp","digest_secondary","dilution","circularize","cleanup_circles","forward","reverse","polymerase","ipcr","ngs_prep","prepare_library","library")


def render_page():
    return assemble_page(protocol_graph(),(head(TITLE),'<div class="wrap">',
        '<h1>'+TITLE+'</h1><p class="research-notes"><a href="01_tla.html">Research notes</a></p>',
        info('<a href="https://doi.org/10.1038/nbt.2959">de Vree et al., Nature Biotechnology (2014)</a>'),
        info("INFERRED — generic molecular ordering and primer/library boundaries are schematic. No study-specific anchor pair, target map or clone is selected; exact cleavage/ligation and kit geometry remain unassigned. Crosslinked proximity does not establish a covalent DNA junction before ligation."),
        info("NlaIII recognizes CATG. The NspI recognition class RCATGY retains that core and therefore re-cuts only a subset of primary sites. TLA keeps captured fragment interiors, then fragments its larger PCR products for sequencing; it does not use the 4C adapter-tailed inverse-PCR endpoint."),
        WorkflowSection("Proximity ligation and inverse-PCR library construction"),
        PlanLimitationSection("library"),'</div>'))


def selftest():
    from checks import Check, Source, have
    from storage import sources_root
    from chemistry_validation import Status, check_graph
    check=Check()
    check("no contradictory molecular claims",[c for c in check_graph(protocol_graph()) if c.status is Status.CONTRADICTORY],[])
    directory=sources_root()/Path(__file__).resolve().parents[1].name
    main=Source(directory/MAIN,"https://www.hubrecht.eu/app/uploads/2017/11/DeLaat_Key_2014_DeVree_Targeted-sequencing-by-proximity-ligation-for-comprehensive-variant-deletion-and-local-haplotyping.pdf")
    supp=Source(directory/SUPPLEMENT,"Published Supplementary Text and Figures linked from https://doi.org/10.1038/nbt.2959")
    if have(check,main,label="published main Online Methods"):
        text=' '.join(main.path.read_text().split())
        primary,secondary=restriction_hierarchy()
        check("recognition motifs transcribed independently from published original", all(motif in text for motif in (primary,secondary)),True)
        check("original Methods identify paired-end acquisition", "we perform paired-end sequencing" in text,True)
        check("original distinguishes independent-end mapping from acquisition", "Paired ends are therefore treated separately" in text,True)
    if have(check,supp,label="published generic Supplementary Protocol"):
        text=' '.join(supp.path.read_text().split())
        check("generic source protocol independently names the PCR polymerase", "Phire Polymerase" in text,True)
        check("source distinguishes TLA primers from adapter-tailed 4C primers", "TLA primers do not contain dangling sequences" in text,True)
    check.report()
