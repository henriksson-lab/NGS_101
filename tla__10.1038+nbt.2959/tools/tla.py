"""Protocol-owned, source-scoped molecular graph; unresolved coordinates stay explicit."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
from chemdraw import Construct, Row, Scene, complement_segments
from circular import CircularizationGeometry
from pcr import PCRGeometry
from page import head, info, table
from protocol_graph import Disposition, InputRole as R, MaterialKind as K, ProtocolGraph, StepKind as S, SequencingPlan, Terminal
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


PCR_COUNT = 8
PCR_ALIQUOT_UL = 25


def _qc(g, key, stock, amount):
    analytical = g.product(key + "_analytical", amount + " source diagnostic aliquot", K.QC, state=stock.state)
    continuing = g.product(key + "_continuing", "Stock after " + amount + " diagnostic allocation", K.MOLECULAR, state=stock.state)
    g.react(key, "Reserve " + amount + " for source diagnostic QC", S.SPLIT, [(stock, R.TEMPLATE)], [analytical, continuing])
    g.mark_terminal(Terminal(analytical.id, Disposition.QC, "Source diagnostic branch; assay reagents and optional failed-control retries are outside the production path."))
    return continuing


def protocol_graph():
    restriction_hierarchy()
    g=ProtocolGraph(TITLE,authority=SOURCE,representation="published_instance")
    cells=g.input("cells","cells; generic locus, no selected study target",K.BIOLOGICAL,
        state=Drawing((Row(chunks=[("Cells containing genomic DNA; no covalent joining between genomic regions is implied",None,False)]),)),primary=True)
    formaldehyde=g.input("formaldehyde","formaldehyde; molecular crosslinking reagent",K.MOLECULAR)
    fixation_conditions=g.input("fixation_conditions","Source fixation: 10minutes at room temperature",K.CONDITION)
    fixed=g.product("fixed","crosslinked genomic material",K.MOLECULAR,
        state=Drawing((Row(chunks=[("Physically proximal genomic DNA held by formaldehyde crosslinks; not a DNA ligation junction",None,False)]),)))
    g.react("crosslink","Crosslink physically proximal genomic regions",S.CROSS_LINKING,[(cells,R.TEMPLATE),(formaldehyde,R.STRUCTURAL),(fixation_conditions,R.CONDITION)],[fixed])
    glycine=g.input("glycine","glycine quench",K.CONDITION)
    quenched=g.product("quenched","quenched fixed cells",K.MOLECULAR,state=fixed.state)
    g.react("quench","Quench crosslinking",S.QUENCHING,[(fixed,R.TEMPLATE),(glycine,R.CONDITION)],[quenched])
    lysis=g.input("lysis_buffer","source detergent lysis buffer; SDS/Triton preparation",K.CONDITION)
    accessible=g.product("accessible","restriction-accessible crosslinked material",K.MOLECULAR,state=fixed.state)
    supernatant=g.product("lysis_supernatant","removed lysis/wash supernatant",K.WASTE)
    g.react("lyse","Lyse and expose chromatin to restriction digestion",S.SOLUBILIZATION,[(quenched,R.TEMPLATE),(lysis,R.CONDITION)],[accessible,supernatant])
    g.discard(supernatant,"Supplementary Protocol fixation/lysis steps9–11 retain the pellet and remove the supernatant. Exact DNA losses are not assigned.")
    accessible=_qc(g,"undigested_qc",accessible,"5µl")
    nla=g.input("nla","NlaIII; primary CATG recognition motif",K.CATALYST)
    fragments=g.product("fragments","crosslinked primary restriction fragments",K.MOLECULAR,
        state=Drawing(tuple(Row(chunks=[(label,None,False)]) for label in (
            "Anchor-containing restriction-fragment class; exact bases, cuts and primer coordinates unassigned",
            "Physically associated captured-fragment classes; separate molecules before proximity ligation"))))
    g.react("digest_primary","Digest crosslinked DNA with NlaIII",S.DIGESTION,[(accessible,R.TEMPLATE),(nla,R.CATALYST)],[fragments])
    fragments=_qc(g,"digested_qc",fragments,"5µl")
    ligase=g.input("ligase","Actual first Roche ligase addition: 20U, SI step24",K.CATALYST)
    contact_state=duplex("INFERRED source-described chimeric class", "anchor-containing primary fragment", "captured primary fragments; count and order unassigned")
    contact=g.product("contact","proximity-ligated chimeric DNA",K.MOLECULAR,state=contact_state)
    g.react("proximity","INFERRED molecule order — ligate proximal primary fragments",S.PROXIMITY_LIGATION,
        [(fragments,R.TEMPLATE),(ligase,R.CATALYST)],[contact])
    contact=_qc(g,"ligated_qc",contact,"10µl")
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
    trimmed=_qc(g,"secondary_digest_qc",trimmed,"5µl")
    secondary_ligase=g.input("secondary_ligase","Actual second Roche ligase addition: 100U, SI step44; not recovered first ligase",K.CATALYST)
    diluted=g.input("dilution","dilute ligation conditions promoting intramolecular closure",K.CONDITION)
    circular=g.product("circular","circular chimeric template classes",K.MOLECULAR,
        state=circle("INFERRED chimeric circle order", "anchor fragment context", "captured fragments"))
    g.react("circularize","INFERRED circle layout — re-ligate secondary fragments",S.CIRCULARIZATION,
        [(trimmed,R.TEMPLATE),(secondary_ligase,R.CATALYST),(diluted,R.CONDITION)],[circular],
        geometry=CircularizationGeometry.unavailable("The original identifies NspI fragment self-ligation. A generic locus does not specify selected cleavage coordinates or exact junction bases; exact end geometry is not assigned merely from the intended circle."))
    purified=g.product("purified","purified circular TLA templates",K.MOLECULAR,state=circular.state)
    cleanup_waste=g.product("cleanup_waste","removed precipitation/column fractions",K.WASTE)
    g.react("cleanup_circles","Precipitate and column-purify circular templates",S.CLEANUP,[(circular,R.TEMPLATE)],[purified,cleanup_waste])
    g.discard(cleanup_waste,"Supplementary Protocol steps46–53 retain the precipitated/column-purified DNA. Exact off-product depletion is unassigned.")
    purified=_qc(g,"qubit_qc",purified,"unassigned Qubit volume")
    selected=g.product("pcr_template", "Actual 800ng circular-template contribution to one source PCR setup", K.MOLECULAR,state=purified.state)
    remaining=g.product("template_stock", "Untransferred purified circle stock; amount unassigned",K.MOLECULAR,state=purified.state)
    g.react("select_pcr_template","Take 800ng for the source PCR setup",S.SPLIT,[(purified,R.TEMPLATE)],[selected,remaining])
    g.mark_terminal(Terminal(remaining.id,Disposition.UNRESOLVED,"Actual template stock outside the selected PCR contribution."))
    forward=g.input("forward","Forward component of source primer set; Table2 member unselected",K.OLIGO)
    reverse=g.input("reverse","Reverse component of the same source primer set; combined dose only",K.OLIGO)
    polymerase=g.input("polymerase","Actual 4µl source Phire Polymerase contribution",K.CATALYST)
    nucleotides=g.input("pcr_nucleotides","Actual 4µl source 10mM dNTP stock contribution",K.MOLECULAR)
    pcr_buffer=g.input("pcr_buffer","Actual 40µl source 5x PCR buffer contribution",K.CONDITION)
    mix=g.product("pcr_mix","One actual 200µl source setup: 800ng template, combined 10µl primer set and source Phire/dNTP/buffer contributions",K.MOLECULAR,state=purified.state)
    g.react("make_pcr_mix","Assemble one source PCR setup before dividing it",S.GENERIC,
        [(selected,R.TEMPLATE),(forward,R.STRUCTURAL),(reverse,R.STRUCTURAL),(polymerase,R.CATALYST),(nucleotides,R.STRUCTURAL),(pcr_buffer,R.CONDITION)],[mix],changes_molecule=False)
    aliquots=[g.product("pcr"+str(i)+"_aliquot","Actual 25µl reaction aliquot"+str(i)+" from the same 200µl setup",K.MOLECULAR,state=purified.state) for i in range(1,PCR_COUNT+1)]
    g.react("divide_pcr_mix","Divide the source 200µl setup into8x25µl reactions",S.SPLIT,[(mix,R.TEMPLATE)],aliquots)
    pcr_state=duplex("INFERRED linear chimeric amplicon class", "anchor primer boundary A", "captured fragments across the circle", "anchor primer boundary B")
    products=[]
    for i,aliquot in enumerate(aliquots,1):
        key="pcr"+str(i)
        # These ports expose copresent components of that same physical aliquot.
        # They are not fresh stocks, separate physical fractions or new doses.
        template=g.product(key+"_template","Copresent template in aliquot"+str(i),K.MOLECULAR,state=purified.state)
        f=g.product(key+"_f","Copresent unselected forward primer in that same aliquot",K.OLIGO)
        r=g.product(key+"_r","Copresent unselected reverse primer in that same aliquot",K.OLIGO)
        enzyme=g.product(key+"_enzyme","Copresent Phire component of that same aliquot",K.CATALYST)
        nt=g.product(key+"_nt","Copresent dNTP component of that same aliquot",K.MOLECULAR)
        buffer=g.product(key+"_buffer","Copresent buffer component of that same aliquot",K.CONDITION)
        g.react(key+"_components","Inventory view of the same 25µl aliquot; no physical fractionation or fresh additions",S.GENERIC,[(aliquot,R.TEMPLATE)],[template,f,r,enzyme,nt,buffer],changes_molecule=False)
        product=g.product(key+"_product","INFERRED TLA PCR product class from actual aliquot"+str(i),K.MOLECULAR,state=pcr_state)
        g.react(key+"_reaction","INFERRED amplicon boundaries — actual inverse-PCR aliquot"+str(i),S.PCR,
            [(template,R.TEMPLATE),(f,R.STRUCTURAL),(r,R.STRUCTURAL),(enzyme,R.CATALYST),(nt,R.STRUCTURAL),(buffer,R.CONDITION)],[product],
            geometry=PCRGeometry.unavailable(duplex("INFERRED opened circle reference","outward anchor primer sites","captured circular DNA"),reason="Generic source preparation: no study-specific Table2 pair or target map selected. Exact primer sites and amplicon boundaries are unassigned; source primers have no dangling NGS adapters."))
        products.append(product)
    pooled=g.product("amplified","Actual pooled population of the eight source TLA PCR aliquots",K.MOLECULAR,state=pcr_state)
    g.react("pool_pcr","Pool the eight actual source PCR products",S.POOL,[(x,R.TEMPLATE)for x in products],[pooled])
    prep=g.input("ngs_prep","manufacturer Illumina library-preparation components; selected kit unassigned",K.MOLECULAR)
    library=g.product("library","Illumina libraries from randomly fragmented TLA amplicons",K.LIBRARY,
        state=duplex("INFERRED library class; adapter identity/order unassigned", "unassigned sequencing-library end", "random fragment of chimeric amplicon", "unassigned sequencing-library end opposite"))
    g.react("prepare_library","INFERRED end arrangement — random fragmentation and Illumina library preparation",S.MAGIC,
        [(pooled,R.TEMPLATE),(prep,R.STRUCTURAL)],[library],
        reason="Supplementary Fig1 specifies random fragmentation followed by manufacturer library preparation; main Methods identifies Illumina platforms but no selected kit or exact adapter/index geometry in the reviewed generic recipe. No Nextera, TruSeq, tagmentation enzyme or ligation recipe is silently substituted.")
    loading=g.product("loading","Selected sequencing-loading portion; amount unassigned",K.LIBRARY,state=library.state)
    library_stock=g.product("library_stock","Untransferred prepared library stock; amount unassigned",K.LIBRARY,state=library.state)
    g.react("select_loading","Select the source sequencing contribution",S.SPLIT,[(library,R.TEMPLATE)],[loading,library_stock])
    g.mark_terminal(Terminal(library_stock.id,Disposition.UNRESOLVED,"No whole-stock consumption or additional pool is invented."))
    g.sequence(loading,SequencingPlan(acquired_roles=("Read 1","Read 2"),limitation="Published Methods identifies paired-end Illumina sequencing on MiSeq, HiSeq2000 and HiSeq2500. Exact kit, sequencing-primer geometry, indexes and cycle counts are unassigned in this generic graph. Paired ends were mapped independently in the general analysis owing to sequence reshuffling; this is computational handling, not physical read removal."))
    g.declare_reviewed_scope(SOURCE_INVENTORY,SOURCE)
    return g.require_publishable()


# Independently curated source obligations and the source-defined eightfold array.
QC_ITEMS=tuple(key+suffix for key in ("undigested_qc","digested_qc","ligated_qc","secondary_digest_qc","qubit_qc") for suffix in ("","_analytical","_continuing"))
PCR_ITEMS=tuple(key+suffix for key in ("pcr1","pcr2","pcr3","pcr4","pcr5","pcr6","pcr7","pcr8") for suffix in ("_aliquot","_components","_template","_f","_r","_enzyme","_nt","_buffer","_reaction","_product"))
SOURCE_INVENTORY=("cells","formaldehyde","fixation_conditions","crosslink","glycine","quench","lysis_buffer","lyse","nla","digest_primary","ligase","proximity","protease","de_crosslink","rnase","purify","nsp","digest_secondary","dilution","circularize","cleanup_circles","forward","reverse","polymerase","pcr_nucleotides","pcr_buffer","pcr_mix","make_pcr_mix","divide_pcr_mix","pool_pcr","pcr_template","template_stock","select_pcr_template","ngs_prep","prepare_library","library","secondary_ligase","loading","library_stock","select_loading") + QC_ITEMS + PCR_ITEMS


def render_page():
    return assemble_page(protocol_graph(),(head(TITLE),'<div class="wrap" data-protocol-complete="true">',
        '<h1>'+TITLE+'</h1><p class="research-notes"><a href="01_tla.html">Research notes</a></p>',
        info('<a href="https://doi.org/10.1038/nbt.2959">de Vree et al., Nature Biotechnology (2014)</a>'),
        info("INFERRED — generic molecular ordering and primer/library boundaries are schematic. No study-specific anchor pair, target map or clone is selected; exact cleavage/ligation and kit geometry remain unassigned. Crosslinked proximity does not establish a covalent DNA junction before ligation."),
        info("NlaIII recognizes CATG. The NspI recognition class RCATGY retains that core and therefore re-cuts only a subset of primary sites. TLA keeps captured fragment interiors, then fragments its larger PCR products for sequencing; it does not use the 4C adapter-tailed inverse-PCR endpoint."),
        "<h2>From spatial proximity to an anchor-spanning library</h2>",
        "<p>Crosslinking preserves proximity before DNA is digested and ligated. NspI trimming and a second ligation produce circular chimeric template classes. Outward primers on an anchor-containing fragment amplify across the circle; the resulting products undergo a separate, unassigned Illumina library preparation.</p>",
        table(("Source operation","Scope"),(("First and second ligation","Distinct Roche ligase additions: 20U then 100U, separated by de-crosslinking and purification."),("PCR array","One 200µl setup with 800ng template is divided into eight 25µl reactions, then pooled. Component views expose the same aliquots, not fresh reagent additions."),("Diagnostics","Undigested/digested, first-ligation and second-digest aliquots plus Qubit measurement remain analytical. Optional failed-control retries are outside this selected continuation."),("Readout","Paired-end acquisition on source-named Illumina platforms; exact cycle counts, indexes and priming geometry remain unassigned."))),
        '<details><summary>Full source preparation and unresolved chemistry</summary>',
        WorkflowSection("Proximity ligation and inverse-PCR library construction"),'</details>',
        PlanLimitationSection("loading"),'</div>'))


def selftest():
    import re
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
        dimensions=re.search(r"divideover(\d+)x(\d+)ulrxns",re.sub(r"\s+","",text))
        check("PCR array dimensions match the independently printed preparation",(PCR_COUNT,PCR_ALIQUOT_UL),tuple(map(int,dimensions.groups())) if dimensions else None)
    check.report()
