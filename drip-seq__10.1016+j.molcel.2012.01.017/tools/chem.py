"""Source-owned DRIP-seq declaration for Ginno et al.2012.

The later2019 DRIP/DRIPc protocol supplies no reaction choices to this graph.
"""
import seqprimers as sp
from chemdraw import Scene
from adapter_ligation import AdapterLigationGeometry
from protocol_graph import (Highlight,InputRole,MaterialKind,ProtocolGraph,
                            SequencingPlan,StepKind)
from protocol_model import role

TITLE="DRIP-seq — DNA:RNA immunoprecipitation sequencing"
NOTES="01_drip_seq.html"
SOURCE=('Defining source: <a href="https://doi.org/10.1016/j.molcel.2012.01.017">'
        'Ginno et al., Molecular Cell (2012)</a>; Experimental Procedures, '
        'Genome-wide R-loop profiling using DRIP-seq and DRIVE-seq.')
SUMMARY=("Gentle nucleic-acid extraction and restriction digestion preserve RNA:DNA hybrids "
         "for S9.6 immunoprecipitation. The captured material is sonicated, size-selected "
         "and ligated to barcoded Illumina adaptors for DNA-fragment sequencing.")
CAVEAT=("The RNA:DNA hybrid is retained through immunoprecipitation. Its subsequent RNA fate "
        "and complete DNA-library end processing, adaptor/index layout and read mask remain "
        "unassigned in this2012 model. Published supplements are available, with only a bounded "
        "scope inspected. The later2019 protocol is separate and supplies no2012 reaction steps.")
FINAL_LIBRARY=None
SEQ_PRIMERS=()
FINAL_CAPTION="2012 barcoded-adaptor DNA library; molecular ends and lengths remain unassigned."
SEQUENCING_INTRO=("The2012 Methods report Illumina GAIIx and HiSeq sequencing. Exact read roles, "
                  "lengths, index acquisition and primer assignments remain unassigned.")
READ_LENGTHS=None
SEQUENCING_ENDING=SEQUENCING_INTRO
READ_SPECS=(sp.read_spec("Sequencing output",None,"hybrid-associated DNA; read and index roles unassigned","include"),)
SEQUENCING_UNAVAILABLE=""  # No fictitious final Construct is supplied to the primer diagram.


def rloop_scene():
    """Three source-described components; no sequence or molecular length is assigned.

    The paired role tokens identify the hybrid relationship rather than claiming
    base-resolved pairing. The displaced DNA is a separate strand in the assembly.
    """
    sc=Scene()
    sc.strand("template DNA",[role("DNA in RNA:DNA hybrid")],label="hybrid template DNA")
    sc.anneal("hybrid RNA",[role("RNA in RNA:DNA hybrid")],to="template DNA",
        pair=("RNA in RNA:DNA hybrid","DNA in RNA:DNA hybrid"),
        unpaired=("RNA in RNA:DNA hybrid",),label="hybrid RNA (bases unassigned)")
    sc.strand("displaced DNA",[role("displaced single-stranded DNA")],label="displaced DNA")
    sc.note("hybrid RNA","S9.6 recognizes the RNA:DNA hybrid; capture does not remove RNA")
    return sc


def sonicated_material():
    sc=Scene()
    sc.strand("DNA population",[role("sonicated hybrid-associated DNA fragments")],
              label="DNA fragments; lengths and ends unassigned")
    sc.strand("RNA component",[role("co-captured RNA; post-sonication hybrid integrity unresolved")],
              label="RNA component retained in model scope; subsequent fate unresolved")
    return sc


def unknown_dna_substrate():
    sc=Scene()
    sc.strand("library DNA",[role("DNA library substrate; physical ends unassigned")],
              label="DNA component for adaptor ligation; RNA fate is an explicit model boundary")
    return sc


def unknown_library():
    sc=Scene()
    sc.strand("sequencing library",[role("hybrid-associated DNA library; barcoded-adaptor layout unassigned")],
              label="2012 sequencing library; no adapter, index or insert lengths assigned")
    return sc


def oligos():
    return ()


def protocol_graph():
    g=ProtocolGraph(TITLE,authority=SOURCE)
    sample=g.input("rloop_chromatin","Ntera2 cells containing genomic R-loop nucleic acids",
        MaterialKind.BIOLOGICAL,state=rloop_scene(),primary=True)
    extraction=g.input("gentle_extraction","SDS/proteinase K37°C; phenol–chloroform; ethanol precipitation",
        MaterialKind.CONDITION)
    extracted=g.product("extracted_nucleic_acid","gently extracted R-loop-containing nucleic acid",
        MaterialKind.MOLECULAR,state=rloop_scene())
    g.react("extract","Gently extract nucleic acid",StepKind.EXTRACTION,
        ((sample,InputRole.TEMPLATE),(extraction,InputRole.CONDITION)),(extracted,))
    enzymes=g.input("restriction_cocktail","HindIII, EcoRI, BsrGI, XbaI and SspI",MaterialKind.CATALYST)
    digest=g.product("restriction_fragments","restriction-fragmented RNA:DNA-hybrid material",
        MaterialKind.MOLECULAR,state=rloop_scene())
    g.react("restriction_digest","Restriction-fragment genomic DNA while preserving hybrid capture substrate",
        StepKind.DIGESTION,((extracted,InputRole.TEMPLATE),(enzymes,InputRole.CATALYST)),(digest,))
    antibody=g.input("s96_antibody","S9.6 antibody",MaterialKind.SELECTION)
    capture_support=g.input("capture_beads","MeDIP-style capture support (genomic IP implementation unassigned)",
        MaterialKind.SELECTION,inferred=True,
        source="2012 Methods cite Weber2005 MeDIP without denaturation; that predecessor remains uninspected")
    captured=g.product("captured_rloops","S9.6-captured RNA:DNA hybrids with displaced DNA",
        MaterialKind.MOLECULAR,state=rloop_scene())
    rejected=g.product("uncaptured_nucleic_acid","nucleic acid not retained by S9.6 immunoprecipitation",
        MaterialKind.WASTE)
    g.react("immunoprecipitate","Immunoprecipitate RNA:DNA hybrids without denaturation",
        StepKind.AFFINITY_CAPTURE,((digest,InputRole.TEMPLATE),(antibody,InputRole.SELECTION),
        (capture_support,InputRole.SELECTION)),(captured,rejected))
    g.discard(rejected,"non-retained material is excluded by the immunoprecipitation separation; detailed genomic-IP washes remain unassigned")
    sonication=g.input("sonication","post-pull-down sonication condition (settings unassigned)",MaterialKind.CONDITION)
    fragmented=g.product("captured_dna_fragments","sonicated captured material; RNA fate unresolved",
        MaterialKind.MOLECULAR,state=sonicated_material())
    g.react("fragment_captured_dna","Sonicate the pulled-down material",StepKind.FRAGMENTATION,
        ((captured,InputRole.TEMPLATE),(sonication,InputRole.CONDITION)),(fragmented,))
    size_rule=g.input("size_selection","physical size selection; numeric window unassigned",MaterialKind.SELECTION)
    selected=g.product("size_selected_dna","size-selected hybrid-associated material",
        MaterialKind.MOLECULAR,state=fragmented.state)
    off_size=g.product("out_of_range_fragments","material outside the physically selected size fraction",MaterialKind.WASTE)
    g.react("select_size","Size-select the sonicated material",StepKind.SELECTION,
        ((fragmented,InputRole.TEMPLATE),(size_rule,InputRole.SELECTION)),(selected,off_size),changes_molecule=False)
    g.discard(off_size,"material outside the recovered size fraction is not retained; numeric cutoffs are unassigned")
    substrate=g.product("ligation_substrate","DNA component entering library-adaptor ligation",
        MaterialKind.MOLECULAR,state=unknown_dna_substrate())
    g.react("library_substrate_boundary","Unresolved RNA fate and DNA-library substrate processing",
        StepKind.MAGIC,((selected,InputRole.TEMPLATE),),(substrate,),
        reason="The model reaches the reported DNA sequencing library from captured hybrid material, but has not reconstructed RNA recovery/removal, end processing or any amplification. No RNA-destruction enzyme, end repair, dA-tailing or PCR is imported from the separate2019 method.")
    adapters=g.input("library_adapters","2012 Illumina barcoded adaptors (physical sequences/layout unassigned)",
        MaterialKind.ADAPTER)
    ligase=g.input("dna_ligase","adaptor-ligation activity (enzyme formulation unassigned)",MaterialKind.CATALYST)
    library=g.product("sequencing_library","2012 DRIP-seq DNA-fragment sequencing library",
        MaterialKind.LIBRARY,state=unknown_library())
    g.react("adapter_ligation","Ligate the source-named barcoded Illumina adaptors",StepKind.LIGATION,
        ((substrate,InputRole.TEMPLATE),(adapters,InputRole.STRUCTURAL),(ligase,InputRole.CATALYST)),
        (library,),geometry=AdapterLigationGeometry.unavailable(
            "The current2012 declaration has not assigned the physical adapter strands, substrate ends or ligation junction. Later2019 TruSeq/dA chemistry does not establish those values for this source."))
    g.sequence(library,SequencingPlan(limitation=SEQUENCING_INTRO))
    g.annotate(Highlight("immunoprecipitate","Control samples are separate libraries",
        "The2012 Methods process RNaseH-pretreated material and input DNA as parallel controls. Their sequencing is not a fate of the captured experimental hybrid population. FigS5A's agaroseA/G experiment is an in-vitro validation, not an independently established genomic-IP bead recipe."))
    g.annotate(Highlight("library_substrate_boundary","RNA fate is an unresolved boundary",
        "Capture preserves the hybrid in this graph. Subsequent RNA fate, DNA-end processing and complete library structure remain unassigned; the2019 protocol is not represented as2012 chemistry."))
    g.declare_reviewed_scope(("rloop_chromatin","gentle_extraction","extract","restriction_cocktail",
        "restriction_digest","s96_antibody","capture_beads","immunoprecipitate","sonication",
        "fragment_captured_dna","size_selection","select_size","library_substrate_boundary",
        "library_adapters","dna_ligase","adapter_ligation","sequencing_library"),SOURCE)
    return g.require_publishable()
