"""Tan et al. 2018 Dip-C chemistry."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
from restriction import MBOI
from chemdraw import Construct, Segment
from single_cell_hic import ContactWorkflow, unresolved_illumina_library
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind
from protocol_model import molecule, role
WORKFLOW = ContactWorkflow("Dip-C", MBOI, "after ligation", "none", "META")
TITLE = "Dip-C — contact ligation followed by META amplification"
NOTES = "01_dip-c.html"
CITATION = 'Defining source: <a href="https://doi.org/10.1126/science.aat5641">Tan et al., <i>Science</i> (2018)</a>; detailed protocol: <a href="https://pmc.ncbi.nlm.nih.gov/articles/PMC8225968/">Tan et al. (2021)</a>.'
SUMMARY = "Dip-C omits biotin capture and uses multiplex end-tagging amplification (META), a transposon-based single-cell WGA route."
STEPS = (("Digest and ligate", "Process fixed chromatin with a compatible restriction workflow (MboI is drawn) and proximity-ligate without biotin fill-in."),
         ("Isolate one nucleus", "Sort a single ligated nucleus, lyse it and reverse crosslinks."),
         ("META", "Tag many molecule ends by transposition and extend/amplify them as whole-genome material."),
         ("Sequencing library", "Convert the amplified material to a paired-end indexed library for contact calling and haplotype imputation."))
READOUT = "Paired reads spanning chimeric products report contacts; dense single-cell coverage supports diploid haplotype reconstruction."
JUNCTION_CAPTION = "Dip-C retains ligation products without a biotin-selected fill junction. ** marks ligation."
LIBRARY_CAPTION = "A representative contact molecule after META. The defining paper identifies 39 transposon-derived bases at each read start but the accessible text does not disclose their sequence or a complete outer adapter structure."
LIBRARY_CAVEAT = "INFERRED — the 39-base META read prefixes are source-defined by length and role only; their bases and the complete sequencing shell remain unresolved."
SEQUENCING_UNAVAILABLE = "Dip-C reports paired-end sequencing and removal of the first 39 transposon-derived bases from both reads, but the defining source available here does not print the primer or transposon oligos needed for base-level placement."
def contact_product(): return WORKFLOW.junction()
META_READ_PREFIX_NT = 39
def final_library():
    tagged = Construct([
        Segment("META transposon prefix", "T"*META_READ_PREFIX_NT, "me", placeholder=True),
        *list(contact_product()),
        Segment("opposite META transposon prefix", "T"*META_READ_PREFIX_NT, "me", placeholder=True),
    ], name="META-tagged contact")
    return unresolved_illumina_library(tagged, "Dip-C META-derived library", indexed=None)


def protocol_graph() -> ProtocolGraph:
    """Dip-C ligates before nucleus isolation, then uses META without biotin capture."""
    g = ProtocolGraph(TITLE, authority=CITATION)
    sample = g.input("sample", "fixed chromatin nuclei", MaterialKind.BIOLOGICAL,
                     state=molecule("fixed chromatin", (role("chromosomal DNA"),)),
                     primary=True)
    mboi = g.input("mboi", "MboI restriction enzyme", MaterialKind.CATALYST)
    cut = g.product("cut", "MboI-cut chromatin", MaterialKind.MOLECULAR,
                    state=sample.state, chemical_state=WORKFLOW.digest())
    g.react("digest", "Digest chromatin with MboI", StepKind.DIGESTION,
            ((sample, InputRole.TEMPLATE), (mboi, InputRole.CATALYST)), (cut,))
    ligase = g.input("contact_ligase", "contact-ligation DNA ligase",
                     MaterialKind.CATALYST)
    contacts = g.product("contacts", "proximity-ligated contacts",
                         MaterialKind.MOLECULAR, state=contact_product())
    g.react("contact_ligation", "Ligate compatible ends without biotin fill-in",
            StepKind.PROXIMITY_LIGATION,
            ((cut, InputRole.TEMPLATE), (ligase, InputRole.CATALYST)), (contacts,))
    isolation = g.input("isolation", "single-nucleus sorting", MaterialKind.CONDITION)
    nucleus = g.product("nucleus", "one ligated nucleus", MaterialKind.MOLECULAR,
                        state=contacts.state)
    excluded = g.product("excluded", "other ligated nuclei", MaterialKind.WASTE)
    g.react("isolate", "Isolate a nucleus after contact ligation", StepKind.SAMPLING,
            ((contacts, InputRole.TEMPLATE), (isolation, InputRole.CONDITION)),
            (nucleus, excluded), changes_molecule=False)
    g.discard(excluded, "outside the retained single-nucleus reaction")
    transposon = g.input("meta_transposon", "META transposon adapter",
                         MaterialKind.ADAPTER, inferred=True,
                         source=LIBRARY_CAVEAT)
    transposase = g.input("meta_transposase", "META transposase",
                          MaterialKind.CATALYST)
    tagged = g.product("tagged", "META-tagged contact molecules",
                       MaterialKind.MOLECULAR, state=nucleus.state)
    offproducts = g.product("meta_offproducts", "nonamplifiable META products",
                            MaterialKind.WASTE)
    g.react("meta_tag", "Tag and amplify contact DNA by META",
            StepKind.TAGMENTATION,
            ((nucleus, InputRole.TEMPLATE), (transposon, InputRole.STRUCTURAL),
             (transposase, InputRole.CATALYST)), (tagged, offproducts))
    g.discard(offproducts, "does not enter the contact sequencing library")
    final = final_library()
    library = g.product("library", final.name, MaterialKind.LIBRARY,
                        state=final, inferred=True, source=LIBRARY_CAVEAT)
    g.react("meta_library", "Obtain the paired-end META sequencing library",
            StepKind.MAGIC, ((tagged, InputRole.TEMPLATE),), (library,),
            reason="The defining source and accessible detailed protocol do not print the complete transposon oligos or outer adapter chemistry; this transition cannot be reconstructed base by base.")
    g.sequence(library, SequencingPlan(acquired_roles=("Read 1", "Read 2"),
                                       limitation=SEQUENCING_UNAVAILABLE))
    g.declare_reviewed_scope(tuple(g.steps) + tuple(g.materials), CITATION)
    return g.require_publishable()
