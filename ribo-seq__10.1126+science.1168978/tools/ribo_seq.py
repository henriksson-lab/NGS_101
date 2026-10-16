"""The poly(A)-tail library route of the defining Ingolia 2009 Ribo-seq paper.

The detailed Ingolia 2012 linker route is a later version, not evidence for
its oligo sequences or indexed read geometry in the 2009 experiment.
"""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
from chemdraw import Construct
from pcr import PCRGeometry
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind
from protocol_model import role

TITLE = "Ribo-seq — original 2009 poly(A)-tail workflow"
NOTES = "01_ribo-seq.html"
SOURCE = ('<a href="https://doi.org/10.1126/science.1168978">Ingolia et al. (2009)</a>, '
          'Figure 1 and Results; the authors’ later '
          '<a href="https://doi.org/10.1038/nprot.2012.086">2012 protocol</a>, '
          '“Linker ligation,” identifies what changed.')
SUMMARY = ("Nuclease digestion retains ribosome-protected RNA footprints. The original "
           "library route adds a poly(A) tail to each selected footprint before reverse "
           "transcription and DNA library preparation.")
CAVEAT = ("The 2009 supplementary methods are currently inaccessible. The original oligo "
          "sequences, their exact landing sites, and the final run structure remain unresolved; "
          "bracketed regions are roles, not reported bases.")
VERSION_NOTICE = ("<p><b>Source versions:</b> The 2009 Results say that the chosen library "
                  "strategy avoids RNA ligase. The authors’ 2012 Nature Protocols “Linker "
                  "ligation” section explicitly says their earlier method added a poly(A) "
                  "tail, then introduces a preadenylated 3′ linker and changed RT/PCR primers "
                  "for the optimized workflow. The linker and six-base index printed in 2012 "
                  "therefore do not describe the 2009 library on this page.</p>")
CAVEAT += VERSION_NOTICE
FINAL_LIBRARY = Construct((role("5-prime sequencing arm", inferred=True),
                           role("ribosome footprint cDNA", length_bp=30),
                           role("poly(A)-derived junction", inferred=True),
                           role("3-prime sequencing arm", inferred=True)),
                          name="2009 ribosome-footprint sequencing library (roles only)")
FINAL_CAPTION = ("Role-level library endpoint. The poly(A)-derived junction is established by "
                 "the 2012 authors’ retrospective account; exact bases and 2009 read sites "
                 "await the original supplement.")
SEQ_PRIMERS = ()
SEQUENCING_ENDING = ("The 2009 study sequenced ribosome footprints, but the accessible main "
                     "article does not give the library’s exact primer sites or indexed run layout.")
SEQUENCING_UNAVAILABLE = "The original 2009 supplementary methods have not been inspected."


def protocol_graph():
    graph = ProtocolGraph(TITLE, authority=SOURCE)
    mrna = graph.input("translated_mrna", "translated mRNA in ribosome complexes",
                       MaterialKind.BIOLOGICAL, primary=True)
    nuclease = graph.input("nuclease", "footprinting nuclease", MaterialKind.CATALYST)
    protected = graph.product("protected", "ribosome-protected fragments", MaterialKind.MOLECULAR)
    digested = graph.product("digested", "unprotected RNA digest", MaterialKind.WASTE)
    graph.react("digest", "Digest RNA outside ribosomes", StepKind.DIGESTION,
                ((mrna, InputRole.TEMPLATE), (nuclease, InputRole.CATALYST)),
                (protected, digested))
    graph.discard(digested, "unprotected RNA is digested")
    selection = graph.input("size_selection", "monosome and footprint-size selection",
                            MaterialKind.SELECTION)
    footprint = graph.product("footprint", "selected ribosome footprint RNA",
                              MaterialKind.MOLECULAR,
                              state=Construct((role("ribosome footprint RNA", length_bp=30),),
                                              name="selected footprint"))
    rejected = graph.product("rejected", "other RNA fragments", MaterialKind.WASTE)
    graph.react("select", "Recover and size-select ribosome footprints", StepKind.SELECTION,
                ((protected, InputRole.TEMPLATE), (selection, InputRole.SELECTION)),
                (footprint, rejected))
    graph.discard(rejected, "not in the selected footprint fraction")
    tailing_enzyme = graph.input("tailing_enzyme", "poly(A)-tailing enzyme",
                                 MaterialKind.CATALYST)
    adenylate = graph.input("adenylate", "adenosine nucleotide substrate",
                            MaterialKind.CONDITION)
    tailed = graph.product("tailed_footprint", "poly(A)-tailed footprint RNA",
                           MaterialKind.MOLECULAR,
                           state=Construct((role("ribosome footprint RNA", length_bp=30),
                                            role("poly(A) tail, length unreported", inferred=True)),
                                           name="poly(A)-tailed footprint"))
    graph.react("polyadenylate", "Add a 3-prime poly(A) tail to selected footprints",
                StepKind.TAILING,
                ((footprint, InputRole.TEMPLATE), (tailing_enzyme, InputRole.CATALYST),
                 (adenylate, InputRole.CONDITION)), (tailed,))
    rt_primer = graph.input("rt_primer", "2009 reverse-transcription primer: sequence unresolved",
                            MaterialKind.OLIGO)
    rtase = graph.input("rtase", "reverse transcriptase", MaterialKind.CATALYST)
    cdna = graph.product("cdna", "footprint cDNA from tailed RNA", MaterialKind.MOLECULAR,
                         state=Construct((role("ribosome footprint cDNA", length_bp=30),
                                          role("poly(A)-derived junction", inferred=True),
                                          role("RT-primer-derived sequence", inferred=True)),
                                         name="first-strand cDNA"))
    graph.react("rt", "Reverse-transcribe poly(A)-tailed footprints",
                StepKind.REVERSE_TRANSCRIPTION,
                ((tailed, InputRole.TEMPLATE), (rt_primer, InputRole.STRUCTURAL),
                 (rtase, InputRole.CATALYST)), (cdna,))
    forward = graph.input("forward_primer", "2009 library forward primer: sequence unresolved",
                          MaterialKind.OLIGO)
    reverse = graph.input("reverse_primer", "2009 library reverse primer: sequence unresolved",
                          MaterialKind.OLIGO)
    polymerase = graph.input("polymerase", "library PCR polymerase", MaterialKind.CATALYST)
    library = graph.product("library", FINAL_LIBRARY.name, MaterialKind.LIBRARY,
                            state=FINAL_LIBRARY)
    graph.react("library_pcr", "Complete and amplify the original DNA library",
                StepKind.PCR,
                ((cdna, InputRole.TEMPLATE), (forward, InputRole.STRUCTURAL),
                 (reverse, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)),
                (library,), geometry=PCRGeometry.unavailable(
                    FINAL_LIBRARY, reason="Original 2009 supplement and primer sites unavailable; "
                                          "the 2012 indexed primers belong to a later version."))
    graph.sequence(library, SequencingPlan(limitation=SEQUENCING_ENDING))
    graph.declare_reviewed_scope(("translated_mrna", "nuclease", "digest", "size_selection",
                                  "select", "tailing_enzyme", "adenylate", "polyadenylate",
                                  "rt_primer", "rtase", "rt", "forward_primer",
                                  "reverse_primer", "polymerase", "library_pcr", "library"),
                                 SOURCE)
    return graph.require_publishable()
