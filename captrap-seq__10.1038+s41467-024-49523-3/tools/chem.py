"""CapTrap-seq dual full-length selection and long-read endpoints."""
from adapter_ligation import AdapterLigationGeometry
from chemdraw import Construct, MolecularState, Row, Segment, Workflow
from protocol_extensions import capped_full_length_cdna
from pcr import PCRGeometry
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph,
                            SequencingPlan, StepKind)
TITLE="CapTrap-seq — cap- and poly(A)-selected full-length cDNA"
NOTES="01_captrap-seq.html"
SOURCE='Defining source: <a href="https://doi.org/10.1038/s41467-024-49523-3">Glinos et al., <i>Nature Communications</i> (2024)</a>.'
SUMMARY="Anchored oligo(dT) selects the RNA 3-prime end, cap trapping selects cDNA that reached the native 5-prime cap, and end-dependent linker ligations retain molecules complete at both ends for long-read sequencing."
CAVEAT="CapTrap-seq cDNA was sequenced on both ONT and PacBio. The endpoint below shows the ONT adapter geometry; PacBio SMRTbell conversion is an alternative platform step, not a different CapTrap construct."
FINAL_LIBRARY=Construct([Segment("motor-loaded ONT adapter","[motor-loaded adapter]","r2",placeholder=True,inferred=True),Segment("5-prime linker","L"*14,placeholder=True),Segment("cap-selected full-length cDNA","X"*42,placeholder=True),Segment("copied poly(A) tail","T"*16),Segment("3-prime linker","R"*14,placeholder=True)],name="ONT-adapted CapTrap-seq cDNA")
FINAL_CAPTION="INFERRED — full-length-selected cDNA after platform adapter attachment; ONT adapter bases are proprietary."
SEQ_PRIMERS=()
SEQUENCING_ENDING="There is no synthesis sequencing primer on the ONT branch. A motor-loaded ligation adapter meters the full-length cDNA strand through the pore. For PacBio, hairpin ligation instead creates a SMRTbell and the platform sequencing primer anneals in the hairpin adapter."
def sections():
    return [("Prime at the poly(A) tail and synthesize first-strand cDNA",[Row(chunks=[("5′ cap—RNA—poly(A) || anchored oligo(dT) → cDNA", "r1", False)])],"Anchored oligo(dT) is the first 3-prime-end selection."),("Biotinylate the cap and capture RNA–cDNA hybrids",[Row(chunks=[("biotin—cap—RNA/cDNA ● streptavidin", "w1", False)])],"Only cDNA still connected to a capped RNA 5-prime end is retained."),("Release full-length cDNA and ligate both end-dependent linkers",capped_full_length_cdna().rows(),"Sequential 5-prime and 3-prime linker ligations select molecules carrying both native-end signatures."),("INFERRED — attach a long-read platform adapter",[Row(chunks=[("5′ linker—full-length cDNA—poly(A/T)—3′ linker → ONT or PacBio", None, False)])],"LA-PCR enriches doubly linked molecules before platform-specific library preparation; proprietary adapter bases are not disclosed."),]

def workflow():
    s=sections();w=Workflow(MolecularState(s[0][0],tuple(s[0][1])))
    w.react(s[1][0],s[1][1],note=s[1][2])
    selection=AdapterLigationGeometry.unavailable("CapTrap-seq establishes sequential 5′- and 3′-linker selection, but strand-resolved linker termini have not yet been verified from the supplementary source; no overhang is asserted.")
    w.ligate_adapter(s[2][0],selection,s[2][1],note=s[2][2])
    g=AdapterLigationGeometry.unavailable("The defining CapTrap-seq paper names the Oxford Nanopore ligation kit but does not disclose the motor-adapter end sequence or overhang geometry.")
    w.ligate_adapter("Attach the ONT ligation adapter",g,s[3][1],note=s[3][2])
    return w

def protocol_graph():
    source="Defining source: https://doi.org/10.1038/s41467-024-49523-3"
    ont_source=source+" (ONT motor-adapter bases are proprietary)"
    graph=ProtocolGraph("CapTrap-seq",authority=source)
    rna=graph.input("rna","capped polyadenylated RNA",MaterialKind.BIOLOGICAL,primary=True)
    oligo_dt=graph.input("anchored_oligo_dt","anchored oligo(dT) primer",MaterialKind.OLIGO)
    rt=graph.input("reverse_transcriptase","reverse transcriptase",MaterialKind.CATALYST)
    hybrid=graph.product("rna_cdna_hybrid","full-length RNA–cDNA hybrids",MaterialKind.MOLECULAR)
    graph.react("reverse_transcription","Prime at poly(A) and synthesize first-strand cDNA",
                StepKind.REVERSE_TRANSCRIPTION,
                ((rna,InputRole.TEMPLATE),(oligo_dt,InputRole.STRUCTURAL),
                 (rt,InputRole.CATALYST)),(hybrid,))
    cap_reagent=graph.input("cap_biotinylation","cap-biotinylation reagent",
                            MaterialKind.CONDITION)
    biotinylated=graph.product("biotinylated_hybrids","biotinylated capped RNA–cDNA hybrids",
                               MaterialKind.MOLECULAR)
    graph.react("cap_biotinylation_step","Biotinylate the native RNA cap",
                StepKind.CHEMICAL_CONVERSION,
                ((hybrid,InputRole.TEMPLATE),(cap_reagent,InputRole.CONDITION)),(biotinylated,))
    beads=graph.input("streptavidin","streptavidin capture beads",MaterialKind.SELECTION)
    captured=graph.product("captured_hybrids","cap-retained full-length hybrids",
                            MaterialKind.MOLECULAR)
    truncated=graph.product("truncated_cdna","uncapped or truncated cDNA",MaterialKind.WASTE)
    graph.react("cap_capture","Capture cap-connected RNA–cDNA hybrids",
                StepKind.AFFINITY_CAPTURE,
                ((biotinylated,InputRole.TEMPLATE),(beads,InputRole.SELECTION)),
                (captured,truncated))
    graph.discard(truncated,"cDNA that did not reach a native cap is removed")
    linkers=graph.input("end_linkers","5′- and 3′-end-dependent CapTrap linkers",
                        MaterialKind.ADAPTER)
    ligase=graph.input("ligase","linker and platform-adapter ligase",MaterialKind.CATALYST)
    linked=graph.product("doubly_linked_cdna","doubly linked full-length cDNA",
                         MaterialKind.MOLECULAR,state=capped_full_length_cdna())
    graph.react("linker_ligation","Release cDNA and ligate both end-dependent linkers",
                StepKind.LIGATION,
                ((captured,InputRole.TEMPLATE),(linkers,InputRole.STRUCTURAL),
                 (ligase,InputRole.CATALYST)),(linked,),
                geometry=AdapterLigationGeometry.unavailable(
                    "Strand-resolved linker termini have not been verified from the supplementary source."))
    primers=graph.input("la_pcr_primers","long-amplicon PCR-primer pair",MaterialKind.OLIGO)
    polymerase=graph.input("la_polymerase","long-amplicon PCR polymerase",MaterialKind.CATALYST)
    enriched=graph.product("enriched_cdna","LA-PCR-enriched doubly linked cDNA",
                           MaterialKind.MOLECULAR)
    graph.react("la_pcr","Enrich molecules carrying both end linkers",StepKind.PCR,
                ((linked,InputRole.TEMPLATE),(primers,InputRole.STRUCTURAL),
                 (polymerase,InputRole.CATALYST)),(enriched,),
                geometry=PCRGeometry.unavailable(
                    FINAL_LIBRARY,reason="The source establishes double-linker selection but the page does not model exact LA-PCR landing objects.",symmetric_species=True))
    ont_adapter=graph.input("ont_adapter","motor-loaded ONT ligation adapter",
                            MaterialKind.ADAPTER,inferred=True,source=ont_source)
    library=graph.product("library",FINAL_LIBRARY.name,MaterialKind.LIBRARY,state=FINAL_LIBRARY,
                          inferred=True,source=ont_source)
    graph.react("ont_ligation","Attach the motor-loaded ONT adapter",StepKind.LIGATION,
                ((enriched,InputRole.TEMPLATE),(ont_adapter,InputRole.STRUCTURAL),
                 (ligase,InputRole.CATALYST)),(library,),
                geometry=AdapterLigationGeometry.unavailable(
                    "The named Oxford Nanopore kit does not disclose motor-adapter end sequence or geometry."),
                source=ont_source)
    graph.sequence(library,SequencingPlan(limitation=SEQUENCING_ENDING))
    graph.declare_reviewed_scope((
        "rna","anchored_oligo_dt","reverse_transcriptase","reverse_transcription",
        "cap_biotinylation","cap_biotinylation_step","streptavidin","cap_capture",
        "end_linkers","ligase","linker_ligation","la_pcr_primers","la_polymerase",
        "la_pcr","ont_adapter","ont_ligation","library"),source)
    return graph.require_publishable()
