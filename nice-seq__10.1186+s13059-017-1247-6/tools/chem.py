"""Protocol-owned source, oligos, workflow and reaction graph."""
from adapter_ligation import AdapterLigationGeometry
from chemdraw import Construct
from pcr import PCRGeometry
from protocol_graph import Disposition, Terminal, InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind
from protocol_model import da_tailed, molecule, role


def nice_library():
    """Source-established endpoint roles, without selecting a historical oligo kit."""
    return Construct([
        role("Illumina-compatible left library end"),
        role("sonicated genomic DNA carrying internal biotin-dA and biotin-dC labels"),
        role("Illumina-compatible right library end"),
    ], name="NicE-seq library"), ()


BEAD_RESUSPENSION_UL = 40
PCR_ALIQUOT_UL = 4
BEAD_REMAINDER_UL = BEAD_RESUSPENSION_UL - PCR_ALIQUOT_UL

FINAL_LIBRARY, SEQ_PRIMERS = nice_library()
TITLE = 'NicE-seq — nick-translation labelling of accessible chromatin'
NOTES = '01_nice_seq.html'
SOURCE = '<a href="https://doi.org/10.1186/s13059-017-1247-6">Ponnaluri et al., Genome Biology (2017)</a>; downstream library reagents from the named NEBNext Ultra DNA E7370 kit.'
SUMMARY = 'Standard-reaction illustration: Nt.CviPII nicks exposed CCD sites inside nuclei and DNA polymerase I converts each accessible nick into a tract containing biotin-dATP and/or biotin-dCTP. After random sonication and adapter ligation, streptavidin selects labelled fragments for PCR. The standard method uses a 4-µL aliquot of a 40-µL bead suspension; the remaining 36 µL has an unresolved fate in this model. The separate 25-cell branch sonicates all recovered DNA and captures it using 10 µL beads; its template handling is not substituted for this standard branch.'
CAVEAT = ("<b>Published article/deposit differences.</b> The article Methods specify "
          "10 U DNA polymerase I, overnight SDS/proteinase K treatment at 65°C, and "
          "0.05% Triton X-100 in B&W washes. The authors' GSE97889 NicE-seq sample "
          "extraction metadata, including GSM2579632, instead specify 50 U, 1 hour, "
          "and 0.005%, respectively. The model retains the shared operations without "
          "assigning universal amounts, duration or detergent concentration. Whether "
          "these differences reflect run variants or transcription errors remains unresolved. "
          "<b>Oligo identity under review.</b> The article names E7370S. The available "
          "E7370 manual is Version 9.1_1/23 and requires separate oligo kits; it does "
          "not establish the study's historical oligo choice. The endpoint therefore "
          "shows molecular roles without selecting index architecture or adapter bases.")
FINAL_CAPTION = 'Role-level endpoint: internally biotin-labelled genomic DNA between Illumina-compatible library ends. Fragment boundaries arise during later random sonication and end preparation. Exact adapter arms, index architecture and primer sites are not assigned by this model; the study-specific oligo identity remains an open source question.'
READ_LENGTHS = None
SEQUENCING_INTRO = ''
SEQUENCING_ENDING = 'The article processes paired-end NicE-seq reads. Its explicit 75-bp NextSeq statement describes the separate ATAC-seq comparison branch. Exact NicE-seq run cycles and study-specific oligos remain under source review; no index width, number of index reads or base-level primer binding site is selected here.'
SEQUENCING_UNAVAILABLE = ''
def oligos():
    return ()

def protocol_graph() -> ProtocolGraph:
    final = FINAL_LIBRARY
    source = ("Ponnaluri et al., Genome Biology (2017); doi:10.1186/s13059-017-1247-6; "
              "downstream library kit: NEBNext Ultra DNA E7370")
    graph = ProtocolGraph("NicE-seq", authority=source)
    cells = graph.input(
        "cells", "native or formaldehyde-fixed cells",
        MaterialKind.BIOLOGICAL,
        state=molecule("native or fixed cells", (
            role("nuclei containing open and nucleosome-protected chromatin"),)),
        primary=True)
    isolation = graph.input("nuclei_isolation", "NP-40 nuclei-isolation condition",
                            MaterialKind.CONDITION)
    nuclei = graph.product(
        "nuclei", "intact nuclei with accessible chromatin exposed",
        MaterialKind.MOLECULAR,
        state=molecule("intact nuclei", (role("open chromatin"),
                                          role("protein-occluded chromatin"))))
    graph.react("isolate_nuclei", "Isolate nuclei with NP-40", StepKind.EXTRACTION,
                ((cells, InputRole.TEMPLATE), (isolation, InputRole.CONDITION)), (nuclei,))
    nickase = graph.input("nickase", "Nt.CviPII nicking endonuclease",
                          MaterialKind.CATALYST)
    nicked = graph.product(
        "nicked_chromatin", "accessible CCD sites carrying single-strand nicks",
        MaterialKind.MOLECULAR,
        state=molecule("nicked accessible chromatin", (
            role("genomic DNA"), role("Nt.CviPII nick at accessible CCD"),
            role("genomic DNA"))))
    graph.react("nick", "Nick exposed CCD sites", StepKind.DIGESTION,
                ((nuclei, InputRole.TEMPLATE), (nickase, InputRole.CATALYST)), (nicked,))
    biotin_dntps = graph.input("biotin_dntps", "dNTPs including biotin-14-dATP and biotin-16-dCTP",
                               MaterialKind.MOLECULAR)
    polymerase_i = graph.input("polymerase_i", "DNA polymerase I",
                               MaterialKind.CATALYST)
    labelled = graph.product(
        "labelled_chromatin", "internally biotin-labelled open chromatin",
        MaterialKind.MOLECULAR,
        state=molecule("internally labelled open chromatin", (
            role("genomic DNA"), role("Nt.CviPII nick at accessible CCD"),
            role("DNA polymerase I replacement tract containing internal biotin-dA and biotin-dC"),
            role("genomic DNA"))))
    graph.react("nick_translation", "Perform biotin nick translation (same source labelling reaction)",
                StepKind.EXTENSION,
                ((nicked, InputRole.TEMPLATE), (biotin_dntps, InputRole.STRUCTURAL),
                 (polymerase_i, InputRole.CATALYST)), (labelled,))
    stop = graph.input("edta", "EDTA labelling stop", MaterialKind.CONDITION)
    rnase = graph.input("rnase", "RNase A", MaterialKind.CATALYST)
    stopped = graph.product("stopped_labelled_dna", "stopped reaction with digested RNA",
                            MaterialKind.MOLECULAR, state=labelled.state)
    graph.react("stop_and_rnase", "Stop with EDTA and digest RNA", StepKind.DIGESTION,
                ((labelled, InputRole.TEMPLATE), (stop, InputRole.CONDITION),
                 (rnase, InputRole.CATALYST)), (stopped,))
    reversal = graph.input("crosslink_reversal", "SDS and proteinase K at 65°C",
                           MaterialKind.CONDITION)
    deproteinized = graph.product("deproteinized_dna", "deproteinized labelled DNA mixture",
                                  MaterialKind.MOLECULAR, state=labelled.state)
    graph.react("deproteinize", "Remove proteins and reverse fixed-cell cross-links",
                StepKind.DEPROTEINIZATION,
                ((stopped, InputRole.TEMPLATE), (reversal, InputRole.CONDITION)),
                (deproteinized,))
    extraction = graph.input("dna_extraction", "phenol/chloroform DNA extraction",
                             MaterialKind.CONDITION)
    purified = graph.product(
        "purified_dna", "purified biotin-marked genomic DNA",
        MaterialKind.MOLECULAR,
        state=molecule("purified biotin-marked genomic DNA", (
            role("long genomic molecule with sparse internal biotin labels"),)))
    extraction_waste = graph.product("extraction_waste", "DNA extraction removed fraction",
                                    MaterialKind.WASTE)
    graph.react("purify", "Extract labelled genomic DNA", StepKind.EXTRACTION,
                ((deproteinized, InputRole.TEMPLATE), (extraction, InputRole.CONDITION)),
                (purified, extraction_waste))
    graph.discard(extraction_waste, "is removed during phenol/chloroform DNA recovery")
    sonication = graph.input("sonication", "approximately 150-bp sonication condition",
                             MaterialKind.CONDITION)
    fragments = graph.product(
        "fragments", "randomly sheared internally labelled fragments",
        MaterialKind.MOLECULAR,
        state=molecule("randomly sheared labelled fragments", (
            role("random shear end"), role("biotin-labelled accessible-chromatin tract"),
            role("random shear end"))))
    graph.react("sonicate", "Sonicate labelled genomic DNA", StepKind.FRAGMENTATION,
                ((purified, InputRole.TEMPLATE), (sonication, InputRole.CONDITION)),
                (fragments,))
    repair = graph.input("end_repair_mix", "NEBNext Ultra end-repair enzyme mix",
                         MaterialKind.CATALYST)
    repaired = graph.product("repaired_fragments", "end-repaired labelled fragments",
                             MaterialKind.MOLECULAR, state=fragments.state)
    graph.react("end_repair", "End-repair labelled fragments", StepKind.END_REPAIR,
                ((fragments, InputRole.TEMPLATE), (repair, InputRole.CATALYST)),
                (repaired,))
    tailing = graph.input("tailing_enzyme", "NEBNext Ultra dA-tailing enzyme",
                          MaterialKind.CATALYST)
    tailed = graph.product(
        "tailed_fragments", "dA-tailed labelled genomic fragments",
        MaterialKind.MOLECULAR,
        state=da_tailed("dA-tailed labelled genomic fragment",
                        "labelled or unlabelled genomic insert"))
    graph.react("a_tail", "Add one 3-prime dA to each end", StepKind.TAILING,
                ((repaired, InputRole.TEMPLATE), (tailing, InputRole.CATALYST)), (tailed,))
    adapter = graph.input("nebnext_adapter", "NEBNext Ultra E7370 adapter",
                          MaterialKind.ADAPTER,
                          state=molecule("study adapter identity under review", (role("NEB Illumina adapter"),)))
    ligase = graph.input("adapter_ligase", "NEBNext Ultra adapter ligase",
                         MaterialKind.CATALYST)
    adapted = graph.product("adapted_pool", "adapter-ligated labelled and unlabelled fragments",
                            MaterialKind.MOLECULAR, state=molecule("adapter-ligated pool", (
                                role("ligated adapter"), role("labelled or unlabelled genomic fragment"),
                                role("ligated adapter"))))
    graph.react("adapter_ligation", "Ligate the NEBNext Ultra adapter",
                StepKind.LIGATION,
                ((tailed, InputRole.TEMPLATE), (adapter, InputRole.STRUCTURAL),
                 (ligase, InputRole.CATALYST)), (adapted,),
                geometry=AdapterLigationGeometry.unavailable(
                    "The study-specific oligo identity and adapter end geometry remain under source review; E7370S alone does not select a complete oligo set."))
    streptavidin = graph.input("streptavidin", "streptavidin affinity beads",
                              MaterialKind.SELECTION)
    bound = graph.product("bound_library", "bead-bound internally biotinylated molecules",
                          MaterialKind.MOLECULAR, state=molecule("streptavidin-bound adapted DNA", (
                              role("ligated adapter"), role("internally biotin-labelled genomic fragment"),
                              role("ligated adapter"))))
    unlabelled = graph.product("unlabelled_library", "unbound capture and wash fraction",
                               MaterialKind.WASTE)
    graph.react("capture", "Capture internally biotinylated molecules directly after ligation",
                StepKind.AFFINITY_CAPTURE,
                ((adapted, InputRole.TEMPLATE), (streptavidin, InputRole.SELECTION)),
                (bound, unlabelled))
    graph.discard(unlabelled, "is removed by stringent streptavidin-bead washes")
    aliquot_rule = graph.input("aliquot_rule", "standard 4-of-40-µL aliquot selection",
                               MaterialKind.CONDITION)
    aliquot = graph.product("pcr_aliquot", f"{PCR_ALIQUOT_UL}-µL bead-suspension PCR aliquot",
                            MaterialKind.MOLECULAR, state=bound.state)
    remainder = graph.product("bead_remainder", f"{BEAD_REMAINDER_UL}-µL bead-suspension remainder",
                              MaterialKind.MOLECULAR, state=bound.state)
    graph.react("sample_beads", f"Resuspend beads in {BEAD_RESUSPENSION_UL} µL; take {PCR_ALIQUOT_UL} µL for PCR",
                StepKind.SAMPLING, ((bound, InputRole.TEMPLATE), (aliquot_rule, InputRole.CONDITION)),
                (aliquot, remainder), changes_molecule=False)
    graph.mark_terminal(Terminal(remainder.id, Disposition.UNRESOLVED,
                                 "The standard source says 4 µL of the 40-µL suspension was used; subsequent handling of the remainder is not assigned."))
    forward = graph.input("pcr_forward", "study library PCR forward primer (identity under review)",
                          MaterialKind.OLIGO)
    reverse = graph.input("pcr_reverse", "study library PCR reverse primer (identity under review)",
                          MaterialKind.OLIGO)
    pcr_polymerase = graph.input("pcr_polymerase", "E7370S library PCR polymerase",
                                 MaterialKind.CATALYST)
    library = graph.product("sequencing_library", "NicE-seq library",
                            MaterialKind.LIBRARY, state=final)
    graph.react("library_pcr", "Amplify the 4-µL bead aliquot (ten cycles usually sufficient)",
                StepKind.PCR,
                ((aliquot, InputRole.TEMPLATE), (forward, InputRole.STRUCTURAL),
                 (reverse, InputRole.STRUCTURAL), (pcr_polymerase, InputRole.CATALYST)),
                (library,), geometry=PCRGeometry.unavailable(
                    final,
                    reason="The study-specific oligos and pre-PCR adapter-ligated template geometry remain unresolved."))
    graph.sequence(library, SequencingPlan(library=final, limitation=SEQUENCING_ENDING))
    graph.declare_reviewed_scope((
        "cells", "nuclei_isolation", "isolate_nuclei", "nickase", "nick",
        "biotin_dntps", "polymerase_i", "nick_translation", "crosslink_reversal",
        "edta", "rnase", "stop_and_rnase", "deproteinize", "dna_extraction",
        "purify", "sonication", "sonicate", "end_repair_mix", "end_repair",
        "tailing_enzyme", "a_tail", "nebnext_adapter", "adapter_ligase",
        "adapter_ligation", "streptavidin", "capture", "sample_beads", "aliquot_rule", "bead_remainder", "pcr_forward", "pcr_reverse",
        "pcr_polymerase", "library_pcr", "sequencing_library",
    ), source)
    return graph.require_publishable()
