"""ORGANIC source-owned native workflow and explicitly separate Orsi library illustration."""
from adapter_ligation import AdapterLigationGeometry
from chemdraw import Construct, Scene, Segment, feature
import illumina as il
import seqprimers as sp
from pcr import PCRGeometry
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            SourceClaim, SourceDiscrepancy, StepKind)
from protocol_model import da_tailed, molecule, role

ORSI_UNIVERSAL_ADAPTER = il.P5[:-4] + il.TRUSEQ_READ1

ORSI_INDEXED_ADAPTER = il.CANONICAL_ARM_READ2 + "N" * 6 + il.P7_RC

ORSI_P5_PCR = il.P5[:21]

ORSI_P7_PCR = il.P7[:22]

SAMPLE_INDEX = feature("orsi_illustrative_i7", "sample_index", "unknown")


def indexed_adapter_segments():
    return [Segment("Orsi Read 2 adapter arm", ORSI_INDEXED_ADAPTER[:33], "r2"),
            Segment("six-position sample index; oligo orientation", "N" * 6, "i7",
                    placeholder=True, feature=SAMPLE_INDEX),
            Segment("Orsi P7 complement", ORSI_INDEXED_ADAPTER[39:], "p7")]


def orsi_adapter_scene():
    # The universal terminal T is a one-base 3′ overhang, not part of the paired core.
    scene = Scene()
    scene.strand("universal", [Segment("Orsi P5/Read1 Y arm", ORSI_UNIVERSAL_ADAPTER[:-13]),
                                Segment("universal paired core", ORSI_UNIVERSAL_ADAPTER[-13:-1]),
                                Segment("adapter terminal dT", ORSI_UNIVERSAL_ADAPTER[-1])],
                 label="Orsi Universal Adapter", mod3="terminal phosphorothioate")
    bottom = [Segment("indexed paired core", ORSI_INDEXED_ADAPTER[:12]),
              Segment("indexed Read2 Y arm", ORSI_INDEXED_ADAPTER[12:33]),
              *indexed_adapter_segments()[1:]]
    scene.anneal("indexed", bottom, to="universal",
                 pair=("indexed paired core", "universal paired core"),
                 unpaired=("indexed Read2 Y arm", "six-position sample index; oligo orientation", "Orsi P7 complement"),
                 label="Orsi Indexed Adapter", mod5="/5Phos/", mod3="terminal phosphorothioate")
    scene.overhang("universal", "adapter terminal dT", polarity="3-prime", owner="adapter",
                   text="unpaired source terminal dT")
    return scene


ORSI_DA_LIGATION = AdapterLigationGeometry.known("cohesive", orsi_adapter_scene(),
    da_tailed("Orsi illustrative dA-tailed insert", "native MNase footprint"),
    statement="Separate Orsi Basic Protocol3 supplies both protected Y-adapter oligos and the dA-tailed insert workflow. The universal adapter's terminal T is unpaired and complements the fragment dA; these exact oligos are an illustrated variant, not a demonstrated 2014 sample-index assignment.")


def organic_library():
    library = Construct([Segment("complete Orsi Universal Adapter; P5/Read1 overlap", ORSI_UNIVERSAL_ADAPTER, "r1"),
                         role("native MNase footprint retained with the affinity-purified factor"),
                         Segment("fragment-derived dA junction", "A"), *indexed_adapter_segments()],
                        name="ORGANIC readout illustrated with the separately documented Orsi adapter variant")
    primers = (sp.TRUSEQ["R1"], sp.TRUSEQ["I1"], sp.TRUSEQ["R2"])
    errors = sp.verify(library, primers, required_roles=tuple(p.role for p in primers))
    if errors:
        raise ValueError("Orsi-derived ORGANIC readout: " + "; ".join(errors))
    return library, primers




def _recover_native_dna(graph, prefix, template, reviewed):
    """Native branch inherited by the deposited protocol from Zentner2013§ChIP."""
    rnase = graph.input(prefix + "_rnase", "RNase A; native-protocol predecessor", MaterialKind.CATALYST)
    rna_digested = graph.product(prefix + "_rna_digested", "RNA-digested sample mixture", MaterialKind.MOLECULAR)
    graph.react(prefix + "_digest_rna", "Digest RNA in the native sample", StepKind.DIGESTION,
                ((template, InputRole.TEMPLATE), (rnase, InputRole.CATALYST)), (rna_digested,))
    proteinase = graph.input(prefix + "_proteinase", "Proteinase K; native-protocol predecessor", MaterialKind.CATALYST)
    sds = graph.input(prefix + "_sds", "Source SDS deproteinization condition", MaterialKind.CONDITION)
    digested = graph.product(prefix + "_protein_digested", "Protein/RNA-digested mixture before purification", MaterialKind.MOLECULAR)
    graph.react(prefix + "_digest_protein", "Digest protein after RNase treatment", StepKind.DEPROTEINIZATION,
                ((rna_digested, InputRole.TEMPLATE), (proteinase, InputRole.CATALYST), (sds, InputRole.CONDITION)), (digested,))
    extraction = graph.input(prefix + "_dna_extraction", "Source phenol/chloroform and ethanol DNA recovery", MaterialKind.CONDITION)
    dna = graph.product(prefix + "_dna", "Purified native DNA", MaterialKind.MOLECULAR,
                        state=molecule("recovered native DNA", (role("MNase-bounded DNA"),)))
    waste = graph.product(prefix + "_recovery_waste", "Material excluded during DNA purification", MaterialKind.WASTE)
    graph.react(prefix + "_purify_dna", "Extract and precipitate native DNA", StepKind.EXTRACTION,
                ((digested, InputRole.TEMPLATE), (extraction, InputRole.CONDITION)), (dna, waste))
    graph.discard(waste, "Physical DNA purification removes other material; nuclease/protease treatment itself is not called a physical discard.")
    reviewed.extend((prefix + "_rnase", prefix + "_digest_rna", prefix + "_proteinase", prefix + "_sds",
                     prefix + "_digest_protein", prefix + "_dna_extraction", prefix + "_purify_dna"))
    return dna

def _orsi_library_branch(graph, prefix, template, reviewed):
    """Later same-author Orsi variant; the defining Supplementary Figure1c gives a common input/IP library workflow."""
    final, primers = FINAL_LIBRARY, SEQ_PRIMERS
    polish = graph.input(prefix + "_repair_polymerases", "Orsi end-polish T4 DNA polymerase and Klenow fragment", MaterialKind.CATALYST)
    kinase = graph.input(prefix + "_kinase", "Orsi T4 PNK in the same end-polish mix", MaterialKind.CATALYST)
    repaired = graph.product(prefix + "_repaired", "End-polished/phosphorylated DNA; illustrated Orsi variant", MaterialKind.MOLECULAR)
    graph.react(prefix + "_repair", "Combined end polishing and phosphorylation (Orsi variant)", StepKind.END_REPAIR,
                ((template, InputRole.TEMPLATE), (polish, InputRole.CATALYST), (kinase, InputRole.CATALYST)), (repaired,))
    organic = graph.input(prefix + "_organic_cleanup", "Orsi phenol/chloroform and S-300 cleanup", MaterialKind.CONDITION)
    clean = graph.product(prefix + "_repair_clean", "Cleaned end-polished DNA", MaterialKind.MOLECULAR)
    repair_waste = graph.product(prefix + "_repair_waste", "Material removed by Orsi end-polish cleanup", MaterialKind.WASTE)
    graph.react(prefix + "_clean_repair", "Clean up after end polishing", StepKind.CLEANUP,
                ((repaired, InputRole.TEMPLATE), (organic, InputRole.CONDITION)), (clean, repair_waste))
    graph.discard(repair_waste, "Reported phenol/chloroform and S-300 cleanup; no exact molecular identities assigned.")
    tailing = graph.input(prefix + "_tailing_enzyme", "Orsi Klenow exo−, source materials M0212M", MaterialKind.CATALYST)
    datp = graph.input(prefix + "_datp", "Orsi dATP tailing condition", MaterialKind.CONDITION)
    tailed = graph.product(prefix + "_tailed", "Orsi dA-tailed ORGANIC fragment", MaterialKind.MOLECULAR,
                           state=da_tailed("dA-tailed DNA", "ORGANIC insert"))
    graph.react(prefix + "_tail", "Add fragment 3′ dA tails (Orsi variant)", StepKind.TAILING,
                ((clean, InputRole.TEMPLATE), (tailing, InputRole.CATALYST), (datp, InputRole.CONDITION)), (tailed,))
    tail_clean = graph.product(prefix + "_tail_clean", "Cleaned dA-tailed DNA", MaterialKind.MOLECULAR, state=tailed.state)
    tail_waste = graph.product(prefix + "_tail_waste", "Material removed by Orsi A-tail cleanup", MaterialKind.WASTE)
    graph.react(prefix + "_clean_tail", "Clean up after A-tailing", StepKind.CLEANUP,
                ((tailed, InputRole.TEMPLATE), (organic, InputRole.CONDITION)), (tail_clean, tail_waste))
    graph.discard(tail_waste, "Reported phenol/chloroform and S-300 cleanup.")
    adapters = graph.input(prefix + "_adapters", "Separately documented protected Orsi Y adapters", MaterialKind.ADAPTER,
                           state=orsi_adapter_scene().rows())
    ligase = graph.input(prefix + "_adapter_ligase", "Orsi Rapid DNA ligase, Enzymatics L6030-HC-L", MaterialKind.CATALYST)
    adapted = graph.product(prefix + "_adapted", "Y-adapter-ligated DNA; physical Y arms not fully constructed", MaterialKind.MOLECULAR)
    graph.react(prefix + "_ligation", "Ligate the illustrated Orsi adapters", StepKind.LIGATION,
                ((tail_clean, InputRole.TEMPLATE), (adapters, InputRole.STRUCTURAL), (ligase, InputRole.CATALYST)),
                (adapted,), geometry=ORSI_DA_LIGATION)
    positive_beads = graph.input(prefix + "_positive_beads", "Orsi positive AMPure cleanup; retain bead-bound DNA", MaterialKind.SELECTION)
    ligation_clean = graph.product(prefix + "_ligation_clean", "Bead-purified adapter-ligated DNA", MaterialKind.MOLECULAR)
    ligation_waste = graph.product(prefix + "_ligation_waste", "Unbound material from post-ligation cleanup", MaterialKind.WASTE)
    graph.react(prefix + "_clean_ligation", "Retain bead-bound ligated DNA", StepKind.CLEANUP,
                ((adapted, InputRole.TEMPLATE), (positive_beads, InputRole.SELECTION)), (ligation_clean, ligation_waste))
    graph.discard(ligation_waste, "Source positive cleanup retains bead-bound DNA and discards supernatant.")
    p5 = graph.input(prefix + "_p5_primer", "Orsi P5 enrichment primer", MaterialKind.OLIGO,
                     state=Construct([Segment("P5 primer", ORSI_P5_PCR)]))
    p7 = graph.input(prefix + "_p7_primer", "Orsi P7 enrichment primer", MaterialKind.OLIGO,
                     state=Construct([Segment("P7 primer", ORSI_P7_PCR)]))
    polymerase = graph.input(prefix + "_pcr_polymerase", "KAPA HiFi HotStart polymerase; later Orsi illustration", MaterialKind.CATALYST)
    amplified = graph.product(prefix + "_amplified", "PCR-amplified illustrated Orsi library", MaterialKind.MOLECULAR, state=final)
    graph.react(prefix + "_pcr", "PCR-enrich the illustrated Orsi library", StepKind.PCR,
                ((ligation_clean, InputRole.TEMPLATE), (p5, InputRole.STRUCTURAL),
                 (p7, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)), (amplified,),
                geometry=PCRGeometry.unavailable(molecule("initial Y-adapter-ligated DNA", (role("Y-adapter-ligated DNA"),)),
                    reason="Both Orsi enrichment primers are source-known, but the complete initial Y-arm strand state and first PCR copies are not constructed here. A fully complementary amplified endpoint is not substituted as the physical pre-PCR template."))
    library = graph.product(prefix + "_library", final.name, MaterialKind.LIBRARY, state=final)
    pcr_waste = graph.product(prefix + "_pcr_waste", "Unbound material from post-PCR cleanup", MaterialKind.WASTE)
    graph.react(prefix + "_clean_pcr", "Retain bead-bound amplified library", StepKind.CLEANUP,
                ((amplified, InputRole.TEMPLATE), (positive_beads, InputRole.SELECTION)), (library, pcr_waste))
    graph.discard(pcr_waste, "Source post-PCR positive AMPure cleanup; no consensus or mapping filter is represented as physical removal.")
    graph.sequence(library, SequencingPlan(final, tuple(primers), tuple(p.role for p in primers), READ_LENGTHS))
    reviewed.extend((prefix + "_repair_polymerases", prefix + "_kinase", prefix + "_repair", prefix + "_organic_cleanup", prefix + "_clean_repair",
        prefix + "_tailing_enzyme", prefix + "_datp", prefix + "_tail", prefix + "_clean_tail", prefix + "_adapters", prefix + "_adapter_ligase",
        prefix + "_ligation", prefix + "_positive_beads", prefix + "_clean_ligation", prefix + "_p5_primer", prefix + "_p7_primer",
        prefix + "_pcr_polymerase", prefix + "_pcr", prefix + "_clean_pcr", prefix + "_library"))

def oligos():
    return (
        f"Orsi Universal Adapter — 5′-{ORSI_UNIVERSAL_ADAPTER[:-1]}*{ORSI_UNIVERSAL_ADAPTER[-1]}-3′",
        f"Orsi Indexed Adapter — 5′-/5Phos/{ORSI_INDEXED_ADAPTER[:-1]}*{ORSI_INDEXED_ADAPTER[-1]}-3′",
        f"P5 library-enrichment primer — 5′-{ORSI_P5_PCR}-3′",
        f"P7 library-enrichment primer — 5′-{ORSI_P7_PCR}-3′",
    )


FINAL_LIBRARY, SEQ_PRIMERS = organic_library()
TITLE = 'ORGANIC profiling — native affinity-purified chromatin footprints'
NOTES = '01_organic_profiling.html'
SOURCE = '<a href="https://doi.org/10.1038/nmeth.2766">Kasinathan et al., Nature Methods (2014)</a>; the fixed oligo example comes from the authors’ later <a href="https://doi.org/10.1002/0471142727.mb2131s110">detailed Current Protocols procedure (2015)</a>. The 2014 article instead cites <a href="https://doi.org/10.1073/pnas.1110731108">Henikoff et al. (2011)</a> for its modified library preparation.'
SUMMARY = 'MNase digests unfixed native chromatin before affinity purification of factor-associated DNA. The graph illustrates the yeast 80 mM salt branch; the original source also describes other salt conditions and a separate fly branch. Its paired endpoints and insert lengths report protected fragments.'
CAVEAT = 'ORGANIC relies on stable native protein–DNA association: weak or rapidly exchanging complexes can be lost, while MNase sequence preference and digestion depth affect boundaries. A titrated MNase ladder, total-digest input and biological replicates are integral controls. The exact adapter arms drawn here are source-backed by a later same-author ORGANIC procedure; the defining Supplementary Figure1c names TruSeq adapters, while the exact 2014 adapter and sample-index assignments remain unreconciled. The cited 2011 predecessor original has been recovered, but its detailed SI remains gated; the local author manuscript also lacks its cited Online Methods. These gaps do not establish that historical details were unpublished. The physical initial Y-adapter PCR geometry remains unconstructed.'
FINAL_CAPTION = 'The insert reports repaired MNase-derived fragment boundaries around DNA retained with the affinity-purified factor. The exact adapter arms and six-position index shown here illustrate the separate later Orsi implementation, whose full Universal Adapter is included once with its P5/Read1 overlap. Paired-end geometry is essential because short factor-associated fragments and longer co-occupied or nucleosomal fragments are distinguished by their full insert span.'
READ_LENGTHS = {'Read 1': 25, 'Index 1 (i7)': 6, 'Read 2': 25}
SEQUENCING_INTRO = 'The defining study used 25-cycle paired-end HiSeq 2000 sequencing. Standard TruSeq primers read inward from the two native MNase boundaries. The six-position Index1 view belongs to the separately illustrated Orsi variant; it is not a recovered 2014 sample-index assignment or instrument mask.'
SEQUENCING_ENDING = ''
SEQUENCING_UNAVAILABLE = ''

def protocol_graph() -> ProtocolGraph:
    final = FINAL_LIBRARY
    primers = SEQ_PRIMERS
    source = ("Kasinathan et al., Nature Methods (2014); doi:10.1038/nmeth.2766; "
              "separate illustrated library variant: Orsi et al., Current Protocols (2015)")
    graph = ProtocolGraph("ORGANIC profiling", authority=source)
    cells = graph.input(
        "native_cells", "unfixed yeast cells containing factor-bound chromatin",
        MaterialKind.BIOLOGICAL,
        state=molecule("unfixed cells", (
            role("native chromatin with factor-bound and unbound DNA"),)), primary=True)
    isolation = graph.input("nuclei_isolation", "source yeast nuclei-isolation condition",
                            MaterialKind.CONDITION)
    nuclei = graph.product(
        "native_nuclei", "native nuclei preserving factor–DNA complexes",
        MaterialKind.MOLECULAR,
        state=molecule("native nuclei", (role("factor–DNA complexes"),)))
    graph.react("isolate_nuclei", "Isolate nuclei without formaldehyde",
                StepKind.EXTRACTION,
                ((cells, InputRole.TEMPLATE), (isolation, InputRole.CONDITION)), (nuclei,))
    calcium = graph.input("calcium", "source calcium activation condition", MaterialKind.CONDITION)
    mnase = graph.input("mnase", "titrated micrococcal nuclease", MaterialKind.CATALYST)
    digested = graph.product(
        "digested_chromatin", "native chromatin trimmed at exposed DNA",
        MaterialKind.MOLECULAR,
        state=molecule("native chromatin trimmed at exposed DNA", (
            role("MNase boundary"), role("protein-protected DNA footprint"),
            role("MNase boundary"), role("nucleosome-protected DNA"))))
    graph.react("mnase_digest", "Titrate MNase to a mostly mononucleosome ladder",
                StepKind.DIGESTION,
                ((nuclei, InputRole.TEMPLATE), (mnase, InputRole.CATALYST), (calcium, InputRole.CONDITION)), (digested,))
    edta = graph.input("edta_stop", "EDTA stop condition; selected yeast branch", MaterialKind.CONDITION)
    stopped = graph.product("stopped_chromatin", "quenched native yeast digest", MaterialKind.MOLECULAR, state=digested.state)
    graph.react("stop_mnase", "Stop yeast MNase digestion with EDTA", StepKind.QUENCHING,
                ((digested, InputRole.TEMPLATE), (edta, InputRole.CONDITION)), (stopped,))
    low_salt = graph.input("low_salt", "selected yeast 80 mM NaCl extraction condition",
                           MaterialKind.CONDITION)
    soluble = graph.product(
        "soluble_complexes", "soluble native protein–DNA complexes",
        MaterialKind.MOLECULAR,
        state=molecule("soluble native protein–DNA complexes", (
            role("epitope-tagged or antibody-recognized factor"),
            role("MNase-bounded protected fragment"))))
    insoluble = graph.product("insoluble_material", "insoluble nuclear material",
                              MaterialKind.WASTE)
    graph.react("solubilize", "Extract and combine source yeast soluble chromatin fractions",
                StepKind.SOLUBILIZATION,
                ((stopped, InputRole.TEMPLATE), (low_salt, InputRole.CONDITION)),
                (soluble, insoluble))
    graph.discard(insoluble, "does not enter the soluble native protein–DNA fraction")
    total_input = graph.product("total_digest_input", "matched total-digest input aliquot",
                                MaterialKind.MOLECULAR, state=soluble.state)
    capture_input = graph.product("capture_input", "complexes assigned to affinity purification",
                                  MaterialKind.MOLECULAR, state=soluble.state)
    graph.react("split_input", "Retain matched total-digest input",
                StepKind.SPLIT, ((soluble, InputRole.TEMPLATE),),
                (total_input, capture_input), changes_molecule=False)
    affinity = graph.input("affinity_matrix", "anti-FLAG affinity matrix; selected yeast branch",
                           MaterialKind.SELECTION)
    immunoprecipitate = graph.product(
        "immunoprecipitate", "factor-associated ORGANIC immunoprecipitate",
        MaterialKind.MOLECULAR,
        state=molecule("ORGANIC immunoprecipitate", (
            role("factor"), role("short protected DNA fragment"))))
    unbound = graph.product("unbound_complexes", "unbound affinity-capture/wash fraction",
                            MaterialKind.WASTE)
    graph.react("affinity_capture", "Affinity-purify factor-bound DNA",
                StepKind.AFFINITY_CAPTURE,
                ((capture_input, InputRole.TEMPLATE), (affinity, InputRole.SELECTION)),
                (immunoprecipitate, unbound))
    graph.discard(unbound, "is removed as the unbound/wash fraction; this does not prove that every discarded complex lacked the selected factor")
    recovery_reviewed = []
    footprints = _recover_native_dna(graph, "native_ip", immunoprecipitate, recovery_reviewed)
    input_dna = _recover_native_dna(graph, "native_input", total_input, recovery_reviewed)
    reviewed = [
        "native_cells", "nuclei_isolation", "isolate_nuclei", "calcium", "mnase", "mnase_digest", "edta_stop", "stop_mnase",
        "low_salt", "solubilize", "split_input", "affinity_matrix", "affinity_capture",
    ]
    reviewed.extend(recovery_reviewed)
    _orsi_library_branch(graph, "ip", footprints, reviewed)
    _orsi_library_branch(graph, "input", input_dna, reviewed)
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()
