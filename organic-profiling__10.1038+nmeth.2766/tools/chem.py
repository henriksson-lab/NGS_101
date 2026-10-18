"""ORGANIC source-owned native workflow and explicitly separate Orsi library illustration."""
from adapter_ligation import AdapterLigationGeometry
from chemdraw import Construct, Scene, Segment, feature, panel
import illumina as il
import seqprimers as sp
from pcr import PCRGeometry
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            SourceClaim, SourceDiscrepancy, StepKind)
from protocol_model import molecule, role
from endprep import EndPreparedDuplex, dA_tailed_scene
from protocol_graph import Disposition, Terminal

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
    dA_tailed_scene(EndPreparedDuplex(Construct([role("reference native MNase footprint")]), five_prime_phosphate=True, three_prime_overhang="A")),
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
    """Defining supplementary FigS1a: nuclease/protease phase then physical recovery."""
    rnase = graph.input(prefix + "_rnase", "RNase A; defining FigS1", MaterialKind.CATALYST)
    rna_digested = graph.product(prefix + "_rna_digested", "RNA-digested sample mixture", MaterialKind.MOLECULAR,state=dna_population("native DNA with RNA hydrolysates",rna=True,protein=True))
    graph.react(prefix + "_digest_rna", "Digest RNA in the native sample", StepKind.DIGESTION,
                ((template, InputRole.TEMPLATE), (rnase, InputRole.CATALYST)), (rna_digested,))
    proteinase = graph.input(prefix + "_proteinase", "Proteinase K; defining FigS1", MaterialKind.CATALYST)
    sds = graph.input(prefix + "_sds", "Native deproteinization condition; detailed formulation unassigned", MaterialKind.CONDITION)
    digested = graph.product(prefix + "_protein_digested", "Protein/RNA-digested mixture before purification", MaterialKind.MOLECULAR,state=dna_population("native DNA after protein digestion; products remain until recovery",rna=True,protein=True))
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

def _native_library_branch(graph, prefix, template, reviewed):
    """Defining Supplementary Figure1c workflow; historical enzymes/oligos unassigned."""
    repair=graph.input(prefix+"_repair_enzyme","native end-repair reagents; identities unassigned",MaterialKind.CATALYST)
    repaired=graph.product(prefix+"_repaired","end-repaired native fragments",MaterialKind.MOLECULAR,state=dna_population("repaired native DNA"))
    graph.react(prefix+"_repair","End repair (defining FigS1c)",StepKind.END_REPAIR,((template,InputRole.TEMPLATE),(repair,InputRole.CATALYST)),(repaired,))
    def recover(key, before, label):
        reagent=graph.input(key+"_reagent",label,MaterialKind.CONDITION)
        kept=graph.product(key+"_kept","recovered "+before.label,MaterialKind.MOLECULAR,state=before.state)
        waste=graph.product(key+"_waste","physically removed recovery fraction",MaterialKind.WASTE)
        graph.react(key,label,StepKind.CLEANUP,((before,InputRole.TEMPLATE),(reagent,InputRole.CONDITION)),(kept,waste),changes_molecule=False)
        graph.discard(waste,"Physical fraction removed; class-specific depletion unassigned.")
        reviewed.extend((key+"_reagent",key));return kept
    clean=recover(prefix+"_repair_phenol",repaired,"Phenol/chloroform after end repair")
    clean=recover(prefix+"_repair_s300",clean,"S-300 spin-column purification after end repair")
    enzyme=graph.input(prefix+"_tail_enzyme","native dA-tailing reagent; enzyme unassigned",MaterialKind.CATALYST)
    prep=EndPreparedDuplex(Construct([role("native ORGANIC fragment")]),five_prime_phosphate=False,three_prime_overhang="A")
    tailed=graph.product(prefix+"_tailed","opposed native 3′dA ends;5′ chemistry unassigned",MaterialKind.MOLECULAR,state=dA_tailed_scene(prep))
    graph.react(prefix+"_tail","Add3′A overhangs (defining FigS1c)",StepKind.TAILING,((clean,InputRole.TEMPLATE),(enzyme,InputRole.CATALYST)),(tailed,))
    tailed=recover(prefix+"_tail_phenol",tailed,"Phenol/chloroform after A-tailing")
    tailed=recover(prefix+"_tail_s300",tailed,"S-300 spin-column purification after A-tailing")
    adapter=graph.input(prefix+"_adapter","native Illumina TruSeq adapters named by FigS1c; historical member/strands unassigned",MaterialKind.ADAPTER)
    ligase=graph.input(prefix+"_ligase","native adapter ligase; identity unassigned",MaterialKind.CATALYST)
    adapted=graph.product(prefix+"_adapted","native adapter-ligated DNA",MaterialKind.MOLECULAR,state=native_library_state("actual adapter-ligated native precursor"))
    graph.react(prefix+"_ligate","Ligate native adapters (defining FigS1c)",StepKind.LIGATION,((tailed,InputRole.TEMPLATE),(adapter,InputRole.STRUCTURAL),(ligase,InputRole.CATALYST)),(adapted,),geometry=AdapterLigationGeometry.unavailable("Native3′A substrate known; the exact historical TruSeq reagent and its terminal/strand chemistry remain unassigned. Later Orsi oligos are a separate reference."))
    adapted=recover(prefix+"_postligation_beads",adapted,"AMPureXP bead purification before PCR (FigS1c)")
    f=graph.input(prefix+"_forward","native library forward-primer role unassigned",MaterialKind.OLIGO,inferred=True)
    r=graph.input(prefix+"_reverse","native library opposite-primer role unassigned",MaterialKind.OLIGO,inferred=True)
    p=graph.input(prefix+"_polymerase","native library-PCR enzyme; identity unassigned",MaterialKind.CATALYST)
    amplified=graph.product(prefix+"_amplified","INFERRED — native amplified library",MaterialKind.MOLECULAR,state=native_library_state("native copied DNA; end/index scheme unassigned"),inferred=True)
    residual=graph.product(prefix+"_residual","original adapter-ligated molecules after copying",MaterialKind.MOLECULAR,state=adapted.state)
    graph.react(prefix+"_pcr","INFERRED boundaries — native library amplification (FigS1c)",StepKind.PCR,((adapted,InputRole.TEMPLATE),(f,InputRole.STRUCTURAL),(r,InputRole.STRUCTURAL),(p,InputRole.CATALYST)),(amplified,residual),geometry=PCRGeometry.unavailable(adapted.state,reason="Actual recovered native adapter-ligated molecules; native primer identities and first-cycle boundaries remain unassigned, not inherited from the later reference."))
    graph.mark_terminal(Terminal(residual.id,Disposition.UNRESOLVED,"PCR copying does not physically remove original templates; their post-PCR bead recovery is unassigned."))
    final=recover(prefix+"_postpcr_beads",amplified,"AMPureXP bead purification after PCR (FigS1c)")
    library=graph.product(prefix+"_library","native "+prefix+" library",MaterialKind.LIBRARY,state=final.state,inferred=True)
    graph.react(prefix+"_retain","Retain recovered library",StepKind.GENERIC,((final,InputRole.TEMPLATE),),(library,),changes_molecule=False)
    graph.sequence(library,SequencingPlan(acquired_roles=("Read 1","Read 2"),read_lengths=READ_LENGTHS,limitation=SEQUENCING_ENDING))
    reviewed.extend((prefix+"_repair_enzyme",prefix+"_repair",prefix+"_tail_enzyme",prefix+"_tail",prefix+"_adapter",prefix+"_ligase",prefix+"_ligate",prefix+"_forward",prefix+"_reverse",prefix+"_polymerase",prefix+"_pcr",prefix+"_retain",prefix+"_library"))

def oligos():
    return (
        f"Orsi Universal Adapter — 5′-{ORSI_UNIVERSAL_ADAPTER[:-1]}*{ORSI_UNIVERSAL_ADAPTER[-1]}-3′",
        f"Orsi Indexed Adapter — 5′-/5Phos/{ORSI_INDEXED_ADAPTER[:-1]}*{ORSI_INDEXED_ADAPTER[-1]}-3′",
        f"P5 library-enrichment primer — 5′-{ORSI_P5_PCR}-3′",
        f"P7 library-enrichment primer — 5′-{ORSI_P7_PCR}-3′",
    )


ORSI_REFERENCE, ORSI_REFERENCE_PRIMERS = organic_library()
FINAL_LIBRARY = None
SEQ_PRIMERS = ()
TOPOLOGY_PANELS = (("Separate later Orsi adapter reference", panel(orsi_adapter_scene().rows(), caption="2015 Basic Protocol3 reagent reference; not a study-specific adapter assignment.")),)
TITLE = 'ORGANIC profiling — native affinity-purified chromatin footprints'
NOTES = '01_organic_profiling.html'
SOURCE = '<a href="https://doi.org/10.1038/nmeth.2766">Kasinathan et al., Nature Methods (2014)</a>; the fixed oligo example comes from the authors’ later <a href="https://doi.org/10.1002/0471142727.mb2131s110">detailed Current Protocols procedure (2015)</a>. The 2014 article instead cites <a href="https://doi.org/10.1073/pnas.1110731108">Henikoff et al. (2011)</a> for its modified library preparation.'
SUMMARY = 'MNase digests unfixed native chromatin before affinity purification of factor-associated DNA. The graph illustrates the yeast 80 mM salt branch; the original source also describes other salt conditions and a separate fly branch. Its paired endpoints and insert lengths report protected fragments.'
CAVEAT = 'ORGANIC relies on stable native protein–DNA association: weak or rapidly exchanging complexes can be lost, while MNase sequence preference and digestion depth affect boundaries. A titrated MNase ladder, total-digest input and biological replicates are integral controls. The exact adapter arms drawn here are source-backed by a later same-author ORGANIC procedure; the defining Supplementary Figure1c names TruSeq adapters, while the exact 2014 adapter and sample-index assignments remain unreconciled. The cited 2011 predecessor original has been recovered, but its detailed SI remains gated; the local author manuscript also lacks its cited Online Methods. These gaps do not establish that historical details were unpublished. The physical initial Y-adapter PCR geometry remains unconstructed.'
FINAL_CAPTION = 'The insert reports repaired MNase-derived fragment boundaries around DNA retained with the affinity-purified factor. The exact adapter arms and six-position index shown here illustrate the separate later Orsi implementation, whose full Universal Adapter is included once with its P5/Read1 overlap. Paired-end geometry is essential because short factor-associated fragments and longer co-occupied or nucleosomal fragments are distinguished by their full insert span.'
READ_LENGTHS = {'Read 1': 25, 'Read 2': 25}
SEQUENCING_INTRO = 'The defining study used 25-cycle paired-end HiSeq 2000 sequencing. Standard TruSeq primers read inward from the two native MNase boundaries. The six-position Index1 view belongs to the separately illustrated Orsi variant; it is not a recovered 2014 sample-index assignment or instrument mask.'
SEQUENCING_ENDING = ''
SEQUENCING_UNAVAILABLE = ''

def dna_population(label, *, rna=False, protein=False, carrier=False):
    sc=Scene();sc.strand("dna",[role(label)],label="native DNA population; strand pairing/ends unassigned")
    if rna:sc.strand("rna",[role("independent native RNA/hydrolysates; stage-dependent")],label="separate RNA-derived population")
    if protein:sc.strand("protein",[role("associated protein/proteolysis products; stage-dependent")],label="separate non-DNA material")
    if carrier:sc.strand("carrier",[role("GlycoBlue carrier; component-specific recovery unassigned")],label="independent recovery carrier")
    return sc

def native_library_state(label, *, carrier=False):
    sc=Scene();parts=[role("native library end A; strand/boundary unassigned",inferred=True),role(label,inferred=True),
         role("native indexed end(s); scheme/width unassigned",inferred=True,molecular_feature=feature("native_sample_index","sample_index","unknown"))]
    sc.strand("native",parts,label="INFERRED boundaries — native library molecule")
    sc.strand("opposite",[role("native opposite strand class; copied/paired extent unassigned",inferred=True)],label="native opposite DNA strand")
    if carrier:sc.strand("carrier",[role("possible co-carried GlycoBlue; recovery unassigned")],label="independent carrier")
    return sc

def protocol_graph() -> ProtocolGraph:
    final = FINAL_LIBRARY
    primers = SEQ_PRIMERS
    source = ("Kasinathan et al., Nature Methods (2014); doi:10.1038/nmeth.2766; "
              "separate illustrated library variant: Orsi et al., Current Protocols (2015)")
    graph = ProtocolGraph("ORGANIC profiling", authority=source)
    cells = graph.input(
        "native_cells", "unfixed yeast cells containing factor-bound chromatin",
        MaterialKind.BIOLOGICAL,
        state=dna_population("unfixed native chromatin DNA",rna=True,protein=True), primary=True)
    isolation = graph.input("nuclei_isolation", "source yeast nuclei-isolation condition",
                            MaterialKind.CONDITION)
    nuclei = graph.product(
        "native_nuclei", "native nuclei preserving factor–DNA complexes",
        MaterialKind.MOLECULAR,
        state=dna_population("native nuclear DNA complexes",rna=True,protein=True))
    graph.react("isolate_nuclei", "Isolate nuclei without formaldehyde",
                StepKind.EXTRACTION,
                ((cells, InputRole.TEMPLATE), (isolation, InputRole.CONDITION)), (nuclei,))
    calcium = graph.input("calcium", "source calcium activation condition", MaterialKind.CONDITION)
    mnase = graph.input("mnase", "titrated micrococcal nuclease", MaterialKind.CATALYST)
    digested = graph.product(
        "digested_chromatin", "native chromatin trimmed at exposed DNA",
        MaterialKind.MOLECULAR,
        state=dna_population("native MNase-cut DNA fragment populations",rna=True,protein=True))
    graph.react("mnase_digest", "Titrate MNase to a mostly mononucleosome ladder",
                StepKind.DIGESTION,
                ((nuclei, InputRole.TEMPLATE), (mnase, InputRole.CATALYST), (calcium, InputRole.CONDITION)), (digested,))
    edta = graph.input("edta_stop", "EDTA stop condition; selected yeast branch", MaterialKind.CONDITION)
    stopped = graph.product("stopped_chromatin", "quenched native yeast digest", MaterialKind.MOLECULAR, state=digested.state)
    graph.react("stop_mnase", "Stop yeast MNase digestion with EDTA", StepKind.QUENCHING,
                ((digested, InputRole.TEMPLATE), (edta, InputRole.CONDITION)), (stopped,))
    needle=graph.input("needle_condition","defining yeast FigS1a: sequential20/26gauge needle extraction",MaterialKind.CONDITION)
    s1=graph.product("s1","source soluble fractionS1",MaterialKind.MOLECULAR,state=dna_population("native S1 DNA/protein complexes",rna=True,protein=True))
    pellet=graph.product("first_pellet","first needle-extraction pellet for salt extraction",MaterialKind.MOLECULAR,state=stopped.state)
    graph.react("needle_extract","Needle-extract chromatin; centrifuge and retainS1",StepKind.EXTRACTION,((stopped,InputRole.TEMPLATE),(needle,InputRole.CONDITION)),(s1,pellet))
    low_salt=graph.input("low_salt","selected80mMNaCl extraction of firstpellet; source other150/600mM branches separate",MaterialKind.CONDITION)
    s2=graph.product("s2","source soluble fractionS2",MaterialKind.MOLECULAR,state=dna_population("native S2 DNA/protein complexes",rna=True,protein=True))
    insoluble=graph.product("insoluble_material","pellet after second salt extraction",MaterialKind.WASTE)
    graph.react("salt_extract","Extract firstpellet and retainS2",StepKind.SOLUBILIZATION,((pellet,InputRole.TEMPLATE),(low_salt,InputRole.CONDITION)),(s2,insoluble))
    graph.discard(insoluble,"Defining yeast FigS1a discards pellet after second extraction.")
    adjust=graph.input("adjust_salt","adjustS1 salt concentration to matchS2",MaterialKind.CONDITION)
    soluble=graph.product("soluble_complexes","combined native solubleS1/S2 complexes",MaterialKind.MOLECULAR,state=s1.state)
    graph.react("solubilize","AdjustS1 salt and combineS1/S2",StepKind.POOL,((s1,InputRole.TEMPLATE),(s2,InputRole.TEMPLATE),(adjust,InputRole.CONDITION)),(soluble,),changes_molecule=False)
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
        state=dna_population("affinity-bound native DNA fragment population",rna=True,protein=True))
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
        "needle_condition","needle_extract","low_salt","salt_extract","adjust_salt", "solubilize", "split_input", "affinity_matrix", "affinity_capture",
    ]
    reviewed.extend(recovery_reviewed)
    _native_library_branch(graph, "ip", footprints, reviewed)
    _native_library_branch(graph, "input", input_dna, reviewed)
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()

# The separate 2015 ordering reference above is not the defining experiment's endpoint.
SHOW_READ_CONTENT=False
READ_SPECS=(sp.read_spec("Read 1",25,"Defining source paired-end genomic boundary read","include"),
            sp.read_spec("Read 2",25,"Defining source paired-end opposite boundary read","include"),
            sp.read_spec("Index 1 (i7)",None,"Native historical index acquisition/width unassigned","undetermined"),
            sp.read_spec("Index 2 (i5)",None,"Native historical index scheme/acquisition unassigned","undetermined"))
SEQUENCING_ENDING="Defining source25-cycle paired-end sequencing. Native adapter/index identity, acquired index cycles and run-primer landing remain unassigned in this capped pass."
SEQUENCING_INTRO=""
CAVEAT += " INFERRED — native library ends and opposing primer roles are unassigned. The exact Orsi oligos are a separate2015 reference and do not assign historical native index geometry or sequencing cycles."
