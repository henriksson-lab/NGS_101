"""SPLASH source-specific page specification and psoralen choices."""
from adapter_ligation import AdapterLigationGeometry
from chemdraw import Construct, Scene, Segment, feature
from circular import CircularizationGeometry
from interaction_capture import ReversibleDuplexCrosslink
from pcr import PCRGeometry
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            SourceClaim, SourceDiscrepancy, StepKind, Highlight, Disposition, Terminal)
from protocol_specs import Definition
from protocol_model import molecule, role
from splash_oligos import (INDEX1_PCR_PRIMER, RNA_LINKER,
                           UNIVERSAL_PCR_PRIMER, RT_LEFT, RT_MIDDLE, RT_RIGHT, oligos)

DEFINITION = Definition(
    "SPLASH — sequencing psoralen-crosslinked RNA hybrids",
    "10.1016/j.molcel.2016.04.028", "Aw et al., Molecular Cell (2016)",
    "Biotinylated psoralen photocrosslinks RNA duplexes for streptavidin selection, proximity ligation and sequencing.",
    "custom", ("RNA A-derived cDNA", "RNA B-derived cDNA"),
    "Aw et al. 2016 Experimental Procedures sequence the libraries on a NextSeq 500 after 9–12 PCR cycles with NEB Primers Set 1; Figure S2G specifies 2×150-base paired-end reads, as does the 2017 JoVE protocol. The official JoVE Materials table names NEB E7335L for library PCR. Exact circle-to-amplicon and read-primer geometry remain unresolved.",
    "RNA",
    "The official 2017 JoVE Materials table prints the custom adapter and RT primer, including 5rApp, 3ddC and two internal iSp18 modifications, and names the NEB E7335L PCR primer kit. The available NEB manual is revision 8.0 (June 2024), so its printed primer family is a later kit reference, not proof of the 2017 lot's selected index. Full circle-derived read geometry is unresolved. The recovered 2016 Molecular Cell supplement (Figure S2F and Supplemental Experimental Procedures) independently prints the same adapter and two-spacer RT primer, labels the cDNA 5′-phosphate, and shows circularization/PCR. Its 3′ adapter string does not annotate the activation/blocking groups supplied by the later author protocol. The defining 2016 main article (Experimental Procedures, read 2026-10-17) gives the reaction order and conditions the graph now uses where it differs from the 2017 protocol (Q5 PCR, 200–220 nt cDNA band); those differences are shown as source notices. Base-resolved circle-to-amplicon/read geometry is still a model gap.",
)

# Workflow reagent, workflow selection, graph reagent and graph selection.
PSORALEN = ("biotinylated psoralen",
            "streptavidin selects biotin-psoralen RNA duplexes",
            "biotinylated psoralen", "streptavidin capture")

# The author materials list names NEB E7335L; the available NEB Rev 8.0
# manual prints these Set 1 primer bases with a phosphorothioate bond to
# terminal T. Historical 2017 lot identity remains unverified.
# The reaction's exact circle-derived PCR template remains unavailable.
NEB_SET1_SOURCE = "https://www.neb.com/en/-/media/nebus/files/manuals/manuale7335_e7500_-e7710_e7730.pdf"

TITLE=DEFINITION.title
NOTES='01_splash.html'
SOURCE='<a href="https://doi.org/10.1016/j.molcel.2016.04.028">Aw et al., Molecular Cell2016</a>, Experimental Procedures and original supplement; JoVE2017 identified explicitly as later detailed reference.'
SUMMARY=DEFINITION.summary
CAVEAT=DEFINITION.caveat+' Native2016 PrimersSet1 members/widths and exact circle-to-amplicon geometry are unassigned; the later NEBRev8 primer reference is not assigned to the executed2016 lot. Actual modified ss-cDNA/circle states enter native PCR.'
INDEX=feature('splash_sample_index','sample_index','unknown',note='Native selected member/physical width unassigned')
FINAL_LIBRARY=Construct([role('native forward PCR end; full copied primer architecture unassigned',inferred=True),role('RNA-A/RNA-B-derived chimera cDNA'),role('native indexed opposite PCR end; member/width unassigned',inferred=True,molecular_feature=INDEX)],name='SPLASH native endpoint boundary')
SEQ_PRIMERS=()
FINAL_CAPTION='INFERRED native endpoint — no default completed library shell or selected later index is substituted.'
SEQUENCING_UNAVAILABLE='Source custom RNA linker and two-spacer RT primer are explicit; complete native circle junction/first-copy/PCR geometry and acquired indexes remain unassigned.'
SEQUENCING_ENDING='Defining2016 article: NextSeq500 sequencing after9–12PCR cycles. Existing FigureS2G review records PE150; later2017JoVE likewise states PE150, a reference rather than proof of native index acquisition.'
READ_SPECS=(__import__('seqprimers').read_spec('Read 1',150,'chimera-derived cDNA; native primer landing unassigned','include'),__import__('seqprimers').read_spec('Read 2',150,'opposite chimera-derived cDNA end','include'),__import__('seqprimers').read_spec('Index 1 (i7)','unassigned','Native index identity/acquisition unassigned','undetermined'))

def rt_segments():
    return [Segment('RT primer left arm',RT_LEFT),role('internal iSp18 spacer1'),Segment('central tract',RT_MIDDLE),role('internal iSp18 spacer2'),Segment('RT primer right arm',RT_RIGHT)]

def rt_state(*,extended=False,rna=False,circle=False):
    sc=Scene();parts=rt_segments()
    if extended:parts.append(role('antisense copy of linked RNA chimera; native extent unassigned'))
    sc.strand('cDNA' if extended else 'primer',parts,label='native circular ss-cDNA; exact closure sequence unassigned' if circle else 'modified ss-cDNA' if extended else 'source RT/circularization primer',mod5='' if circle else 'p',mod3='' if circle else 'OH')
    if circle:sc.note('cDNA','Covalent circle returns the copied3′end to the RT-primer5′end; exact base-resolved junction/firstPCR phase unassigned. 5′P consumed in sealing.')
    if rna:sc.strand('RNA',[role('independent linker-ligated RNA template; retained through RT')],label='RNA template, not covalently joined to cDNA')
    return sc


def recover(g,key,template,label,state=None):
    reagent=g.input(key+'_reagent',label,MaterialKind.SELECTION)
    out=g.product(key+'_out','source recovered '+template.label,MaterialKind.LIBRARY if template.kind==MaterialKind.LIBRARY else MaterialKind.MOLECULAR,state=template.state if state is None else state)
    waste=g.product(key+'_waste','physically removed selection supernatant/washes or gel fraction',MaterialKind.WASTE)
    g.react(key,label,StepKind.CLEANUP,((template,InputRole.TEMPLATE),(reagent,InputRole.SELECTION)),(out,waste),changes_molecule=False)
    g.discard(waste,'Source physical separation; no presumed complete depletion of every off-target class.')
    return out

SOURCE_OLIGO_EXCEPTIONS = {
    "rna_linker": ("The retrieved JoVE author PDF contains the reaction protocol but "
                   "omits the web article's materials-list sequences; the author's "
                   "HTML materials list prints 5rApp-CTGTAGGCACCATCAAT-3ddC."),
}


def linker_input():
    sc=Scene();sc.strand('linker',[Segment('activated RNA linker',RNA_LINKER)],mod5='rApp',mod3='ddC');return sc


def rna_linker_geometry() -> AdapterLigationGeometry:
    """The author protocol's activated linker meets a repaired RNA 3′-OH."""
    sc = Scene()
    sc.strand("RNA chimera", [Segment("RNA 3-prime end", "X" * 18,
                                      placeholder=True)], mod3="OH")
    sc.ligation_end("RNA chimera", "RNA 3-prime end", kind="single-strand",
                    polarity="3-prime", owner="substrate", text="RNA 3′-OH")
    sc.strand("RNA linker", [Segment("activated RNA linker", RNA_LINKER)],
              mod5="rApp", mod3="ddC")
    sc.ligation_end("RNA linker", "activated RNA linker", kind="single-strand",
                    polarity="5-prime", owner="adapter", text="linker 5′-rApp")
    return AdapterLigationGeometry.known(
        "single-strand", sc,
        statement="Aw et al. JoVE 2017 section 5 repairs RNA 3′ ends to OH; section 7 "
                  "ligates the printed 5′-rApp, 3′-ddC RNA linker with T4 Rnl2(tr).")


def protocol_graph():
    source = f"{DEFINITION.source_label}; doi:{DEFINITION.doi}"
    graph = ProtocolGraph(DEFINITION.title, authority=source)
    reagent, selection_label = PSORALEN[2:]
    topology = ReversibleDuplexCrosslink(reagent, selection_label)
    rna = graph.input("rna", "cellular RNA containing paired regions",
                      MaterialKind.BIOLOGICAL,
                      state=molecule("cellular RNA", (role("RNA strand classes; paired regions selected subsequently"),)), primary=True)
    crosslinker = graph.input("crosslinker", reagent + " and long-wave UV condition",
                              MaterialKind.CONDITION)
    crosslinked = graph.product("crosslinked", reagent + "-crosslinked RNA duplexes",
                                MaterialKind.MOLECULAR, state=topology.duplex())
    graph.react("crosslink", "Photocrosslink paired RNA with " + reagent,
                StepKind.CROSS_LINKING,
                ((rna, InputRole.TEMPLATE), (crosslinker, InputRole.CONDITION)),
                (crosslinked,))
    fragmentation = graph.input("fragmentation_condition", "Mg-based RNA fragmentation:9mM MgCl2/225mM KCl/150mM TrisHCl pH8.3;95°C5min",
                                MaterialKind.CONDITION)
    fragments = graph.product("fragments", "fragmented crosslinked RNA duplexes",
                              MaterialKind.MOLECULAR, state=crosslinked.state)
    graph.react("fragment", "Fragment the crosslinked RNA", StepKind.FRAGMENTATION,
                ((crosslinked, InputRole.TEMPLATE),
                 (fragmentation, InputRole.CONDITION)), (fragments,))
    fragments=recover(graph,"fragment_size",fragments,"6%TBE8Murea PAGE:90–110nt RNA fraction, defining2016")
    selector = graph.input("selector", selection_label, MaterialKind.SELECTION)
    selected = graph.product("selected", "enriched crosslinked RNA duplexes",
                             MaterialKind.MOLECULAR, state=fragments.state)
    rejected = graph.product("rejected", "nonselected RNA fragments", MaterialKind.WASTE)
    graph.react("select_duplexes", "Enrich crosslinked RNA duplexes",
                StepKind.SELECTION,
                ((fragments, InputRole.TEMPLATE), (selector, InputRole.SELECTION)),
                (selected, rejected))
    graph.discard(rejected, "does not retain the selected crosslinked duplex population")
    phosphatase = graph.input("end_phosphatase", "T4 PNK 3′-cyclic-phosphate repair",
                              MaterialKind.CATALYST)
    three_oh = graph.product("three_oh", "crosslinked RNA duplexes with 3′-OH ends",
                             MaterialKind.MOLECULAR, state=selected.state)
    graph.react("repair_three_prime", "Convert RNA 3′ cyclic phosphates to 3′-OH",
                StepKind.DEPHOSPHORYLATION,
                ((selected, InputRole.TEMPLATE), (phosphatase, InputRole.CATALYST)),
                (three_oh,))
    kinase = graph.input("end_kinase", "T4 PNK plus ATP for RNA 5′ phosphorylation",
                         MaterialKind.CATALYST)
    five_phosphate = graph.product("five_phosphate", "ligatable RNA duplex ends",
                                   MaterialKind.MOLECULAR, state=three_oh.state)
    graph.react("phosphorylate_five_prime", "Convert RNA 5′-OH ends to 5′ phosphate",
                StepKind.PHOSPHORYLATION,
                ((three_oh, InputRole.TEMPLATE), (kinase, InputRole.CATALYST)),
                (five_phosphate,))
    ligase = graph.input("proximity_ligase", "T4 RNA ligase 1 (NEB M0204L), 2.5 U/µL, 16 °C overnight",
                         MaterialKind.CATALYST)
    chimera = graph.product("chimera", "proximity-ligated RNA chimera",
                            MaterialKind.MOLECULAR, state=topology.chimera())
    graph.react("proximity_ligation", "Proximity-ligate the paired RNA fragments",
                StepKind.PROXIMITY_LIGATION,
                ((five_phosphate, InputRole.TEMPLATE),
                 (ligase, InputRole.CATALYST)), (chimera,))
    chimera=recover(graph,"release_chimeras",chimera,"Elute from streptavidin95°C10min inPKbuffer; TRIzol extraction and RNeasy cleanup")
    reversal = graph.input("photoreversal", "254 nm UV, 5 min on ice (2016 Experimental Procedures; JoVE section 6.1)",
                           MaterialKind.CONDITION)
    reversed_state = Scene()
    reversed_state.strand("chimera", [role("RNA partner A fragment"), role("RNA partner B fragment")], label="photoreversed RNA proximity-ligation product")
    reversed_state.junction("chimera", "RNA partner A fragment", "RNA partner B fragment", "retained RNA proximity-ligation junction")
    reversed_state.note("chimera", "Psoralen crosslink photoreversed; the covalent RNA ligation junction remains.")
    reversed_chimera = graph.product("reversed_chimera", "photoreversed RNA chimera",
                                     MaterialKind.MOLECULAR, state=reversed_state)
    graph.react("photoreverse", "Photoreverse the psoralen crosslink",
                StepKind.CHEMICAL_CONVERSION,
                ((chimera, InputRole.TEMPLATE), (reversal, InputRole.CONDITION)),
                (reversed_chimera,))
    linker = graph.input("rna_linker", "published activated 3′ RNA linker",
                         MaterialKind.ADAPTER, state=linker_input(),
                         source="Aw et al. JoVE 2017, materials list")
    linker_ligase = graph.input("linker_ligase", "T4 RNA ligase 2, truncated KQ (NEB M0373L), 25 °C, 2.5 h",
                                MaterialKind.CATALYST)
    linked_state = Scene()
    linked_state.strand("linked RNA", [
        Segment("RNA chimera", "X" * 36, placeholder=True),
        Segment("3-prime RNA linker", RNA_LINKER)], label="3′-linker-ligated RNA", mod5="p", mod3="ddC")
    linked_state.junction("linked RNA", "RNA chimera", "3-prime RNA linker",
                          "RNA 3′-OH to linker 5′-rApp")
    linked_rna = graph.product("linked_rna", "3′-linker-ligated RNA chimera",
                               MaterialKind.MOLECULAR, state=linked_state)
    graph.react("ligate_rna_linker", "Ligate the activated 3′ RNA linker",
                StepKind.LIGATION,
                ((reversed_chimera, InputRole.TEMPLATE),
                 (linker, InputRole.STRUCTURAL),
                 (linker_ligase, InputRole.CATALYST)),
                (linked_rna,), geometry=rna_linker_geometry())
    linked_rna=recover(graph,"linked_rna_size",linked_rna,"6%TBE8Murea PAGE:110–130nt linked RNA, defining2016")
    rt_primer = graph.input("library_rt_primer", "published 96-nt RT primer",
                            MaterialKind.OLIGO, state=rt_state(),
                            source="Original2016 supplement two-iSp18 primer; retained5′P interpretation from prior FigureS2F review; laterJoVE materials agrees on sequence")
    rt_enzyme = graph.input("library_rt_enzyme",
                            "SuperScript III, 208 nM RT primer (50 °C, 30 min)",
                            MaterialKind.CATALYST)
    cdna = graph.product("library_cdna", "RT-primer-extended SPLASH cDNA",
                         MaterialKind.MOLECULAR,
                         state=rt_state(extended=True,rna=True))
    graph.react("library_rt", "Reverse-transcribe the linker-ligated RNA",
                StepKind.REVERSE_TRANSCRIPTION,
                ((linked_rna, InputRole.TEMPLATE),
                 (rt_primer, InputRole.STRUCTURAL),
                 (rt_enzyme, InputRole.CATALYST)), (cdna,))
    hydrolysis = graph.input("rna_hydrolysis", "100 mM NaOH, 98 °C for 20 min (2016 Experimental Procedures; JoVE 7.13 adds 1 N stock)",
                             MaterialKind.CONDITION)
    free_cdna = graph.product("rna_free_cdna", "single-stranded SPLASH cDNA without RNA template",
                              MaterialKind.MOLECULAR, state=rt_state(extended=True))
    hydrolysed = graph.product("hydrolysed_rna", "alkali-hydrolysed RNA template",
                               MaterialKind.MOLECULAR)
    graph.react("hydrolyse_rna", "Hydrolyse the RNA template with NaOH",
                StepKind.CHEMICAL_CONVERSION,
                ((cdna, InputRole.TEMPLATE), (hydrolysis, InputRole.CONDITION)),
                (free_cdna, hydrolysed))
    graph.mark_terminal(Terminal(hydrolysed.id,Disposition.UNRESOLVED,"Chemical hydrolysis precedes physical PAGE. Exact hydrolysate recovery/clearance is not quantified; no complete removal assumed."))
    cdna_gel = graph.input("cdna_gel", "6% TBE-urea gel, excise 200–220 nt (2016; JoVE 7.17: 200–240 nt)",
                           MaterialKind.SELECTION)
    sized_cdna = graph.product("sized_cdna", "gel-purified 200–220 nt SPLASH cDNA",
                               MaterialKind.MOLECULAR, state=rt_state(extended=True))
    off_size = graph.product("off_size_cdna", "cDNA and primer outside the excised band",
                             MaterialKind.WASTE)
    graph.react("select_cdna", "Gel-select cDNA 96 nt longer than the ligated RNA",
                StepKind.SELECTION,
                ((free_cdna, InputRole.TEMPLATE), (cdna_gel, InputRole.SELECTION)),
                (sized_cdna, off_size))
    graph.discard(off_size, "left in the urea gel outside the excised slice")
    sized_cdna=recover(graph,"recover_cdna",sized_cdna,"Ethanol precipitation of eluted cDNA, defining2016")
    circle_ligase = graph.input("circle_ligase", "CircLigase II ssDNA ligase",
                                MaterialKind.CATALYST)
    circle = graph.product("circular_cdna", "circularized SPLASH cDNA",
                           MaterialKind.MOLECULAR,state=rt_state(extended=True,circle=True),
                           source="Aw et al. Molecular Cell2016 circularization; source5′P previously read in FigureS2F")
    graph.react("circularize_cdna", "Circularize single-stranded cDNA",
                StepKind.CIRCULARIZATION,
                ((sized_cdna, InputRole.TEMPLATE),
                 (circle_ligase, InputRole.CATALYST)), (circle,),
                geometry=CircularizationGeometry.known(
                    "single-strand intramolecular ligation",
                    statement="Aw et al. 2016 original Figure S2F labels first-strand cDNA 5′P and shows its circularization; "
                              "the 2017 detailed protocol section 7/Table 6 names CircLigase II. "
                              "A base-resolved circular junction is not assembled in this model."))
    circle=recover(graph,"circle_cleanup",circle,"DNA Clean&Concentrator-5 column after CircLigaseII, defining2016")
    reviewed = ["rna", "crosslinker", "crosslink", "fragmentation_condition", "fragment",
                "selector", "select_duplexes", "end_phosphatase", "repair_three_prime",
                "end_kinase", "phosphorylate_five_prime", "proximity_ligase",
                "proximity_ligation", "photoreversal", "photoreverse", "rna_linker",
                "linker_ligase", "ligate_rna_linker", "library_rt_primer",
                "library_rt_enzyme", "library_rt", "rna_hydrolysis", "hydrolyse_rna",
                "cdna_gel", "select_cdna", "circle_ligase", "circularize_cdna"]
    forward = graph.input("library_forward", "Native PrimersSet1 universal primer; historical strand/lot unassigned",
                          MaterialKind.OLIGO, state=None,
                          source="Native2016 Set1 strand/lot unassigned; later reference only: " + NEB_SET1_SOURCE)
    reverse = graph.input("library_reverse", "Native PrimersSet1 index primer; selected member/width unassigned",
                          MaterialKind.OLIGO, state=None,
                          source="Native2016 Set1 member/width unassigned; later reference only: " + NEB_SET1_SOURCE)
    polymerase = graph.input("library_polymerase",
                             "Q5 DNA polymerase, 9–12 cycles (2016 Experimental Procedures)",
                             MaterialKind.CATALYST)
    library = graph.product("library_library", "SPLASH sequencing library",
                            MaterialKind.LIBRARY, state=FINAL_LIBRARY)
    pcr_residual=graph.product('pcr_residual','native circular template strands and other PCR species; recovery unassigned',MaterialKind.MOLECULAR,state=circle.state)
    graph.react("library_pcr", "Amplify and index the SPLASH cDNA circle",
                StepKind.PCR,
                ((circle, InputRole.TEMPLATE),
                 (forward, InputRole.STRUCTURAL),
                 (reverse, InputRole.STRUCTURAL),
                 (polymerase, InputRole.CATALYST)),
                (library,pcr_residual), geometry=PCRGeometry.unavailable(
                    circle.state,
                    reason="Actual modified circular ss-cDNA is the input; exact native junction/first-copy primer geometry and historical kit member remain unassigned. A completed library is not its template."))
    graph.mark_terminal(Terminal(pcr_residual.id,Disposition.UNRESOLVED,'PCR does not establish destruction/removal of the modified originalcircle or off-target products. Native class recovery at final gel selection remains unassigned.'))
    graph.annotate(SourceDiscrepancy(
        "photoreverse", "Photoreversal wavelength in the 2017 author protocol",
        (SourceClaim("Aw et al. JoVE 2017, section 6.1",
                     "irradiate the eluted RNA under UV 254 nm for 5 min on ice"),
         SourceClaim("Aw et al. JoVE 2017, Figure 1 legend",
                     "reverse crosslinking under UV light at 265 nm")),
        model_choice="The graph uses 254 nm, the value in the executable procedure step.",
        evidence="The numbered procedure gives the irradiation condition, duration and "
                 "temperature; the figure legend is a workflow summary.",
        uncertainty="The defining 2016 Molecular Cell Experimental Procedures also "
                    "irradiate at UV 254 nm for 5 min on ice, supporting the model "
                    "choice; the 265 nm legend value remains unexplained."))
    graph.annotate(Highlight('linked_rna_size','Defining2016 route and later2017 reference','The model selects2016 linkedRNA110–130nt, cDNA200–220nt and Q5 9–12PCR. Later2017JoVE uses110–140nt,200–240nt and Phusion/testPCR. These are version-scoped choices; ordinary later method changes do not establish errors in either publication.'))
    library=recover(graph,'final_size',library,'3%GTGNusieve agarose:200–300bp PCR products; QIAGEN gel extraction, defining2016')
    graph.sequence(library, SequencingPlan(library=FINAL_LIBRARY,acquired_roles=("Read 1","Read 2"),read_lengths={"Read 1":150,"Read 2":150},limitation=SEQUENCING_UNAVAILABLE))
    reviewed.extend(("library_forward", "library_reverse", "library_polymerase",
                     "library_pcr", "library_library"))
    graph.declare_reviewed_scope(REVIEWED_TARGETS, source)
    return graph.require_publishable()


REVIEWED_TARGETS=('rna','crosslinker','crosslink','fragmentation_condition','fragment','fragment_size_reagent','fragment_size','selector','select_duplexes','end_phosphatase','repair_three_prime','end_kinase','phosphorylate_five_prime','proximity_ligase','proximity_ligation','release_chimeras_reagent','release_chimeras','photoreversal','photoreverse','rna_linker','linker_ligase','ligate_rna_linker','linked_rna_size_reagent','linked_rna_size','library_rt_primer','library_rt_enzyme','library_rt','rna_hydrolysis','hydrolyse_rna','cdna_gel','select_cdna','recover_cdna_reagent','recover_cdna','circle_ligase','circularize_cdna','circle_cleanup_reagent','circle_cleanup','library_forward','library_reverse','library_polymerase','library_pcr','final_size_reagent','final_size','final_size_out')
