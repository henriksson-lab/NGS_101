"""PARIS2016 adapter route with explicit physical RT primer, ss circle and source boundaries."""
from types import SimpleNamespace
from chemdraw import Construct,Scene,Segment,feature,oligo,revcomp
from protocol_model import molecule,role
from adapter_ligation import AdapterLigationGeometry
from circular import circularize_ssdna,CircularizationGeometry
from pcr import PCRGeometry
from protocol_graph import ProtocolGraph,InputRole,MaterialKind,StepKind,SequencingPlan,Terminal,Disposition
import seqprimers as sp
from paris_2018_oligos import oligos as later_oligos
TITLE='PARIS — psoralen analysis of RNA interactions and structures'
NOTES='01_paris.html'
SOURCE='Lu et al., Cell2016, <a href="https://doi.org/10.1016/j.cell.2016.04.028">doi:10.1016/j.cell.2016.04.028</a>; original mmc1 Steps1–3 and expressly inherited Flynn2016 icSHAPE Steps26–61, PMC4896316.'
SUMMARY='AMT selects base-paired RNA partners; proximity ligation joins the selected fragments before photoreversal. The selected2016 HEK293/mES route then ligates an adapter, reverse-transcribes, circularizes ss-cDNA and amplifies.'
CAVEAT='Three distinct source scopes: original2016 HeLa random-hexamer RT is a separate route; selected2016 HEK293/mES inherits Flynn2016 adapter priming;2018 comparison reagents remain outside actual graph molecules/read plans. The source does not assign which inherited3′-biotin versus3′-ddC linker variant PARIS used. PostRT variant-specific capture/RNA handling remains an explicit unresolved recovery boundary. Native ss-circle PCR copying/amplicon and read acquisition are not reconstructed from later oligos.'
ADAPTER='AGATCGGAAGAGCGGTTCAG'
RT_BODY='AGATCGGAAGAGCGTCGTGGA'
RT_BIND='TACTGAACCGC'
P5='AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT'
P3='CAAGCAGAAGACGGCATACGAGATCGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT'
INDEX=feature('paris_2016_sample','sample_index','unknown',note='4-base sample barcode; no historical member assigned')
UMI1=feature('paris_2016_tag_part1','umi','random',group='paris_2016_tag',part=1)
UMI2=feature('paris_2016_tag_part2','umi','random',group='paris_2016_tag',part=2)
def rt_parts():return [Segment('source DDDNN tag','DDDNN',placeholder=True,length_bp=5,feature=UMI1),Segment('source4nt sample barcode','NNNN',placeholder=True,length_bp=4,feature=INDEX),Segment('source N4 tag','NNNN',placeholder=True,length_bp=4,feature=UMI2),Segment('2016 RT body',RT_BODY),role('non-nucleotide iSp18 spacer1'),Segment('spacer-flanked BamHI tract','GGATCC'),role('non-nucleotide iSp18 spacer2'),Segment('adapter-binding RT arm',RT_BIND)]
# Compute the longest source-compatible3′ overlap instead of copying the blocked adapter suffix.
BIND_NT=max(n for n in range(1,len(RT_BIND)+1) if ADAPTER.endswith(revcomp(RT_BIND[-n:])))
def cdna_parts():return [*rt_parts(),Segment('copied upstream adapter sequence',revcomp(ADAPTER[:-BIND_NT])),role('RNA B-derived cDNA interval'),role('RNA A-derived cDNA interval')]
CDNA=Construct(cdna_parts(),name='PARIS actual2016 adapter-route first-strand cDNA')
CLOSED=circularize_ssdna(CDNA,five_prime_phosphate=True)
class CircleScene(Scene):
 def __init__(self):
  super().__init__();self.strand('circle',list(CDNA),label='actual spacer-containing ssDNA precursor');self.chemical_state=CLOSED
 def rows(self):return CLOSED.rows()
def cdna_scene():
 sc=Scene();sc.strand('cDNA',list(CDNA),label='actual first-strand target',mod5='P',mod3='OH');return sc

def duplex_scene():
 sc=Scene();sc.strand('RNA A',[role('RNA partnerA paired interval; extent unassigned')],label='RNA A');sc.strand('RNA B',[role('RNA partnerB paired interval; extent unassigned')],label='independent RNA B')
 sc.note('RNA A','AMT crosslinked base-paired duplex; native sequences/paired coordinates are unassigned, not fixed24nt fragments.')
 return sc

def chimera_scene(*,adapted=False):
 parts=[role('RNA A fragment'),role('RNA B fragment')]
 if adapted:parts.extend([Segment('inherited adapter core',ADAPTER),role('terminal blocker variant; biotin or appended-ddC unassigned')])
 sc=Scene();sc.strand('chimera',parts,label='single covalent RNA chimera'+(' with adapter' if adapted else ''),mod3='blocked;3′biotin/ddC variant unassigned' if adapted else 'OH')
 sc.junction('chimera',parts[0].name,parts[1].name,'RNA proximity ligation')
 if adapted:sc.junction('chimera',parts[1].name,parts[2].name,'3′adapter ligation')
 return sc
def adapter_scene():
 sc=Scene();sc.strand('adapter',[Segment('adapter core',ADAPTER),role('3′blocker variant; appendedddC versusbiotin unassigned')],label='source preadenylylated/blocked adapter family',mod5='rApp',mod3='blocked; variant unassigned');return sc

SEQ_PRIMERS=()
READ_LENGTHS=None
SEQUENCING_ENDING='Original PARIS2016 supplement specifies multiplexed IlluminaNextSeq; read lengths, acquired index cycles and exact custom read-entry geometry remain unassigned.'
SEQUENCING_INTRO=SEQUENCING_ENDING
SEQUENCING_UNAVAILABLE='Circular-template first-copy/amplicon proof is deferred model work, not a claim that publication data are absent.'
READ_SPECS=tuple(sp.read_spec(x,None,'Read/index acquisition unassigned for selected2016 route; no2018 or2015 run mask substituted','undetermined') for x in ('Read 1','Read 2','Index 1 (i7)','Index 2 (i5)'))
FINAL_LIBRARY=Construct([Segment('source P5 Solexa PCR end',P5),role('PCR-selected RNA-chimera/RT-derived body; tag/index geometry unassigned',molecular_feature=INDEX),Segment('source P3 Solexa reverse complement',revcomp(P3))],name='PARIS2016 source-limited expected PCR endpoint')
FINAL_CAPTION='Expected primer-bounded endpoint only; exact internal circle-selected geometry and modified-spacer bypass boundaries are unassigned. Actual PCR input remains a spacer-containing ss circle.'
def oligos():return ('<b>Actual2016 selected adapter route; four-base index member unassigned</b>',oligo('2016 inherited RT-primer family',rt_parts(),five='5′-/5Phos/-'),oligo('2016 inherited adapter core;3′block variant unassigned',[Segment('adapter core',ADAPTER)],five='5′-/5rApp/-',three='-/3′block unassigned/-3′'),oligo('2016 inherited P5-Solexa',[Segment('P5',P5)]),oligo('2016 inherited P3-Solexa',[Segment('P3',P3)]),'<b>Separate2016 HeLa random-hexamer route</b> —5′-/5phos/DDDNNXXXXNNNNNNTACCCTTCGCTTCACACACAAG/iSp18/GGATCC/iSp18/TACTGAACCGCNNNNNN-3′; XXXX=4nt sample barcode; firstN6 duplicate tag, terminalN6 RNA priming','<b>Later2018 comparator inventory only; not used in actual2016 graph/read plans</b>',*later_oligos())
def protocol_graph():
    """Build a psoralen contact graph from source-specific reagent choices."""
    source = SOURCE
    reagent="AMT psoralen"
    selection_label="12%native then20%urea2Dgel;excise crosslinkedRNA above diagonal"
    graph = ProtocolGraph(TITLE, authority=source)
    rna = graph.input("rna", "cellular RNA containing paired regions",
                      MaterialKind.BIOLOGICAL,
                      state=duplex_scene(), primary=True)
    crosslinker = graph.input("crosslinker", reagent + " 0.5 mg/ml in PBS; UV365 nm, 30 min on ice (mmc1 Step 1)",
                              MaterialKind.CONDITION)
    crosslinked = graph.product("crosslinked", reagent + "-crosslinked RNA duplexes",
                                MaterialKind.MOLECULAR,
                                state=duplex_scene())
    graph.react("crosslink", "Photocrosslink paired RNA with " + reagent,
                StepKind.CROSS_LINKING,
                ((rna, InputRole.TEMPLATE), (crosslinker, InputRole.CONDITION)),
                (crosslinked,))
    s1 = graph.input("s1_nuclease", "S1 nuclease in urea/SDS lysate, room temperature 10 min (mmc1 Step 1)",
                     MaterialKind.CATALYST)
    extracted = graph.product("extracted", "S1-digested, proteinase K/TRIzol-purified crosslinked RNA",
                              MaterialKind.MOLECULAR, state=crosslinked.state)
    graph.react("s1_extract", "S1 digestion, proteinase K (50 °C 30 min) and TRIzol extraction",
                StepKind.DIGESTION,
                ((crosslinked, InputRole.TEMPLATE), (s1, InputRole.CATALYST)), (extracted,))
    fragmentation = graph.input("fragmentation_condition",
                                "ShortCut RNase III with MnCl2, 37 °C 20 min (mmc1 Step 2)",
                                MaterialKind.CATALYST)
    fragments = graph.product("fragments", "fragmented crosslinked RNA duplexes",
                              MaterialKind.MOLECULAR, state=crosslinked.state)
    graph.react("fragment", "Fragment the crosslinked RNA with RNase III", StepKind.FRAGMENTATION,
                ((extracted, InputRole.TEMPLATE),
                 (fragmentation, InputRole.CATALYST)), (fragments,))
    selector = graph.input("selector", selection_label, MaterialKind.SELECTION)
    selected = graph.product("selected", "enriched crosslinked RNA duplexes",
                             MaterialKind.MOLECULAR, state=fragments.state)
    rejected = graph.product("rejected", "nonselected RNA fragments", MaterialKind.WASTE)
    graph.react("select_duplexes", "Enrich crosslinked RNA duplexes",
                StepKind.SELECTION,
                ((fragments, InputRole.TEMPLATE), (selector, InputRole.SELECTION)),
                (selected, rejected))
    graph.discard(rejected, "does not retain the selected crosslinked duplex population")
    ligase = graph.input("proximity_ligase", "T4 RNA ligase 1 (NEB M0437M), room temperature 16–20 h; boiled 2 min to stop (mmc1 Step 3)",
                         MaterialKind.CATALYST)
    chimera = graph.product("chimera", "proximity-ligated RNA chimera",
                            MaterialKind.MOLECULAR, state=chimera_scene())
    graph.react("proximity_ligation", "Proximity-ligate the paired RNA fragments",
                StepKind.PROXIMITY_LIGATION,
                ((selected, InputRole.TEMPLATE), (ligase, InputRole.CATALYST)), (chimera,))
    reversal = graph.input("photoreversal", "UV254 nm, 15 min on ice (mmc1 Step 3)",
                           MaterialKind.CONDITION)
    reversed_chimera = graph.product("reversed_chimera", "photoreversed RNA chimera",
                                     MaterialKind.MOLECULAR, state=chimera.state)
    graph.react("photoreverse", "Photoreverse the psoralen crosslink",
                StepKind.CHEMICAL_CONVERSION,
                ((chimera, InputRole.TEMPLATE), (reversal, InputRole.CONDITION)),
                (reversed_chimera,))
    adapter=graph.input("rna_adapter","inherited5′App adapter core;3′biotin/ddC variant unassigned",MaterialKind.ADAPTER,state=adapter_scene())
    adapter_ligase=graph.input("adapter_ligase","FlynnStep26 T4RNAligase1/highconcentration;15%PEG8000,25°C3h",MaterialKind.CATALYST)
    tagged=graph.product("tagged_rna","adapter-tagged RNA chimera",MaterialKind.MOLECULAR,state=chimera_scene(adapted=True))
    graph.react("ligate_adapter","Ligate source adapter on selected2016 HEK293/mES route",StepKind.LIGATION,((reversed_chimera,InputRole.TEMPLATE),(adapter,InputRole.STRUCTURAL),(adapter_ligase,InputRole.CATALYST)),(tagged,),geometry=AdapterLigationGeometry.unavailable("5′App/3′blocked core is source-known; exact3′biotin versus appended-ddC adapter variant remains unassigned for PARIS, so complete physical substrate is not fixed from Flynn's icSHAPE choices."))
    rna_gel=graph.input("rna_gel","denaturingPAGE;PARIS retains>40nt inserts;FlynnSteps27–39",MaterialKind.SELECTION)
    sized=graph.product("sized_tagged_rna","gel-recovered adapter-tagged chimera",MaterialKind.MOLECULAR,state=tagged.state)
    free=graph.product("free_linker","gel-excluded linker/short RNA products",MaterialKind.WASTE)
    graph.react("rna_size_select","Gel-purify adapter-ligatedRNA",StepKind.SELECTION,((tagged,InputRole.TEMPLATE),(rna_gel,InputRole.SELECTION)),(sized,free));graph.discard(free,"Outside physically excised source gel band.")
    primer=graph.input("library_rt_primer","Flynn2016 RT-primer family;AACC example not assigned",MaterialKind.OLIGO,state=Construct(rt_parts(),name="inherited source RT-primer family"))
    rt=graph.input("library_rt_enzyme","SuperScriptIII;25°C3min,42°C5min,52°C30min;FlynnStep42",MaterialKind.CATALYST)
    cdna=graph.product("library_cdna","actual source first-strand cDNA target",MaterialKind.MOLECULAR,state=cdna_scene())
    residual=graph.product("rt_template_rna","adapter-tagged RNA after copying",MaterialKind.MOLECULAR,state=tagged.state)
    graph.react("library_rt","Reverse-transcribe adapted RNA using inherited primer",StepKind.REVERSE_TRANSCRIPTION,((sized,InputRole.TEMPLATE),(primer,InputRole.STRUCTURAL),(rt,InputRole.CATALYST)),(cdna,residual))
    handling=graph.input("post_rt_handling","inheritedFlynn43–50 variant-specificcapture/RNAelution;PARISbranchunassigned",MaterialKind.CONDITION)
    recovered=graph.product("recovered_cdna","cDNA target entering source denaturingPAGE",MaterialKind.MOLECULAR,state=cdna.state)
    rna_fate=graph.product("unresolved_post_rt_rna","RNA templates/products after variant-dependent handling",MaterialKind.MOLECULAR)
    graph.react("post_rt_recovery","PostRT recovery boundary before source cDNA gel",StepKind.MAGIC,((cdna,InputRole.TEMPLATE),(residual,InputRole.TEMPLATE),(handling,InputRole.CONDITION)),(recovered,rna_fate),reason="FlynnSteps43–50 capture/biotin/RNase elution depends on linker variant; PARIS inherits this range but does not assign3′linker variant. No universal capture, RNA digestion or physical discard is invented.")
    graph.mark_terminal(Terminal(rna_fate.id,Disposition.UNRESOLVED,"RNA templates or hydrolysis products remain molecular; exact variant-specific treatment/recovery not reconstructed."))
    gel=graph.input("cdna_gel","FlynnSteps52–57 denaturingPAGE;above unextendedprimer;PARIS>40ntinserts",MaterialKind.SELECTION)
    selected=graph.product("sized_cdna","gel-purified ssDNA target",MaterialKind.MOLECULAR,state=cdna.state)
    primer_waste=graph.product("free_rt_primer","unextended RT primer/excluded gel products",MaterialKind.WASTE)
    graph.react("cdna_size_select","Separate actual cDNA from unextended primer",StepKind.SELECTION,((recovered,InputRole.TEMPLATE),(gel,InputRole.SELECTION)),(selected,primer_waste));graph.discard(primer_waste,"Excluded by source physical gel excision.")
    circligase=graph.input("circle_ligase","CircLigaseII/MnCl2;60°C2h;FlynnStep58",MaterialKind.CATALYST)
    circle=graph.product("cdna_circle","actual spacer-containing ssDNA circle",MaterialKind.MOLECULAR,state=CircleScene())
    graph.react("circularize_cdna","Seal cDNA3′OH to source5′phosphate",StepKind.CIRCULARIZATION,((selected,InputRole.TEMPLATE),(circligase,InputRole.CATALYST)),(circle,),geometry=CircularizationGeometry.known("single-strand intramolecular ligation",statement="Source inherited5′P primer andCircLigaseII close actualss cDNA. Physical iSp18 features remain in circle; this does not establish PCR copying through spacers."))
    column=graph.input("circle_column","FlynnStep59 ZymoDNAclean/concentrator5;elute~20µL",MaterialKind.SELECTION)
    clean=graph.product("clean_circle","purified ssDNA circle target",MaterialKind.MOLECULAR,state=circle.state)
    wash=graph.product("circle_wash","column flowthrough and wash fractions",MaterialKind.WASTE)
    graph.react("purify_circle","Purify circles before qPCR-monitored amplification",StepKind.CLEANUP,((circle,InputRole.TEMPLATE),(column,InputRole.SELECTION)),(clean,wash));graph.discard(wash,"Physically separated source column fractions.")
    f=graph.input("library_forward","Flynn2016 P5Solexa primer",MaterialKind.OLIGO,state=P5)
    r=graph.input("library_reverse","Flynn2016 P3Solexa primer",MaterialKind.OLIGO,state=P3)
    pol=graph.input("library_polymerase","PhusionHF/SYBRGreen;source qPCR-monitoredcycles",MaterialKind.CATALYST)
    library=graph.product("library_library","source-limited expected2016 PCR library",MaterialKind.LIBRARY,state=FINAL_LIBRARY)
    graph.react("library_pcr","Amplify actual ss-circle target with inherited primer pair",StepKind.PCR,((clean,InputRole.TEMPLATE),(f,InputRole.STRUCTURAL),(r,InputRole.STRUCTURAL),(pol,InputRole.CATALYST)),(library,),geometry=PCRGeometry.unavailable(clean.state,reason="Actual spacer-containing ss-circle is the template. Native first-copy/recopy/primer-bounded arc remains unconstructed model work; source primers do not justify a completed duplex input or copying through non-nucleotide spacers."))
    graph.sequence(library,SequencingPlan(library=FINAL_LIBRARY,limitation=SEQUENCING_ENDING+' '+SEQUENCING_UNAVAILABLE))
    graph.declare_reviewed_scope(('rna','crosslinker','crosslink','s1_nuclease','s1_extract','fragmentation_condition','fragment','selector','select_duplexes','proximity_ligase','proximity_ligation','photoreversal','photoreverse','rna_adapter','adapter_ligase','ligate_adapter','rna_gel','rna_size_select','library_rt_primer','library_rt_enzyme','library_rt','post_rt_handling','post_rt_recovery','cdna_gel','cdna_size_select','circle_ligase','circularize_cdna','circle_column','purify_circle','library_forward','library_reverse','library_polymerase','library_pcr','library_library'),SOURCE)
    return graph.require_publishable()

SPEC=SimpleNamespace(**{k:v for k,v in globals().copy().items() if k.isupper()},oligos=oligos,protocol_graph=protocol_graph)
