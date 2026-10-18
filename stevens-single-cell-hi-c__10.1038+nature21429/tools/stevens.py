"""Stevens et al. 2017 single-cell Hi-C chemistry."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
import illumina as il
import seqprimers as sp
from restriction import MBOI
from chemdraw import Construct, Scene, Segment, revcomp, complement_segments, feature
from endprep import dA_tailed_scene, EndPreparedDuplex
from pcr import PCRGeometry
from protocol_graph import (Disposition, InputRole, MaterialKind, MolecularRelation, ProtocolGraph,
                            SequencingPlan, StepKind, Terminal)
from adapter_ligation import AdapterLigationGeometry
from protocol_model import molecule, role
from restriction import Digest
HISTORICAL_PE_READ2=il.PE_READ2_PRIMER
DIGEST=Digest(MBOI,"X"*16,"X"*16)
FILLED=DIGEST.fill_in(biotin_base="A")
TITLE = "Stevens single-cell Hi-C — single nuclei before contact ligation"
NOTES = "01_stevens-single-cell-hi-c.html"
CITATION = 'Defining source: <a href="https://doi.org/10.1038/nature21429">Stevens et al., <i>Nature</i> (2017)</a>.'
SUMMARY = "Individual nuclei are isolated first; MboI contact ends are biotin-filled and proximity-ligated. AluI is used later to fragment captured DNA for library preparation."
STEPS = (("Isolate nuclei", "Lyse fixed cells and distribute individual nuclei before library construction."),
         ("Contact digest", "MboI cuts GATC in each isolated nucleus."),
         ("Mark and ligate", "Fill MboI ends with biotin-14-dATP, then ligate spatial neighbours."),
         ("Capture and library prep", "Capture contacts, digest captured DNA with AluI, A-tail, ligate a 3-bp-tagged Illumina adapter and PCR. Nextera XT was a published alternative."))
READOUT = "Paired reads report the two MboI fragments joined within one nucleus; the first three bases identify the library adapter."
JUNCTION_CAPTION = "Biotin-filled MboI ends form the contact boundary. AluI acts only after contact capture. ** marks ligation."
LIBRARY_CAPTION = "The source-printed adapter-ligation route, with the CAA 3-bp identification tag shown. The paper also reports a dual-index Nextera XT alternative."

ADAPTER_F=il.TRUSEQ_READ1+"CAA"+"T"
ADAPTER_R=revcomp(ADAPTER_F[-18:-1])+"GTTCAGCAGGAATGCCGAG"

INLINE_A=feature("cell_inline_a","cell_barcode","whitelist",whitelist="published12 identification-tag set",group="cell_id",part="Read1 end")
INLINE_B=feature("cell_inline_b","cell_barcode","whitelist",whitelist="published12 identification-tag set",group="cell_id",part="Read2 end")
def restriction_scene():
 sc=Scene()
 for k in ("left","right"):
  body=role(k+" restriction fragment; locus extent unassigned")
  sc.strand(k,[Segment(k+"5′GATC",DIGEST.enzyme.overhang),body],label="source restriction5′cohesive end",mod5="p")
  opp=complement_segments([body]);sc.anneal(k+" opposite",opp,to=k,pair=(opp[0].name,body.name),unpaired=(opp[0].name,),label="recessed3′OH; outer end unspecified")
  sc.overhang(k,k+"5′GATC",polarity="5-prime",owner="substrate",text="source restriction overhang")
 return sc

def filled_scene():
 sc=Scene()
 for key,parts in (("left",[role("left fragment body"),Segment("left filled end",DIGEST.enzyme.site[:DIGEST.enzyme.bottom_cut])]),("right",[Segment("right filled end",DIGEST.enzyme.site[DIGEST.enzyme.top_cut:]),role("right fragment body")])):
  sc.strand(key,parts,label="blunt filled end; biotin-dA in newly synthesized GATC",mod5="p")
  opposite=complement_segments(parts);sc.anneal(key+" opposite",opposite,to=key,pair=(opposite[0].name,parts[-1].name),unpaired=tuple(p.name for p in opposite if p.is_role_token()),label="opposite fragment strand; native outer end unassigned",mod5="p")
 return sc

def contact_product():
 return Construct([role("locus A; extent unassigned"),Segment("contact junction",FILLED.junction),role("locus B; extent unassigned")],name="source filled-end contact reference")
def contact_scene(*,modified=True):
 sc=Scene();parts=list(contact_product());sc.strand("original top",parts,label="original contact DNA; genomic methylation retained"+("; junction biotin-dA at positions "+str(FILLED.junction_biotin_top) if modified else ""))
 bottom=complement_segments(parts);sc.anneal("original bottom",bottom,to="original top",pair=(bottom[0].name,parts[0].name),unpaired=tuple(p.name for p in bottom if p.is_role_token()),label="opposite strand"+("; junction biotin-dA at drawn positions "+str(FILLED.junction_biotin_bottom) if modified else "; normal daughter DNA; original chemical labels not copied"))
 return sc
ALUI_CONTACT=contact_product()
ALUI_TAILED=EndPreparedDuplex(ALUI_CONTACT,five_prime_phosphate=True,three_prime_overhang="A")
def native_ligated_scene():
 # Separate covalent strands with source Y forks; full native contact lengths remain roles.
 sc=Scene();parts=[Segment("left Y arm",il.TRUSEQ_READ1[:19]),Segment("left paired core",il.TRUSEQ_READ1[19:]),Segment("inline tag","CAA",feature=INLINE_A),Segment("left adapter dT","T",note="Original adapter retains the source A*T phosphorothioate linkage upstream of this ligating T. The backbone modification is not copied by polymerase."),*contact_product(),Segment("right insert dA","A"),Segment("right reverse tag","TTG",feature=INLINE_B),Segment("right R core",revcomp(il.TRUSEQ_READ1[19:])),Segment("right Y arm","GTTCAGCAGGAATGCCGAG")]
 sc.strand("ligated strand",parts,label="actual covalent native Y-adapter/contact strand; original biotin/methylation retained")
 # The opposite strand is an independent molecule with the same source fork structure.
 opposite=[*parts[:4],*complement_segments(list(contact_product())),*parts[-4:]]
 sc.strand("opposite ligated strand",opposite,label="opposite native5′→3′ strand; full fork alignment not reconstructed")
 return sc

def amplified_scene():
 sc=Scene.duplex(list(final_library()),label="source-primer-derived daughter library",unpaired=tuple(p.name+"'" for p in final_library() if p.is_role_token()))
 sc.strand("residual originals",[role("biotinylated original contact templates/other products; physical recovery unresolved")],label="independent original-template lineage; labels not copied")
 return sc

def physical_wash(g,key,x,label):
 reagent=g.input(key+"_selection",label,MaterialKind.SELECTION);out=g.product(key+"_out",label+" retained fraction",MaterialKind.MOLECULAR,state=x.state);w=g.product(key+"_waste",label+" removed fraction",MaterialKind.WASTE)
 g.react(key,label,StepKind.CLEANUP,((x,InputRole.TEMPLATE),(reagent,InputRole.SELECTION)),(out,w));g.discard(w,"source physical wash; class-specific purity/depletion unassigned");return out

def source_adapter_scene():
 sc=Scene();stem=il.TRUSEQ_READ1[19:]
 sc.strand('F',[Segment('F Y arm',il.TRUSEQ_READ1[:19]),Segment('paired stem',stem),Segment('adapter identification tag','CAA',feature=INLINE_A),Segment('source3′dT','T')],label='source CAA adapter F',mod3='terminal phosphorothioate')
 sc.anneal('R',[Segment('opposite identification tag','TTG',feature=INLINE_B),Segment('opposite paired stem',revcomp(stem)),Segment('R Y arm','GTTCAGCAGGAATGCCGAG')],to='F',pair=('opposite paired stem','paired stem'),unpaired=('R Y arm',),label='source adapter R',mod5='p')
 sc.ligation_end('F','source3′dT',kind='cohesive',polarity='3-prime',owner='adapter',text='source3′dT pairs with insert3′dA')
 return sc

def adapter_scene():
    return source_adapter_scene()


def adapter_geometry():
    return AdapterLigationGeometry.known("cohesive", dA_tailed_scene(ALUI_TAILED), adapter_scene(),
        statement="Stevens2017 SI prints the CAA-tagged adapter pair, including its 17 bp core, 5-prime phosphate and terminal phosphorothioate.")

SEQ_PRIMERS = (
    sp.custom("Read 1", "historical Illumina Read 1", il.TRUSEQ_READ1,
              "Stevens et al. 2017 Supplementary Methods"),
    sp.custom("Read 2", "historical paired-end Read 2", HISTORICAL_PE_READ2,
              "Stevens et al. 2017 Supplementary Methods"),
)
REQUIRED_ROLES = ("Read 1", "Read 2")
READ_LENGTHS = {"Read 1": 75, "Read 2": 75}
PCR_FORWARD = il.P5 + il.TRUSEQ_READ1[4:]
PCR_REVERSE = il.P7 + HISTORICAL_PE_READ2


def primer_scene(name, sequence):
    sc = Scene()
    sc.strand(name, [Segment(name, sequence)], mod3="terminal phosphorothioate")
    return sc

def final_library():
    lib=Construct([Segment("P5",il.P5,"p5"),Segment("historical Read1 arm",il.TRUSEQ_READ1[4:],"r1",note="Forward daughter strand retains the source PCR primer terminal C*T phosphorothioate bond. Complementary copies do not inherit this backbone modification."),
        Segment("inline identification tag","CAA","cbc",feature=INLINE_A),Segment("source adapter dT junction","T"),
        *contact_product(),Segment("dA junction","A"),Segment("opposite inline tag","TTG","cbc",feature=INLINE_B),
        Segment("historical Read2 site",revcomp(HISTORICAL_PE_READ2),"r2",note="Opposite daughter strand retains the reverse PCR primer terminal phosphorothioate bond; the top-strand sequence copy does not copy that backbone chemistry."),Segment("P7 reverse complement",il.P7_RC,"p7")],name=TITLE+" source-derived adapter-route library")
    problems=sp.verify(lib,SEQ_PRIMERS,required_roles=REQUIRED_ROLES)
    if problems:raise ValueError("; ".join(problems))
    return lib


def protocol_graph() -> ProtocolGraph:
    """Stevens' single-nucleus MboI route with AluI only after contact capture."""
    g = ProtocolGraph(TITLE, authority=CITATION)
    sample = g.input("sample", "fixed chromatin nuclei", MaterialKind.BIOLOGICAL,
                     state=molecule("fixed chromatin", (role("chromosomal DNA"),)),
                     primary=True)
    isolation = g.input("isolation", "single-nucleus distribution", MaterialKind.CONDITION)
    nucleus = g.product("nucleus", "one isolated nucleus", MaterialKind.MOLECULAR,
                        state=sample.state)
    excluded = g.product("excluded", "other nuclei", MaterialKind.MOLECULAR)
    g.react("isolate", "Distribute individual nuclei before library construction",
            StepKind.SAMPLING,
            ((sample, InputRole.TEMPLATE), (isolation, InputRole.CONDITION)),
            (nucleus, excluded), changes_molecule=False)
    g.mark_terminal(Terminal(excluded.id, Disposition.UNRESOLVED, "Other nuclei are outside this retained single-nucleus lineage; their disposition is unrepresented."))
    agarose = g.input("agarose", "1% low-melting-point agarose overlay after imaging", MaterialKind.CONDITION)
    mboi = g.input("mboi", "MboI", MaterialKind.CATALYST)
    cut = g.product("cut", "MboI-cut chromatin", MaterialKind.MOLECULAR,
                    state=restriction_scene(), chemical_state=DIGEST)
    g.react("digest", "Digest single-nucleus chromatin with MboI",
            StepKind.DIGESTION,
            ((nucleus, InputRole.TEMPLATE), (mboi, InputRole.CATALYST), (agarose, InputRole.CONDITION)), (cut,))
    biotin = g.input("biotin_datp", "biotin-14-dATP", MaterialKind.MOLECULAR)
    fill_polymerase = g.input("fill_polymerase", "DNA Polymerase I Large (Klenow) Fragment",
                             MaterialKind.CATALYST)
    marked = g.product("marked", "biotin-filled MboI ends", MaterialKind.MOLECULAR,
                       state=filled_scene(), chemical_state=FILLED)
    g.react("fill", "Fill MboI overhangs with biotin-dATP", StepKind.EXTENSION,
            ((cut, InputRole.TEMPLATE), (biotin, InputRole.STRUCTURAL),
             (fill_polymerase, InputRole.CATALYST)), (marked,),
            relations=(MolecularRelation.RESTRICTION_FILL,))
    ligase = g.input("contact_ligase", "T4 DNA ligase",
                     MaterialKind.CATALYST)
    contacts = g.product("contacts", "single-nucleus contact molecules",
                         MaterialKind.MOLECULAR, state=contact_scene())
    g.react("contact_ligation", "Ligate biotin-filled MboI ends",
            StepKind.PROXIMITY_LIGATION,
            ((marked, InputRole.TEMPLATE), (ligase, InputRole.CATALYST)),
            (contacts,), relations=(MolecularRelation.CONTACT_JUNCTION,))
    contacts=physical_wash(g,"post_ligation_wash",contacts,"PBS wash after in-well ligation")
    heat = g.input("crosslink_reversal", "65 C overnight crosslink reversal", MaterialKind.CONDITION)
    reversed_dna = g.product("reversed_contacts", "crosslink-reversed single-cell contacts", MaterialKind.MOLECULAR, state=contacts.state)
    g.react("reverse_crosslinks", "Reverse crosslinks before streptavidin capture", StepKind.DEPROTEINIZATION,
            ((contacts, InputRole.TEMPLATE), (heat, InputRole.CONDITION)), (reversed_dna,))
    beads = g.input("streptavidin", "streptavidin beads", MaterialKind.SELECTION)
    captured = g.product("captured", "captured contact DNA", MaterialKind.MOLECULAR,
                         state=contacts.state)
    unbound = g.product("unbound", "contact-poor unbound DNA", MaterialKind.WASTE)
    g.react("capture", "Capture biotin-marked contact junctions",
            StepKind.AFFINITY_CAPTURE,
            ((reversed_dna, InputRole.TEMPLATE), (beads, InputRole.SELECTION)),
            (captured, unbound))
    g.discard(unbound, "physically removed agarose-containing binding supernatant; no class purity assigned")
    captured=physical_wash(g,"capture_wash",captured,"Source BW/BWT and reaction-buffer washes after capture")
    alui = g.input("alui", "AluI fragmentation enzyme", MaterialKind.CATALYST)
    fragments = g.product("fragments", "AluI-fragmented captured contacts",
                          MaterialKind.MOLECULAR, state=contact_scene())
    g.react("fragment", "Digest captured contacts with AluI", StepKind.DIGESTION,
            ((captured, InputRole.TEMPLATE), (alui, InputRole.CATALYST)),
            (fragments,))
    fragments=physical_wash(g,"alui_wash",fragments,"Source washes after AluI digestion")
    prep = g.input("end_prep", "Klenow Fragment (3-prime to 5-prime exo-minus)",
                   MaterialKind.CATALYST)
    tailed = g.product("tailed", "dA-tailed AluI fragments", MaterialKind.MOLECULAR,
                       state=ALUI_TAILED)
    g.react("tail", "dA-tail AluI fragments with Klenow exo-minus", StepKind.TAILING,
            ((fragments, InputRole.TEMPLATE), (prep, InputRole.CATALYST)), (tailed,))
    tailed=physical_wash(g,"tail_wash",tailed,"Source washes after dA tailing")
    adapter = g.input("adapter", "source-printed CAA-tagged Illumina adapter",
                      MaterialKind.ADAPTER, state=adapter_scene())
    adapter_ligase = g.input("adapter_ligase", "T4 DNA ligase for adapter ligation",
                             MaterialKind.CATALYST)
    ligated = g.product("ligated", "adapter-ligated contact fragments",
                        MaterialKind.MOLECULAR,state=native_ligated_scene())
    g.react("adapter_ligation", "Ligate the 3-bp-tagged adapter",
            StepKind.LIGATION,
            ((tailed, InputRole.TEMPLATE), (adapter, InputRole.STRUCTURAL),
             (adapter_ligase, InputRole.CATALYST)), (ligated,),
            geometry=adapter_geometry())
    ligated=physical_wash(g,"adapter_wash",ligated,"Source washes after adapter ligation and into Pfx buffer")
    forward = g.input("pcr_forward", "source-printed PCR forward primer", MaterialKind.OLIGO, state=primer_scene("PCR forward", PCR_FORWARD))
    reverse = g.input("pcr_reverse", "source-printed PCR reverse primer", MaterialKind.OLIGO, state=primer_scene("PCR reverse", PCR_REVERSE))
    pcr = g.input("pcr_polymerase", "Platinum Pfx DNA polymerase", MaterialKind.CATALYST)
    final = final_library()
    library = g.product("amplified", final.name, MaterialKind.MOLECULAR, state=amplified_scene())
    g.react("library_pcr", "Amplify adapter-ligated contacts for 25 cycles with Platinum Pfx", StepKind.PCR,
            ((ligated, InputRole.TEMPLATE), (forward, InputRole.STRUCTURAL),
             (reverse, InputRole.STRUCTURAL), (pcr, InputRole.CATALYST)),
            (library,), geometry=PCRGeometry.unavailable(
                ligated.state, reason="Actual native Y-adapter/contact strands; complete first-cycle alignment bounded."))
    ampure = g.input("ampure", "AMPure XP beads, 0.7x post-PCR recovery", MaterialKind.SELECTION)
    cleaned = g.product("cleaned", "AMPure-purified single-cell library", MaterialKind.MOLECULAR, state=amplified_scene())
    waste = g.product("cleanup_waste", "material not recovered in post-PCR cleanup", MaterialKind.WASTE)
    g.react("post_pcr_cleanup", "Purify adapter-route library with 0.7x AMPure XP", StepKind.CLEANUP,
            ((library, InputRole.TEMPLATE), (ampure, InputRole.SELECTION)), (cleaned, waste))
    g.discard(waste, "not retained in library eluate; exact composition unrepresented")
    selection = g.input("size_selection", "concentrated AMPure double selection, 300–700 bp", MaterialKind.SELECTION)
    selected = g.product("library", final.name, MaterialKind.LIBRARY, state=amplified_scene())
    offsize = g.product("offsize", "products outside pooled library-size interval", MaterialKind.WASTE)
    g.react("recover_size_window", "Recover pooled 300–700 bp libraries by double AMPure selection", StepKind.SELECTION,
            ((cleaned, InputRole.TEMPLATE), (selection, InputRole.SELECTION)), (selected, offsize))
    g.discard(offsize, "outside recovered size interval")
    g.sequence(selected, SequencingPlan(library=final, primers=SEQ_PRIMERS,
                                       acquired_roles=REQUIRED_ROLES, read_lengths=READ_LENGTHS))
    g.declare_reviewed_scope(('sample', 'isolation', 'isolate', 'mboi', 'digest', 'biotin_datp', 'fill_polymerase', 'fill', 'contact_ligase', 'contact_ligation', 'post_ligation_wash', 'post_ligation_wash_selection', 'crosslink_reversal', 'reverse_crosslinks', 'streptavidin', 'capture', 'capture_wash', 'capture_wash_selection', 'alui', 'fragment', 'alui_wash', 'alui_wash_selection', 'end_prep', 'agarose', 'tail', 'tail_wash', 'tail_wash_selection', 'adapter', 'adapter_ligase', 'adapter_ligation', 'adapter_wash', 'adapter_wash_selection', 'pcr_forward', 'pcr_reverse', 'pcr_polymerase', 'library_pcr', 'ampure', 'post_pcr_cleanup', 'size_selection', 'recover_size_window', 'library'), CITATION)
    return g.require_publishable()

SOURCE=CITATION
CAVEAT="Source-owned adapter-ligation route. The CAA member is a source example;3nt inline cell identification tags are not an assigned index-read cycle mask. Original biotin/methylation and adapter backbone modifications are not copied into daughter DNA; physically incorporated PCR primers retain their own source terminal phosphorothioate bonds. Native Y forks and first-cycle PCR alignment remain bounded despite exact source component strands. Read4 includes the adapter dT/A junction before genomic sequence."
FINAL_LIBRARY=final_library()
FINAL_CAPTION=LIBRARY_CAPTION
SEQUENCING_INTRO=READOUT
SEQUENCING_ENDING=READOUT
SHOW_READ_CONTENT=False
READ_SPECS=(sp.read_spec('Read 1',READ_LENGTHS['Read 1'],'first3nt inline identification tag, dT junction, genomic contact arm','extract_then_include'),sp.read_spec('Read 2',READ_LENGTHS['Read 2'],'first3nt inline identification tag, dT junction, opposite contact arm','extract_then_include'),sp.read_spec('Index 1 (i7)','unassigned','Inline tags do not assign separate index acquisition','undetermined'),sp.read_spec('Index 2 (i5)','unassigned','Index acquisition unassigned in selected historical route','undetermined'))

def oligos():
 return ["CAA adapter F — 5′-"+il.TRUSEQ_READ1+"CAA*T-3′", "CAA adapter R — 5′-pTTGAGATCGGAAGAGCGGTTCAGCAGGAATGCCGAG-3′", "PCR forward — 5′-"+PCR_FORWARD[:-1]+"*"+PCR_FORWARD[-1]+"-3′", "PCR reverse — 5′-"+PCR_REVERSE[:-1]+"*"+PCR_REVERSE[-1]+"-3′"]

XT_INDEX5=feature('xt_i5','sample_index','unknown',note='source named i5; member/width unassigned')
XT_INDEX7=feature('xt_i7','sample_index','unknown',note='source named i7; member/width unassigned')
XT_LIBRARY=Construct([role('native XT P5/i5 end',molecular_feature=XT_INDEX5),*contact_product(),role('native XT P7/i7 end',molecular_feature=XT_INDEX7)],name='Stevens published Nextera XT alternative')
def xt_state(stage,*,indexed=True):
 sc=Scene();parts=[role(stage+' native XT/contact geometry; extent unassigned')]
 if indexed:parts.extend([role('installed i5 lineage',molecular_feature=XT_INDEX5),role('installed i7 lineage',molecular_feature=XT_INDEX7)])
 sc.strand('XT material',parts,label='native XT tagged/copy lineage; no CAA inline adapter assigned')
 sc.strand('originals',[role('possible original biotin contact molecules, offproducts and prior copies; recovery unresolved')],label='independent original/background lineages')
 return sc

def nextera_graph():
 g=ProtocolGraph(TITLE+' — published Nextera XT alternative',authority=CITATION)
 start=g.input('xt_captured','Alternative starting point: captured single-nucleus contacts from common capture stage; not an aliquot of the adapter route',MaterialKind.BIOLOGICAL,state=contact_scene(),primary=True)
 x=physical_wash(g,'xt_pre_wash',start,'Two BWT and two Tris washes')
 enzyme=g.input('xt_alui','AluI1U,37°C1h',MaterialKind.CATALYST);cut=g.product('xt_cut','AluI-cut captured contacts',MaterialKind.MOLECULAR,state=contact_scene())
 g.react('xt_fragment','AluI digest for XT route',StepKind.DIGESTION,((x,InputRole.TEMPLATE),(enzyme,InputRole.CATALYST)),(cut,))
 x=physical_wash(g,'xt_post_alui_wash',cut,'Two BWT/two Tris washes after AluI')
 t=g.input('xt_transposase','Nextera XT ATM1/100 dilution plus TD,55°C5min',MaterialKind.CATALYST);tagged=g.product('xt_tagged','native XT-tagmented precursors',MaterialKind.MOLECULAR,state=xt_state('tagmented',indexed=False))
 loaded=g.input('xt_loaded_adapter','ATM-loaded native transposon ends; exact source strand unassigned',MaterialKind.ADAPTER)
 g.react('xt_tag','Tagment on captured DNA',StepKind.TAGMENTATION,((x,InputRole.TEMPLATE),(loaded,InputRole.STRUCTURAL),(t,InputRole.CATALYST)),(tagged,))
 nt=g.input('xt_nt','NTbuffer RT5min',MaterialKind.CONDITION);neutral=g.product('xt_neutral','neutralized native tagmentation mixture',MaterialKind.MOLECULAR,state=tagged.state)
 g.react('xt_neutralize','Neutralize tagmentation',StepKind.CHEMICAL_CONVERSION,((tagged,InputRole.TEMPLATE),(nt,InputRole.CONDITION)),(neutral,))
 x=physical_wash(g,'xt_tag_wash',neutral,'Two BWT/two Tris washes; transfer beads to PCR tube')
 i5=g.input('xt_i5','appropriate native i5 primer; member/width unassigned',MaterialKind.OLIGO);i7=g.input('xt_i7','appropriate native i7 primer; member/width unassigned',MaterialKind.OLIGO);pol=g.input('xt_kapa','KAPA HiFi polymerase mix',MaterialKind.CATALYST)
 for n,cycles in ((1,9),(2,16)):
  amp=g.product('xt_amplified'+str(n),'XT round'+str(n)+' amplified mixture',MaterialKind.MOLECULAR,state=xt_state('round'+str(n)+' amplified'))
  g.react('xt_pcr'+str(n),'72°C3min;95°C30s;'+str(cycles)+'cycles95/55/72;72°C5min',StepKind.PCR,((x,InputRole.TEMPLATE),(i5,InputRole.STRUCTURAL),(i7,InputRole.STRUCTURAL),(pol,InputRole.CATALYST)),(amp,),geometry=PCRGeometry.unavailable(x.state,reason='Actual native XT-tagmented/copy precursor; selected indexed primer members and exact strand extents unassigned.'))
  x=physical_wash(g,'xt_cleanup'+str(n),amp,'Remove streptavidin beads/1×AMPure recovery' if n==1 else '1×AMPure recovery after second PCR')
 pool=g.input('xt_other_libraries','other independently prepared libraries for pooling',MaterialKind.MOLECULAR);pooled=g.product('xt_pool','pooled XT libraries',MaterialKind.MOLECULAR,state=x.state)
 g.react('xt_pooling','Pool independently indexed libraries',StepKind.POOL,((x,InputRole.TEMPLATE),(pool,InputRole.TEMPLATE)),(pooled,))
 selection=g.input('xt_size_selection','concentratedAMPure double0.5–0.7× selection300–700bp',MaterialKind.SELECTION);lib=g.product('xt_library',XT_LIBRARY.name,MaterialKind.LIBRARY,state=xt_state('size-selected amplified'));w=g.product('xt_size_waste','outside physically recovered fraction',MaterialKind.WASTE)
 g.react('xt_size','Select pooled size interval',StepKind.SELECTION,((pooled,InputRole.TEMPLATE),(selection,InputRole.SELECTION)),(lib,w));g.discard(w,'Physical double-selection fractions; no complete class purge assumed.')
 g.sequence(lib,SequencingPlan(library=XT_LIBRARY,primers=(),acquired_roles=('Read 1','Read 2'),read_lengths={'Read 1':75,'Read 2':75},limitation='Source MiSeqPE75; native i5/i7 index acquisition/width remains unassigned.'))
 g.declare_reviewed_scope(('xt_captured','xt_pre_wash','xt_pre_wash_selection','xt_alui','xt_fragment','xt_post_alui_wash','xt_post_alui_wash_selection','xt_transposase','xt_loaded_adapter','xt_tag','xt_nt','xt_neutralize','xt_tag_wash','xt_tag_wash_selection','xt_i5','xt_i7','xt_kapa','xt_pcr1','xt_cleanup1','xt_cleanup1_selection','xt_pcr2','xt_cleanup2','xt_cleanup2_selection','xt_other_libraries','xt_pooling','xt_size_selection','xt_size','xt_library'),CITATION)
 return g.require_publishable()
