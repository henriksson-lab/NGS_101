"""Nagano et al. 2013 single-cell Hi-C chemistry."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
from restriction import BGLII
import illumina as il
import seqprimers as sp
from chemdraw import Construct, Scene, Segment, bases_pair, revcomp, complement_segments, feature
from adapter_ligation import AdapterLigationGeometry
from endprep import dA_tailed_scene, EndPreparedDuplex
from pcr import PCRGeometry
from protocol_graph import (Disposition, InputRole, MaterialKind, MolecularRelation, ProtocolGraph,
                            SequencingPlan, StepKind, Terminal)
from protocol_model import molecule, role

from restriction import Digest
HISTORICAL_PE_READ2=il.PE_READ2_PRIMER
DIGEST=Digest(BGLII,"X"*16,"X"*16)
FILLED=DIGEST.fill_in(biotin_base="A")
TITLE = "Nagano single-cell Hi-C — ligate first, isolate nuclei second"
NOTES = "01_nagano-single-cell-hi-c.html"
CITATION = 'Defining source: <a href="https://doi.org/10.1038/nature12593">Nagano et al., <i>Nature</i> (2013)</a>; <a href="https://www.nature.com/nature/journal/v502/n7469/extref/nature12593-s1.pdf">original supplementary methods</a>; detailed protocol: <a href="https://doi.org/10.1038/nprot.2015.127">Nagano et al. (2015)</a>.'
SUMMARY = "BglII contact ends are biotin-filled and ligated in fixed nuclei; individual nuclei are isolated only after ligation."
STEPS = (("Crosslink and digest", "Digest fixed chromatin in intact nuclei with BglII."),
         ("Fill and ligate", "Fill the 5′ overhang with biotin-14-dATP, then proximity-ligate."),
         ("Isolate one nucleus", "Pick individual ligated nuclei under a microscope with mouth-controlled glass capillaries before reversal of crosslinks."),
         ("Capture and amplify", "Capture biotin-bearing junctions, cut with AluI, add a 3′ dA, ligate the printed Y-adapter, and PCR-amplify."))
READOUT = "Each paired read starts with a 3-bp inline identification tag, then reads one genomic fragment flanking the captured contact junction."
JUNCTION_CAPTION = "BglII-compatible ends form the contact boundary; the experimental junction is selected through the biotin-bearing fill-in. ** marks ligation."
LIBRARY_CAPTION = "PCR-completed paired-end library with the paper's 3-bp identification tag at each read end. The CAA member of the twelve-adapter set is shown."
JUNCTIONS = (("locus A; extent unassigned", "contact junction", "proximity ligation"),)
SEQ_PRIMERS = (
    sp.custom("Read 1", "historical Illumina Read 1", il.TRUSEQ_READ1,
              "Nagano et al. 2013 Supplementary Methods"),
    sp.custom("Read 2", "historical paired-end Read 2", HISTORICAL_PE_READ2,
              "Nagano et al. 2013 Supplementary Methods"),
)
REQUIRED_ROLES = ("Read 1", "Read 2")
READ_LENGTHS = {"Read 1": 40, "Read 2": 40}
PCR_FORWARD = il.P5 + il.TRUSEQ_READ1[4:]
PCR_REVERSE = il.P7 + HISTORICAL_PE_READ2


def pcr_primer_scene(name, sequence):
    scene = Scene()
    scene.strand(name, [Segment(name, sequence)], mod3="terminal phosphorothioate")
    return scene

# Reuse only the endpoint: original 2013 SI pp3–5 independently establishes these
# CAA-tagged arms and both historical PCR primers; no later reaction settings are inherited.
def final_library():
    lib=Construct([Segment("P5",il.P5,"p5"),Segment("historical Read1 arm",il.TRUSEQ_READ1[4:],"r1",note="Forward daughter strand retains the source PCR primer terminal C*T phosphorothioate bond. Complementary copies do not inherit this backbone modification."),
        Segment("inline identification tag","CAA","cbc",feature=INLINE_A),Segment("source adapter dT junction","T"),
        *contact_product(),Segment("dA junction","A"),Segment("opposite inline tag","TTG","cbc",feature=INLINE_B),
        Segment("historical Read2 site",revcomp(HISTORICAL_PE_READ2),"r2",note="Opposite daughter strand retains the reverse PCR primer terminal phosphorothioate bond; the top-strand sequence copy does not copy that backbone chemistry."),Segment("P7 reverse complement",il.P7_RC,"p7")],name=TITLE+" source-derived adapter-route library")
    problems=sp.verify(lib,SEQ_PRIMERS,required_roles=REQUIRED_ROLES)
    if problems:raise ValueError("; ".join(problems))
    return lib





def adapter_geometry() -> AdapterLigationGeometry:
    """The CAA member of the source-written adapter, before ligation."""
    if not bases_pair("T", ALUI_TAILED.three_prime_overhang):
        raise ValueError("published adapter 3′ dT cannot pair with the prepared contact end")
    sc=source_adapter_scene()
    return AdapterLigationGeometry.known(
        "cohesive", dA_tailed_scene(ALUI_TAILED, label="AluI-cut contact fragment"), sc,
        statement="2013 Supplementary Information pp. 2–5 prints dA tailing and both CAA adapter strands, including the 3′ T overhang, 5′ phosphate and terminal phosphorothioate.")


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

LIBRARY_ADAPTER_GEOMETRY = adapter_geometry()
LIBRARY_ADAPTER_CAPTION = ("Source-printed CAA Y-adapter before ligation: the 17-bp core pairs, "
                           "the sequencing arms do not pair, and the terminal 3′ dT meets the fragment's 3′ dA.")


def protocol_graph() -> ProtocolGraph:
    """Bulk contact ligation precedes the single-nucleus sampling boundary."""
    g = ProtocolGraph(TITLE, authority=CITATION)
    sample = g.input("sample", "fixed intact nuclei", MaterialKind.BIOLOGICAL,
                     state=molecule("fixed chromatin", (role("chromosomal DNA"),)),
                     primary=True)
    enzyme = g.input("bgl2", "BglII", MaterialKind.CATALYST)
    cut = g.product("cut", "BglII-cut chromatin", MaterialKind.MOLECULAR,
                    state=restriction_scene(), chemical_state=DIGEST)
    g.react("digest", "Digest fixed chromatin with BglII", StepKind.DIGESTION,
            ((sample, InputRole.TEMPLATE), (enzyme, InputRole.CATALYST)), (cut,))
    biotin = g.input("biotin_datp", "biotin-14-dATP", MaterialKind.MOLECULAR)
    polymerase = g.input("fill_polymerase", "DNA Polymerase I Large (Klenow) Fragment",
                         MaterialKind.CATALYST)
    marked = g.product("marked", "biotin-filled BglII ends", MaterialKind.MOLECULAR,
                       state=filled_scene(), chemical_state=FILLED)
    g.react("fill", "Fill BglII overhangs with biotin-dATP", StepKind.EXTENSION,
            ((cut, InputRole.TEMPLATE), (biotin, InputRole.STRUCTURAL),
             (polymerase, InputRole.CATALYST)), (marked,),
            relations=(MolecularRelation.RESTRICTION_FILL,))
    ligase = g.input("contact_ligase", "T4 DNA ligase",
                     MaterialKind.CATALYST)
    contacts = g.product("contacts", "bulk proximity-ligated contacts",
                         MaterialKind.MOLECULAR, state=contact_scene())
    g.react("contact_ligation", "Ligate contacts before isolating individual nuclei",
            StepKind.PROXIMITY_LIGATION,
            ((marked, InputRole.TEMPLATE), (ligase, InputRole.CATALYST)),
            (contacts,), relations=(MolecularRelation.CONTACT_JUNCTION,))
    contacts=physical_wash(g,"post_ligation_wash",contacts,"PBS wash after bulk ligation")
    isolation = g.input("isolation", "microscope picking with mouth-controlled glass capillaries",
                        MaterialKind.CONDITION)
    nucleus = g.product("nucleus", "one ligated nucleus", MaterialKind.MOLECULAR,
                        state=contacts.state)
    excluded = g.product("excluded", "other ligated nuclei", MaterialKind.MOLECULAR)
    g.react("isolate", "Isolate one nucleus after contact ligation", StepKind.SAMPLING,
            ((contacts, InputRole.TEMPLATE), (isolation, InputRole.CONDITION)),
            (nucleus, excluded), changes_molecule=False)
    g.mark_terminal(Terminal(excluded.id, Disposition.UNRESOLVED,
        "Other nuclei lie outside this picked-nucleus lineage; their experimental fate is not represented."))
    reversal = g.input("reversal", "crosslink reversal at 65 °C", MaterialKind.CONDITION)
    fragments = g.product("fragments", "crosslink-reversed contact DNA",
                          MaterialKind.MOLECULAR, state=nucleus.state)
    g.react("reverse_crosslinks", "Reverse crosslinks after single-nucleus isolation",
            StepKind.DEPROTEINIZATION,
            ((nucleus, InputRole.TEMPLATE), (reversal, InputRole.CONDITION)), (fragments,))
    beads = g.input("streptavidin", "streptavidin beads", MaterialKind.SELECTION)
    captured = g.product("captured", "biotin-selected contacts", MaterialKind.MOLECULAR,
                         state=fragments.state)
    unbound = g.product("unbound", "fragments without a captured biotin contact",
                        MaterialKind.WASTE)
    g.react("capture", "Enrich biotin-bearing contact junctions",
            StepKind.AFFINITY_CAPTURE,
            ((fragments, InputRole.TEMPLATE), (beads, InputRole.SELECTION)),
            (captured, unbound))
    g.discard(unbound, "physically removed binding supernatant; no complete class-specific depletion assumed")
    captured=physical_wash(g,"capture_wash",captured,"Source BW/BWT and reaction-buffer washes after capture")
    alui = g.input("alui", "AluI", MaterialKind.CATALYST)
    cut_contacts = g.product("alui_cut", "AluI-cut captured contacts",
                             MaterialKind.MOLECULAR, state=contact_scene())
    g.react("alui_digest", "Cut captured contacts with AluI", StepKind.DIGESTION,
            ((captured, InputRole.TEMPLATE), (alui, InputRole.CATALYST)),
            (cut_contacts,))
    cut_contacts=physical_wash(g,"alui_wash",cut_contacts,"Source washes after AluI digestion")
    datp = g.input("datp", "dATP", MaterialKind.MOLECULAR)
    klenow = g.input("klenow_exo", "Klenow fragment (3′→5′ exo−)", MaterialKind.CATALYST)
    tailed = g.product("tailed", "3′ dA-tailed contact fragments",
                       MaterialKind.MOLECULAR, state=ALUI_TAILED)
    g.react("da_tail", "Add a 3′ dA to AluI-cut contact fragments", StepKind.TAILING,
            ((cut_contacts, InputRole.TEMPLATE), (datp, InputRole.STRUCTURAL),
             (klenow, InputRole.CATALYST)), (tailed,))
    tailed=physical_wash(g,"tail_wash",tailed,"Source washes after dA tailing")
    adapter = g.input("adapter", "3-bp-tagged historical Illumina Y-adapter",
                      MaterialKind.ADAPTER,state=source_adapter_scene())
    adapter_ligase = g.input("adapter_ligase", "T4 DNA ligase for adapter ligation",
                             MaterialKind.CATALYST)
    ligated = g.product("ligated", "adapter-ligated contact fragments",
                        MaterialKind.MOLECULAR, state=native_ligated_scene())
    g.react("adapter_ligation", "Install historical paired-end adapters",
            StepKind.LIGATION,
            ((tailed, InputRole.TEMPLATE), (adapter, InputRole.STRUCTURAL),
             (adapter_ligase, InputRole.CATALYST)), (ligated,),
            geometry=LIBRARY_ADAPTER_GEOMETRY)
    ligated=physical_wash(g,"adapter_wash",ligated,"Source washes after adapter ligation and into Pfx buffer")
    forward = g.input("pcr_forward", "source-printed library PCR forward primer", MaterialKind.OLIGO,
                      state=pcr_primer_scene("PCR forward", PCR_FORWARD))
    reverse = g.input("pcr_reverse", "source-printed library PCR reverse primer", MaterialKind.OLIGO,
                      state=pcr_primer_scene("PCR reverse", PCR_REVERSE))
    pcr = g.input("pcr_polymerase", "Platinum Pfx DNA polymerase", MaterialKind.CATALYST)
    final = final_library()
    library = g.product("amplified", final.name, MaterialKind.MOLECULAR, state=amplified_scene())
    g.react("library_pcr", "Amplify for 25 cycles with Platinum Pfx",
            StepKind.PCR,
            ((ligated, InputRole.TEMPLATE), (forward, InputRole.STRUCTURAL),
             (reverse, InputRole.STRUCTURAL), (pcr, InputRole.CATALYST)),
            (library,), geometry=PCRGeometry.unavailable(
                ligated.state, reason="Actual heterogeneous native Y-adapter/contact strands and source primers; complete first-cycle alignment bounded."))
    ampure = g.input("ampure", "AMPure XP magnetic beads after removal of streptavidin beads", MaterialKind.SELECTION)
    cleaned = g.product("cleaned", "AMPure-purified amplified contact library", MaterialKind.MOLECULAR, state=amplified_scene())
    cleanup_waste = g.product("cleanup_waste", "material not recovered from post-PCR AMPure purification", MaterialKind.WASTE)
    g.react("post_pcr_cleanup", "Remove streptavidin beads and purify PCR library with AMPure XP",
            StepKind.CLEANUP, ((library, InputRole.TEMPLATE), (ampure, InputRole.SELECTION)),
            (cleaned, cleanup_waste))
    g.discard(cleanup_waste, "not recovered in the library eluate; exact molecular composition not assigned")
    gel = g.input("size_gel", "2% agarose gel, 300–700 bp recovery window", MaterialKind.SELECTION)
    selected = g.product("library", final.name, MaterialKind.LIBRARY, state=amplified_scene())
    offsize = g.product("offsize", "products outside the recovered gel interval", MaterialKind.WASTE)
    g.react("size_selection", "Recover 300–700 bp libraries from 2% agarose gel",
            StepKind.SELECTION, ((cleaned, InputRole.TEMPLATE), (gel, InputRole.SELECTION)),
            (selected, offsize))
    g.discard(offsize, "outside the recovered 300–700 bp gel interval")
    g.sequence(selected, SequencingPlan(library=final, primers=SEQ_PRIMERS,
                                       acquired_roles=REQUIRED_ROLES, read_lengths=READ_LENGTHS))
    g.declare_reviewed_scope(('sample', 'isolation', 'isolate', 'bgl2', 'digest', 'biotin_datp', 'fill_polymerase', 'fill', 'contact_ligase', 'contact_ligation', 'post_ligation_wash', 'post_ligation_wash_selection', 'reversal', 'reverse_crosslinks', 'streptavidin', 'capture', 'capture_wash', 'capture_wash_selection', 'alui', 'alui_digest', 'alui_wash', 'alui_wash_selection', 'datp', 'klenow_exo', 'da_tail', 'tail_wash', 'tail_wash_selection', 'adapter', 'adapter_ligase', 'adapter_ligation', 'adapter_wash', 'adapter_wash_selection', 'pcr_forward', 'pcr_reverse', 'pcr_polymerase', 'library_pcr', 'ampure', 'post_pcr_cleanup', 'size_gel', 'size_selection', 'library'), CITATION)
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
