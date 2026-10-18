"""INVADEseq: targeted bacterial 16S capture beside 10x 5' host RNA capture."""
from __future__ import annotations

import illumina as il
import seqprimers as sp
from chemdraw import Construct, Scene, Segment, feature, revcomp, oligo, complement_segments
from pcr import PCRGeometry, PCRPrimer
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind, Disposition, Terminal, Highlight)
from protocol_model import molecule,role
from adapter_ligation import AdapterLigationGeometry

TITLE = "INVADEseq — invasion–adhesion-directed expression sequencing"
NOTES = "01_invadeseq.html"
SOURCE = ('Defining protocol: <a href="https://doi.org/10.1038/s41596-023-00888-7">'
          'Galeano Niño et al., <i>Nature Protocols</i> (2023)</a>.')
SUMMARY = ("A bacterial 16S reverse-transcription primer is spiked into a 10x Chromium "
           "5′ single-cell reaction beside the ordinary poly(dT) primer. Host and bacterial "
           "cDNAs therefore inherit the same cell barcode; an aliquot is then amplified "
           "twice with nested 16S primers to make a separate bacterial library.")
CAVEAT=("This graph selects the defining2023 protocol, whose described baseline is5′v1 with explicit v2 substitutions and a Materials/Step85 kit-label inconsistency: Materials identifies PN1000263 as v2 reagent kit and PN1000002 as5′library kit, while Step85 calls PN1000002 a5′v2 kit. The model follows Materials role identities; intended Step85 kit combination remains unresolved. Cited CG000086RevJ and actual selected TT index/primer species remain unassigned; a laterRevK bead does not establish historical complete geometry. Source-known cell-barcode/UMI roles are retained without borrowing exact bead bases or widths. The2023 Equipment setup gives GEX paired75bp; Step84 allows user-chosen platform settings. The2022 application's26/10/10/90 layout is a separate historical run, not substituted here. Bacterial PE300 is positive; index acquisition/cycles remain unassigned. Enrichment retains a host-dominated mixed cDNA background, as the source explicitly states. Native first-cycle and complete TT product assembly remain model boundaries.")
RT_16S = "GGGTTGCGCTCGTTG"
ENRICH_FORWARD = "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTC"
ENRICH_REVERSE = "GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTTCACGRCACGAGCTGACGAC"
ENRICH_FORWARD_BINDING = ENRICH_FORWARD[len(il.P5):]
ENRICH_REVERSE_BINDING = ENRICH_REVERSE[len(il.TRUSEQ_READ2):]

CELL_BARCODE=feature("cell_barcode","cell_barcode","whitelist",whitelist="historical10x5′whitelist; member/width unassigned")
UMI=feature("umi","umi","random")
I5=feature("sample_index_i5","sample_index","unknown")
I7=feature("sample_index_i7","sample_index","unknown")

def identities():
    return [role("10x cell barcode; historical width unassigned",molecular_feature=CELL_BARCODE),role("molecular UMI; historical width unassigned",molecular_feature=UMI)]

def input_scene():
    sc=Scene();sc.strand("host",[role("host poly(A) RNA")],label="independent host RNA population")
    sc.strand("bacterial",[role("associated bacterial16S RNA")],label="independent bacterial RNA population")
    return sc

def source_da_scene():
    body=[role("mixed GEX cDNA fragment; sequence/length variable")]
    bottom=complement_segments(body)
    sc=Scene();sc.strand("top",[*body,Segment("top3′A","A")],label="source dA-tailed fragment")
    sc.anneal("bottom",[*bottom,Segment("bottom3′A","A")],to="top",pair=(bottom[-1].name,body[0].name),unpaired=tuple(x.name for x in bottom),label="source opposed strand;5′chemistry unassigned")
    sc.overhang("top","top3′A",polarity="3-prime",owner="substrate",text="source3′dA")
    sc.overhang("bottom","bottom3′A",polarity="3-prime",owner="substrate",text="source opposed3′dA")
    return sc

def copied_scene(switched=False,duplex=False):
    sc=Scene()
    for kind in ("host","bacterial"):
        parts=([role("source kit poly(dT) RT primer; exact historical species unassigned")] if kind=="host" else [Segment("16S1100R primer",RT_16S)])+[role(kind+" antisense first-strand cDNA")]
        if switched:parts += [role("copied TSO; exact historical end unassigned"),*reversed(identities()),role("copied bead PCR handle; full sequence unassigned")]
        sc.strand(kind,parts,label=kind+(" amplified strand class; complete duplex geometry unassigned" if duplex else " single-strand cDNA"))
    return sc

ENRICHED_AMPLICON=Construct([Segment("source enrichment forward",ENRICH_FORWARD),role("remaining historical bead arm"),*identities(),role("TSO/cDNA16S interval; native extent unassigned"),Segment("nested1061R footprint complement",revcomp(ENRICH_REVERSE_BINDING),placeholder=True,length_bp=20),Segment("reverse-primer-installed Read2 arm",revcomp(il.TRUSEQ_READ2))],name="source-primer-derived16S target component")

def enriched_population():
    sc=Scene();sc.strand("16S target",list(ENRICHED_AMPLICON),label="enriched16S target component; complete product assembly bounded")
    sc.strand("host background",[role("host-dominated cDNA background retained through enrichment")],label="independent host/background molecules")
    return sc

def endpoint(name):
    return Construct([role("TT indexed5′end; selected primer/width and complete junction unassigned",inferred=True,molecular_feature=I5),*identities(),role(name+" cDNA population; complete native adapter assembly unassigned",inferred=True),role("TT indexed opposite end; selected primer/width unassigned",inferred=True,molecular_feature=I7)],name=name)
GEX_LIBRARY=endpoint("INVADEseq GEX (host and off-target bacterial reads)")
BACTERIAL_LIBRARY=endpoint("INVADEseq bacterial-enriched mixed-background library")
PRIMERS=()
FINAL_LIBRARIES=(("Final host gene-expression library",GEX_LIBRARY,(),"Defining2023 GEX endpoint; native historical bead/TT geometry remains unassigned.","Published Equipment setup gives paired75bp; Step84 allows alternative user-selected settings."),("Final bacterial16S enrichment library",BACTERIAL_LIBRARY,(),"Bacterial enrichment within host-dominated background; no modality purification inferred.","Published MiSeqV3 PE300; repeated955–1215bp selection after15cycle TT PCR."))
SEQUENCING_UNAVAILABLE="Historical physical read/index primer geometry and selected TT index members/widths remain unassigned; no later commercial revision supplies them."
READ_SPECS={name:(sp.read_spec("Read 1",cycles,"10x identifiers and cDNA; native priming/identifier widths unassigned","extract_then_include"),sp.read_spec("Read 2",cycles,"cDNA from GEX fragments" if cycles==75 else "nested-primer-side cDNA from enriched mixed library","include"),sp.read_spec("Index 1 (i7)","unassigned","physical TT sample index; acquisition/cycles unassigned","undetermined"),sp.read_spec("Index 2 (i5)","unassigned","physical TT sample index; acquisition/cycles unassigned","undetermined")) for name,cycles in ((FINAL_LIBRARIES[0][0],75),(FINAL_LIBRARIES[1][0],300))}
READ_SPECS_NOTE="Positive2023 read-length statements retained separately from unassigned index acquisition. The2022 application run is not transferred to this protocol."

def oligos():
    return tuple(oligo(name,[Segment(name,seq,placeholder="R" in seq)]) for name,seq in (("16S1100R",RT_16S),("16SEnrichForward",ENRICH_FORWARD),("1061REnrichNestedReverse",ENRICH_REVERSE)))

def clean(g,key,template,recipe,kind=StepKind.CLEANUP,library=False):
    rule=g.input(key+"_rule",recipe,MaterialKind.SELECTION)
    kept=g.product(key+"_clean","physically recovered "+template.label,MaterialKind.LIBRARY if library else MaterialKind.MOLECULAR,state=template.state)
    waste=g.product(key+"_waste","source removed bead/size/wash fraction",MaterialKind.WASTE)
    g.react(key,"Recover source-selected fraction: "+recipe,kind,((template,InputRole.TEMPLATE),(rule,InputRole.SELECTION)),(kept,waste),changes_molecule=False);g.discard(waste,"Physical source selection/wash transfer, not nonamplification or computational QC")
    return kept

def protocol_graph():
    g=ProtocolGraph(TITLE,authority=SOURCE)
    cells=g.input("cells","host-cell partitions with associated bacterial RNA",MaterialKind.BIOLOGICAL,state=input_scene(),primary=True)
    gem_rule=g.input("gem_rule","source ChromiumGEM generation",MaterialKind.CONDITION)
    gem=g.product("gem","GEM mixed host/bacterial RNA populations",MaterialKind.MOLECULAR,state=input_scene())
    other=g.product("other_partitions","empty/colliding/other GEM outcomes",MaterialKind.MOLECULAR)
    g.react("partition","Partition host cells with associated bacterialRNA",StepKind.PARTITION,((cells,InputRole.TEMPLATE),(gem_rule,InputRole.CONDITION)),(gem,other),changes_molecule=False)
    g.mark_terminal(Terminal(other.id,Disposition.UNRESOLVED,"Partition interpretability does not establish physical removal of empty/multiplet material."))
    poly=g.input("poly_dt","source poly(dT)RT primer PN2000007; historical species unassigned",MaterialKind.OLIGO)
    targeted=g.input("primer_16s","source16S1100R",MaterialKind.OLIGO,state=RT_16S)
    rt=g.input("rt","RTenzymeC PN2000085;53°C45min",MaterialKind.CATALYST)
    first=g.product("first","host and bacterial actual first-strand cDNA",MaterialKind.MOLECULAR,state=copied_scene())
    g.react("reverse_transcribe","Prime host poly(A) and bacterial16S RNA",StepKind.REVERSE_TRANSCRIPTION,((gem,InputRole.TEMPLATE),(poly,InputRole.STRUCTURAL),(targeted,InputRole.STRUCTURAL),(rt,InputRole.CATALYST)),(first,))
    bead=g.input("bead_tso","barcoded/UMI10x TSO; exact citedRevJ geometry unassigned",MaterialKind.OLIGO)
    switched=g.product("switched","template-switched single-strand cDNA populations",MaterialKind.MOLECULAR,state=copied_scene(True),inferred=True)
    g.react("template_switch","Copy barcode/UMI-bearing beadTSO during RT",StepKind.TEMPLATE_SWITCHING,((first,InputRole.TEMPLATE),(bead,InputRole.STRUCTURAL),(rt,InputRole.CATALYST)),(switched,))
    recovered=clean(g,"post_gem",switched,"source emulsion recovery/silane purification;35µl eluate")
    f=g.input("cdna_f","cDNA primer mix PN2000089 first binding role; species unassigned",MaterialKind.OLIGO)
    r=g.input("cdna_r","same commercial primer mix opposite binding role",MaterialKind.OLIGO)
    pol=g.input("cdna_pol","amplification master mix PN2000103;16cycles98/63/72°C",MaterialKind.CATALYST)
    amp=g.product("amp","amplified mixed host/bacterial cDNA",MaterialKind.MOLECULAR,state=copied_scene(True,True))
    g.react("cdna_pcr","Amplify source mixed cDNA pool16cycles",StepKind.PCR,((recovered,InputRole.TEMPLATE),(f,InputRole.STRUCTURAL),(r,InputRole.STRUCTURAL),(pol,InputRole.CATALYST)),(amp,),geometry=PCRGeometry.unavailable(recovered.state,reason="Actual incoming material is switched single-strand cDNA, not an enriched/final duplex. Historical commercial primer-pair geometry remains unassigned."))
    amp=clean(g,"cdna_cleanup",amp,"0.6×SPRIselect;45µlEB")
    host=g.product("host_aliquot","GEX aliquot containing host and bacterial cDNA",MaterialKind.MOLECULAR,state=amp.state)
    bacterial=g.product("bacterial_aliquot","2–10µl enrichment aliquot with same mixed cDNA",MaterialKind.MOLECULAR,state=amp.state)
    g.react("split","Sample GEX and enrichment inputs from same mixed pool",StepKind.SPLIT,((amp,InputRole.TEMPLATE),),(host,bacterial),changes_molecule=False)
    frag=g.input("fragmentation","fragmentation mix PN2000090/91;32°C5min then65°C30min endrepair/A-tail",MaterialKind.CATALYST)
    repaired=g.product("gex_fragments","fragmented/end-repaired/dA-tailed GEX cDNA population",MaterialKind.MOLECULAR,state=source_da_scene())
    g.react("gex_fragment","Fragment and endprepare50ngGEX input",StepKind.FRAGMENTATION,((host,InputRole.TEMPLATE),(frag,InputRole.CATALYST)),(repaired,))
    repaired=clean(g,"gex_fragment_cleanup",repaired,"source30µl then10µl double-sidedSPRI;50µlEB")
    adapter=g.input("gex_adapter","adapter oligos PN2000094; native end geometry unassigned",MaterialKind.ADAPTER)
    ligase=g.input("gex_ligase","DNA ligase PN220110/ligationbuffer PN2000092;20°C15min",MaterialKind.CATALYST)
    ligated=g.product("gex_ligated","actual pre-index adapter-ligated GEX fragments",MaterialKind.MOLECULAR,state=molecule("pre-index GEX",(role("adapter-ligated mixed cDNA; native species/junction geometry unassigned"),)))
    g.react("gex_ligate","Ligate source-named GEX adapters",StepKind.LIGATION,((repaired,InputRole.TEMPLATE),(adapter,InputRole.STRUCTURAL),(ligase,InputRole.CATALYST)),(ligated,),geometry=AdapterLigationGeometry.unavailable("Source adapter part number is positive; complete source-revision strand/end assembly remains unassigned."))
    ligated=clean(g,"gex_ligation_cleanup",ligated,"0.8×SPRI;30µlEB")
    ef=g.input("enrich_f","published16SEnrichForward",MaterialKind.OLIGO,state=ENRICH_FORWARD)
    er=g.input("enrich_r","published1061R nestedreverse",MaterialKind.OLIGO,state=ENRICH_REVERSE)
    ep=g.input("enrich_pol","amplification mix PN220125;PN220067 additive or explicit v2 MasterAmp substitution",MaterialKind.CATALYST)
    current=bacterial
    for number,cycles in ((1,35),(2,20)):
        enriched=g.product(f"enrichment_{number}",f"round{number}16S-enriched host-background cDNA",MaterialKind.MOLECULAR,state=enriched_population(),inferred=True)
        g.react(f"enrich_{number}",f"Enrich16S targets{cycles}cycles98/67/72°C",StepKind.PCR,((current,InputRole.TEMPLATE),(ef,InputRole.STRUCTURAL),(er,InputRole.STRUCTURAL),(ep,InputRole.CATALYST)),(enriched,),geometry=PCRGeometry.unavailable(current.state,reason="The actual input remains mixed source cDNA. Source primer-derived target component does not establish complete native first-cycle assembly or host-background removal."))
        current=clean(g,f"enrich_cleanup_{number}",enriched,"source single-sided0.8×SPRI;35µlEB")
        if number==1:current=clean(g,"enrich_size",current,"BluePippin955–1215bp",StepKind.SELECTION)
    for key,template,final,cycles in (("gex",ligated,GEX_LIBRARY,"14–16"),("bacterial",current,BACTERIAL_LIBRARY,"15")):
        fi=g.input(key+"_index_f","DualIndexTTSetA PN1000215 first primer; species/width unassigned",MaterialKind.OLIGO)
        ri=g.input(key+"_index_r","DualIndexTTSetA opposite primer; species/width unassigned",MaterialKind.OLIGO)
        ip=g.input(key+"_index_pol","amplification master mix PN2000103;"+cycles+"cycles",MaterialKind.CATALYST)
        indexed=g.product(key+"_indexed","source-indexed "+key+" product",MaterialKind.MOLECULAR,state=final,inferred=True)
        g.react(key+"_index_pcr","Index "+key+" branch "+cycles+"cycles",StepKind.PCR,((template,InputRole.TEMPLATE),(fi,InputRole.STRUCTURAL),(ri,InputRole.STRUCTURAL),(ip,InputRole.CATALYST)),(indexed,),geometry=PCRGeometry.unavailable(template.state,reason="Actual incoming pre-index branch retained; selected historical TT primer pair and full endpoint assembly unassigned."))
        indexed=clean(g,key+"_final_cleanup",indexed,"0.6×/0.8×double-sidedSPRI" if key=="gex" else "single-sided0.8×SPRI",library=key=="gex")
        if key=="bacterial":indexed=clean(g,"bacterial_final_size",indexed,"repeatBluePippin955–1215bp",StepKind.SELECTION,library=True)
        g.sequence(indexed,SequencingPlan(library=final,acquisition="unspecified",limitation=SEQUENCING_UNAVAILABLE))
    g.annotate(Highlight("reverse_transcribe","RNA-template fate unassigned","Native RNA persists as RT template; its later physical recovery/degradation is not assigned or treated as discarded."))
    g.declare_reviewed_scope(REVIEWED_TARGETS,SOURCE)
    return g.require_publishable()

REVIEWED_TARGETS=("cells","gem_rule","partition","poly_dt","primer_16s","rt","reverse_transcribe","bead_tso","template_switch","post_gem_rule","post_gem","cdna_f","cdna_r","cdna_pol","cdna_pcr","cdna_cleanup_rule","cdna_cleanup","split","fragmentation","gex_fragment","gex_fragment_cleanup_rule","gex_fragment_cleanup","gex_adapter","gex_ligase","gex_ligate","gex_ligation_cleanup_rule","gex_ligation_cleanup","enrich_f","enrich_r","enrich_pol","enrich_1","enrich_cleanup_1_rule","enrich_cleanup_1","enrich_size_rule","enrich_size","enrich_2","enrich_cleanup_2_rule","enrich_cleanup_2","gex_index_f","gex_index_r","gex_index_pol","gex_index_pcr","gex_final_cleanup_rule","gex_final_cleanup","gex_final_cleanup_clean","bacterial_index_f","bacterial_index_r","bacterial_index_pol","bacterial_index_pcr","bacterial_final_cleanup_rule","bacterial_final_cleanup","bacterial_final_size_rule","bacterial_final_size","bacterial_final_size_clean")
