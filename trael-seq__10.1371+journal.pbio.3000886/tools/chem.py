"""Kara2021 TrAEL-seq: native modified hairpins, USER precursor and source-qualified kit ends."""
from chemdraw import Construct,Scene,Segment,feature,revcomp
from protocol_model import role
from protocol_graph import ProtocolGraph,MaterialKind,InputRole,StepKind,SequencingPlan,Terminal,Disposition
from adapter_ligation import AdapterLigationGeometry
from endprep import EndPreparedDuplex,dA_tailed_scene
from pcr import PCRGeometry
from topology import HairpinTopology
import seqprimers as sp
TITLE='TrAEL-seq — native DNA3′-end sequencing'
NOTES='01_trael-seq.html'
SOURCE='Kara et al., PLOS Biology2021, <a href="https://doi.org/10.1371/journal.pbio.3000886">doi:10.1371/journal.pbio.3000886</a>; published Methods199–225/S2File; NEB E7645 manual v6.1_5/20.'
SUMMARY='TdT adds a variable1–4rA tail; activated biotin hairpin1 captures the native3′OH. Hairpin extension, fragmentation, streptavidin capture, on-bead end preparation/hairpin2 ligation and USER elution precede sample-index PCR.'
CAVEAT='The paper calls UltraII EndPrep blunting/polishing; the same named kit manual v6.1_5/20 explicitly describes end repair,5′phosphorylation and dA-tailing. The historical guide revision/lot used in each experiment is not identified, so the TA-compatible junction is a revision-qualified interpretation, not a publication-error claim. The dA drawing remains inferred. Exact library index set/member/width and full native PCR/USER end geometry remain unassigned; E7335S is only a source example. Chemical rA range1–4 and computational removal of up to3T are different statements.'
UMI_FEATURE=feature('umi','umi','random')
INDEX_FEATURE=feature('trael_sample','sample_index','unknown',note='Native NEBNext index member/set/width unassigned; E7335S only source example.')
ADAPTER2 = "GATCGGAAGAGCACACGTCTGAACTCCAGTCUUUUGACTGGAGTTCAGACGTGTGCTCTTCCGATCT"
ADAPTER1 = ("NNNNNNNNAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTU"
            "GCGCAGGCCATTGGCCTGCGCUACACTCTTTCCCTACACGACGCT")
ADAPTER1_TOPOLOGY = HairpinTopology(
    "TrAEL adapter 1", ADAPTER1, first_stem_start=18,
    second_stem_start=64, stem_length=23, colour="var(--c-r1)")
ADAPTER2_TOPOLOGY = HairpinTopology(
    "TrAEL adapter 2", ADAPTER2, first_stem_start=0,
    second_stem_start=35, stem_length=31, colour="var(--c-r2)")
def oligos():
    return ["<b>TrAEL adapter 1</b> — 5′-[rApp]NNNNNNNNAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTUGCGCAGGCCATTGGCC[BtndT]GCGCUACACTCTTTCCCTACACGACGCT-3′",f"<b>Adapter 2</b> — 5′-[Phos]{ADAPTER2[:-1]}*{ADAPTER2[-1]}-3′"]

FINAL_LIBRARY=Construct([role('native adapter1/NEBNext primer-derived end'),Segment('copied UMI','N'*8,'umi',placeholder=True,feature=UMI_FEATURE),role('T tract copied from variable1–4rA tail; length unassigned'),role('genomic complement of captured3′ end'),role('native adapter2/indexed end; exact USER/PCR boundary unassigned',molecular_feature=INDEX_FEATURE)],name='TrAEL native source-limited endpoint')
SEQ_PRIMERS=()
FINAL_CAPTION='Source-directed read strand: UMI8, variable T tract, then genomic complement. Outer native ends and index width remain unassigned.'
SEQUENCING_INTRO='MiSeq50-base single-read or NextSeq500 high-output75-base single-end. Read1 reports UMI8/T tract/antisense genomic sequence; index acquisition allocation remains unassigned.'
SEQUENCING_ENDING=SEQUENCING_INTRO
SHOW_READ_CONTENT=False
READ_SPECS=(sp.read_spec('Read 1','50 MiSeq;75 NextSeq500','UMI8 followed by variable T tract and genomic complement; computationally extract UMI and trim up to3T','extract_then_include'),sp.read_spec('Index 1 (i7)','unassigned','Native sample index acquisition/extent unassigned','undetermined'),sp.read_spec('Index 2 (i5)','unassigned','Index acquisition unassigned; no molecular i5 window invented','undetermined'))

def adapter1_parts():
 biotin=ADAPTER1.index('UGCGCAGGCCATTGGCC')+len('UGCGCAGGCCATTGGCC')
 if ADAPTER1[biotin]!='T':raise ValueError('Source biotin-dT position must retain its T base')
 return [Segment('UMI','N'*8,'umi',placeholder=True,feature=UMI_FEATURE),Segment('adapter1 before biotin-dT',ADAPTER1[8:biotin],note='U is deoxyuridine'),Segment('internal biotin-dT',ADAPTER1[biotin],note='source [BtndT]'),Segment('adapter1 after biotin-dT',ADAPTER1[biotin+1:],note='U is deoxyuridine')]

def adapter1_scene():
 sc=Scene();sc.strand('adapter1',adapter1_parts(),label='source hairpin1; two dU flank biotin-bearing loop',mod5='rApp');sc.ligation_end('adapter1','UMI',kind='single-strand',polarity='5-prime',owner='adapter',text='activated5′App');return sc

def adapter2_scene():
 sc=Scene();sc.strand('adapter2',[Segment('stem1',ADAPTER2[:31]),Segment('dU4 loop',ADAPTER2[31:35],note='four deoxyuridines'),Segment('stem2',ADAPTER2[35:66]),Segment('terminal3′T',ADAPTER2[66:],note='preceded by source phosphorothioate bond *')],label='source hairpin2,31bp stem and terminal*T',mod5='p');sc.ligation_end('adapter2','terminal3′T',kind='cohesive',polarity='3-prime',owner='adapter',text='single3′T');return sc

def end_state(*,tailed=False,tagged=False):
 sc=Scene();parts=[role('genomic DNA upstream of native3′OH')]
 if tailed:parts.append(role('TdT-added1–4rA tail; exact count unassigned'))
 if tagged:parts.extend(adapter1_parts())
 sc.strand('captured strand',parts,label='native capture strand',mod3='OH' if not tagged else '')
 if tagged:sc.junction('captured strand','TdT-added1–4rA tail; exact count unassigned','UMI','RNA-ligase capture')
 return sc

def extended_state():
 sc=Scene()
 parts=[role('genomic capture strand'),role('variable1–4rA tail'),*adapter1_parts(),Segment('copied adapter1 prefix',revcomp(ADAPTER1[8:18])),Segment('copied UMI','N'*8,'umi',placeholder=True,feature=UMI_FEATURE),role('T copy of variable rA tail'),role('new genomic complementary copy')]
 sc.strand('hairpin-connected DNA',parts,label='one covalent original–hairpin1–Bst-copy polymer; full duplex fold not reconstructed')
 return sc

def endprep_state():
 product=EndPreparedDuplex(Construct([role('biotin-hairpin1-tagged genomic fragment'),Segment('retained UMI lineage','N'*8,'umi',placeholder=True,feature=UMI_FEATURE)]),five_prime_phosphate=False,three_prime_overhang='A',overhang_inferred=True)
 sc=dA_tailed_scene(product,label='INFERRED—TA-compatible kit EndPrep; historical revision unassigned')
 sc.strand('carrier',[role('independent salmon-sperm carrier DNA; end geometry/partition unassigned')],label='carrier molecule population, not genomic capture strand')
 return sc

def dual_state():
 sc=Scene()
 sc.strand('capture construct',[role('hairpin1-connected original/copy body with hairpin2 at opposite end; full native nick/seal/closure topology unassigned'),Segment('UMI lineage','N'*8,'umi',placeholder=True,feature=UMI_FEATURE)],label='one two-hairpin molecular object; UMI region shown, full strand order unassigned')
 sc.strand('carrier',[role('independent carrier/ligation background; recovery unassigned')],label='carrier DNA population')
 return sc

def opened_state():
 sc=Scene();sc.strand('released target',[role('USER-opened adapter1 end; exact retained boundary unassigned'),Segment('UMI lineage','N'*8,'umi',placeholder=True,feature=UMI_FEATURE),role('variable rA/T tract and paired captured genomic/copy body'),role('USER-opened adapter2 end; native footprint unassigned')],label='actual USER-eluted PCR substrate; no completed P5/P7 shell')
 sc.strand('possible carrier',[role('possible co-eluted carrier/background DNA; partition/integrity unassigned')],label='independent possible background lineage')
 return sc

def post_pcr_state():
 sc=Scene.duplex(list(FINAL_LIBRARY),label='expected amplified capture library',unpaired=tuple(p.name+"'" for p in FINAL_LIBRARY if p.is_role_token()))
 sc.strand('source templates',[role('possible co-carried USER-opened original/copy templates; integrity/partition unassigned')],label='independent residual template lineage')
 sc.strand('carrier background',[role('possible carrier and USER byproducts; integrity/partition unassigned')],label='independent possible residual background')
 return sc

POST_PCR=post_pcr_state()

def protocol_graph():
 g=ProtocolGraph(TITLE,authority=SOURCE)
 dna=g.input('native_dna','agarose-embedded DNA with native3′OH ends',MaterialKind.BIOLOGICAL,state=end_state(),primary=True)
 tdt=g.input('tdt','TdT plus ATP',MaterialKind.CATALYST)
 tailed=g.product('tailed','variable1–4rA-tailed3′ ends',MaterialKind.MOLECULAR,state=end_state(tailed=True))
 g.react('tail','TdT/ATP37°C2h in agarose plugs',StepKind.TAILING,((dna,InputRole.TEMPLATE),(tdt,InputRole.CATALYST)),(tailed,))
 a1=g.input('adapter1','source5′adenylated internal-biotin hairpin1',MaterialKind.ADAPTER,state=adapter1_scene())
 rl=g.input('rna_ligase','T4 RNA ligase2 truncatedKQ/PEG',MaterialKind.CATALYST)
 tagged=g.product('tagged','adapter1-linked native capture strand',MaterialKind.MOLECULAR,state=end_state(tailed=True,tagged=True))
 g.react('ligate1','Ligate activated hairpin1 at25°C overnight',StepKind.LIGATION,((tailed,InputRole.TEMPLATE),(a1,InputRole.STRUCTURAL),(rl,InputRole.CATALYST)),(tagged,),geometry=AdapterLigationGeometry.known('single-strand',adapter1_scene(),statement='Source rA-tailed3′OH joins activated5′App hairpin1.'))
 wash=g.input('plug_wash','three Tris washes plus overnight plug wash',MaterialKind.SELECTION)
 washed=g.product('washed','washed capture-tagged plug DNA',MaterialKind.MOLECULAR,state=tagged.state);waste=g.product('free_a1','unligated adapter1/plug-wash material',MaterialKind.WASTE)
 g.react('wash_plug','Wash plug extensively after capture ligation',StepKind.CLEANUP,((tagged,InputRole.TEMPLATE),(wash,InputRole.SELECTION)),(washed,waste));g.discard(waste,'Source plug-wash fraction; no complete depletion claim.')
 agarase=g.input('agarase','beta-agarase after65°C melting/42°C incubation',MaterialKind.CATALYST)
 agarose_conditions=g.input('agarose_conditions','source65°C melt/42°C agarase-buffer conditions',MaterialKind.CONDITION)
 released=g.product('released','agarose-free adapter1-linked DNA',MaterialKind.MOLECULAR,state=tagged.state)
 g.react('remove_agarose','Melt plug and beta-agarase digest42°C1h',StepKind.EXTRACTION,((washed,InputRole.TEMPLATE),(agarase,InputRole.CATALYST),(agarose_conditions,InputRole.CONDITION)),(released,))
 ethanol=g.input('ethanol','NH4OAc/ethanol precipitation conditions',MaterialKind.CONDITION)
 recovered=g.product('recovered','precipitated adapter1-linked DNA',MaterialKind.MOLECULAR,state=tagged.state);w=g.product('precip_waste','pre-Bst precipitation supernatant',MaterialKind.WASTE)
 g.react('precipitate1','Ethanol precipitate released DNA',StepKind.CLEANUP,((released,InputRole.TEMPLATE),(ethanol,InputRole.CONDITION)),(recovered,w));g.discard(w,'Physical precipitation fraction.')
 bst=g.input('bst','Bst2.0 WarmStart',MaterialKind.CATALYST)
 filled=g.product('filled','hairpin-extended capture DNA',MaterialKind.MOLECULAR,state=extended_state())
 g.react('extend_hairpin','Bst extension65°C30min across UMI/rA/genomic template',StepKind.EXTENSION,((recovered,InputRole.TEMPLATE),(bst,InputRole.CATALYST)),(filled,))
 rec2=g.product('recovered2','precipitated extended DNA in130µlTE',MaterialKind.MOLECULAR,state=filled.state);w2=g.product('precip2_waste','post-Bst precipitation supernatant',MaterialKind.WASTE)
 g.react('precipitate2','Ethanol precipitate after Bst extension',StepKind.CLEANUP,((filled,InputRole.TEMPLATE),(ethanol,InputRole.CONDITION)),(rec2,w2));g.discard(w2,'Physical precipitation fraction.')
 sonication=g.input('sonication','CovarisE220 source sonication settings',MaterialKind.CONDITION)
 fragmented=g.product('fragmented','sonicated biotin-tagged fragments and other DNA',MaterialKind.MOLECULAR,state=filled.state)
 g.react('sonicate','Fragment extended capture DNA',StepKind.FRAGMENTATION,((rec2,InputRole.TEMPLATE),(sonication,InputRole.CONDITION)),(fragmented,))
 beads=g.input('beads','MyOne streptavidinC1',MaterialKind.SELECTION)
 captured=g.product('captured','washed bead-bound biotin capture fragments',MaterialKind.MOLECULAR,state=filled.state);uncaptured=g.product('uncaptured','unbound DNA and capture-wash material',MaterialKind.WASTE)
 g.react('capture','Bind biotin fragments and wash with salt/TE',StepKind.AFFINITY_CAPTURE,((fragmented,InputRole.TEMPLATE),(beads,InputRole.SELECTION)),(captured,uncaptured));g.discard(uncaptured,'Source binding supernatant/washes; exact species partition unassigned.')
 prep=g.input('endprep_mix','NEBNext UltraII EndPrep named in paper; v6.1_5/20 kit interpretation',MaterialKind.CATALYST)
 carrier=g.input('carrier','1ng sonicated salmon-sperm carrier DNA',MaterialKind.MOLECULAR)
 prepared=g.product('prepared','INFERRED—TA-compatible on-bead end-prepared DNA plus carrier',MaterialKind.MOLECULAR,state=endprep_state())
 g.react('endprep','On-bead EndPrep20°C30min/65°C30min; historical kit revision unassigned',StepKind.END_REPAIR,((captured,InputRole.TEMPLATE),(prep,InputRole.CATALYST),(carrier,InputRole.SPIKE_IN)),(prepared,))
 a2=g.input('adapter2','source5′P hairpin2 with dU4 and terminal*T',MaterialKind.ADAPTER,state=adapter2_scene());ligase=g.input('ligase2','UltraII ligation mix/enhancer',MaterialKind.CATALYST)
 dual=g.product('dual','bead-bound two-hairpin capture DNA plus ligation mixture',MaterialKind.MOLECULAR,state=dual_state())
 g.react('ligate2','INFERRED—TA-compatible adapter2 ligation20°C30min',StepKind.LIGATION,((prepared,InputRole.TEMPLATE),(a2,InputRole.STRUCTURAL),(ligase,InputRole.CATALYST)),(dual,),geometry=AdapterLigationGeometry.known('cohesive',endprep_state(),adapter2_scene(),statement='Source adapter2 has one3′T outside its31bp stem. Source-named kit EndPrep dA activity described by May2020 manual; exact historical kit revision remains unassigned.'))
 salt=g.input('bead_wash','source salt and0.1×TE bead washes',MaterialKind.SELECTION)
 clean=g.product('clean','washed bead-bound two-adapter material',MaterialKind.MOLECULAR,state=dual.state);dimers=g.product('ligation_waste','ligation supernatant/adapter2 homodimer-containing wash fraction',MaterialKind.WASTE)
 g.react('wash_beads','Discard ligation mix and wash beads',StepKind.CLEANUP,((dual,InputRole.TEMPLATE),(salt,InputRole.SELECTION)),(clean,dimers));g.discard(dimers,'Source ligation supernatant and bead washes; no total carrier depletion assumed.')
 user=g.input('user','USER enzyme',MaterialKind.CATALYST)
 opened=g.product('opened','combined two USER eluates',MaterialKind.MOLECULAR,state=opened_state())
 remnants=g.product('user_remnants','biotin-loop/bead-associated or hydrolysed adapter remnants; exact fate unassigned',MaterialKind.MOLECULAR)
 g.react('user_elute','USER cleave/elute twice15min37°C; combine eluates',StepKind.DIGESTION,((clean,InputRole.TEMPLATE),(user,InputRole.CATALYST)),(opened,remnants));g.mark_terminal(Terminal(remnants.id,Disposition.UNRESOLVED,reason='Cleavage and target elution do not establish exact remnant physical disposal.'))
 aliquot=g.input('aliquot_amounts','source1.25µl testPCR and21µl mainPCR aliquots',MaterialKind.CONDITION)
 test_input=g.product('test_input','testPCR aliquot',MaterialKind.MOLECULAR,state=opened.state);main_input=g.product('main_input','mainPCR aliquot',MaterialKind.MOLECULAR,state=opened.state)
 g.react('aliquot','Separate test and main USER-eluate inputs',StepKind.PARTITION,((opened,InputRole.TEMPLATE),(aliquot,InputRole.CONDITION)),(test_input,main_input))
 f=g.input('forward','NEBNext Universal primer; native strand unassigned',MaterialKind.OLIGO);r=g.input('reverse','chosen NEBNext Index primer; actual set/member unassigned',MaterialKind.OLIGO);q5=g.input('q5','UltraII Q5 PCR mix',MaterialKind.CATALYST)
 test_amp=g.product('test_amp','18-cycle test amplification for Bioanalyzer cycle choice',MaterialKind.QC,state=POST_PCR)
 g.react('test_pcr','TestPCR18cycles; cleanup/Bioanalyzer QC outside selected sequencing output',StepKind.PCR,((test_input,InputRole.TEMPLATE),(f,InputRole.STRUCTURAL),(r,InputRole.STRUCTURAL),(q5,InputRole.CATALYST)),(test_amp,),geometry=PCRGeometry.unavailable(test_input.state,reason='Actual USER-opened precursor; selected Index strand and full native retained end boundaries unassigned.'))
 g.mark_terminal(Terminal(test_amp.id,Disposition.QC,reason='Source testPCR/Bioanalyzer branch determines main-cycle choice; not a sequencing endpoint.'))
 choice=g.input('cycle_choice','main ncycles selected using testPCR/Bioanalyzer',MaterialKind.CONDITION)
 amp=g.product('amplified','main amplified capture-library mixture',MaterialKind.LIBRARY,state=POST_PCR)
 g.react('indexed_pcr','MainPCR:98°C30s; n×98°C10s/65°C75s;65°C5min',StepKind.PCR,((main_input,InputRole.TEMPLATE),(f,InputRole.STRUCTURAL),(r,InputRole.STRUCTURAL),(q5,InputRole.CATALYST),(choice,InputRole.CONDITION)),(amp,),geometry=PCRGeometry.unavailable(main_input.state,reason='Actual USER-opened mixture; no canonical final library or fixed six-base index primer is used as the substrate.'))
 ampure=g.input('ampure','AMPureXP',MaterialKind.SELECTION)
 clean1=g.product('clean1','first40µl-bead cleanup eluate26µl',MaterialKind.LIBRARY,state=POST_PCR);w1=g.product('waste1','firstAMPure non-retained fraction',MaterialKind.WASTE)
 g.react('cleanup1','First40µlAMPure cleanup; elute26µl',StepKind.CLEANUP,((amp,InputRole.TEMPLATE),(ampure,InputRole.SELECTION)),(clean1,w1));g.discard(w1,'Physical first cleanup fraction.')
 library=g.product('library',FINAL_LIBRARY.name,MaterialKind.LIBRARY,state=POST_PCR);w3=g.product('waste2','secondAMPure non-retained fraction',MaterialKind.WASTE)
 g.react('cleanup2','Take25µl;20µlAMPure cleanup, elute11µl and QC',StepKind.CLEANUP,((clean1,InputRole.TEMPLATE),(ampure,InputRole.SELECTION)),(library,w3));g.discard(w3,'Physical second cleanup fraction.')
 g.sequence(library,SequencingPlan(library=FINAL_LIBRARY,primers=(),limitation='Source MiSeq50SR/NextSeq75SE; UMI8 and antisense readout known, native index allocation and primer/USER end geometry unassigned.'))
 g.declare_reviewed_scope(('native_dna','tdt','tail','tailed','adapter1','rna_ligase','ligate1','tagged','plug_wash','wash_plug','washed','agarase','agarose_conditions','remove_agarose','released','ethanol','precipitate1','recovered','bst','extend_hairpin','filled','precipitate2','recovered2','sonication','sonicate','fragmented','beads','capture','captured','endprep_mix','endprep','carrier','prepared','adapter2','ligase2','ligate2','dual','bead_wash','wash_beads','clean','user','user_elute','opened','aliquot_amounts','aliquot','test_input','main_input','forward','reverse','q5','test_pcr','test_amp','cycle_choice','indexed_pcr','amplified','ampure','cleanup1','clean1','cleanup2','library'),SOURCE)
 return g.require_publishable()
