"""Published Subassembly SI metagenomic route, with independent tag-pairing output."""
from dataclasses import dataclass
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
import illumina as il
from chemdraw import Scene,Construct,Segment,feature,revcomp,panel,workflow_from_protocol_graph,workflow_panel
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,Terminal,Disposition,SourceItem,SequencingPlan
from pcr import PCRGeometry
from adapter_ligation import AdapterLigationGeometry
from page import head,info
PAGE_FILE='chem.html'
TITLE='Subassembly — published metagenomic preparation'
SOURCE='Published NatureMethods DOI10.1038/nmeth.1416 SI MOESM439,Table4 printedp11 and Protocols1–3;mainPMC2848820 is an accepted manuscript,not substituted for published chemistry'
SI='nmeth.1416_41592_2010_BFnmeth1416_MOESM439_ESM.pdf.txt'
LIMIT='Private incomplete model. Native fragment-end tags are not synthetic random UMIs. Complete heterogeneous source-DNA,concatemer,breakpoint and final run-primer geometry remain unassembled. Published SI owns this preparation;accepted-manuscript overview is version-qualified separately.'
AD1='TCGCAATACAGAGTTTACCGCATT'
AD2='CTCTTCCGCATCTCACAACCTACT'
BP='CTCGGCATTCCTGCTGAACCGCTCTTCCGATCT'
SEQ1='CAATGGAGC'+AD1
SEQ2='ATCGAGAGC'+AD2
SEQBP='CGGT'+BP
F1=il.P5+SEQ1
F2=il.P5+SEQ2
RBP=il.P7+SEQBP
MERGE_R=il.P7+SEQ2

@dataclass(frozen=True)
class SourceOligo:
    name:str
    sequence:str
    phosphate5:bool=False
    ps_after:tuple[int,...]=()
    def ordering(self):
        return ('/5Phos/' if self.phosphate5 else '')+''.join(b+('*' if i in self.ps_after else '') for i,b in enumerate(self.sequence))
    def as_scene(self):
        sc=Scene();sc.strand('oligo',[Segment(self.name,self.sequence)],mod5='phosphate' if self.phosphate5 else '')
        if self.ps_after:sc.note('oligo','Source physical phosphorothioate bonds: '+self.ordering())
        sc.source_physical_oligo=self
        return sc
    def rows(self):return self.as_scene().rows()

OLIGOS={
 'ad1amp':SourceOligo('Ad1_amp',AD1,True),'ad2amp':SourceOligo('Ad2_amp',AD2,True),
 'f1':SourceOligo('Illum_amp_f_Ad1',F1),'f2':SourceOligo('Illum_amp_f_Ad2',F2),
 'rbp':SourceOligo('Illum_amp_r',RBP),'merge_r':SourceOligo('Illum_amp_r_Ad2',MERGE_R)}

@dataclass(frozen=True)
class AdapterMixture:
    components:tuple[SourceOligo,...]
    scenes:tuple[Scene,...]
    def as_scene(self):
        sc=Scene()
        for top,bottom in zip(self.components[::2],self.components[1::2]):adapter_scene(top.sequence,top.name,bool(top.ps_after),sc)
        sc.source_adapter_mixture=self
        return sc
    def rows(self):return self.as_scene().rows()

def adapter_scene(body,name,ps=False,scene=None):
    # Table4 source pairs have a 3′T extension on the top and5′P bottom.
    sc=scene if scene is not None else Scene();top=name+'_top';bottom=name+'_bottom'
    sc.strand(top,[Segment(name+' core',body[:-1]),Segment(name+' 3primeT',body[-1])])
    sc.anneal(bottom,[Segment(name+' bottom',revcomp(body[:-1]))],to=top,pair=(name+' bottom',name+' core'),mod5='phosphate')
    sc.ligation_end(top,name+' 3primeT',kind='cohesive',polarity='3-prime',owner='adapter')
    if ps:sc.note(top,'Source phosphorothioate bond immediately before terminal T;this modification is not assigned to polymerase-made copies.')
    return sc

def native_tag_scene():
    sc=Scene();sc.strand('tag',[Segment('native genomic end tag','N'*20,placeholder=True,feature=feature('subassembly.native_tag','linked_read_barcode','unknown',note='Native fragment-end sequence;not a synthesized random barcode.'))])
    sc.note('tag','Identifier component only;no arbitrary genomic sequence or full final library assembled.');return sc

def unknown_geometry(reason):
    sc=Scene();sc.strand('unknown',[Segment('unassembled strands','[unassembled heterogeneous molecular geometry]',placeholder=True)])
    return AdapterLigationGeometry('unavailable','unavailable',(sc,),reason)

def protocol_graph():
    g=ProtocolGraph(TITLE,authority=SOURCE)
    def reagent(key,label,kind=K.CATALYST,state=None):return g.input(key,label,kind,state=state.as_scene() if hasattr(state,'as_scene') else state)
    def op(key,x,label,kind,reagents=(),geometry=None):
        y=g.product(key+'_product',label+';reaction species remain until workup',K.MOLECULAR)
        kw={} if geometry is None else {'geometry':geometry}
        g.react(key,label,kind,((x,R.TEMPLATE),)+tuple(reagents),(y,),**kw);return y
    def clean(key,x,label):
        sel=reagent(key+'_selection',label,K.SELECTION);y=g.product(key+'_kept','Retained DNA after '+label,K.MOLECULAR);w=g.product(key+'_waste','Removed physical workup fractions',K.WASTE)
        g.react(key,label,S.CLEANUP,((x,R.TEMPLATE),(sel,R.SELECTION)),(y,w));g.discard(w,'Physical column/gel workup;no quantitative purity/recovery assumed');return y
    def pcr(key,x,f,r,label):
        enzyme=reagent(key+'_polymerase','Source iProofHF/SYBR real-time PCR' if key=='merging_pcr' else 'Source Phusion Hot-Start/SYBR real-time PCR')
        return op(key,x,label,S.PCR,((f,R.STRUCTURAL),(r,R.STRUCTURAL),(enzyme,R.CATALYST)),PCRGeometry.unavailable(Scene(),reason='Source primers known;full target/concatemer and heterogeneous product geometry unassembled.'))
    def split(key,x,labels):
        ys=tuple(g.product(key+'_'+str(i),label,K.MOLECULAR) for i,label in enumerate(labels))
        g.react(key,'Take independent physical aliquots',S.SPLIT,((x,R.TEMPLATE),),ys);return ys
    dna=g.input('dna','Metagenomic source DNA;source extraction/enrichment outside graph',K.MOLECULAR,primary=True)
    shear=reagent('bioruptor','Source Bioruptor sonication');x=op('shear',dna,'Fragment source DNA',S.FRAGMENTATION,((shear,R.CATALYST),))
    x=clean('shear_cleanup',x,'QIAquick purification')
    repair=reagent('endit','Epicentre End-It DNA End Repair kit,source45min roomtemperature');x=op('repair1',x,'End repair',S.END_REPAIR,((repair,R.CATALYST),));x=clean('repair1_cleanup',x,'QIAquick purification')
    x=clean('input_size',x,'Metagenomic400–500bp6%TBE selection/ethanol recovery')
    taq=reagent('taq','Source Taq/dATP70C20min');x=op('atail1',x,'Add terminal3primeA',S.TAILING,((taq,R.CATALYST),));x=clean('atail1_cleanup',x,'QIAquick purification')
    adapters=reagent('tag_adapters','Annealed Ad1/Ad1_rc andAd2/Ad2_rc;two source T-tailed adapters',K.ADAPTER,state=AdapterMixture((SourceOligo('Ad1',AD1),SourceOligo('Ad1_rc',revcomp(AD1[:-1]),True),SourceOligo('Ad2',AD2),SourceOligo('Ad2_rc',revcomp(AD2[:-1]),True)),(adapter_scene(AD1,'Ad1'),adapter_scene(AD2,'Ad2'))))
    ligase=reagent('quick_ligase','NEB Quick Ligase,roomtemperature15min')
    geom=AdapterLigationGeometry.known('cohesive',adapter_scene(AD1,'Ad1'),adapter_scene(AD2,'Ad2'),statement='Published Table4 exact paired cores,top3primeT andbottom5primeP;Protocol2step6 ligates to A-tailed DNA. Full source insert/product geometry unassembled.')
    x=op('tag_ligation',x,'Ligate tag-adjacent adaptors',S.LIGATION,((adapters,R.STRUCTURAL),(ligase,R.CATALYST)),geom)
    x=clean('ligated_size',x,'400–800bp6%TBE selection/ethanol recovery')
    chosen,unused=split('bottleneck',x,('Dilution chosen for complexity bottleneck','Untransferred ligated material'))
    g.mark_terminal(Terminal(unused.id,Disposition.UNRESOLVED,'Untransferred material not assumed disposed;alternative dilution trials are not physically pooled.'))
    a=reagent('ad1_amp','Physical5primeP Ad1_amp',K.OLIGO,OLIGOS['ad1amp']);b=reagent('ad2_amp','Physical5primeP Ad2_amp',K.OLIGO,OLIGOS['ad2amp'])
    x=pcr('bottleneck_pcr',chosen,a,b,'Step8a log-phase-monitored bottleneck PCR');x=clean('bottleneck_cleanup',x,'QIAquick cleanup')
    x=clean('metagenomic_resizing',x,'Metagenomic-only500–600bp size selection')
    x=pcr('after_resizing_pcr',x,a,b,'Repeat8a after metagenomic resizing');x=clean('resized_pcr_cleanup',x,'QIAquick cleanup')
    x=pcr('bulk_reamp',x,a,b,'Step8c amplification in eight reactions,pool products');x=clean('bulk_cleanup',x,'Pool andQIAquick purify')
    main,pair,rem=split('readout_split',x,('Main concatemer preparation aliquot','Optional metagenomic TDRG merging aliquot','Unallocated amplified source material'))
    g.mark_terminal(Terminal(rem.id,Disposition.UNRESOLVED,'Source aliquot remainder not assigned a disposal fate.'))
    f1=reagent('f1','Source Illum_amp_f_Ad1',K.OLIGO,OLIGOS['f1']);f2=reagent('f2','Source Illum_amp_f_Ad2',K.OLIGO,OLIGOS['f2']);rm=reagent('merge_r','Source Illum_amp_r_Ad2',K.OLIGO,OLIGOS['merge_r']);rbp=reagent('rbp','Source Illum_amp_r',K.OLIGO,OLIGOS['rbp'])
    pair=pcr('merging_pcr',pair,f1,rm,'Optional8d iProofHF merging-library PCR;Ad1_seq/Ad2_seq36×36readout')
    g.sequence(pair,SequencingPlan(acquired_roles=('Read 1','Read 2'),read_lengths={'Read 1':36,'Read 2':36},limitation='Published SI8d separate36×36merging-library sequencing,Ad1_seq/Ad2_seq;exactfullproduct/primeracquisition geometry and workup not assembled. No physical merge with main molecular library implied.'))
    # Concatemerization is intermolecular joining, not a source circularization recipe.
    main=op('concatemerize',main,'Step9 blunt-ligate phosphorylatedPCR products into high-molecular-weight concatemers',S.LIGATION,((main,R.STRUCTURAL),(ligase,R.CATALYST)),unknown_geometry('The existing PCR-product population supplies both ligating partners;no new primer or adapter is added. Heterogeneous full concatemer end/strand population unassembled.'))
    main=op('reshear',main,'Step10 Bioruptor fragmentation of concatemers',S.FRAGMENTATION,((shear,R.CATALYST),));main=clean('reshear_cleanup',main,'QIAquick purification')
    main=op('repair2',main,'Step11 End-It end repair',S.END_REPAIR,((repair,R.CATALYST),));main=clean('repair2_cleanup',main,'QIAquick purification')
    main=op('atail2',main,'Step12 Taq3primeA tailing',S.TAILING,((taq,R.CATALYST),));main=clean('atail2_cleanup',main,'QIAquick purification')
    bp=reagent('breakpoint_adapter','Illum_rev/Illum_rev_rc;source terminalPS bond andbottom5primeP',K.ADAPTER,state=AdapterMixture((SourceOligo('Illum_rev',BP,ps_after=(len(BP)-2,)),SourceOligo('Illum_rev_rc',revcomp(BP[:-1]),True)),(adapter_scene(BP,'breakpoint',True),)))
    main=op('breakpoint_ligation',main,'Step13 breakpoint-adjacent adaptor ligation',S.LIGATION,((bp,R.STRUCTURAL),(ligase,R.CATALYST)),AdapterLigationGeometry.known('cohesive',adapter_scene(BP,'breakpoint',True),statement='Published source calls T-tailed adapter;Table4 two oligos pair with a3primeT overhang. Full insert/junction unassembled.'))
    main=clean('breakpoint_cleanup',main,'QIAquick purification')
    left,right,unused=split('end_pcr_split',main,('Ad1-targeted PCR aliquot','Ad2-targeted PCR aliquot','Untransferred breakpoint-ligation preparation'))
    g.mark_terminal(Terminal(unused.id,Disposition.UNRESOLVED,'Untransferred source eluate not presumed discarded.'))
    for key,y,fp in (('ad1',left,f1),('ad2',right,f2)):
        y=pcr(key+'_readout_pcr',y,fp,rbp,'Step14 separate endpoint-specific log-phasePCR')
        y=clean(key+'_pcr_cleanup',y,'QIAquick purification')
        y=clean(key+'_readout_size',y,'Step15 metagenomic300–450bp size selection')
        y=pcr(key+'_final_pcr',y,fp,rbp,'Step15 finalPCR for adequate sequencing material');y=clean(key+'_final_cleanup',y,'QIAquick purification')
        g.sequence(y,SequencingPlan(acquired_roles=('Read 1','Read 2'),read_lengths={'Read 1':20,'Read 2':76},limitation='Published SI16:20bp native tagRead1/76bp breakpointRead2;Ad1_seq orAd2_seq plusIllum_seq_r. Full final run geometry unassembled.'))
    g.declare_reviewed_scope(tuple(g.materials)+tuple(g.steps),SOURCE)
    g.add_source_item(SourceItem('unassembled_source','Full molecular/read geometry and remaining original scope',LIMIT));return g

def render_page():
    g=protocol_graph()
    return '\n'.join([head(TITLE),'<main class="wrap"><h1>'+TITLE+'</h1>',info(LIMIT),'<p><a href="01_subassembly.html">Research notes</a></p>',
        ''.join(workflow_panel(workflow_from_protocol_graph(g,e)) for e in ('ad1_final_cleanup_kept','ad2_final_cleanup_kept','merging_pcr_product')),panel(native_tag_scene().rows(),caption='Native-fragment identifier component,not an added random barcode'),'</main>'])

def selftest():
    from checks import Check
    from storage import sources_root
    c=Check();g=protocol_graph();c('Private incomplete source model',g.publishable,False);c('No false completion marker','data-protocol-complete=' in render_page(),False)
    c('Merging branch derives from preconcatemer aliquot',g.steps['merging_pcr'].inputs[0].material,'readout_split_1')
    c('ReadoutPCRs physically separate',[i.material for i in g.steps['end_pcr_split'].inputs],['breakpoint_cleanup_kept'])
    c('Native tag not marked random synthesis',native_tag_scene().strands['tag'].segs[0].feature.encoding,'unknown')
    p=sources_root()/Path(__file__).resolve().parents[1].name/SI
    if p.exists():
        source=''.join(p.read_text().split()).upper()
        for o in OLIGOS.values():c('PublishedTable4 '+o.name,o.ordering().upper() in source,True)
        c('Source breakpoint phosphorothioate',BP[:-1]+'*'+BP[-1] in source,True)
    else:c.skip('PublishedTable4 transcriptions','External original unavailable')
    c.report()
