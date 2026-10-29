"""hiCLIP: local source-qualified adapter-block gating, not a finished library."""
from dataclasses import dataclass,replace
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Construct,Segment,Scene,feature,revcomp
from adapter_ligation import AdapterLigationGeometry,reported_adapter_end
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,Terminal,Disposition,SourceItem
import illumina as il
SOURCE='Sugimoto et al.DOI10.1038/nature14280; accepted manuscript PMC4376666 Methods104–140,160–174; published-version applicability unverified'
A_RNA=revcomp(il.PE_READ2_ARM)[:20].replace('T','U')
# Source-specific bridging linker; no inherited recipe.
B_RNA='CUGUAGGCACCAUACAAUG'
RT_BODY=il.INDEX2_PRIMER_RC[:14]+'TCGTGGATCCTGAACCGC'

@dataclass(frozen=True)
class BlockedNucleotide:
    base: str='C'
    sugar: str='2,3-dideoxyribose'
    oh3: bool=False
    def __post_init__(self):
        if (self.base,self.sugar,self.oh3)!=('C','2,3-dideoxyribose',False):raise ValueError('A terminates in blocked cytidine nucleotide ddC')

@dataclass(frozen=True)
class Adapter:
    name: str
    activation5: str='phosphate'
    block3: str='source_default'
    terminal_nucleotide: BlockedNucleotide|None=None
    def __post_init__(self):
        if self.name not in {'A','B'} or self.activation5 not in {'phosphate','rApp','ligated'}:raise ValueError('invalid source adapter/5′state')
        if self.block3=='source_default':object.__setattr__(self,'block3','ddC' if self.name=='A' else 'phosphate')
        if self.name=='A':
            if self.block3!='ddC':raise ValueError('PNK cannot restore3′OH to ddC')
            if self.terminal_nucleotide is None:object.__setattr__(self,'terminal_nucleotide',BlockedNucleotide())
        elif self.block3 not in {'phosphate','OH','ligated'} or self.terminal_nucleotide is not None:raise ValueError('B block is a removable phosphate, not an extra nucleotide')
    @property
    def bases(self):return A_RNA+'C' if self.name=='A' else B_RNA
    def adenylate(self):
        if self.activation5!='phosphate':raise ValueError('adenylation requires source5′phosphate')
        return replace(self,activation5='rApp')
    def attach_to_rna3(self):
        if self.activation5!='rApp':raise ValueError('first ligation requires preadenylated adapter')
        return replace(self,activation5='ligated')
    def pnk(self):return self if self.name=='A' else replace(self,block3='OH')
    def bridge(self):
        if (self.name,self.activation5,self.block3)!=('B','ligated','OH'):raise ValueError('B bridge must be attached at5′ and deblocked at3′')
        return replace(self,block3='ligated')
    def construct(self):
        parts=[Segment('adapter '+self.name+' RNA body',(A_RNA if self.name=='A' else B_RNA).replace('U','T'))]
        if self.terminal_nucleotide:parts.append(Segment('terminal ddC nucleotide','C'))
        return Construct(parts)

@dataclass(frozen=True)
class SelectedBridge:
    upstream_rna: str
    downstream_rna: str
    a: Adapter
    b: Adapter
    downstream5: str='unassigned'
    joined: bool=False
    def __post_init__(self):
        if not self.upstream_rna or not self.downstream_rna or set(self.upstream_rna+self.downstream_rna)-set('ACGU'):raise ValueError('two explicit RNA arms required')
        if self.a.name!='A' or self.b.name!='B' or self.a.activation5!='ligated' or self.b.activation5!='ligated':raise ValueError('selected arms must carry source adapters')
        if self.joined and (self.downstream5!='consumed' or self.b.block3!='ligated'):raise ValueError('joined junction must consume both selected ends')
    def pnk(self):
        if self.joined:raise ValueError('PNK represented before bridge formation')
        return replace(self,a=self.a.pnk(),b=self.b.pnk(),downstream5='phosphate')
    def ligate(self):
        if self.joined or self.downstream5!='phosphate':raise ValueError('downstreamRNA5′phosphate required')
        return replace(self,b=self.b.bridge(),downstream5='consumed',joined=True)
    @property
    def ordered_bases(self):
        if not self.joined:raise ValueError('two arms not covalently joined yet')
        return self.upstream_rna+self.b.bases+self.downstream_rna+self.a.bases

def rt_primer(barcode='TATT'):
    if barcode not in {'TATT','AACC','CGCC','GTGG'}:raise ValueError('this local inventory covers source hiCLIP samples only')
    return Construct([Segment('UMI part1','NN',placeholder=True,feature=feature('hiclip_umi_2','umi','random',group='hiclip_composite_umi',part='first')),
      Segment('sample barcode',barcode,feature=feature('hiclip_sample','inline_barcode','fixed')),
      Segment('UMI part2','NNN',placeholder=True,feature=feature('hiclip_umi_3','umi','random',group='hiclip_composite_umi',part='second')),
      Segment('RT body',RT_BODY)])

def primer_scene():
    binding=revcomp(RT_BODY[-9:]);a=Adapter('A');body=a.bases.replace('U','T');start=body.index(binding)
    sc=Scene();sc.strand('A',[Segment('prefix',body[:start]),Segment('binding',binding),Segment('terminal ddC','C')],mod3='ddC: no3′OH')
    sc.anneal('RT',[*list(rt_primer())[:-1],Segment('RT body upstream',RT_BODY[:-9]),Segment('RT3′binding',RT_BODY[-9:])],to='A',pair=('RT3′binding','binding'),unpaired=('RT body upstream',))
    return sc

def protocol_graph():
    g=ProtocolGraph('hiCLIP — selected adapter-gating prefix',authority=SOURCE,representation='published_instance')
    rna=g.input('bound_duplex','STAU1-associated RNA duplex on beads after source3′dephosphorylation; pairing/arm sequences variable',K.MOLECULAR,primary=True)
    a=g.input('adapter_a','Source PAGE-purified preadenylated A, terminal ddC retained',K.ADAPTER,state=Adapter('A').adenylate().construct(),chemical_state=Adapter('A').adenylate())
    b=g.input('adapter_b','Source PAGE-purified preadenylated B,3′phosphate blocked',K.ADAPTER,state=Adapter('B').adenylate().construct(),chemical_state=Adapter('B').adenylate())
    ligase=g.input('truncated_ligase','T4RNA ligase2 truncatedK227Q,16°C overnight',K.CATALYST)
    attached=g.product('attached','Selected one-A/one-B duplex class; AA/BB/unligated populations not assigned proportions',K.MOLECULAR,chemical_state={'A':Adapter('A').adenylate().attach_to_rna3(),'B':Adapter('B').adenylate().attach_to_rna3(),'arm5':'unassigned'})
    ends=[reported_adapter_end('source'+name+'5′rApp',kind='single-strand',polarity='5-prime',terminal_sequence=seq[0],body_sequence=seq[1:].replace('U','T'),statement='Source preadenylated adapter donor; exact transcript arm sequence variable') for name,seq in [('A',A_RNA+'C'),('B',B_RNA)]]
    geo=AdapterLigationGeometry.known('single-strand',*(e.scenes[0] for e in ends),statement='Source A/B5′rApp donor nucleotides attach to prepared RNA3′OH; typed adapter state preserves terminal blocks')
    g.react('first_ligation','Source equimolar A/B; displayed class is conditional, not guaranteed50% recovery',S.LIGATION,((rna,R.TEMPLATE),(a,R.STRUCTURAL),(b,R.STRUCTURAL),(ligase,R.CATALYST)),(attached,),geometry=geo)
    def clean(identifier,prior,label):
        keep=g.product(identifier,label,K.MOLECULAR,chemical_state=prior.chemical_state)
        waste=g.product(identifier+'_waste','Physical bead wash rejects; free adapter/reagent fractions unquantified',K.WASTE)
        g.react(identifier+'_wash','Source bead washing',S.CLEANUP,((prior,R.TEMPLATE),),(keep,waste));g.discard(waste,'Physical wash separation');return keep
    washed=clean('washed_first',attached,'Adapter-bearing duplexes after first-ligation washes')
    pnk=g.input('pnk','T4PNK plusATP37°C30min',K.CATALYST)
    activated=g.product('deblocked','Selected A/B class: RNA5′phosphorylated; B3′OH restored; A remains ddC-blocked',K.MOLECULAR,chemical_state={'A':attached.chemical_state['A'].pnk(),'B':attached.chemical_state['B'].pnk(),'arm5':'phosphate'})
    g.react('pnk_deblock','Source simultaneous RNA5′phosphorylation and adapterB3′dephosphorylation',S.END_REPAIR,((washed,R.TEMPLATE),(pnk,R.CATALYST)),(activated,))
    washed2=clean('washed_pnk',activated,'Activated duplexes after source PNK washes')
    ligase1=g.input('ligase1','T4RNA ligase1 withATP,16°C overnight',K.CATALYST)
    joined=g.product('hybrid_rna','Conditional RNAarm1–B–RNAarm2–A product; terminal A ddC remains; alternative/unjoined products unassigned',K.MOLECULAR,chemical_state={'order':('RNAarm1','B','RNAarm2','A'),'A':activated.chemical_state['A'],'B':activated.chemical_state['B'].bridge()})
    g.react('second_ligation','Selected intraduplex bridge by B3′OH→otherRNA5′p; no yield prediction',S.PROXIMITY_LIGATION,((washed2,R.TEMPLATE),(ligase1,R.CATALYST)),(joined,),geometry=AdapterLigationGeometry.unavailable('Conditional terminal compatibility/order modeled locally; actual duplex3D arm placement and ligation geometry are unassembled.'))
    washed3=clean('washed_hybrid',joined,'Source washed hybrid/unjoined RNA population before stringent purification')
    g.mark_terminal(Terminal(washed3.id,Disposition.UNRESOLVED,'Accepted-manuscript procedure applicability to final VOR unverified. Upstream capture, stringent purification/proteolysis, RT truncation/RNA hydrolysis, cDNA circle/BamHI/PCR and sequencing geometry remain unassembled.'))
    g.declare_reviewed_scope(['bound_duplex','adapter_a','adapter_b','truncated_ligase','first_ligation','washed_first_wash','pnk','pnk_deblock','washed_pnk_wash','ligase1','second_ligation','washed_hybrid_wash'],SOURCE)
    g.add_source_item(SourceItem('unreviewed_final_vor','Final published procedure applicability and remaining preparation/read scope','https://doi.org/10.1038/nature14280'))
    return g

def selftest():
    from checks import Check,Source,have
    from storage import sources_root
    c=Check();a=Adapter('A').adenylate().attach_to_rna3();b=Adapter('B').adenylate().attach_to_rna3()
    p=SelectedBridge('ACGU','UGCA',a,b)
    try:p.ligate();bad=False
    except ValueError:bad=True
    c('blockedB cannot bridge beforePNK',bad,True)
    q=p.pnk().ligate();c('bridge preserves both arms and full adapter nucleotide inventory',q.ordered_bases,'ACGU'+B_RNA+'UGCA'+A_RNA+'C')
    c('PNK removesB phosphate without converting terminalddC toOH',(q.a.terminal_nucleotide.oh3,q.b.block3),(False,'ligated'))
    primer_scene();c('RT3′binding occurs insideA body, before terminalddC',revcomp(RT_BODY[-9:]) in A_RNA.replace('U','T'),True)
    parts=[s.feature for s in rt_primer() if s.feature and s.feature.role=='umi'];c('compositeUMI uses unique partIDs and shared group',(len({x.id for x in parts}),len({x.group for x in parts})),(2,1))
    original=Source(sources_root()/'hiclip__10.1038+nature14280'/'nature14280_PMC4376666.xml.txt','accepted-manuscript Methods104/164')
    if have(c,original,label='source oligos'):
        text=''.join(original.path.read_text().split())
        written='/5Phos/'+''.join('r'+x for x in A_RNA)+'/3ddC/'
        c('sourceA body and terminalddC transcribed',written in text,True)
        c('sourceB terminalphosphate transcribed','/5Phos/'+''.join('r'+x for x in B_RNA)+'/3Phos/' in text,True)
        c('sourceRT sample/UMI/body transcribed','/5Phos/NNTATTNNN'+RT_BODY in text,True)
    c.report()
if __name__=='__main__':selftest()
