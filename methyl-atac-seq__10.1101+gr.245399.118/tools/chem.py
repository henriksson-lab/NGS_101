"""Protocol-owned selected chemical classes; incomplete source-qualified graph."""
from dataclasses import dataclass
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Scene, Segment, revcomp, feature
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,SourceItem,SequencingPlan,Terminal,Disposition
from pcr import PCRGeometry
from adapter_ligation import AdapterLigationGeometry

def scene(name,text):
    sc=Scene();sc.strand(name,[Segment(name,text,placeholder=True)]);return sc

def clean(g,key,x,label):
    condition=g.input(key+'_condition',label,K.CONDITION)
    kept=g.product(key+'_kept','Source recovered material; no perfect purity or quantitative partition assumed',K.MOLECULAR)
    waste=g.product(key+'_waste','Actually separated cleanup fractions with unquantified molecular carryover',K.WASTE)
    g.react(key,label,S.CLEANUP,((x,R.TEMPLATE),(condition,R.CONDITION)),(kept,waste));g.discard(waste,'Source physical separation, not perfect molecular classification.');return kept

def op(g,key,x,label,kind=S.MAGIC,*,reagent='source mixture',catalyst=False,product=None,state=None,chemical_state=None):
    stock=g.input(key+'_reagent',reagent,K.CATALYST if catalyst else K.CONDITION)
    y=g.product(key+'_product',product or label+'; all physical reaction companions retained until separation',K.MOLECULAR,state=state,chemical_state=chemical_state)
    kwargs={}
    if kind==S.MAGIC:kwargs['reason']='Source operation established; exact internal chemistry/geometry not assembled from unreviewed historical originals.'
    if kind==S.PCR:kwargs['geometry']=PCRGeometry.unavailable(Scene(),reason='Exact source template/primer endpoint geometry unassembled, no ordinary complete library asserted.')
    if kind in (S.LIGATION,S.PROXIMITY_LIGATION):kwargs['geometry']=AdapterLigationGeometry('unavailable','unavailable',(scene('unassembled','[source-known ligation; full substrate/end geometry unassembled]'),),'Unreviewed geometry; not a claim of unavailable disclosure.')
    inputs=[(x,R.TEMPLATE),(stock,R.CATALYST if catalyst else R.CONDITION)]
    if kind in {S.DIGESTION,S.END_REPAIR,S.EXTENSION,S.TAILING,S.PROXIMITY_LIGATION} and not catalyst:
        enzyme=g.input(key+'_enzyme','Explicit catalytic component of source mixture: '+reagent,K.CATALYST);inputs.append((enzyme,R.CATALYST))
    if kind==S.LIGATION:
        adapter=g.input(key+'_adapter','Source ligating molecular adapter; exact historical oligo/end inventory unassembled',K.ADAPTER,state=scene(key+'_adapter','[source adapter geometry unassembled]'));inputs.append((adapter,R.STRUCTURAL))
    if kind==S.PCR:
        for arm in ('f','r'):
            primer=g.input(key+'_primer_'+arm,'Source PCR primer family; exact oligo/endpoint unassembled',K.OLIGO,state=scene(key+'_primer_'+arm,'[source primer body/index geometry unassembled]'));inputs.append((primer,R.STRUCTURAL))
    g.react(key,label,kind,inputs,(y,),**kwargs);return y

def qc_split(g,key,x,label):
    assay=g.product(key+'_assay','Physically separate source analytical aliquot',K.QC)
    kept=g.product(key+'_retained','Unassayed retained preparation for downstream library work',K.MOLECULAR)
    g.react(key,label,S.SPLIT,((x,R.TEMPLATE),),(assay,kept));g.mark_terminal(Terminal(assay.id,Disposition.QC,'Source analytical aliquot never feeds retained library.'));return kept
from nextera import ADAPTOR_S5,ADAPTOR_S7,ME_RC,S5,S7
from illumina import P5,P7
SOURCE='Spektor et al.2019 DOI10.1101/gr.245399.118 published PMC6581052 Methods102–116; published Supplemental TableS3 (original rendered) and FigS1–S6 captions.'
# Tokens transcribed from rendered TableS3, not inferred from all cytosines.
WRITTEN_A='T/iMe-dC/GT/iMe-dC/GG/iMe-dC/AG/iMe-dC/GT/iMe-dC/AGATGTGTATAAGAGA/iMe-dC/AG'
WRITTEN_B='GT/iMe-dC/T/iMe-dC/GTGGG/iMe-dC/T/iMe-dC/GGAGATGTGTATAAGAGA/iMe-dC/AG'
@dataclass(frozen=True)
class ProtectedDNA:
    sequence:str
    methyl_c:tuple=()
    def __post_init__(self):
        if not self.sequence or set(self.sequence)-set('ACGT'):raise ValueError('selected exact DNA sequence required')
        if len(set(self.methyl_c))!=len(self.methyl_c) or any(type(i)!=int or not 0<=i<len(self.sequence) or self.sequence[i]!='C' for i in self.methyl_c):raise ValueError('5mC positions must identify cytidines')
    @classmethod
    def from_written(cls,written):
        seq='';marked=[]
        while written:
            if written.startswith('/iMe-dC/'):
                marked.append(len(seq));seq+='C';written=written[8:]
            else:seq+=written[0];written=written[1:]
        return cls(seq,tuple(marked))
    def converted(self):
        return ''.join('U' if b=='C' and i not in self.methyl_c else b for i,b in enumerate(self.sequence))
    def copied_reference(self):return revcomp(self.converted().replace('U','T'))
    def scene(self,name):
        sc=Scene();sc.strand(name,[Segment(name,self.sequence,note='5mC positions0based:'+str(self.methyl_c))],label='selected methylated DNA oligo; methylation is a covalent base modification');return sc

def adapter(which):
    if which not in {'a','b'}:raise ValueError('source adapter a/b required')
    return ProtectedDNA.from_written(WRITTEN_A if which=='a' else WRITTEN_B)

def index_primer(which):
    if which not in {'i5','i7'}:raise ValueError('source index type required')
    sc=Scene();sc.strand(which,[Segment('platform',P5 if which=='i5' else P7),Segment('index','[source '+which+' member; length unassigned]',placeholder=True,feature=feature('matac.'+which,'sample_index','unknown')),Segment('entry',S5 if which=='i5' else S7)]);return sc

def protocol_graph():
    g=ProtocolGraph('methyl-ATAC-seq',authority=SOURCE,representation='published_instance')
    x=g.input('cells','Source250000 cells; parental/DKO comparison is an experimental branch',K.MOLECULAR,primary=True)
    x=op(g,'nuclei',x,'Source modified-Omni lysis/wash nuclei, ice3min and500g recovery',S.EXTRACTION,reagent='SourceNP40/Tween/digitonin lysis andRBS wash')
    x=clean(g,'nuclei_recovery',x,'Source500g10min4C nuclei recovery')
    for key in ('a','b'):
        a=adapter(key);g.input('me_'+key,'Published methylated ME-'+key.upper()+' oligo; source tokens retained',K.OLIGO,state=a.scene('me_'+key),chemical_state=a)
    sc=Scene();sc.strand('me_rev',[Segment('me_rev',ME_RC)],mod5='phosphate');g.input('me_rev','Published Tn5MErev /5Phos/ oligo',K.OLIGO,state=sc)
    assembly=g.product('transposome','Source methylated A/B/Tn5 complex; occupancy/distribution unassembled',K.MOLECULAR)
    enzyme=g.input('tn5','Source Picelli2014-produced Tn5, production original unreviewed',K.CATALYST)
    g.react('load_tn5','Source adapter anneal94C2min→25C0.1C/s andTn5 loading',S.BINDING,tuple((g.materials[k],R.TEMPLATE if k=='me_a' else R.STRUCTURAL) for k in ('me_a','me_b','me_rev'))+((enzyme,R.CATALYST),),(assembly,))
    spike=g.input('lambda','Source0.01ng unmethylated lambda DNA conversion control',K.MOLECULAR)
    y=g.product('tagged','Tagmented chromatin/lambda with methylated adapter strands, other end combinations and physical Tn5 companions; no endpoint composition predicted',K.MOLECULAR)
    g.react('tagment','Source7ul methylated transposome37C30min',S.TAGMENTATION,((x,R.TEMPLATE),(assembly,R.STRUCTURAL),(enzyme,R.CATALYST),(spike,R.SPIKE_IN)),(y,))
    x=op(g,'proteinase',y,'SourceSDS/100ugPK55C30min inactivation',S.DEPROTEINIZATION,reagent='SourceSDS/proteinaseK mixture')
    x=clean(g,'precipitate',x,'Sourceisopropanol/GlycoBlue DNA recovery')
    x=clean(g,'size_select',x,'Source0.5x large-fragment bead removal then1.8x finalSPRI recovery')
    x=op(g,'fill',x,'SourceKlenowexo-minus37C30min end repair with5-methyl-dCTP',S.END_REPAIR,reagent='SourceKlenowexo-minus/dATP/dGTP/dTTP/5-mdCTP mixture',product='End-repaired DNA with synthetic5mC in newly copied C positions; original biological marks and companion material retained, exact gaps unassembled')
    x=clean(g,'fill_cleanup',x,'Source1.8xAMPure cleanup')
    x=qc_split(g,'saved_ten_percent',x,'Source ten percent saved for QC PCR, remaining preparation for bisulfite')
    x=op(g,'bisulfite',x,'SourceEZ DNA Methylation-Lightning D5030T conversion',reagent='Source historicalZymoD5030T, kit internals unreviewed',product='Bisulfite-converted recovered DNA population: selected unprotectedC→U, protectedC retained; other base chemistry/kit losses unassembled')
    for key in ('i5','i7'):g.input(key,'Published '+key+' PCR primer family; exact study member unassigned',K.OLIGO,state=index_primer(key))
    enzyme=g.input('pfu','SourcePfuTurboCx',K.CATALYST);amp=g.product('amp','PCR DNA descendants plus source converted templates/unextended primers/remnants; no synthetic mark inheritance to unmodified PCR dCTP asserted',K.MOLECULAR)
    g.react('pcr','Source immediate13-cyclePfuTurboCx PCR after conversion',S.PCR,((x,R.TEMPLATE),(g.materials['i5'],R.STRUCTURAL),(g.materials['i7'],R.STRUCTURAL),(enzyme,R.CATALYST)),(amp,),geometry=PCRGeometry.unavailable(Scene(),reason='Exact converted-end/product geometry and index members unassembled; source primer bodies known.'))
    x=clean(g,'final_cleanup',amp,'Source1.8xAMPure cleanup')
    x=qc_split(g,'qubit',x,'SourceQubit measurement and retained sequencing preparation')
    g.sequence(x,SequencingPlan(acquisition='unspecified',limitation='SourceNextSeq500single-end75bp; exact converted primer/cycle geometry unassembled. Computational exclusion of9nt atbothends does not physically trim DNA.'))
    g.add_source_item(SourceItem('pending','Picelli2014Tn5 production; Omni2017 protocol; historicalZymoD5030T andPfuTurboCx instructions','Direct named originals unreviewed within cap; full end/primer/acquisition assembly remains open.'))
    return g

def selftest():
    from checks import Check,Source,have
    from storage import sources_root
    c=Check();a=adapter('a');b=adapter('b')
    c('source methyl tokens retain canonical A/B nucleotide bodies',(a.sequence,b.sequence),(ADAPTOR_S5,ADAPTOR_S7))
    d=ProtectedDNA('ACCGC',(2,));c('selected5mC retained while unprotected C physically becomesU',d.converted(),'AUCGU')
    c('copying creates complementaryDNA not uracil DNA output',d.copied_reference(),'ACGAT')
    try:ProtectedDNA('AT',(1,));bad=False
    except ValueError:bad=True
    c('constructor rejects methyl mark on nonC',bad,True)
    g=protocol_graph();c('saved QC does not feed converted preparation','saved_ten_percent_assay' not in tuple(i.material for i in g.steps['bisulfite'].inputs),True)
    c('Qubit assayed fraction is QC',g.materials['qubit_assay'].kind,K.QC)
    p=Source(sources_root()/'methyl-atac-seq__10.1101+gr.245399.118'/'gr.245399.118_supp_supp_gr.245399.118_Supplemental_Table_S1-S3.pdf.txt','get_sources --slugs methyl-atac-seq__10.1101+gr.245399.118; originalS3render inspected')
    if have(c,p,label='sourceS3'):
        t=p.path.read_text();c('source explicitphosphate andmethyl-dC tokens','/5Phos/'+ME_RC in t and '/iMe-dC/' in t,True)
    c.report()
if __name__=='__main__':selftest()
