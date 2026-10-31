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
    kept=g.product(key+'_kept','Source recovered material; no perfect purity or quantitative partition assumed',K.LIBRARY if x.kind is K.LIBRARY else K.MOLECULAR)
    waste=g.product(key+'_waste','Actually separated cleanup fractions with unquantified molecular carryover',K.WASTE)
    g.react(key,label,S.CLEANUP,((x,R.TEMPLATE),(condition,R.CONDITION)),(kept,waste));g.discard(waste,'Source physical separation, not perfect molecular classification.');return kept

def op(g,key,x,label,kind=S.MAGIC,*,reagent='source mixture',catalyst=False,product=None,state=None,chemical_state=None,additional=()):
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
    g.react(key,label,kind,tuple(inputs)+tuple(additional),(y,),**kwargs);return y

def qc_split(g,key,x,label):
    assay=g.product(key+'_assay','Physically separate source analytical aliquot',K.QC)
    kept=g.product(key+'_retained','Unassayed retained preparation for downstream library work',K.LIBRARY if x.kind is K.LIBRARY else K.MOLECULAR)
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
    reverse_stock=g.materials['me_rev']
    rev_a=g.product('me_rev_a','Actual MErev aliquot for independent A duplex preparation; amount unassigned',K.OLIGO,state=sc)
    rev_b=g.product('me_rev_b','Actual MErev aliquot for independent B duplex preparation; amount unassigned',K.OLIGO,state=sc)
    g.react('reverse_allocations','Allocate source reverse-oligo contributions before separate A/B annealing',S.SPLIT,((reverse_stock,R.TEMPLATE),),(rev_a,rev_b))
    duplexes=[]
    for key,reverse in (('a',rev_a),('b',rev_b)):
        known=adapter(key);duplex_scene=Scene()
        duplex_scene.strand('top',[Segment('prefix',known.sequence[:-len(ME_RC)]),Segment('ME',known.sequence[-len(ME_RC):])])
        duplex_scene.anneal('bottom',[Segment('MErev',ME_RC)],to='top',pair=('MErev','ME'),mod5='phosphate')
        duplex=g.product('duplex_'+key,'Actual source methylated '+key.upper()+' duplex stock; paired nucleotide relation known, source methyl tokens preserved separately',K.ADAPTER,state=duplex_scene,chemical_state=known)
        anneal=g.input('anneal_'+key,'Actual94C2min then0.1C/sec ramp to25C',K.CONDITION)
        g.react('anneal_'+key+'_step','Independent source '+key.upper()+'/MErev duplex annealing',S.HYBRIDIZATION,((g.materials['me_'+key],R.TEMPLATE),(reverse,R.STRUCTURAL),(anneal,R.CONDITION)),(duplex,))
        duplexes.append(duplex)
    glycerol=g.input('glycerol','Actual50ul glycerol contribution to25ul20uM A plus25ul20uM B duplex stocks',K.MOLECULAR)
    assembly=g.product('transposome','Actual loaded-source stock:25ulA+25ulB+50ulglycerol mixed with equal-volume Tn5; occupied stoichiometry/free companions unassembled',K.MOLECULAR)
    enzyme=g.input('tn5','Actual source Picelli-produced Tn5 stock atA2803.0; production original unreviewed',K.CATALYST)
    g.react('load_tn5','Load actual mixed source A/B duplex/glycerol preparation with Tn5; no exact occupancy asserted',S.BINDING,((duplexes[0],R.TEMPLATE),(duplexes[1],R.STRUCTURAL),(glycerol,R.STRUCTURAL),(enzyme,R.CATALYST)),(assembly,))
    selected=g.product('loaded_7ul','Actual7ul source loaded-transposome contribution',K.MOLECULAR)
    remainder=g.product('loaded_stock','Untransferred source loaded stock; later use/actual amount unassigned',K.MOLECULAR)
    g.react('loaded_allocation','Actual7ul library-reaction contribution selected from prepared loaded stock',S.SPLIT,((assembly,R.TEMPLATE),),(selected,remainder))
    g.mark_terminal(Terminal(remainder.id,Disposition.UNRESOLVED,'Prepared factory stock retained separately; no assumed disposal.'))
    loaded_enzyme=g.product('loaded_enzyme','Catalytic Tn5 class co-present in actual7ul loaded stock',K.CATALYST)
    loaded_adapters=g.product('loaded_adapters','Methylated A/B donor classes co-present in same7ul; whole occupied complex geometry unassembled',K.ADAPTER)
    loaded_other=g.product('loaded_other','Other actual7ul loaded-stock material, including unreacted/buffer companions',K.MOLECULAR)
    g.react('loaded_components','Expose same actual loaded-stock components, not new free Tn5/adaptor doses',S.GENERIC,((selected,R.STRUCTURAL),),(loaded_enzyme,loaded_adapters,loaded_other),changes_molecule=False)
    spike=g.input('lambda','Source0.01ng unmethylated lambda DNA conversion control',K.MOLECULAR)
    y=g.product('tagged','Tagmented chromatin/lambda with methylated adapter strands, other end combinations and physical Tn5 companions; no endpoint composition predicted',K.MOLECULAR)
    g.react('tagment','Source7ul methylated transposome37C30min',S.TAGMENTATION,((x,R.TEMPLATE),(loaded_adapters,R.STRUCTURAL),(loaded_enzyme,R.CATALYST),(loaded_other,R.STRUCTURAL),(spike,R.SPIKE_IN)),(y,))
    sds=g.input('sds_buffer','Actual5volumesSDSlysisbuffer(100mMTris7.4/50mMNaCl/10mMEDTA/.5percentSDS), source totalvolume unassigned',K.MOLECULAR)
    x=op(g,'proteinase',y,'SourceSDS/100ugPK55C30min inactivation',S.DEPROTEINIZATION,reagent='Actual100ugProteinaseK',catalyst=True,additional=((sds,R.STRUCTURAL),))
    x=clean(g,'precipitate',x,'Sourceisopropanol/GlycoBlue DNA recovery')
    x=clean(g,'size_select',x,'Source0.5x large-fragment bead removal then1.8x finalSPRI recovery')
    fill_nt=g.input('fill_nt','Actual0.5mM each dATP,dGTP,dTTP,5-methyl-dCTP molecular contribution',K.MOLECULAR)
    fill_buffer=g.input('fill_buffer','Actual1xNEBBuffer2',K.CONDITION)
    x=op(g,'fill',x,'SourceKlenowexo-minus37C30min end repair with5-methyl-dCTP',S.END_REPAIR,reagent='Actual5UKlenowexo-minus',catalyst=True,additional=((fill_nt,R.STRUCTURAL),(fill_buffer,R.CONDITION)),product='End-repaired DNA with synthetic5mC in newly copied C positions; original biological marks and companion material retained, exact gaps unassembled')
    x=clean(g,'fill_cleanup',x,'Source1.8xAMPure cleanup')
    x=qc_split(g,'saved_ten_percent',x,'Source ten percent saved for QC PCR, remaining preparation for bisulfite')
    x=op(g,'bisulfite',x,'SourceEZ DNA Methylation-Lightning D5030T conversion',reagent='Source historicalZymoD5030T, kit internals unreviewed',product='Bisulfite-converted recovered DNA population: selected unprotectedC→U, protectedC retained; other base chemistry/kit losses unassembled')
    for key in ('i5','i7'):g.input(key,'Published '+key+' PCR primer family; exact study member unassigned',K.OLIGO,state=index_primer(key))
    pcr_nt=g.input('pcr_nt','Actual PCR nucleotide functional class; composition unreviewed, no repair5mdCTP inherited',K.MOLECULAR)
    pcr_buffer=g.input('pcr_buffer','Actual PfuTurboCx reaction-buffer class; recipe unreviewed',K.CONDITION)
    enzyme=g.input('pfu','SourcePfuTurboCx',K.CATALYST);amp=g.product('amp','PCR DNA descendants plus source converted templates/unextended primers/remnants; no synthetic mark inheritance to PCR copies asserted',K.LIBRARY)
    g.react('pcr','Source immediate13-cyclePfuTurboCx PCR after conversion',S.PCR,((x,R.TEMPLATE),(g.materials['i5'],R.STRUCTURAL),(g.materials['i7'],R.STRUCTURAL),(enzyme,R.CATALYST),(pcr_nt,R.STRUCTURAL),(pcr_buffer,R.CONDITION)),(amp,),geometry=PCRGeometry.unavailable(Scene(),reason='Exact converted-end/product geometry and index members unassembled; source primer bodies known.'))
    x=clean(g,'final_cleanup',amp,'Source1.8xAMPure cleanup')
    x=qc_split(g,'qubit',x,'SourceQubit measurement and retained sequencing preparation')
    g.sequence(x,SequencingPlan(acquired_roles=('Read 1',),read_lengths={'Read 1':75},limitation='SourceNextSeq500single-end75bp; exact converted primer/cycle geometry unassembled. Computational exclusion of9nt atbothends does not physically trim DNA.'))
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
    c('A and B receive distinct pre-anneal reverse allocations',tuple(g.steps['reverse_allocations'].outputs),('me_rev_a','me_rev_b'))
    c('Actual7ul comes from prepared loaded factory stock',g.steps['loaded_allocation'].inputs[0].material,'transposome')
    c('Tagmentation catalyst derives from same7ul, not fresh Tn5',g.producer_of('loaded_enzyme').id,'loaded_components')
    c('Original free Tn5 not supplied again to tagmentation',any(i.material=='tn5' for i in g.steps['tagment'].inputs),False)
    c('Fill uses actual protected nucleotide contribution',any(i.material=='fill_nt' for i in g.steps['fill'].inputs),True)
    p=Source(sources_root()/'methyl-atac-seq__10.1101+gr.245399.118'/'gr.245399.118_supp_supp_gr.245399.118_Supplemental_Table_S1-S3.pdf.txt','get_sources --slugs methyl-atac-seq__10.1101+gr.245399.118; originalS3render inspected')
    if have(c,p,label='sourceS3'):
        t=p.path.read_text();c('source explicitphosphate andmethyl-dC tokens','/5Phos/'+ME_RC in t and '/iMe-dC/' in t,True)
    c.report()
if __name__=='__main__':selftest()
