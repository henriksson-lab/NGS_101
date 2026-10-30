"""BAsE-Seq source-owned workflow; whole genomic/product geometry unassembled."""
from dataclasses import dataclass
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent
sys.path[:0]=[str(HERE),str(HERE.parents[1]/'lib')]
import illumina as il
import nextera as nx
from chemdraw import Construct,Scene,Segment,feature,revcomp,panel,workflow_from_protocol_graph,workflow_panel
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,SequencingPlan,SourceItem
from pcr import PCRGeometry
from circular import CircularizationGeometry
from page import head,info
PAGE_FILE='base-seq.html'
SOURCE='Hong et al.2014 DOI10.1186/s13059-014-0517-9; published PMC4269956 Methods356; Additionalfile2,2July2014,Steps1–9 and oligo tablep3'
HBV_F='ACTGTTCAAGCCTCCAAGCTG'
HBV_R='AAAAGTTGCATGGTGCTGGTGA'
UNIA='ATGCTGCCTGCAGGATGTCG'
UNIB='GTTCGGTAACTGGAGCTGAC'
SBFI='CCTGCAGG'
CUSTOM_R2_HEAD='GAACCGCTCTTCCGATCT'
CUSTOM_R2=CUSTOM_R2_HEAD+revcomp(HBV_F)
FINAL_F=il.P5+'TAGATCGC'+nx.S5
FINAL_R=il.P7+CUSTOM_R2
FEATURE=feature('base_seq.genome_barcode','linked_read_barcode','degenerate')

def assignment_primer(*,internal_standard=False):
    return Construct([Segment('UniA/SbfI',UNIA),*([Segment('internal-standard insertion','AT')] if internal_standard else []),Segment('genome barcode','N'*20,'umi',placeholder=True,feature=FEATURE),Segment('HBV forward priming site',HBV_F)])

@dataclass(frozen=True)
class ModifiedPrimer:
    """Source physical modifications, not modifications of polymerase-made copies."""
    name:str
    sequence:str
    biotin_at:int|None=None
    ps_after:tuple[int,...]=()
    def __post_init__(self):
        if self.biotin_at is not None and (not 0<=self.biotin_at<len(self.sequence) or self.sequence[self.biotin_at]!='T'):
            raise ValueError('Source internal biotin must occupy a real dT')
        if len(set(self.ps_after))!=len(self.ps_after) or any(not 0<=i<len(self.sequence)-1 for i in self.ps_after):
            raise ValueError('Phosphorothioates join actual adjacent bases')
    def ordering(self):
        return ''.join(('/iBiodT/' if i==self.biotin_at else b)+('*' if i in self.ps_after else '') for i,b in enumerate(self.sequence))
    def rows(self):
        sc=Scene()
        sc.strand('primer',[Segment(self.name,self.sequence)],label=self.name)
        sc.note('primer','Physical ordering: '+self.ordering())
        return sc.rows()

UNI_A_BIO=ModifiedPrimer('Uni-A-SbfI-Bio',UNIA,UNIA.index('T',15),tuple(range(len(UNIA)-6,len(UNIA)-1)))
R2_PRIMER=ModifiedPrimer('Illumina_PE2 SVS.SeqPrimer',CUSTOM_R2,ps_after=(len(CUSTOM_R2)-3,len(CUSTOM_R2)-2))

def binding_scene():
    sc=Scene();sc.strand('forward',[Segment('HBV forward arm',HBV_F)],label='source HBV forward arm; not full library')
    sc.anneal('r2',[Segment('custom sequencing tail',CUSTOM_R2_HEAD),Segment('custom HBV-binding arm',revcomp(HBV_F))],to='forward',pair=('custom HBV-binding arm','HBV forward arm'),label='R2 primer body; two terminal PS bonds')
    sc.note('r2','This component binding relation does not establish all final library strands or the complete Nextera end geometry.')
    return sc

def protocol_graph():
    g=ProtocolGraph('BAsE-Seq — published HBV route',authority=SOURCE,representation='published_instance')
    scope=[]
    def inp(key,label,kind,**kw):
        scope.append(key);return g.input(key,label,kind,**kw)
    def step(key,label,kind,inputs,outlabel,**kw):
        out=g.product(key+'_mixture',outlabel,K.MOLECULAR)
        g.react(key,label,kind,inputs,(out,),**kw);scope.append(key);return out
    def cleanup(key,material,procedure):
        reagent=inp(key+'_procedure',procedure,K.SELECTION)
        kept=g.product(key,'Source recovered preparation; carryover/yield unquantified',K.MOLECULAR)
        waste=g.product(key+'_reject','Physically removed workup fractions; composition unassigned',K.WASTE)
        g.react(key+'_step',procedure,S.CLEANUP,((material,R.TEMPLATE),(reagent,R.SELECTION)),(kept,waste));scope.append(key+'_step');g.discard(waste,'Actual source physical separation; no perfect recovery/purity asserted.');return kept
    def pcr(key,material,f,r,label):
        poly=inp(key+'_polymerase','Source Long PCR Enzyme Mix' if key!='final_pcr' else 'Historical NexteraXT NPM',K.CATALYST)
        return step(key,label,S.PCR,((material,R.TEMPLATE),(f,R.STRUCTURAL),(r,R.STRUCTURAL),(poly,R.CATALYST)),'PCR product preparation plus residual template, primers, nucleotides and reaction products',geometry=PCRGeometry.unavailable(Scene(),reason='Full HBV template, heterogeneous deletion products and historical transposome/end assembly not yet constructed; source chemistry is available, not claimed absent.'))
    dna=inp('hbv','Source full-length HBV genomes; upstream extraction/qPCR outside graph',K.MOLECULAR,primary=True)
    f=inp('barcode_f','A-SbfI-ID20-HBV_For;20random genome-barcode bases',K.OLIGO,state=assignment_primer())
    r=inp('barcode_r','B-HBV_Rev',K.OLIGO,state=Construct([Segment('UniB',UNIB),Segment('HBV reverse arm',HBV_R)]))
    cur=pcr('assignment',dna,f,r,'Two-cycle barcode assignment; not a collision-free or one-genome-one-strand identity guarantee')
    exo=inp('assignment_exoi','ExonucleaseI37°C1h then98°C5min',K.CATALYST)
    cur=step('primer_digest','Digest excess assignment primers; products remain until cleanup',S.DIGESTION,((cur,R.TEMPLATE),(exo,R.CATALYST)),'Tagged DNA plus ExoI-resistant/unreacted material and soluble digestion products')
    cur=cleanup('assignment_recovered',cur,'Source1xAMPure cleanup')
    a=inp('unia_bio','Physical Uni-A-SbfI-Bio with internal dT-biotin and five PS bonds',K.OLIGO,state=UNI_A_BIO)
    b=inp('unib','Uni-B universal primer',K.OLIGO,state=Construct([Segment('UniB',UNIB)]))
    cur=pcr('pcr1',cur,a,b,'Source PCR1 with qPCR-selected log-linear cycle count; no error/chimera-free copying implied')
    cur=cleanup('pcr1_recovered',cur,'Source0.7xAMPure cleanup')
    cur=pcr('pcr2',cur,a,b,'Source PCR2 across eight wells; pool products; qPCR-selected cycles')
    cur=cleanup('pcr2_recovered',cur,'Pool eight reactions and0.7xAMPure cleanup')
    sbfi=inp('sbfi','SbfI-HF37°C1h then80°C20min',K.CATALYST)
    cur=step('sbfi_digest','Create barcode-proximal4nt3′overhang; cleavage fragments remain in mixture',S.DIGESTION,((cur,R.TEMPLATE),(sbfi,R.CATALYST)),'SbfI-cleaved tagged DNA plus released end fragments; source overhang inhibits ExoIII at proximal end')
    cur=cleanup('sbfi_recovered',cur,'Source1.8xAMPure cleanup')
    exo3=inp('exo3','ExonucleaseIII30°C; aliquots1/3/5/7/9min',K.CATALYST)
    cur=step('nested_deletion','Barcode-distal ExoIII deletion time course; cleavage products persist',S.DIGESTION,((cur,R.TEMPLATE),(exo3,R.CATALYST)),'Heterogeneous barcode-retaining deletion molecules and digestion products; damage/specificity unquantified')
    edta=inp('edta','Separate aliquots into EDTA;80°C15min;pool',K.CONDITION)
    cur=step('stop_exo3','Quench and heat-inactivate ExoIII aliquots then pool; not purification',S.QUENCHING,((cur,R.TEMPLATE),(edta,R.CONDITION)),'Pooled quenched deletion reaction, retaining soluble products')
    cur=cleanup('exo3_recovered',cur,'Source1.8xAMPure cleanup before S1')
    s1=inp('s1','S1nuclease roomtemperature30min',K.CATALYST)
    cur=step('s1_digest','Remove exposed single-stranded regions',S.DIGESTION,((cur,R.TEMPLATE),(s1,R.CATALYST)),'Source S1-treated deletion DNA plus digestion products; exact ends unassembled')
    stop=inp('s1_edta','SourceEDTA quench',K.CONDITION)
    cur=step('stop_s1','Quench S1; do not delete reaction products',S.QUENCHING,((cur,R.TEMPLATE),(stop,R.CONDITION)),'Quenched S1 mixture')
    cur=cleanup('s1_recovered',cur,'Source1.8xAMPure cleanup')
    beads=inp('streptavidin','Source streptavidin Dynabeads,roomtemperature15h',K.SELECTION)
    bound=g.product('bound_biotin','Magnet-retained biotinylated DNA/bead complex; exact capture purity unquantified',K.MOLECULAR)
    unbound=g.product('unbound_biotin','Source unbound fraction and washes',K.WASTE)
    g.react('capture_biotin','Bind source biotinylated DNA;magnet and washes',S.AFFINITY_CAPTURE,((cur,R.TEMPLATE),(beads,R.SELECTION)),(bound,unbound));scope.append('capture_biotin');g.discard(unbound,'Source physical magnetic separation and washes')
    cur=cleanup('bio_recovered',bound,'TwoPBI50°C30min dissociations;magnet each time;combine supernatants on MinElute column and elute')
    # Capture/release is a source physical workup, not copying a biotin into every DNA strand.
    repair=inp('end_repair','T4DNApolymerase+T4PNK+dNTPs20°C20min then75°C20min',K.CATALYST)
    cur=step('repair','Source blunt ending and5′phosphorylation',S.END_REPAIR,((cur,R.TEMPLATE),(repair,R.CATALYST)),'End-repaired DNA preparation; exact genomic end coordinates unassembled')
    cur=cleanup('repair_recovered',cur,'Source0.7xAMPure cleanup')
    ligase=inp('t4ligase','T4DNAligase10°C16h,≤100ng/45ml',K.CATALYST)
    closure=CircularizationGeometry.known('double-strand intramolecular ligation',statement='Published detailed Step7 uses T4DNA ligase on blunt5′phosphorylated dsDNA. Full heterogeneous genomic closure coordinates and strand/product population remain unassembled.')
    cur=step('circularize','Source dilute intramolecular dsDNA ligation; not CircLigase ssDNA closure',S.CIRCULARIZATION,((cur,R.TEMPLATE),(ligase,R.CATALYST)),'Intended dsDNA circles plus uncircularized DNA and possible intermolecular products; closure geometry/yield unassembled',geometry=closure)
    cur=cleanup('circles_precipitated',cur,'Split reaction;GlycoBlue/sodiumacetate/isopropanol precipitation;ethanol washes;combine pellets')
    enz=inp('circle_exonucleases','Lambdaexonuclease+ExonucleaseI37°C30min then80°C20min',K.CATALYST)
    cur=step('remove_linear','Digest accessible uncircularized DNA; no perfect exclusion inferred',S.DIGESTION,((cur,R.TEMPLATE),(enz,R.CATALYST)),'Circle-enriched mixture with residual DNA and digestion products')
    cur=cleanup('circle_recovered',cur,'Source1.8xAMPure cleanup thenSpeedVac concentration')
    xt=inp('nextera_xt','Historical NexteraXT transposome;source-linked kit guide/version/end chemistry not inspected',K.CATALYST)
    cur=step('tagment','Fragment and tag source circular preparation;NexteraXT',S.FRAGMENTATION,((cur,R.TEMPLATE),(xt,R.CATALYST)),'Tagmented circle-derived DNA fragments; exact duplicated target/end distribution unassembled')
    pf=inp('final_f','NexteraXT-PE1',K.OLIGO,state=Construct([Segment('sourcePCRforward',FINAL_F)]))
    pr=inp('final_r','Illumina-PE2_HBV-For.RevComp',K.OLIGO,state=Construct([Segment('sourcePCRreverse',FINAL_R)]))
    cur=pcr('final_pcr',cur,pf,pr,'Source NPM custom-primer PCR72°C3min fill stage;14cycles unless low input; source55°C annealing')
    cur=cleanup('library',cur,'Source0.65xAMPure size selection;KAPAqPCR quantification downstream')
    g.sequence(cur,SequencingPlan(acquired_roles=('Read 1','Read 2'),read_lengths={'Read 1':150,'Read 2':150},limitation='Published MiSeq PE150; source customR2 starts at20ntgenome barcode. Full physical run-primer/product/index geometry remains unassembled. No extra index or sample barcode inferred.'))
    g.declare_reviewed_scope(scope,SOURCE)
    g.add_source_item(SourceItem('pending_geometry','Full HBV product, source enzyme/end chemistry, historical NexteraXT guide and sequencing-primer declarations','Actual linked historical kit original needs version identification; all exact full-strand/substrate geometry not yet assembled. No source nondisclosure claim.'))
    return g

def render_page():
    g=protocol_graph()
    return head('BAsE-Seq — source-limited HBV model')+'<main class="wrap"><h1>BAsE-Seq</h1><p><a href="01_base-seq.html">Research notes</a></p>'+info('Private incomplete source-qualified model. Nested-deletion and circle populations are heterogeneous; no complete barcode specificity, perfect cleanup, chimera-free copying or full genomic geometry is claimed.')+panel(binding_scene().rows(),caption='Printed custom R2 binding subrelation')+workflow_panel(workflow_from_protocol_graph(g,'library'))+'</main>'

def selftest():
    from checks import Check
    from storage import sources_root
    c=Check();p=sources_root()/HERE.parent.name/'s13059-014-0517-9_supp_13059_2014_517_MOESM2_ESM.pdf.txt'
    if p.exists():
        txt=''.join(p.read_text().split()).upper()
        for name,written in [('assignmentF',assignment_primer().top()),('standardF',assignment_primer(internal_standard=True).top()),('assignmentR',UNIB+HBV_R),('UniAmodified',UNI_A_BIO.ordering()),('UniB',UNIB),('finalF',FINAL_F),('finalR',FINAL_R),('customR2modified',R2_PRIMER.ordering())]:
            c('Detailed original transcription '+name,written.upper() in txt,True)
    else:c.skip('Original transcription','External source unavailable')
    binding_scene();g=protocol_graph();c('Incomplete physical/source scope stays unpublished',g.publishable,False)
    c('No finished-page marker','data-protocol-complete=' in render_page(),False)
    c.report()
