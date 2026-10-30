from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Construct, Scene, Segment, feature, workflow_from_protocol_graph, workflow_panel
from protocol_graph import ProtocolGraph, MaterialKind as K, InputRole as R, StepKind as S, SourceItem, SequencingPlan
from page import head,info
PAGE_FILE='chem.html'

def operations(g):
    def op(key,x,label,reagent=None,state=None,kind=S.GENERIC,reagent_kind=K.CONDITION,reagent_role=R.CONDITION):
        inputs=[(x,R.TEMPLATE)]
        if reagent:
            z=g.input(key+'_reagent',reagent,reagent_kind);inputs.append((z,reagent_role))
        y=g.product(key+'_product',label+'; preparation includes unremoved reaction components',K.MOLECULAR,state=state)
        g.react(key,label,kind,inputs,(y,));return y
    def clean(key,x,label):
        z=g.input(key+'_selection',label,K.SELECTION)
        y=g.product(key+'_kept','Physically retained source fraction: '+label,K.MOLECULAR)
        w=g.product(key+'_removed','Physically removed fractions: '+label,K.WASTE)
        g.react(key,label,S.CLEANUP,((x,R.TEMPLATE),(z,R.SELECTION)),(y,w));g.discard(w,'Source physical workup only; no perfect specificity or purity assertion.');return y
    return op,clean

def end(g,x,label):g.sequence(x,SequencingPlan(acquisition='unspecified',limitation=label))
def incomplete(g):
    g.declare_reviewed_scope(tuple(g.materials)+tuple(g.steps),SOURCE)
    g.add_source_item(SourceItem('unassembled','Unreviewed linked originals and full molecular geometry',LIMIT));return g

from dataclasses import dataclass
import nextera as nx
TITLE='Apt-seq — original polyadenylated affinity aptamers and transcripts'
SOURCE='Published PMC5811598 Methods85–99 and Fig1; published SI MOESM1 Tables1–2 (original pp3–4)'
LIMIT='Private partial model. Source DNA aptamer copying differs from RNA reverse transcription. Published B-before-A-tail notation remains a constrained nucleotide class; no biotin modification inferred. Drop-seq linked original/protocol and complete bead,product and read-primer geometry remain unassembled.'
APT_HANDLE='AAGGAGCAGCGTGGAGGATA'
TD05_CORE='AACACCGTGGAGGATAGTTCGGTGGCTGTTCAGGGTCTCCTCCCGGTG'
APT_RIGHT='TTAGGGTGTGTCGTCGTGGT'
SMART='AAGCAGTGGTATCAACGCAGAGT'
# Source-local printed extension is not replaced by a platform-wide recipe.
APT_EXTENSION=nx.S7
@dataclass(frozen=True)
class AptamerPopulation:
    folded:bool
    copied:bool=False
    def scene(self):
        sc=Scene()
        sc.strand('aptamer',[Segment('source left flank',APT_HANDLE),Segment('TD05 affinity sequence',TD05_CORE,feature=feature('aptseq.td05','feature_barcode','fixed')),Segment('source right flank',APT_RIGHT),Segment('printed B before A tail','B',placeholder=True,note='Source IUPAC B is C/G/T;not biotin.'),Segment('source polyA tail','A'*21)],label='DNA TD05 example;source folded population' if self.folded else 'DNA TD05 source strand;fold state no longer asserted')
        sc.note('aptamer','One published panel member;other source aptamers differ. No binding specificity/yield or predicted fold supplied.')
        if self.copied:sc.strand('copy',[Segment('copied aptamer component','[cell-labelled DNA copy;full strand unassembled]',placeholder=True)],label='DNA copy made during Maxima reaction;not RNA→DNA conversion')
        sc.aptamer_population=self;return sc

def protocol_graph():
    g=ProtocolGraph(TITLE,authority=SOURCE);op,clean=operations(g)
    apt=g.input('aptamers','Source synthetic DNA aptamer panel;TD05 sequence example',K.OLIGO,state=AptamerPopulation(False).scene())
    folded=op('fold',apt,'0.5uM aptamers inPBS/Mg/glucose:95C5min thenice10min','Source folding conditions',state=AptamerPopulation(True).scene())
    cells=g.input('cells','Source mixed Ramos/3T3 cell suspension',K.BIOLOGICAL,primary=True)
    stained=g.product('stained','Cells with affinity-bound aptamers plus unbound aptamers and buffer',K.MOLECULAR,state=AptamerPopulation(True).scene())
    g.react('stain','31.25nM each aptamer;30min onice',S.GENERIC,((cells,R.TEMPLATE),(folded,R.STRUCTURAL)),(stained,))
    x=clean('five_washes',stained,'Five cold binding-buffer centrifugation/resuspension washes;retain cell pellet')
    x=op('dilute',x,'Dilute inicePBS/0.1%BSA tosource cell concentration')
    bead=g.input('bead_primer','Historical Drop-seq bead/polyT/cell identifier/UMI components;whole strand unassembled',K.OLIGO)
    drops=g.product('droplets','Source cell/bead/lysis droplets;DNA aptamers andcellRNA coexist',K.MOLECULAR)
    g.react('encapsulate','Source Drop-seq encapsulation andcelllysis',S.GENERIC,((x,R.TEMPLATE),(bead,R.STRUCTURAL)),(drops,))
    x=op('capture',drops,'PolyA-bearing DNA aptamers andRNA bindbeadpolyT;source emulsion/bead workup detailsunreviewed')
    tso=g.input('tso','Published mixed DNA/RNA template-switch oligo',K.OLIGO,state=tso_inventory())
    rt=g.input('maxima','MaximaHMinus RT;RNA-dependentandDNA-dependent copying classes distinct',K.CATALYST)
    copied=g.product('copied','Source barcoded cDNA–RNA hybrids AND double-strandedaptamerDNA;freeoligos/reagents remain',K.MOLECULAR,state=AptamerPopulation(False,True).scene())
    g.react('copy','Source Maxima reaction withtemplate switching efficientforcappedmRNA,rare onaptamers',S.REVERSE_TRANSCRIPTION,((x,R.TEMPLATE),(bead,R.STRUCTURAL),(tso,R.STRUCTURAL),(rt,R.CATALYST)),(copied,))
    primer=g.input('aptamer_primer','Published TC/TD/TE Nextera-tailed aptamer primer',K.OLIGO,state=Construct([Segment('source Nextera overhang',APT_EXTENSION),Segment('source aptamer binding flank',APT_HANDLE)]))
    smart=g.input('smart','Published DropSeq smartPCRprimer',K.OLIGO,state=Construct([Segment('smartPCR',SMART)]))
    product=g.product('mixed_amp','Co-amplified transcript andaptamerDNA population;stillone physical mixture',K.MOLECULAR)
    g.react('preamplify','13cycle jointPCR;aptamer-specificprimers1:20 relative tosmartPCRprimer',S.GENERIC,((copied,R.TEMPLATE),(primer,R.STRUCTURAL),(smart,R.STRUCTURAL)),(product,))
    x=op('nextera',product,'Source Nextera tagmentation ofmixedamplifiedsample;transcripts~500bp,sourceexpectsrare aptameraction because short','Historical Nextera reagents;no perfect aptamer exclusion asserted')
    x=op('downstream',x,'Continue source-referenced Drop-seq library processing;exact later recipe unreviewed','Historical referenced Drop-seq indexing/cleanup components')
    end(g,x,'Published MiSeqReagentKitv3 paired-end;aptamersandtranscripts sequencedjointly,computationallyseparated. Exactindex/runlength/readprimerlandingsunassembled.')
    return incomplete(g)
@dataclass(frozen=True)
class MixedSugarTSO:
    dna:str=SMART+'GAAT'
    rna:str='GGG'
    def scene(self):
        sc=Scene();sc.strand('tso',[Segment('DNA template-switch part',self.dna),Segment('three terminal riboguanosines',self.rna,note='Source aswrittenrGrGrG;these areRNA sugars,notdeoxyG')],label=self.dna+'rGrGrG (mixedDNA/RNA oligo)')
        sc.physical_tso=self;return sc
def tso_inventory():return MixedSugarTSO().scene()

def render_page():return head(TITLE)+'<main class="wrap"><h1>'+TITLE+'</h1>'+info(LIMIT)+'<p><a href="01_apt.html">Research notes</a></p>'+workflow_panel(workflow_from_protocol_graph(protocol_graph(),'downstream_product'))+'</main>'
def selftest():
    from checks import Check
    from storage import sources_root
    c=Check();g=protocol_graph()
    c('Source aptamer is DNA before Maxima copying',g.materials['aptamers'].state.aptamer_population.copied,False)
    c('Fold assertion ends aftercopying',g.materials['copied'].state.aptamer_population.folded,False)
    c('Copied DNA aptamer retained as distinct class',g.materials['copied'].state.aptamer_population.copied,True)
    c('Actual TSO state retains terminal ribonucleotide polymer',g.materials['tso'].state.physical_tso.rna,'GGG')
    c('Joint PCR precedes mixedsample tagmentation',g.steps['nextera'].inputs[0].material,'mixed_amp')
    c('No invented physical aptamer/transcript split','split' in g.steps,False)
    c('Source-limited graph private',g.publishable,False)
    p=sources_root()/Path(__file__).resolve().parents[1].name/'s41598-018-21153-y_supp_41598_2018_21153_MOESM1_ESM.pdf.txt'
    if p.exists():
        t=p.read_text();c('Published TD05core andNextera-tailedprimer',TD05_CORE in t and APT_EXTENSION+APT_HANDLE in t,True);c('Mixed-sugarTSOnotation',SMART+'GAATrGrGrG' in t,True)
    else:c.skip('PublishedSI transcription','Externaloriginalabsent')
    c.report()
