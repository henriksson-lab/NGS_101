"""Protocol-owned source-qualified incomplete molecular graph."""
from dataclasses import dataclass, replace
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Scene, Segment, revcomp, feature
from protocol_graph import ProtocolGraph, MaterialKind as K, InputRole as R, StepKind as S, SourceItem, SequencingPlan
from pcr import PCRGeometry

def clean(g,key,x,label):
    condition=g.input(key+'_condition',label,K.CONDITION)
    kept=g.product(key+'_kept','Source recovered preparation; exact carryover unquantified',K.MOLECULAR)
    rest=g.product(key+'_rest','Actual physically separated workup fractions; no perfect molecular partition inferred',K.WASTE)
    g.react(key,label,S.CLEANUP,((x,R.TEMPLATE),(condition,R.CONDITION)),(kept,rest));g.discard(rest,'Actual source workup separation.');return kept

def op(g,key,x,label,kind=S.MAGIC,*,reagent='source reagents; precise geometry unassembled',catalyst=False,product=None,reason='Source operation known; internal molecular geometry and remnant fates unassembled.'):
    stock=g.input(key+'_reagent',reagent,K.CATALYST if catalyst else K.CONDITION)
    out=g.product(key+'_product',product or label+'; accompanying physical reaction remnants retained',K.MOLECULAR)
    g.react(key,label,kind,((x,R.TEMPLATE),(stock,R.CATALYST if catalyst else R.CONDITION)),(out,),reason=reason);return out

def oligo_scene(name,sequence,*,mod5='',mod3=''):
    sc=Scene();sc.strand(name,[Segment(name,sequence)],mod5=mod5,mod3=mod3);return sc

from protocol_graph import Terminal,Disposition
SOURCE='Brannan et al. DOI10.1038/s41592-021-01128-0; associated accepted manuscript PMC8148648 Methods124–202, final-VOR procedural applicability unverified; published oligoTableS11 workbook identity confirmed by manifestMD5.'
@dataclass(frozen=True)
class EditedRNA:
    original:str
    edited_positions:frozenset[int]
    def __post_init__(self):
        if not self.original or set(self.original)-set('ACGU') or any(type(i)is not int or not 0<=i<len(self.original) or self.original[i]!='C' for i in self.edited_positions):raise ValueError('APOBEC-selected editing needs supplied actual C-origin RNA positions')
    @property
    def physical_rna(self):return ''.join('U' if i in self.edited_positions else b for i,b in enumerate(self.original))
    def selected_first_dna(self):return ''.join({'A':'T','C':'G','G':'C','U':'A'}[b] for b in self.physical_rna[::-1])

def handoff(g,key,x,label):return op(g,key,x,label,reagent='source named historical kit/procedure; internal molecular recipe and remnant fates unassembled')
def acquire(g,x,limitation):g.sequence(x,SequencingPlan(acquisition='unspecified',limitation=limitation+' Complete molecular end/primer/index/run geometry and final-VOR applicability unassembled.'))

def graph_for_branch(branch='bulk_kapa'):
    if branch not in {'bulk_kapa','bulk_truseq','nanopore_cdna','pacbio_cdna','single_cell'}:raise ValueError('source-defined library branch required')
    g=ProtocolGraph('STAMP — '+branch,authority=SOURCE,representation='published_instance')
    x=g.input('cells','Source cells expressing selected RBP/ribosomal-protein–APOBEC1 fusion; editing sites supplied per molecular example, not predicted from binding',K.MOLECULAR,primary=True)
    x=op(g,'editing',x,'Source in-vivo APOBEC1-mediated RNA C-to-U editing; expression duration and fusion choice branch-dependent',S.CHEMICAL_CONVERSION,reagent='source expressed RBP/ribosome–APOBEC1 and cellular substrates',product='Selected U-edited and unedited RNA plus original cell material; RNA editing is distinct from reference-oriented T calls')
    if branch=='single_cell':
        x=handoff(g,'cells_prepare',x,'Source72h1µg/ml doxycycline, TrypLE, count/resuspend1000cells/µl in0.04%BSA/PBS')
        x=handoff(g,'gems',x,'Source Chromium3primev3 GEM partition, cell lysis and barcoded reverse transcription')
        x=handoff(g,'amplify_cdna',x,'Source cDNA amplification; mixedRBFOX2:TIA1 includes FeaturecDNAprimers1')
        expression=g.product('expression_preparation','Source gene-expression preparation from cDNA workup; exact molecular partition unassembled',K.MOLECULAR)
        featureprep=g.product('saved_feature_supernatant','Source saved Ampure cleanup supernatant for3primeFeatureBarcode libraries; composition unassembled',K.MOLECULAR)
        companion=g.product('cleanup_companion','Other physical cleanup material/remnants, exact disposition unassembled',K.MOLECULAR)
        condition=g.input('ampure','Source post-cDNA Ampure cleanup; exact historical geometry unreviewed',K.CONDITION)
        g.react('branch_cleanup','Source gene-expression continuation and SAVED feature supernatant',S.CLEANUP,((x,R.TEMPLATE),(condition,R.CONDITION)),(expression,featureprep,companion));g.mark_terminal(Terminal(companion.id,Disposition.UNRESOLVED,'No ideal recovery or source-unsupported disposal assigned.'))
        expression=op(g,'fragment',expression,'Source gene-expression fragmentation',S.FRAGMENTATION,reagent='source Chromium3primev3 preparation')
        expression=op(g,'end_repair',expression,'Source gene-expression end repair',S.END_REPAIR,catalyst=True,reagent='source historical end-repair enzyme mixture')
        expression=op(g,'a_tail',expression,'Source gene-expression A-tailing',S.TAILING,catalyst=True,reagent='source historical tailing enzymes/nucleotides')
        expression=handoff(g,'expression_adapt',expression,'Source5prime-adapter and sample-index attachment; exact ligation geometry unassembled')
        featureprep=handoff(g,'feature_amplify',featureprep,'Source saved-supernatant amplification using Feature and TemplateSwitchOligo primers')
        featureprep=handoff(g,'feature_index',featureprep,'Source feature-library indexing')
        acquire(g,expression,'AssociatedAMMethods202:3primepolyA library sequenced on NovaSeq6000.')
        acquire(g,featureprep,'AssociatedAMMethods202:3primefeature library sequenced on NovaSeq6000.')
    else:
        x=handoff(g,'extraction',x,'Source TRIzol/DirectZol RNA extraction'+(' withDNase treatment in stated long-read preparation' if branch in {'nanopore_cdna','pacbio_cdna'} else ' following manufacturer procedure'))
        if branch in {'bulk_kapa','bulk_truseq'}:
            x=handoff(g,'library',x,'Source500ngRNA '+('KAPARNAHyperPrep withRiboErase' if branch=='bulk_kapa' else 'TruSeqStrandedmRNA')+' library; source calls totalRNAseq but kit names kept distinct')
            acquire(g,x,'AssociatedAMMethods144:bulkRNAseq single-end100nt; platform not borrowed from separate eCLIP experiment.')
        elif branch=='nanopore_cdna':
            x=handoff(g,'polya',x,'Source20µgRNA NEBE7490L polyA magnetic-resin selection')
            x=handoff(g,'cdna',x,'Source100ngpolyA RNA SQK-DCS109 DIRECT cDNA library, not nativeRNA sequencing')
            acquire(g,x,'AssociatedAMMethods176–178:PromethION FLO-PRO002 flowcells approximately48h, real-time high-accuracy Guppy.')
        else:
            x=handoff(g,'isoseq',x,'Source300ngtotalRNA IsoSeqExpress cDNA preparation')
            f=g.input('forward','Source IsoSeq forward primer',K.OLIGO,state=oligo_scene('forward','GGCAATGAAGTCGCAGGGTTG'))
            r=g.input('reverse','Source IsoSeq reverse primer',K.OLIGO,state=oligo_scene('reverse','AAGCAGTGGTATCAACGCAGAG'))
            cat=g.input('pcr_enzyme','Source IsoSeqExpress amplification reagent, edition-specific internal details unreviewed',K.CATALYST)
            amp=g.product('amplified','Source13–15cycle amplified dsDNA plus original cDNA/reagent remnants',K.MOLECULAR)
            g.react('pcr','Source IsoSeq13–15cycle PCR',S.PCR,((x,R.TEMPLATE),(f,R.STRUCTURAL),(r,R.STRUCTURAL),(cat,R.CATALYST)),(amp,),geometry=PCRGeometry.unavailable(Scene(),reason='Exact full IsoSeq primer/template/strand geometry unassembled.'))
            x=handoff(g,'smrtbell',amp,'Source SMRTbellExpressTemplatePrep2.0 and separate BarcodedOverhangAdapterKit, per-sample barcoding')
            x=handoff(g,'pool',x,'Source equimolar pooling; initial barcode counts used to rebalance later pooled sequencing')
            acquire(g,x,'AssociatedAMMethods188–192:SequelII, SMRTcell8M,2.0chemistry/2.1polymerase,2hpre-extension/30hmovie; circular consensus is DNA sequencing, not directRNA.')
    g.add_source_item(SourceItem('pending','Actual historical library originals and main final-VOR reconciliation','Named TruSeq,KAPARiboErase,Chromiumv3,SQK-DCS109,IsoSeqExpress,SMRTbell2.0 vendor procedures unreviewed; eCLIP comparator not imported into STAMP.'))
    return g

def protocol_graph():return graph_for_branch()
def selftest():
    from checks import Check
    c=Check();r=EditedRNA('ACCG',frozenset({1}));first=r.selected_first_dna()
    c('selected C-to-U changes physical RNA',r.physical_rna,'AUCG')
    c('complementary firstDNA distinct from RNA-oriented report',(first,revcomp(first)),('CGAT','ATCG'))
    c('RNA editing retains original transcript identity',r.original,'ACCG')
    try:EditedRNA('ACCG',frozenset({0}));bad=False
    except ValueError:bad=True
    c('cytidine deamination cannot silently edit adenosine',bad,True)
    c('all five source branches retain explicit publication limits',all(graph_for_branch(b).publication_problems() for b in ('bulk_kapa','bulk_truseq','nanopore_cdna','pacbio_cdna','single_cell')),True)
    g=graph_for_branch('single_cell');c('feature-library path consumes saved cleanup supernatant',g.steps['feature_amplify'].inputs[0].material,'saved_feature_supernatant')
    from checks import Source,have
    from storage import sources_root
    source=Source(sources_root()/'stamp__10.1038+s41592-021-01128-0'/'s41592-021-01128-0_PMC8148648.xml.txt',"python3 tools/get_sources.py --slugs stamp__10.1038+s41592-021-01128-0; inspected original selected chemistry sections")
    if have(c,source,label="local original transcription guard"):
        c('associated main transcribes IsoSeq forward','GGCAATGAAGTCGCAGGGTTG' in source.path.read_text(),True)
        c('associated main transcribes IsoSeq reverse','AAGCAGTGGTATCAACGCAGAG' in source.path.read_text(),True)
    c.report()
if __name__=='__main__':selftest()
