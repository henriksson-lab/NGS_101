"""Protocol-owned source-qualified incomplete molecular graph."""
from dataclasses import dataclass, replace
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Scene, Segment, revcomp, feature
from protocol_graph import ProtocolGraph, MaterialKind as K, InputRole as R, StepKind as S, SourceItem, SequencingPlan
from pcr import PCRGeometry
from protocol_graph import Terminal,Disposition

def clean(g,key,x,label):
    condition=g.input(key+'_condition',label,K.CONDITION)
    kept=g.product(key+'_kept','Source recovered preparation; exact carryover unquantified',K.LIBRARY if x.kind is K.LIBRARY else K.MOLECULAR)
    rest=g.product(key+'_rest','Actual physically separated workup fractions; no perfect molecular partition inferred',K.WASTE)
    g.react(key,label,S.CLEANUP,((x,R.TEMPLATE),(condition,R.CONDITION)),(kept,rest));g.discard(rest,'Actual source workup separation.');return kept

def op(g,key,x,label,kind=S.MAGIC,*,reagent='source reagents; precise geometry unassembled',catalyst=False,product=None,reason='Source operation known; internal molecular geometry and remnant fates unassembled.',molecular_reagent=False,additional=(),changes_molecule=True):
    stock=g.input(key+'_reagent',reagent,K.CATALYST if catalyst else (K.MOLECULAR if molecular_reagent else K.CONDITION))
    out=g.product(key+'_product',product or label+'; accompanying physical reaction remnants retained',K.MOLECULAR)
    g.react(key,label,kind,((x,R.TEMPLATE),(stock,R.CATALYST if catalyst else (R.STRUCTURAL if molecular_reagent else R.CONDITION)))+tuple(additional),(out,),reason=reason,changes_molecule=changes_molecule);return out

def oligo_scene(name,sequence,*,mod5='',mod3=''):
    sc=Scene();sc.strand(name,[Segment(name,sequence)],mod5=mod5,mod3=mod3);return sc


SOURCE='Christy et al.2021 DOI10.1021/acs.biochem.1c00270; associated accepted manuscript PMC8256721 Methods126–142 and Fig1; final-VOR applicability and gated SI unreviewed.'
@dataclass(frozen=True)
class ProbedRNA:
    sequence:str
    crosslinks:tuple=()
    monoadducts:frozenset[int]=frozenset()
    def __post_init__(self):
        if not self.sequence or set(self.sequence)-set('ACGU'):raise ValueError('explicit native RNA required')
        for pair in self.crosslinks:
            if not isinstance(pair,tuple) or len(pair)!=2 or any(type(i)is not int or not 0<=i<len(self.sequence) for i in pair) or pair[0]>=pair[1]:raise ValueError('crosslinks connect two distinct ordered RNA2OH sites')
        if any(type(i)is not int or not 0<=i<len(self.sequence) for i in self.monoadducts):raise ValueError('monoadduct positions outside RNA')
    def selected_copy(self,path):
        # Caller supplies a conditional template path; no exact crosslink->deletion
        # endpoints, enzyme specificity or non-templated additions are predicted.
        if not path or any(type(i)is not int or not 0<=i<len(self.sequence) for i in path) or any(a<=b for a,b in zip(path,path[1:])):raise ValueError('explicit unique3to5template-coordinate path required')
        return ''.join({'A':'T','C':'G','G':'C','U':'A'}[self.sequence[i]] for i in path)

def graph_for_control(control='tbia',target='rnase_p'):
    if control not in {'tbia','ia','none'} or target not in {'rnase_p','p546','mbox','vs','groupii'}:raise ValueError('source target/control branch required')
    g=ProtocolGraph('SHAPE-JuMP — '+target+'/'+control,authority=SOURCE,representation='published_instance')
    x=g.input('template','Source synthetic PCR-generated DNA with5prime/3prime structure cassettes; exact sequenceTableS1 unreviewed/gated',K.MOLECULAR,primary=True)
    ntp=g.input('transcription_ntp','Actual2.5mM each NTP molecular contribution',K.MOLECULAR)
    tx_buffer=g.input('transcription_buffer','Actual40mMTris8/25mMMgCl2/2.5mMspermidine/.01percentTriton/10mMDTT buffer contribution; stock packaging unassembled',K.CONDITION)
    x=op(g,'transcribe',x,'Source T7 in-vitro transcription37°C4h',S.TRANSCRIPTION,catalyst=True,reagent='Actual source T7RNApolymerase catalyst',additional=((ntp,R.STRUCTURAL),(tx_buffer,R.CONDITION)),product='RNA transcripts plus original DNA template and transcription reagents')
    x=op(g,'dnase',x,'Source TURBO DNase treatment after transcription',S.DIGESTION,catalyst=True,reagent='source TURBO DNase',product='RNA plus DNA digestion products/enzymes; no physical removal before cleanup')
    x=clean(g,'rna_cleanup',x,'Source RNACleanXP purification/exchange into0.5xTE')
    x=op(g,'refold',x,'Source95°C1min,ice5min,37°C30min in target-specific foldingbuffer',S.THERMAL_TREATMENT,reagent='source target-specific HEPES/NaCl orKCl/Mg conditions, no universal fold predicted')
    selected=g.product('probe_5pmol','Actual5pmol foldedRNA input for selected chemical-control reaction',K.MOLECULAR)
    rest=g.product('folded_stock','Untransferred foldedRNA, if any; quantity and nextuse unassembled',K.MOLECULAR)
    g.react('probe_allocation','Source5pmol RNA allocation; no common-stock pooling of controls',S.SPLIT,((x,R.TEMPLATE),),(selected,rest))
    g.mark_terminal(Terminal(rest.id,Disposition.UNRESOLVED,'No source-unsupported amount or later fate assigned.'))
    x=selected
    probe_environment=g.input('probe_environment','Source37C15minTBIA or30minIA/DMSO, selected separately',K.CONDITION)
    x=op(g,'probe',x,'Source '+('TBIA37°C15min' if control=='tbia' else 'IA37°C30min' if control=='ia' else 'DMSO/no-reagent37°C30min')+';RNA5pmol/reagent10:1vol:vol',S.CHEMICAL_CONVERSION if control!='none' else S.GENERIC,reagent='source TBIA or IA40mM stock inDMSO' if control!='none' else 'source neatDMSO',molecular_reagent=True,additional=((probe_environment,R.CONDITION),),changes_molecule=control!='none',product=({'tbia':'TBIA-treated RNA population: selected2OH bridges, monoadducts and unreactedRNA plus reagent remnants','ia':'IA-treated RNA population: monoadducts and unreactedRNA plus reagent remnants; no TBIA bridge installed','none':'DMSO/no-reagent RNA control: no probe adduct installed; originalRNA and buffer retained'}[control]+'; background cDNA deletions remain possible and do not prove an RNA crosslink'))
    x=clean(g,'desalt',x,'Source G50 spin-column desalting')
    primer=g.input('rt_primer','Actual10pmol3prime-structure-cassette primer; original sequence not yet assembled',K.OLIGO)
    primed=g.product('primed_rna','Source RNA/primer annealing mixture; full duplex/occupancy unassembled',K.MOLECULAR)
    anneal_condition=g.input('primer_anneal_conditions','Actual95C1min thenice before RT components',K.CONDITION)
    g.react('primer_anneal','Source10pmol cassetteprimer annealing to selected RNA',S.HYBRIDIZATION,((x,R.TEMPLATE),(primer,R.STRUCTURAL),(anneal_condition,R.CONDITION)),(primed,))
    x=primed
    rt_nt=g.input('rt_nt','Actual200uM dNTP contribution for RT-C8',K.MOLECULAR)
    rt_buffer=g.input('rt_buffer','Actual20mMTris8.8/10mM(NH4)2SO4/10mMKCl/2mMMgCl2/.1percentTriton RT buffer; source packaging unassembled',K.CONDITION)
    rt=g.input('rt_c8','Source engineered RT-C8,65°C4h; source Mg2+ conditions, not RT521K/Mn2+ screening variant',K.CATALYST)
    y=g.product('copy','RNA templates with intact nucleotide inventory plus '+('conditional crosslink-jump/background-deletion cDNA' if control=='tbia' else 'control cDNA with possible background deletions/mutations')+' and reaction remnants',K.MOLECULAR)
    g.react('rt_jump','Source RT-C8 reverse transcription of '+('TBIA-treated RNA; selected crosslink-jump and background copying classes' if control=='tbia' else control+' control RNA; background deletion/mutation copying class'),S.REVERSE_TRANSCRIPTION,((x,R.TEMPLATE),(rt,R.CATALYST),(rt_nt,R.STRUCTURAL),(rt_buffer,R.CONDITION)),(y,))
    x=clean(g,'post_rt',y,'Source post-RT G50 desalting; no RNA destruction imported from separate Marathon screen')
    for stage,cycles in ((1,25),(2,10)):
        if stage==2:x=op(g,'dilute',x,'Source dilute firstPCR to0.5ng/µl',S.GENERIC,reagent='source dilution conditions',changes_molecule=False)
        used=g.product('pcr'+str(stage)+'_aliquot','Source'+('12µl of50µlRTreaction' if stage==1 else '1ng dilutedfirstPCRproduct')+' assignedtoPCR',K.MOLECULAR)
        remaining=g.product('pcr'+str(stage)+'_remaining','Source unallocated reaction preparation; exactlaterfateunassembled',K.MOLECULAR)
        g.react('allocate'+str(stage),'Source aliquot PCRinput; preserve physicalremainder',S.SPLIT,((x,R.TEMPLATE),),(used,remaining));g.mark_terminal(Terminal(remaining.id,Disposition.UNRESOLVED,'Actual unallocatedsource preparation retained; not assumedconsumed/discarded.'));x=used
        f=g.input('pcr'+str(stage)+'_f','Source'+str(stage)+' forward primer'+(' withIllumina adapter' if stage==1 else ' adapter-binding/indexed family'),K.OLIGO)
        r=g.input('pcr'+str(stage)+'_r','Source'+str(stage)+' reverse primer'+(' withIllumina adapter' if stage==1 else ' withTruSeq barcode; actual member unassigned'),K.OLIGO)
        nts=g.input('pcr'+str(stage)+'_nt','Actual PCR nucleotide functional class; composition/dose unreviewed, no RT or dUTP recipe inherited',K.MOLECULAR)
        buffer=g.input('pcr'+str(stage)+'_buffer','Actual PCR buffer functional class; composition/packaging unreviewed',K.CONDITION)
        cat=g.input('q5_'+str(stage),'Source Q5HotStartDNApolymerase',K.CATALYST);amp=g.product('pcr'+str(stage)+'_mix','Amplified cDNA population plus original templates and PCR remnants',K.LIBRARY if stage==2 else K.MOLECULAR)
        g.react('pcr'+str(stage),'Source'+str(cycles)+'cyclePCR'+('12µl of50µl RTreaction into50µlPCR' if stage==1 else ''),S.PCR,((x,R.TEMPLATE),(f,R.STRUCTURAL),(r,R.STRUCTURAL),(cat,R.CATALYST),(nts,R.STRUCTURAL),(buffer,R.CONDITION)),(amp,),geometry=PCRGeometry.unavailable(Scene(),reason='Actual cassettes/primer/index/end geometry in gatedSI and cited predecessor unassembled.'))
        x=clean(g,'pcr'+str(stage)+'_clean',amp,'Source1:1Ampure cleanup removes'+(' excessprimers' if stage==1 else ' shortproducts'))
    g.sequence(x,SequencingPlan(acquisition='unspecified',limitation='AssociatedAMMethods142 MiSeq:300ntkitP546/MBox/VS,500ntRNaseP,600ntgroupII. Kit sizes not silently assigned acquiredR1/R2 cycles; exact primer/run layout and final-VOR applicability unassembled.'))
    g.add_source_item(SourceItem('pending','Published main, actual gatedSI and library/cassette predecessors','OriginalSupplementaryTableS1/figures, Methods126/136 refs5/52; RT-C8 refs14/15/51 unreviewed.'))
    return g

def protocol_graph():return graph_for_control()
def selftest():
    from checks import Check,Source,have
    from storage import sources_root
    c=Check();r=ProbedRNA('ACGUAC',((1,4),));copy=r.selected_copy((5,4,1,0))
    c('selected jumped DNA omits caller-supplied interval only',copy,'GTGT')
    c('RT jump does not physically delete interveningRNA',r.sequence,'ACGUAC')
    c('physical2OHcrosslink survives copiedDNA observation',r.crosslinks,((1,4),))
    c('background control can yield same explicitly selected deletion',ProbedRNA('ACGUAC').selected_copy((5,4,1,0)),copy)
    try:r.selected_copy((5,1,4));bad=False
    except ValueError:bad=True
    c('copying path cannot reverse or duplicate template traversal',bad,True)
    g=protocol_graph()
    c('Actual NTPs are not hidden in T7 catalyst',g.materials['transcription_ntp'].kind,K.MOLECULAR)
    c('RT consumes already primed preparation with no fresh primer addition',tuple(p.material for p in g.steps['rt_jump'].inputs),('primed_rna','rt_c8','rt_nt','rt_buffer'))
    c('Probe enters selected5pmol RNA preparation',g.steps['probe'].inputs[0].material,'probe_5pmol')
    c('No-reagent control does not claim covalent transformation',graph_for_control('none').steps['probe'].changes_molecule,False)
    c('PCR dilution does not claim new covalent DNA',g.steps['dilute'].changes_molecule,False)
    c('all three actual chemical controls remain incomplete',all(graph_for_control(b).publication_problems() for b in ('tbia','ia','none')),True)
    p=Source(sources_root()/'shape-jump__10.1021+acs.biochem.1c00270'/'acs.biochem.1c00270_PMC8256721.xml.txt','get_sources --slugs shape-jump__10.1021+acs.biochem.1c00270;associatedAMMethods136')
    if have(c,p,label='selectedsource operation guard'):c('source prints4h RT-C8 rather than borrowedMarathon hydrolysis','at 65 °C for 4 h' in p.path.read_text(),True)
    c.report()
if __name__=='__main__':selftest()
