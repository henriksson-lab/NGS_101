from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Scene,Segment,revcomp,feature
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,SourceItem,SequencingPlan,Terminal,Disposition
from adapter_ligation import AdapterLigationGeometry
from pcr import PCRGeometry
class B:
    def __init__(self,title,authority,start,kind=K.BIOLOGICAL):
        self.g=ProtocolGraph(title,authority=authority,representation='published_instance');self.x=self.g.input('source',start,kind,primary=True)
    def stock(self,key,label,kind=K.CATALYST,role=R.CATALYST,**kw):return(self.g.input(key,label,kind,**kw),role)
    def cond(self,key,label):return self.stock(key,label,K.CONDITION,R.CONDITION)
    def op(self,key,label,kind,inputs=(),geometry=None,product='Source mixture including incomplete/unreacted/off-target classes; occupancy/recovery unassigned',reason=None):
        y=self.g.product(key+'_product',product,K.MOLECULAR);kw={}
        if kind==S.MAGIC:kw['reason']=reason or 'Positive source operation; complete occupied molecular inputs/products and historical formulation unassembled, not a missing-information conclusion.'
        if kind==S.GENERIC:kw['changes_molecule']=False
        if geometry is not None:kw['geometry']=geometry
        self.g.react(key,label,kind,((self.x,R.TEMPLATE),)+tuple(inputs),(y,),**kw);self.x=y;return y
    def clean(self,key,label,extra=()):
        y=self.g.product(key+'_kept','Actual source retained physical fraction; yield/carryover unassigned',K.MOLECULAR);w=self.g.product(key+'_waste','Actual separated source workup fraction',K.WASTE)
        self.g.react(key,label,S.CLEANUP,((self.x,R.TEMPLATE),self.cond(key+'_condition',label))+tuple(extra),(y,w));self.g.discard(w,'Actual workup fraction excluded; no perfect specificity or molecular yield claim.');self.x=y
    def split(self,key,label,selected,remainder):
        y=self.g.product(key+'_selected',selected,K.MOLECULAR);r=self.g.product(key+'_remainder',remainder,K.MOLECULAR)
        self.g.react(key,label,S.SPLIT,((self.x,R.TEMPLATE),),(y,r));self.g.mark_terminal(Terminal(r.id,Disposition.UNRESOLVED,remainder));self.x=y;return y

def unknown_scene():
    s=Scene();s.strand('unassembled',[Segment('source products unassembled','[unassembled]',placeholder=True)]);return s

AUTHORITY='Associated author manuscript PMC6475595 Methods107–129; no final-version equivalence or inherited predecessor/vendor whole recipe.'
def protocol_graph():
    b=B('CUT&RUN.ChIP native targeted release followed by physical histone immunoprecipitation',AUTHORITY,'One actual thawed native yeast nuclei aliquot ~5e8cell equivalents; ten-aliquot nuclei factory outside selected boundary');g=b.g
    b.op('con_a','Actual~200ulConA beads capture nuclei RT>=5min',S.MAGIC,(b.stock('cona','ActualConA beads',K.SELECTION,R.SELECTION),b.cond('conabuffer','20mMHEPES7.9/10mMKCl/1mMCaCl2/1mMMnCl2')))
    b.clean('nuclei_wash','Source magnetic decant then DigWash+2mMEDTA')
    b.op('flag_ab','Source mouseantiFLAG2h–overnight4C; extracted literal1:2 00 retained, ratio unassigned',S.BINDING,(b.stock('flag_stock','ActualmouseantiFLAG antibody',K.MOLECULAR,R.STRUCTURAL),b.cond('flag_cond','DigWash+EDTA')))
    b.clean('flag_wash','2xDigWash, no EDTA from this point')
    b.op('secondary','RabbitantiMouseIgG1:1001h4C',S.BINDING,(b.stock('secondary_stock','ActualrabbitantiMouseIgG',K.MOLECULAR,R.STRUCTURAL),))
    b.clean('secondary_wash','2xDigWash')
    b.op('pamn','SourcepA-MN360ugml,1:400 descriptor1h4C',S.BINDING,(b.stock('pamn_stock','ActualpA-MN molecular fusion stock; carried catalyst',K.MOLECULAR,R.STRUCTURAL),))
    b.clean('pamn_wash','3xDigWash; resuspend1.2ml/equilibrate0C')
    b.op('digest','2mMCaCl2 on ice5min; carried antibody-bound pA-MN performs digestion',S.MAGIC,(b.cond('ca','Actual2mMfinalCaCl2; full cut/end contexts unassembled'),))
    b.op('stop','Add1.2ml2xSTOP and first actual referenceDNA',S.QUENCHING,(b.stock('stop_rnase','Actual50ugmlRNaseA source2xSTOP descriptor'),b.stock('stop_glycogen','Actual40ugmlglycogen source2xSTOP descriptor',K.MOLECULAR,R.STRUCTURAL),b.stock('first_reference','Actual5–50pgml mostlymononucleosome DNA from FA-crosslinked MNase Drosophila chromatin; first reference molecular stock',K.MOLECULAR,R.SPIKE_IN),b.cond('stopbuffer','Source2xSTOP150mMNaCl/20mMEDTA/4mMEGTA; descriptor concentration convention preserved')))
    b.op('release','Source37C20min then13ugmlFLAGpeptide',S.MAGIC,(b.stock('flag_peptide','Actual13ugmlFLAGpeptide',K.MOLECULAR,R.STRUCTURAL),))
    sup=g.product('released_supernatant','Actual released supernatant; referenceDNA/capture carryover unassigned',K.MOLECULAR);res=g.product('nuclear_residue','Actual retained ConA/nuclei/bead residue',K.WASTE)
    g.react('release_separate','5min13000rpm4C then remove supernatant on magnet',S.SPLIT,((b.x,R.TEMPLATE),),(sup,res));g.discard(res,'Source residue excluded from serial ChIP production')
    assay=g.product('histone_aliquot','Actual selected H3K4me3 ChIP aliquot',K.MOLECULAR);inp=g.product('input_aliquot','Actual saved input at4C',K.MOLECULAR);other=g.product('other_aliquots','Actual other ChIP aliquots including IgG; source says7–8 and lists6 targets, allocation unassembled',K.MOLECULAR)
    g.react('split_released','Split released supernatant for physical ChIPs and saved input',S.SPLIT,((sup,R.TEMPLATE),),(assay,inp,other));g.mark_terminal(Terminal(other.id,Disposition.UNRESOLVED,'Other source ChIPs not fabricated as same vessel or included in selectedH3K4me3/input branches.'))
    b.x=assay;b.op('histone_ab','Actual H3K4me3 antibody1:100 overnight4C',S.BINDING,(b.stock('histone_stock','ActualH3K4me3 antibody',K.MOLECULAR,R.STRUCTURAL),))
    bound=g.product('histone_bound','Actual serial histone ChIP retained fraction',K.MOLECULAR);flow=g.product('histone_flow','Actual serial ChIP nonretained fraction',K.WASTE)
    g.react('histone_capture','Source20ulProteinA for rabbitpolyclonal; ProteinG mouseH3K36me3 alternative not merged',S.AFFINITY_CAPTURE,((b.x,R.TEMPLATE),b.stock('histone_beads','Actual20ulProteinA Dynabeads',K.SELECTION,R.SELECTION),b.cond('histone_capture_cond','4C15min')),(bound,flow));g.discard(flow,'Actual nonretained serialChIP fraction');b.x=bound;b.clean('histone_wash','OnceWash(noDigitonin)+.05percentTween20');selected=b.x
    for name,start in (('histone',selected),('input',inp)):
        b.x=start
        extra=[b.stock(name+'_pk','Actual.2mgmlProteinaseK'),b.stock(name+'_second_reference','Actual5–50pgmlS.pombe mostlymononucleosomeDNA after MNase; second reference stock',K.MOLECULAR,R.SPIKE_IN),b.cond(name+'_extraction_buffer','Source.1percentSDS/.2mgmlPK matched input; ChIP150mMNaCl/10mMEDTA/2mMEGTA')]
        if name=='histone':extra.append(b.stock('histone_extract_glycogen','Actual40ugmlglycogen extractionbuffer',K.MOLECULAR,R.STRUCTURAL))
        b.op(name+'_extract','Source50C20min DNA extraction; input bypasses antibody/bead capture',S.MAGIC,tuple(extra))
        b.clean(name+'_pci','Sourcebuffered25:24:1PCI then phase-lockCHCl3 re-extraction')
        b.clean(name+'_ethanol','Fresh2ul2mgmlglycogen/2–2.5volEtOH overnight-20C,13000rpm10min4C/100percentEtOHwash/airdry',extra=(b.stock(name+'_precip_carrier','Actual fresh4ug glycogen',K.MOLECULAR,R.STRUCTURAL),))
        b.op(name+'_elute','Source25ul.1xTE resuspend',S.GENERIC,(b.cond(name+'_te','1mMTris8/.1mMEDTA'),))
        b.op(name+'_kit','KAPA library prep/TruSeq adapters WITHOUT size-selection; AMPureXP removes unligated adapters',S.MAGIC,(b.stock(name+'_kit_stocks','Actual historical kit/adapters/primers; whole modified ends/nucleotide/catalyst compositions unassembled',K.MOLECULAR,R.STRUCTURAL),),reason='Positive cited historical library procedure; no invented genericKAPA or ordinaryCUTRUN adapter/ligation/PCR recipe.')
        lib=g.product(name+'_library','Actual '+name+' library for separate sequencing readout; no physical pool asserted',K.LIBRARY);g.react(name+'_load','SourceHiSeq2500 paired-end25cycles',S.GENERIC,((b.x,R.TEMPLATE),),(lib,),changes_molecule=False);g.sequence(lib,SequencingPlan(read_lengths={'R1':25,'R2':25},limitation='Whole adapters/primer landing/index mask unassembled. Actual two referenceDNA additions; cross-species read allocation computational.'))
    g.add_source_item(SourceItem('native_serial_scope','Main Methods123 native targetedMNase release followedbyphysicalhistoneChIP, with savedinput. FA applies to first referenceDNA source, not experimental yeast nuclei.','MouseantiFLAG extracted1:2 00 not silently normalized; source pAMN360ugml1:400 descriptor retained.'))
    return g

def selftest():
    from checks import Check,Source,have
    from storage import sources_root
    c=Check();g=protocol_graph();c('Serial ChIP uses actual released aliquot',g.steps['histone_ab'].inputs[0].material,'histone_aliquot');c('Input bypasses second antibody',g.steps['input_extract'].inputs[0].material,'input_aliquot');c('Nuclear residue excluded',g.terminals['nuclear_residue'].disposition,Disposition.DISCARD);c('Two physical reference additions',all(x in g.materials for x in ('first_reference','histone_second_reference','input_second_reference')),True);c('No size selection step fabricated',not any('size_selection' in x for x in g.steps),True);c('Read25paired',g.sequencing['histone_library'].read_lengths['R1'],25)
    p=sources_root()/Path(__file__).resolve().parents[1].name/'j.molcel.2018.10.046_PMC6475595.xml.txt'
    if have(c,Source(p,'https://doi.org/10.1016/j.molcel.2018.10.046')):
        t=p.read_text();c('Original serial supernatant split','divided into 7–8 aliquots for ChIP' in t,True);c('Original no size selection','without size-selection' in t,True);c('Original literal antiFLAG preserved','1:2 00' in t,True)
    c('Privatepartial unready',bool(g.publication_problems()),True);c.report()
if __name__=='__main__':selftest()
