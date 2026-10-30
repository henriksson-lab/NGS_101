"""Private source-qualified partial preparation; no complete molecular/absence claim."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Scene,Segment,revcomp,feature
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,SourceItem,SequencingPlan,Terminal,Disposition
from pcr import PCRGeometry
from illumina import P5,P7,TRUSEQ_READ1,TRUSEQ_READ2
class B:
    def __init__(self,title,authority,label,kind=K.MOLECULAR):
        self.g=ProtocolGraph(title,authority=authority,representation='published_instance');self.x=self.g.input('source',label,kind,primary=True)
    def stock(self,key,label,kind=K.MOLECULAR,role=R.STRUCTURAL,**kw):return(self.g.input(key,label,kind,**kw),role)
    def cat(self,key,label):return self.stock(key,label,K.CATALYST,R.CATALYST)
    def cond(self,key,label):return self.stock(key,label,K.CONDITION,R.CONDITION)
    def op(self,key,label,kind=S.MAGIC,inputs=(),product='Actual source preparation; whole occupied substrates/products/side products and conversion unassembled',**kw):
        y=self.g.product(key+'_product',product,K.MOLECULAR)
        if kind==S.MAGIC:kw.setdefault('reason','Positive source operation; whole occupied molecules, end distributions, detailed molecular products and conversion unassembled, not a source-absence claim.')
        if kind==S.GENERIC:kw['changes_molecule']=False
        self.g.react(key,label,kind,((self.x,R.TEMPLATE),)+tuple(inputs),(y,),**kw);self.x=y;return y
    def clean(self,key,label):
        y=self.g.product(key+'_kept','Actual retained physical source fraction, carryover/recovery unassigned',K.MOLECULAR);w=self.g.product(key+'_waste','Actual excluded workup fraction',K.WASTE)
        self.g.react(key,label,S.CLEANUP,((self.x,R.TEMPLATE),self.cond(key+'_condition',label)),(y,w));self.g.discard(w,'Actual source workup fraction excluded; no quantitative purification claim.');self.x=y;return y
    def aliquot(self,key,label,kept='Selected source aliquot',rest='Untransferred stock, amount/fate unassigned',qc=False):
        y=self.g.product(key+'_selected',kept,K.MOLECULAR);r=self.g.product(key+'_remainder',rest,K.MOLECULAR)
        self.g.react(key,label,S.SPLIT,((self.x,R.TEMPLATE),),(y,r))
        if qc:self.g.mark_terminal(Terminal(y.id,Disposition.QC,kept));self.x=r;return r
        self.g.mark_terminal(Terminal(r.id,Disposition.UNRESOLVED,rest));self.x=y;return y
    def pcr(self,key,label,forward,reverse,**kw):
        # Source one bulk preparation: produced components are copresent, not physical fractions.
        chemicals=kw.pop('primer_chemicals',{})
        stock=self.g.input(key+'_stock','Actual source prepared PCR reagent stock; quantitative internal recipe specified in step, components same tube',K.MOLECULAR)
        specs=[('f',K.OLIGO,forward,R.STRUCTURAL),('r',K.OLIGO,reverse,R.STRUCTURAL),('enzyme',K.CATALYST,None,R.CATALYST),('nt',K.MOLECULAR,None,R.STRUCTURAL),('buffer',K.CONDITION,None,R.CONDITION)]
        out=[];ins=[]
        for suffix,kind,state,role in specs:
            z=self.g.product(key+'_'+suffix,'Same prepared PCR stock '+suffix,kind,state=state,chemical_state=chemicals.get(suffix));out.append(z);ins.append((z,role))
        self.g.react(key+'_components','Copresent components of ONE source PCR stock, not extra additions',S.GENERIC,((stock,R.TEMPLATE),),tuple(out),changes_molecule=False)
        return self.op(key,label,S.PCR,ins,geometry=PCRGeometry.unavailable(unknown_scene(),reason='Whole occupied primer landing and heterogeneous amplicon products unassembled; source exact primers retained where assigned.'),**kw)
    def finish(self,key,label,plan):
        y=self.g.product(key,label,K.LIBRARY);self.g.react(key+'_load',label,S.GENERIC,((self.x,R.TEMPLATE),),(y,),changes_molecule=False);self.g.sequence(y,SequencingPlan(limitation=plan));return y

def unknown_scene():
    s=Scene();s.strand('unassembled',[Segment('whole source preparation unassembled','[unassembled]',placeholder=True)]);return s

def oligo(seq,name='whole source oligo',mod5=None,mod3=None):
    import re
    parts=[Segment(name+' '+str(i),part,placeholder=bool(re.search('[NRYKMSWBDHV]',part))) for i,part in enumerate(re.split('([NRYKMSWBDHV]+)',seq)) if part]
    s=Scene();s.strand('oligo',parts,mod5=mod5,mod3=mod3);return s

PROBE='GA'*10
SOURCE_OLIGO_EXCEPTIONS={'dig_probe':'Published main writes the whole probe as (GA)10 with /3Dig_N/, rather than a continuous20base string. The local source fixture independently extracts the printed repeat count and verifies expansion plus terminalDIG; no base repair or source-absence claim.'}
AUTHORITY='PublishedElsevier-sponsoredPMC8791277Methods606–632 pluspublishedTableS1mmc2PNCTR(NR)sourcecontext. Selectedmain-explicitPNCTRrepeatprobeRNAseqbranch, notMS ornewercompanionrecipe.'
def probe_scene():return oligo(PROBE,mod3='DIG')
def local_repeat_pair():
    s=Scene();s.strand('probe',[Segment('sourceGA10probe',PROBE)],mod3='DIG')
    s.anneal('rna',[Segment('conditionalRNArepeatwindow','UC'*10,inferred=True,note='Conditional20ntwindowclasswithinlargerPNCTRUCnrepeat; fullRNAregister/occupancyunassembled')],to='probe',pair=('conditionalRNArepeatwindow','sourceGA10probe'));return s

def protocol_graph():
    b=B('HyPro-seq selectedPNCTRrepeatRNA proximity preparation',AUTHORITY,'Sourcecellsin10cmdish; upstreamculture/geneticconditionsoutsidegraph',K.BIOLOGICAL);g=b.g
    b.op('dsp_fix','Source.5mgmlDSP/PBS30minRT',S.CROSS_LINKING,(b.stock('dsp','ActualsourceDSPcrosslinkermolecularstock'),b.cond('dsp_conditions','SourcePBS30minRT')))
    b.clean('fix_wash','Source3xPBS+20mMTris8/5mineach')
    b.op('permeabilize','Source70percentethanolRT1h',S.PERMEABILIZATION,(b.cond('ethanol','Source70percentethanol1hRT'),))
    b.op('equilibrate','Source2xSSC/10percentformamide5min',S.GENERIC,(b.cond('equil_buffer','ActualsourceSSC/formamide'),))
    probe=g.input('dig_probe','Source50nM(GA)10/3Dig_N_PNCTRrepeatprobe',K.OLIGO,state=probe_scene(),chemical_state={'source_ordering':'(GA)10/3Dig_N/','three_prime':'sourceDIGnon-nucleotideconjugate','purpose':'RNAhybridizationtargeting, notbarcodeUMI'})
    b.op('hybridize','SourceDIGprobehybridize37Covernight',S.HYBRIDIZATION,((probe,R.STRUCTURAL),b.cond('hyb_conditions','Source2xSSC/10percentformamide/10percentdextransulfate; fulloccupiedRNAtargetunassembled')))
    b.clean('hyb_wash','Source2xSSC10percentformamide37C30min/1xSSC15minRT')
    b.op('block','Source.8percentBSA/4xSSC/100UmlRNaseinhibitor30minRT',S.GENERIC,(b.stock('block_bsa','Actual.8percentBSAmolecularcarrier'),b.stock('block_inhib','ActualmurineRNaseinhibitorprotein'),b.cond('block_buffer','Source4xSSC30minRT')))
    b.op('hypro_bind','Source2.7ugmlHyProAPEX2-DIG10.3binding1hRT',S.BINDING,(b.stock('hypro','ActualpurifiedHyProfusionproteinstock, expression/purificationfactoryunassembled'),b.stock('bind_bsa','Actual.8percentBSAcomponent of freshHyProblockingbuffer'),b.cond('hypro_conditions','Source4xSSC/1hRT')))
    b.clean('hypro_wash','Source4xSSC/4xSSC.1percentTriton/4xSSC10mineach thenPBS5min; retainsprobe-associatedpreparation')
    b.op('label','Sourceequalvolume1mMbiotinphenol/.2mMH2O2addedsolution;1minRT',inputs=(b.stock('biotin_phenol','Actualsourcebiotinphenolmolecularsubstrate, addedsolution1mM notassertedfinalconcentration'),b.stock('peroxide','Actualsourcehydrogenperoxidemolecularsubstrate, addedsolution.2mM'),b.cond('label_conditions','SourceequalvolumePBSaddition/gentleagitation1min; originalretainedHyProenzyme, nofreshaddition')),
         product='Sourceproximitylabeled/unreactedRNA/protein/othermolecularclasses; noexactradius/RNAspecificity/labelingyieldclaim')
    b.op('quench','Source3quickrinses5mMTrolox/10mMsodiumascorbatePBS',S.QUENCHING,(b.stock('trolox','Actual5mMTroloxquencher'),b.stock('ascorbate','Actual10mMsodiumascorbatequencher'),b.cond('quench_rinses','Source3rinsesinPBS')))
    b.op('rna_lysis','SourcehighSDSRIPA1mlper10cmdish/ice10minscrape/ice10min',S.EXTRACTION,(b.stock('rna_lysis_inhib','Actual100UmlmurineRNaseinhibitor'),b.cond('rna_lysis_buffer','Source150mMNaCl/1mMEDTA/50mMTris8/1percentNP40/.5percentdeoxycholate/.5percentSDS/10mMascorbate/5mMTrolox/50mMDTT')))
    b.op('sonicate','SourceBioruptor5–10cycles30sON/OFFHIGHonice',S.FRAGMENTATION,(b.cond('sonicate_program','ActualsourceRNAbranchsonication'),))
    b.op('reverse_dsp','Source37C30mincrosslinkreversal',S.THERMAL_TREATMENT,(b.cond('reverse_conditions','Actual37C30min, sourceDTTalreadyinlysisbuffer'),))
    b.op('pre_capture_pk','Source400ugPK50C1hBEFOREtotalRNAextraction/capture',S.DEPROTEINIZATION,(b.cat('pk400','Actual400ugProteinaseK'),b.cond('pk_pre_program','Source50C1h')))
    b.clean('total_rna','Source3volTRIzolLS/chloroform; retainRNAaqueousphase/PureLinkminiprep')
    b.aliquot('rna_capture_input','Source150ulRNAaliquot equal to150ulC1suspension; actualmassunassigned','Actual150ulRNAcapturealiquot','UntransferredpurifiedtotalRNAstock')
    bound=g.product('rna_bound','SourceC1-retainedRNApreparation, free/offtarget/carryoverclassesunassigned',K.MOLECULAR);flow=g.product('rna_flow','ActualRNAcapturenonretainedfraction',K.MOLECULAR)
    beads=g.input('c1','ActualpreparedC1stock150ulsuspension; beads10ulper20ugRNA scaling, sourcewashfactoryunassembled',K.SELECTION,chemical_state={'source_scale':'10ulbeadsper20ugRNA','source_bead_preparation':'B&W3x/NaOHNaCl/NaOHNaClTween/NaCl/NaClTrisEDTATween washes','actual_RNA_mass':'unassigned'})
    g.react('rna_capture','Sourceequal150ulRNA+C1suspension4C2hrotation',S.AFFINITY_CAPTURE,((b.x,R.TEMPLATE),(beads,R.SELECTION),b.cond('capture_buffer','Source100mMNaCl/10mMTris7.5/1mMEDTA/.2percentTween')),(bound,flow));g.mark_terminal(Terminal(flow.id,Disposition.UNRESOLVED,'NonretainedRNAfractionnotimplicitlysequenced.'));b.x=bound
    b.clean('rna_bead_wash','SourceB&Wbufferthreewashes')
    b.op('onbead_pk','Source100ulRNAbeadreleasebuffer/.80ugulPK50C45min',S.DEPROTEINIZATION,(b.cat('pk80','Actualsource.80ugulPKin100ulbuffer'),b.cond('rna_release_buffer','Source100mMNaCl/50mMTris7.5/1mMEDTA/1percentSDS')))
    b.clean('captured_rna_extract','Source300ulTRIzolLS+80ulchloroform/aqueousRNA/ZymoRNAcleanconcentratorR1015')
    b.op('deplete','SourceNEBNextE6350maturerRNAdepletion',inputs=(b.stock('depletion_stock','ActualhistoricalNEBNextE6350reagentinventory, wholeformulation/internalfractionsunassembled'),))
    b.op('library_construct','SourceNEBNextUltraIIDirectionalRNAE7765/barcodedprimers',inputs=(b.stock('library_stock','ActualhistoricalE7765molecularconstructionstock; RT/NT/adapter/index/PCR/endformulationsunassembled'),),product='SourcebarcodedstrandedRNAsequencinglibrarypreparation; wholeproducts/endgeometry/readlandingunassembled')
    b.aliquot('quant','SourceNGSBIOPB71.15-01qPCRquantification analyticalamountunassigned',qc=True)
    own=b.aliquot('pool_input','Actualselectedindexedlibrarypoolcontribution, amountunassigned')
    others=g.input('other_libraries','Actualotherqualifiedbarcodedlibraryaliquots; sourceorigins/stockremaindersunassembled',K.MOLECULAR)
    pool=g.product('pool','ActualsourcebarcodedlibraryNGSpool',K.MOLECULAR);g.react('pool_libraries','Sourcephysicalpoolofbarcodedlibrariesbeforesequencing',S.POOL,((own,R.TEMPLATE),(others,R.TEMPLATE)),(pool,));b.x=pool
    b.aliquot('loading','ActualNextSeqloadingallocation, amountunassigned')
    b.finish('library','HyPro-seqselectedPNCTRrepeatRNAseqpool','PublishedNextSeq500v2.5/75cycleSE; actualbarcodedkitprimers/index/runlandingunassembled. HyProMSprotein60ulC1/3ml3aliquotroute isnotthisRNApurification or aNGSoutput.')
    # A true separate control preparation; its downstream RNAseq products are not silently merged.
    prior=b.x;b.x=g.input('infusion_cells','ActualseparateDSPfixed/permeabilizedpreparedcells forHyProinfusioncontrol; fullpreparationboundaryunassembled',K.BIOLOGICAL,primary=True)
    b.op('infusion_control','Sourcecontrol5.4ugmlHyProinPBS5minthenproximitylabeling',inputs=(b.stock('infusion_hypro','Actualseparate5.4ugmlHyProfusioncontrolstock'),b.stock('infusion_bp','Actualcontrolbiotinphenoladdedsolution1mM'),b.stock('infusion_h2o2','Actualcontrolperoxideaddedsolution.2mM'),b.cond('infusion_program','Source5minPBSinfusionthen1minbiotinylation; controlallocation/quenchandRNAseqworkupnotinstantiated')))
    g.mark_terminal(Terminal(b.x.id,Disposition.UNRESOLVED,'ActualindependentHyProinfusioncontrolstock; subsequentRNAseq/referenceanalysisoutsidepartialgraph.'));b.x=prior
    return g

def selftest():
    from checks import Check,Source,have
    from storage import sources_root
    c=Check();g=protocol_graph();local_repeat_pair()
    c('Probeiswhole20ntrepeat',PROBE,'GAGAGAGAGAGAGAGAGAGA')
    c('PKbeforecaptureactualorder',g.steps['total_rna'].inputs[0].material,'pre_capture_pk_product')
    c('NoMSproteinbeadamountsubstitution',g.materials['c1'].chemical_state['source_scale'],'10ulbeadsper20ugRNA')
    c('OnbeadPKusesactualwashedC1preparation',g.steps['onbead_pk'].inputs[0].material,'rna_bead_wash_kept')
    c('Actualotherindexedlibrarycontribution',len(g.steps['pool_libraries'].inputs),2)
    c('Controlstocknotreturnedtoproduction',g.terminals['infusion_control_product'].disposition,Disposition.UNRESOLVED)
    p=sources_root()/Path(__file__).resolve().parents[1].name/'j.molcel.2021.10.009_PMC8791277.xml.txt'
    if have(c,Source(p,'https://doi.org/10.1016/j.molcel.2021.10.009')):
        import re
        t=p.read_text();n=''.join(t.split()).replace('µ','μ');m=re.search(r'The\s+\(GA\)([0-9]+)\s+DNA oligonucleotide probe.*?/3Dig_N/',t,re.S);c('Originalwholemodifiedrepeatdescriptor',m is not None,True)
        if m:c('Wholeoligoexpansionfromoriginalrepeatcount',PROBE,'GA'*int(m.group(1)))
        c('Wholemodifiedsourceprobeclass','(GA)10' in t and '/3Dig_N/'in t,True)
        c('Sourceequalvolumeaddedsolution','equalvolumeof1×PBScontaining1mMbiotin-phenol' in n,True)
        c('RNAcapturesecondPKdose','0.80μg/μlproteinaseK'in n,True)
    c('Partialnotpublishedready',bool(g.publication_problems()),True);c.report()
if __name__=='__main__':selftest()
