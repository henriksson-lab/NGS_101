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

from adapter_ligation import AdapterLigationGeometry
A='GGCTTCTTCTTGCTCTTAGGTAGTAGGTTC'
STEM=A[:20]
B_STEM=revcomp(STEM)
B_END='AAAA'+'NN'
SOURCE_OLIGO_EXCEPTIONS={'oligo_a':'Whole5P A is printed in originalSI1FigS1a image, not extracted as text. The whole30base ordering/5P and independent20bpstem fixture are local; originalimage peer verification required, not sourceabsence.'}
AUTHORITY='PublishedPMC12775064mainMethods204/Results84–94 plusactualSI1FigS1ap1. ThirdMOESM8SourceDataheadersnonchemistryconsumedadminslot; noTable4/ONTrevisioninheritance.'
def a_scene():return oligo(A,mod5='P')
def b_parts():
    return [Segment('B5tailunassembled','[B5tail]',placeholder=True,note='Remaining5tailbody/lengthunassembled, notsourceabsence'),Segment('sourceBstem',B_STEM),Segment('sourceB3primingA4','AAAA'),Segment('sourceBterminalN2','NN',placeholder=True,note='Sourceprimingclass, notUMI/index')]
def b_scene():
    s=Scene();s.strand('b',b_parts(),label='PartialsourceBstrand; known20bpstem/6nt3end, 5tailunassembled');return s

def adapter_scene(with_rna=False):
    s=Scene();s.strand('a',[Segment('sourceA5stem',STEM),Segment('sourceA3tail',A[20:])],mod5='P')
    s.anneal('b',b_parts(),to='a',pair=('sourceBstem','sourceA5stem'),unpaired=('B5tailunassembled',))
    s.ligation_end('a','sourceA5stem',kind='splinted',polarity='5-prime',owner='adapter',text='Sourceknown5P A localadapterend; wholeRNA/B5tail/sealedproductunassembled')
    if with_rna:
        s.anneal('rna',[Segment('conditionalRNAend','NNUUUU',placeholder=True,inferred=True,note='INFERRED conditionalpolyU-RNA3OHclass; completeRNA/ends/occupationunassembled')],to='b',pair=('conditionalRNAend','sourceBterminalN2'),mod3='OH')
    return s

def protocol_graph():
    b=B('DRAP3R selectedtotalRNApolyUdirectRNApreparation',AUTHORITY,'Actual20–25ugtotalRNA; cellculture/organism/perturbationoutsideboundary');g=b.g
    polya=g.product('polya_bound','ActualfirstDynabeadsoligodTretainedpolyAfraction, separateDRScomparatorpossible',K.MOLECULAR);non=g.product('non_polya','ActualfirstpulldownSUPERNATANTnonpolyAfraction',K.MOLECULAR)
    g.react('polya_subtract','SourceDynabeadsmRNApurification retainfirstsupernatant',S.AFFINITY_CAPTURE,((b.x,R.TEMPLATE),b.stock('polya_beads','ActualsourceDynabeadsmRNApurificationstock',K.SELECTION,R.SELECTION)),(polya,non));g.mark_terminal(Terminal(polya.id,Disposition.UNRESOLVED,'RetainedpolyAstocknotdiscardedbydefault; separate300ngDRSrouteoutsidegraph.'));b.x=non
    b.clean('non_polya_precip','SourceethanolprecipitateactualnonpolyAsupernatantRNA')
    b.aliquot('ribominus_input','Source10ugselectedRNAforRiboMinus','Actual10ugnonpolyARNA','UntransferrednonpolyARN Astock')
    b.op('ribominus','SourceRiboMinusTranscriptomeIsolationKitdepletion',inputs=(b.stock('ribominus_stock','ActualRiboMinusRNAdepletionpreparedstocks; internalcapture/otherfractionstatesunassembled'),))
    b.clean('depleted_precip','SourceethanolprecipitateafterrRNAdepletion')
    b.aliquot('qubit','ActualQubitanalyticalquantification amountunassigned',qc=True)
    rna=b.aliquot('library_input','Source500ngselecteddepletedRNAforDRAP3R','Actual500ngRNA','UntransferredpolyA/rRNAdepletedRNAstock')
    a=g.input('oligo_a','Sourcewhole5P OligoA',K.OLIGO,state=a_scene(),chemical_state={'source_ordering':'/5PHOS/'+A})
    bb=g.input('oligo_b','SourceOligoB-TNknownpartialordering; whole5tailunassembled',K.OLIGO,state=b_scene(),chemical_state={'known_stem':B_STEM,'three_prime_terminal':'AAAANN','source_reverse_reading':'3primeNNAAAA5prime','remaining_5tail':'unassembled'})
    pair=g.product('adapter_stock','SourcepreannealedcustomRTAstock100ul, remainingBtailunassembled',K.ADAPTER,state=adapter_scene())
    g.react('anneal','Source1.4uMeachA/Bin100ul95C2min/.1Csec to4C',S.HYBRIDIZATION,((a,R.TEMPLATE),(bb,R.STRUCTURAL),b.cond('anneal_buffer','Source10mMTris7.5/50mMNaCl')),(pair,))
    prior=b.x;b.x=pair;adapter=b.aliquot('adapter_input','ActualcustomRTAligationaliquot, volumeunassigned','SelectedRTAaliquot','Untransferred100ulRTApreparationstock');b.x=prior
    ligstock=g.input('first_ligation_stock','Aggregate historical RNA004 first-ligation reagent inventory; actual packaging, enzyme identity and complete formulation unassembled',K.MOLECULAR);enz=g.product('first_ligase','First-ligation inventory enzyme component; identity and packaging unassembled',K.CATALYST);other=g.product('first_lig_other','Other first-ligation inventory components; formulation and packaging unassembled',K.MOLECULAR)
    g.react('first_stock_components','Inventory decomposition for bookkeeping; actual first-ligation packaging/addition schedule unassembled',S.GENERIC,((ligstock,R.TEMPLATE),),(enz,other),changes_molecule=False)
    b.op('first_ligation','SourcecustomRTAfirstligation; RCSomitted',S.LIGATION,((adapter,R.STRUCTURAL),(enz,R.CATALYST),(other,R.STRUCTURAL)),geometry=AdapterLigationGeometry.known('splinted',adapter_scene(with_rna=True),statement='SourceknownRTA5PA/20bpstem/3AAAANNBinterface; inferredconditionalRNA3OHpolyUclassonly, wholeRNA/B5tail/sealedproductunassembled.'),product='SourceRNA-customRTA/unreacted/offtargetclasses; completeBtail/substrate/sealedjunctionunassembled, noRCSaddition')
    rtstock=g.input('rt_stock','Aggregate source SSIII reaction reagent inventory; actual packaging, nucleotide formulation and conditions unassembled',K.MOLECULAR);rte=g.product('rt_enzyme','SSIII inventory enzyme component; packaging unassigned',K.CATALYST);rtnt=g.product('rt_nt','RT inventory nucleotide component; packaging unassigned',K.MOLECULAR);rtbuf=g.product('rt_buffer','RT inventory conditions; packaging unassigned',K.CONDITION)
    g.react('rt_stock_components','RT inventory decomposition for bookkeeping, not evidence of one prepared master mix; no duplicate supply',S.GENERIC,((rtstock,R.TEMPLATE),),(rte,rtnt,rtbuf),changes_molecule=False)
    tmp=g.product('rt_template','SamepreparedRNA/RTAtemplatecomponent, wholeoccupiedmoleculesunassembled',K.MOLECULAR);pr=g.product('retained_b_primer','INFERRED copresentRTA-B3OHprimerassignmentfromnamedRTA/localpair, notnewfreeprimer',K.OLIGO,state=b_scene(),inferred=True)
    g.react('rt_components','CopresentRNA/RTAcomponents, notphysicalfractions orfreshBaddition',S.GENERIC,((b.x,R.TEMPLATE),),(tmp,pr),changes_molecule=False);b.x=tmp
    b.op('rt','SourceSSIIIprotectivecDNAsynthesis, programmeunassembled',S.REVERSE_TRANSCRIPTION,((pr,R.STRUCTURAL),(rte,R.CATALYST),(rtnt,R.STRUCTURAL),(rtbuf,R.CONDITION)),product='SourceRNA/cDNAhybrids plusunprotectedRNA/unreactedRTA; noRNAsequencing-to-cDNAreadoutsubstitution')
    b.op('rnase_a','SourceRNaseA43ugml/.4MNaCl37C15min afterRT',S.DIGESTION,(b.cat('rnase_a_stock','Actual43ugmlRNaseA'),b.cond('rnase_conditions','Source.4MNaCl37C15min')),product='SourceprotectedRNA/hybrids andotherdigestionclasses; unadaptedssRNAremovalnotassumedcomplete')
    b.op('pk','SourcePK37C15min/70C15minRNaseinactivation',S.DEPROTEINIZATION,(b.cat('pk_stock','ActualsourceProteinaseK, doseunassembled'),b.cond('pk_program','Source37C15min then70C15min')))
    b.op('rla','SourceRLAadapterligation/standardintermediateworkupboundary, detailedscheduleunassembled',inputs=(b.stock('rla_stock','ActualsourceRLA/motorligationpreparedmolecularstock; fullends/enzyme/NTformulationunassembled'),),product='SourceRNA-adaptedmotorlibrarypreparation; wholeRLAendgeometry/productsunassembled')
    b.clean('final_cleanup','SourceAMPureXPpurifications2xconcentration; finalrecoveryboundarywithintermediatecleanupplacementsunassembled')
    b.aliquot('loading','ActualRAflowcellloadingallocation, mass/volumeunassigned')
    b.finish('library','DRAP3RselectedpolyURNA nanoporelibrary','PublishedRNA-specificRAflowcell/MinIONMk1b24h. MinKNOWqscorefilterdisabled/minreadretention20nt arecomputationalreadretention, notphysicalRNA20ntsizeselection. WholeRNA/RTA/RLAreadlandingunassembled.')
    return g

def selftest():
    from checks import Check,Source,have
    from storage import sources_root
    c=Check();g=protocol_graph();adapter_scene(True)
    c('SourcewholeA30fixture',A,'GGCTTCTTCTTGCTCTTAGGTAGTAGGTTC')
    c('Source20bpstemwrittenbottom3to5',revcomp(B_STEM),STEM)
    c('SourceB3endnotUMI',B_END,'AAAANN')
    c('Firstsupernatantisproduction',g.steps['non_polya_precip'].inputs[0].material,'non_polya')
    c('PolyAstockretainedseparately',g.terminals['polya_bound'].disposition,Disposition.UNRESOLVED)
    c('RCSnotaddedback',all('rcs'not in x for x in g.materials),True)
    c('No freshBprimeratRT',g.steps['rt'].inputs[1].material,'retained_b_primer')
    c('Bprimerassignmentmarkedderived',g.materials['retained_b_primer'].inferred,True)
    c('RNaseafteractualRTmixture',g.steps['rnase_a'].inputs[0].material,'rt_product')
    p=sources_root()/Path(__file__).resolve().parents[1].name/'s41467-025-68230-1_PMC12775064.xml.txt'
    if have(c,Source(p,'https://doi.org/10.1038/s41467-025-68230-1')):
        t=''.join(p.read_text().split()).replace('µ','μ');c('Literalfirstsupernatant','retainedthesupernatantafterthefirstpulldown'in t,True);c('SourceRNaseproteinandsalt','43μg/mlfinalconcentration'in t and'0.4MNaCl'in t,True)
    c('Partialnotpublishedready',bool(g.publication_problems()),True);c.report()
if __name__=='__main__':selftest()
