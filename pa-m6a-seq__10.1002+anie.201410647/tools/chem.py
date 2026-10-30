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
        state=kw.pop('state',None)
        y=self.g.product(key+'_product',product,K.MOLECULAR,state=state)
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

AUTHORITY='AssociatedEnglishAMNIHMS676610/PMC4396828core140–150/Scheme1text; publishedEnglishmainHTTP403/actualSIblocked. Oneoriginalread, noGermanedition/preprint equivalence.'
def source_rna_end():
    s=Scene();s.strand('rna',[Segment('sourceprocessedRNAclass','[RNA body unassembled]',placeholder=True,note='Sourceintended5P/3OHligation-readyRNAclass; fullsequence/enddistributionandprocessingyieldunassigned')],mod5='P',mod3='OH');return s

def protocol_graph():
    b=B('PA-m6A-seq associatedEnglishsource preparation',AUTHORITY,'ActualHeLacells; upstreamcultureconditionsoutsidegraph',K.BIOLOGICAL);g=b.g
    b.op('label_4su','Source4SUmetabolicincorporationinHeLaRNA',inputs=(b.stock('four_su','Actual4-thiouridinenucleosidemolecularsupply, dose/timeunassembled'),b.cond('label_conditions','Sourcegrowthmediumlabeling; metabolicphosphorylation/incorporationsitesunassembled')),product='SourceRNAcellularpopulations including4SU/m6Apositionsunassigned; no4SU-to-Cchemicalconversion')
    b.op('extract_rna','SourceRNApreparationbeforeoligodTselection, detailedstocksunassembled',S.MAGIC,(b.stock('rna_extract_stock','ActualRNAextractionpreparationmolecularstock; internalworkup/fractionsunassembled'),))
    mrna=g.product('mrna_bound','SourceoligodT-capturedfull-length4SU-mRNApreparation, yield/lengthdistributionunassigned',K.MOLECULAR);flow=g.product('polya_flow','ActualnonretainedRNA/workupfraction',K.MOLECULAR)
    g.react('polya_capture','SourceoligodT-conjugatedmagneticbeadspurify4SU-mRNA',S.AFFINITY_CAPTURE,((b.x,R.TEMPLATE),b.stock('oligodt_beads','ActualsourceoligodTmagneticbeads; wholecapturearm/formulationunassembled',K.SELECTION,R.SELECTION)),(mrna,flow));g.mark_terminal(Terminal(flow.id,Disposition.UNRESOLVED,'NonretainedRNAfractionnotimplicitlysequenced.'));b.x=mrna
    b.op('mrna_recover','Sourcepurifiedfull-lengthmRNAbeforeIP, exactbeadelution/workupunassembled',inputs=(b.stock('mrna_recovery_stock','ActualsourceRNArecoverypreparation, internalstocksunassembled'),))
    bound=g.product('m6a_ip_bound','Sourceanti-m6Aretainedfull-lengthRNA/antibodypreparation beforeUV',K.MOLECULAR);unbound=g.product('m6a_ip_unbound','ActualIPnonretainedRNA/workupfraction',K.MOLECULAR)
    g.react('m6a_ip','Sourcefull-lengthmRNAanti-m6AimmunoprecipitationBEFOREUV',S.AFFINITY_CAPTURE,((b.x,R.TEMPLATE),b.stock('anti_m6a','Actualsourceanti-m6ApolyclonalrabbitIgGpreparedselectionstock; support/dose/formulationunassembled',K.SELECTION,R.SELECTION),b.cond('ip_conditions','ActualIPconditionsunassembled, nofragmentationbeforeIP')),(bound,unbound));g.mark_terminal(Terminal(unbound.id,Disposition.UNRESOLVED,'Nonretainedm6AIPfractionnotsequencinginputbyassumption.'));b.x=bound
    b.op('uv','Source365nmUVcrosslinkAFTERanti-m6AIP',S.CROSS_LINKING,(b.cond('uv_365','Source365nmUV, doseunassembled'),),product='Source4SU/antibodycrosslinkedandunreactedRNAcomplexclasses; exactsites/yieldunassigned')
    b.op('rnase_t1','SourceRNaseT1trimmingaround30ntAFTERcrosslink',S.DIGESTION,(b.cat('t1','ActualsourceRNaseT1, dose/conditionsunassembled'),),product='Sourcearound30ntprotectedandotherRNAfragmentclasses; peak~23ntanalysiswidthnotexactfragmentlength')
    b.op('end_process','SourcefurtherprocessedRNA5P/3OHclass, exactenzymesunassembled',inputs=(b.stock('end_stock','Actualsourceendprocessingmolecularreagentstock; identity/orderinginternalchemistryunassembled, noinventedPNK'),),state=source_rna_end(),product='Sourceend-processedRNAmixture withintended5P/3OHclass; fullsubstrates/residualendsunassembled')
    kept=g.product('released_rna','SourcePK/TRIzolrecoveredRNApreparation; peptidecarryover/enddistributionunassigned',K.MOLECULAR,state=source_rna_end());residue=g.product('workup_other','ActualotherPK/wash/TRIzolworkupfractions, laterfatesunassigned',K.MOLECULAR)
    g.react('pk_rna_workup','SourcePKdigestion/wash/TRIzolRNAextraction; exactwashscheduleunassembled',S.MAGIC,((b.x,R.TEMPLATE),b.cat('pk','ActualsourceProteinaseK, dose/conditionsunassembled'),b.stock('trizol','ActualsourceTRIzolRNAworkuppreparation')),(kept,residue),reason='SourcepositivePKremovesantibody/peptidesbeforelibraryconstruction; exactwashschedule/fractionation/peptidecarryoverunassembled. No perfectdeproteinizationornewcapturestepinferred.');g.mark_terminal(Terminal(residue.id,Disposition.UNRESOLVED,'Otherworkupfractionsnotreturnedtoproduction.'));b.x=kept
    b.op('library_construct','SourceIlluminaTruSeqSmallRNAPrepKitlibraryboundary',inputs=(b.stock('smallrna_kit','ActualhistoricalTruSeqSmallRNAlibrarymolecularpreparedstock; wholeadapters/NTs/RTprimers/index/PCRchemistryunassembled'),),product='SourcesmallRNAsequencinglibrarypreparation; wholeproducts/index/readgeometryunassembled')
    b.aliquot('loading','Actualsequencingallocation, amountunassembled')
    b.finish('library','PA-m6A-seqsourcequalifiedRNAsequencinglibrary','AssociatedEnglishmain statesIlluminaTruSeqSmallRNA/highthroughputsequencing; exactinstrument/runmask/wholeadapter/primerlandingsunassembled.4SUchemicalpositionsremainmarked/unassigned, notconvertedtoCmerelybecauseT-to-Creadsignatureisobserved.')
    prior=b.x;b.x=g.input('igg_control_input','Actualseparate4SU-mRNAparallelIgGvalidationpreparation, originalallocationunassembled',K.MOLECULAR,primary=True)
    b.op('igg_control','Sourceparallel normalrabbitIgGIP/same365UVtreatmentvalidation',inputs=(b.stock('normal_igg','ActualnormalrabbitIgGcontrolselectionstock',K.SELECTION,R.SELECTION),b.cond('igg_conditions','SourceFigure1aparallelvalidation; detectionworkup/amountsnotinstantiated')))
    g.mark_terminal(Terminal(b.x.id,Disposition.QC,'SourceparallelIgGvalidationboundary; nozero-backgroundclaimorproductionlibrarymerge.'));b.x=prior
    g.add_source_item(SourceItem('source_limits','AssociatedEnglishAMcoreorderonly; exactpublishedEnglishWileymainHTTP403 andactualPMC-SIcooldownblocked.','NoGerman-Englishrecipeequivalence, guessedSIpath, kitconstructorinheritance, sourceabsenceorcompletechemicalvalidation.'))
    return g

def selftest():
    from checks import Check,Source,have
    from storage import sources_root
    c=Check();g=protocol_graph()
    c('FullmRNAIPprecedesUV',g.steps['uv'].inputs[0].material,'m6a_ip_bound')
    c('T1trimfollowsactualUVmixture',g.steps['rnase_t1'].inputs[0].material,'uv_product')
    c('EndprocessingbeforePKworkup',g.steps['pk_rna_workup'].inputs[0].material,'end_process_product')
    c('KitusesactualrecoveredRNA',g.steps['library_construct'].inputs[0].material,'released_rna')
    c('NoPNKinventedfromsourceends',all('pnk'not in k for k in g.materials),True)
    c('Controlnotsequencedproduction',g.terminals['igg_control_product'].disposition,Disposition.QC)
    c('UnboundIPfractionnotreturned',g.terminals['m6a_ip_unbound'].disposition,Disposition.UNRESOLVED)
    p=sources_root()/Path(__file__).resolve().parents[1].name/'anie.201410647_PMC4396828.xml.txt'
    if have(c,Source(p,'https://doi.org/10.1002/anie.201410647')):
        t=''.join(p.read_text().split());c('SourceUVafterIP','AftertheIPstep,thesampleisirradiatedby365nmUV'in t,True);c('Sourcefull-lengthratherthanfragmented','full-lengthratherthanfragmentedmRNA'in t,True);c('SourceT1notgenericpreIPfragmentation','around30ntusingRNaseT1'in t,True)
    c('Partialnotpublishedready',bool(g.publication_problems()),True);c.report()
if __name__=='__main__':selftest()
