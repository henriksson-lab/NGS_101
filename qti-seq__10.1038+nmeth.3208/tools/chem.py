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

from dataclasses import dataclass
from illumina import TRUSEQ_SMALL_RNA_RA5
BC={'MCA02':'CAG','LGT03':'GTG','YAG04':'AGG','HTC05':'TCG'}
COMMON_LEFT=revcomp(TRUSEQ_SMALL_RNA_RA5)[1:-4]
RIGHT=P7[:-3]+'T'*20+'VN'
PCR_F=P5[:19]+'ACAG'+TRUSEQ_SMALL_RNA_RA5[:-3]
PCR_R=P7[:18]
READ_PRIMER='CGACAG'+TRUSEQ_SMALL_RNA_RA5
AUTHORITY='AssociatedAMNIHMS642756/PMC4344187Methods110–142(originalPDFp10wholeorderings); publishedSI592Fig1/2captions; publishedNaturemain subscriptionpreview, nofinalMethods equivalence. SelectedLTMculturedcellbranch.'
@dataclass(frozen=True)
class SourceRTPrimer:
    choice:str='MCA02'
    def __post_init__(self):
        if self.choice not in BC:raise ValueError('Select one actual source whole ordering')
    @property
    def left(self):return BC[self.choice]+COMMON_LEFT
    @property
    def right(self):return RIGHT
    @property
    def ordering(self):return 'p'+self.left+'Ø'+self.right
    @property
    def parts(self):return (('DNA',self.left),('source_linker','Ø','chemicalidentity/cleavagebondunassigned'),('DNA',self.right))
    @property
    def known_nucleotide_count(self):return len(self.left)+len(self.right)

def primer_left_scene():
    p=SourceRTPrimer();s=Scene();s.strand('left',[Segment('source leftDNAblock',p.left)],mod5='P',label='PartialleftDNAblockofoneprimer; Ølinker/rightblock retainedinchemical_state, not silentlyjoined');return s

def protocol_graph():
    b=B('QTI-seq selectedLTMculturedcell footprint preparation',AUTHORITY,'Actualtwo10cmdishescells, upstreamgrowth/starvationoutsidegraph',K.BIOLOGICAL);g=b.g
    b.op('lyse','Source400ulicepolysomebuffer+5uMLTM/MatrixD6x20s40sice',S.EXTRACTION,(b.stock('ltm','Actualsource5uMlactimidomycin molecularstock'),b.cond('lysis_buffer','Source10mMHEPES7.4/100mMKCl/5mMMgCl2/ice; noinitialCHXcomparatoraddition')))
    clear=g.product('lysate_supernatant','Actual12000xg10min4Csupernatant',K.MOLECULAR);pel=g.product('lysate_pellet','Actualremovedlysatepellet',K.MOLECULAR)
    g.react('clarify','Source12000xg4C10minretainssupernatant',S.SPLIT,((b.x,R.TEMPLATE),),(clear,pel));g.mark_terminal(Terminal(pel.id,Disposition.UNRESOLVED,'Nonretainedclarificationpellet'));b.x=clear
    b.op('pmy','Source25uMpuromycin15min remodelsphysicalribosomepopulation',inputs=(b.stock('puromycin','Actual25uMpuromycin, molecularnascentchainacceptor'),b.stock('energy_atp','Actual.8mMATP'),b.stock('creatine_phosphate','Actual10mMcreatinephosphate'),b.cat('creatine_kinase','Actual40ugmlcreatinephosphokinase'),b.cond('pmy_buffer','Source16mMHEPES7.4/.1mMspermidine')),product='SourceLTM/PMYtreatedribosome/RNA/peptidepopulations; no perfectinitiatingretention or completeelongatorremoval claim')
    pooled=g.product('gradient_selected','Actualsourcepooledgradientpolysomefractions; exactfractionnumbers/amountsunassembled',K.MOLECULAR);other=g.product('gradient_other','Actualothergradientfractions, exactfateunassigned',K.MOLECULAR)
    g.react('gradient','Source15–45%sucroseSW41/38000rpm100min4C; .375mlmin/1minfractionation',S.SPLIT,((b.x,R.TEMPLATE),b.cond('gradient_buffer','Source10mMHEPES7.4/100mMKCl/5mMMgCl2/CHX100ugml/5mMDTT/20UmlSUPERaseIn')),(pooled,other));g.mark_terminal(Terminal(other.id,Disposition.UNRESOLVED,'Otherphysicalgradientfractions'));b.x=pooled
    b.op('rnase','SourceRNaseI750Uper100A260/4C1h',S.DIGESTION,(b.cat('rnase_i','ActualEcoliRNaseI750Uper100A260'),b.cond('rnase_program','Source4C1h')))
    b.op('rnase_stop','SourceSUPERaseIn50Uper100URNaseI',S.QUENCHING,(b.stock('superase_stop','ActualsourceSUPERaseIn inhibitorprotein'),b.cond('superase_ratio','Source50Uinhibitorper100URNaseI quenchingratio')))
    b.clean('rna_extract','SourceTRIzolRNAextraction; wholeinput/RNAyieldunassigned')
    b.op('dephos','Source15ulPNKreaction20U/10USUPERase37C1h65C20min',S.DEPHOSPHORYLATION,(b.cat('pnk','Actual20UT4PNK'),b.cat('dephos_inhib','Actual10USUPERaseIn'),b.cond('pnk_buffer','Source1xPNKbuffer, nofreshATPintroduced')))
    b.op('rna_gel','Source15percentdenaturingTBEureaPAGE thenSYBRGoldvisualization',S.GENERIC,(b.stock('sybr','ActualSYBRGold molecularstain'),b.cond('gel_buffer','Source2xNovexsamplebuffer/denaturing15percentPAGE')))
    b.clean('rna_band','Source~28ntRNAgelband/crush/overnight300mMNaOAc5.5/1mMEDTA/.1UulSUPERase/SpinX/ethanol')
    b.op('pap_denature','SourceRNA65C30s in10mMTris8',S.THERMAL_TREATMENT,(b.cond('pap_denature_conditions','Actual65C30s'),))
    b.op('pap','Source8ulPAP3U/1mMATP/.75UulSUPERase37C45min',S.TAILING,(b.cat('pap_enzyme','Actual3UEcoliPAP'),b.stock('pap_atp','Actual1mMATP covalentsubstrate'),b.cat('pap_inhib','Source.75UulSUPERaseIn'),b.cond('pap_buffer','Actual1xPAPbuffer')))
    primer=g.input('rt_primer','SourcewholeMCA02representativeRTordering; biologicalsamplechoiceunassigned',K.OLIGO,state=primer_left_scene(),chemical_state={'source_ordering':SourceRTPrimer().ordering,'ordered_parts':SourceRTPrimer().parts,'five_prime':'P','sample_index_role':feature('qti.sample','sample_index','unknown'),'barcode_window':'unassembled; source variableleadingprefixes notmoleculeUMI','whole_primer_context':'PartialleftDNAblockdrawn; Ø/rightDNAblocksremainexplicitorderedparts'})
    b.op('rt_prepare','Source.5mMdNTP/printed2.5mMRTprimer65C5min/ice5min',S.THERMAL_TREATMENT,((primer,R.STRUCTURAL),b.stock('rt_nt','Actualsource.5mMdNTP, perbase/totalconventionunassigned'),b.cond('rt_anneal','Literal2.5mMprimerretained; notsilentlyuMcorrected')))
    tmp=g.product('rt_template','Same-tubeRNAtail/templatecomponent',K.MOLECULAR);pr=g.product('retained_rt_primer','Same-tube retainedexternalMCA02primercomponent',K.OLIGO,state=primer_left_scene(),chemical_state={'source_ordering':SourceRTPrimer().ordering,'ordered_parts':SourceRTPrimer().parts})
    g.react('rt_components','Copresenttemplate andretainedexternalprimer, notphysicalfractions/freshaddition',S.GENERIC,((b.x,R.TEMPLATE),),(tmp,pr),changes_molecule=False);b.x=tmp
    b.op('rt','SourceSSIII200U/RNaseOUT40U, vendorprogrammeunassembled',S.REVERSE_TRANSCRIPTION,((pr,R.STRUCTURAL),b.cat('ssiii','Actual200USSIII'),b.cat('rt_inhib','Actual40URNaseOUT'),b.cond('rt_buffer','Source20mMTris8.4/50mMKCl/5mMMgCl2/10mMDTT')))
    b.clean('cdna_band','Source10percentTBEureaPAGE~100ntextendedcDNA/300mMNaCl1mMEDTArecovery')
    b.op('circle','Source20ulCircLigaseII100U60C1h80C10min',S.MAGIC,(b.cat('circligase','Actual100UCircLigaseII'),b.cond('circle_buffer','Actual1xCircLigasebuffer/2.5mMMnCl2/1Mbetaine')),product='Sourcecircularized/unreactedssDNAclasses; exactwholeproduct/endclosureunassembled')
    b.op('linearize','Source7.5UAPE137C1hrelinearization',S.DIGESTION,(b.cat('ape1','Actual7.5UAPE1'),b.cond('ape_buffer','Actual20mMTrisacetate/50mMKacetate/10mMMgacetate/1mMDTT')),product='SourceAPE1processedlinear/circular/unreactedclasses; Øidentity/cleavagebond/terminalchemistryunassembled')
    b.clean('linear_band','Source10percentdenaturingPAGE~100ntlinearproduct/recovery')
    b.pcr('pcr','SourcePhusion.5U/HF/.2mMdNTP/.5uMeachprimer12cycles98C10/60C20/72C10s',oligo(PCR_F),oligo(PCR_R))
    b.clean('pcr_band','Source8percentnativePAGE~120bpband/recovery')
    b.aliquot('bioanalyzer','SourceBioanalyzerDNA1000analyticalquantification amountunassigned',qc=True)
    own=b.aliquot('pool_input','Actualequal-amountbarcodedpool contribution, quantityunassigned')
    others=g.input('other_barcoded_libraries','Actualotherqualifiedbarcodedlibraryaliquots; wholeindependentpreparations/stockfatesunassembled',K.MOLECULAR)
    pool=g.product('barcoded_pool','Actualsourceequal-amountbarcodedlibrarypool',K.MOLECULAR);g.react('pool','Sourceequalamountqualifiedbarcodedsamplespooled',S.POOL,((own,R.TEMPLATE),(others,R.TEMPLATE)),(pool,));b.x=pool
    b.aliquot('loading','Source3–5pMHiSeqloadingallocation, volumeunassigned')
    b.finish('library','QTIselectedculturedcellbarcodedlibrarypool','AssociatedAMHiSeq/custom32ntsequencingprimer '+READ_PRIMER+'; wholeprimerlanding/physicalrunmaskunassembled. Processing25–35ntreads/Psite12ntoffset is computational, notcompletephysicalreadgeometry or proofexactinitiationselectivity.')
    return g

def selftest():
    from checks import Check,Source,have
    from storage import sources_root
    c=Check();g=protocol_graph()
    c('PMYusesphysicalclarifiedlysate',g.steps['pmy'].inputs[0].material,'lysate_supernatant')
    c('Footprintingusesactualgradientselectedfractions',g.steps['rnase'].inputs[0].material,'gradient_selected')
    c('PAPnotfreshATPinthePNKstep',any(x.material=='pap_atp'for x in g.steps['dephos'].inputs),False)
    c('ExternalRTprimerretainednotreadmitted',g.steps['rt'].inputs[1].material,'retained_rt_primer')
    c('APE1followsactualcirclemixture',g.steps['linearize'].inputs[0].material,'circle_product')
    c('SourceØnotordinarybase',SourceRTPrimer().parts[1][0],'source_linker')
    c('DistinctprimerDNAblocksnotjoined',len(SourceRTPrimer().parts),3)
    c('Actualbarcodedpoolhastwocontributions',len(g.steps['pool'].inputs),2)
    p=sources_root()/Path(__file__).resolve().parents[1].name/'nmeth.3208_PMC4344187.xml.txt'
    if have(c,Source(p,'https://doi.org/10.1038/nmeth.3208')):
        t=''.join(p.read_text().split())
        for name in BC:c('Wholeoriginal5PØprimer '+name,SourceRTPrimer(name).ordering in t,True)
        for name,s in [('qNTI200',PCR_R),('qNTI201',PCR_F),('customread',READ_PRIMER)]:c('Wholeoriginalprimer '+name,s in t,True)
    c('Partialnotpublishedready',bool(g.publication_problems()),True);c.report()
if __name__=='__main__':selftest()
