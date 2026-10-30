from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from dataclasses import dataclass,replace
from chemdraw import Scene,Segment,Construct,feature,workflow_from_protocol_graph,workflow_panel
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,SourceItem,SequencingPlan,Terminal,Disposition
from page import head,info
PAGE_FILE='chem.html'

class Operations:
    def __init__(self,g):self.g=g
    def op(self,key,x,label,reagent=None,state=None,chemical_state=None):
        g=self.g;inputs=[(x,R.TEMPLATE)]
        if reagent:
            y=g.input(key+'_condition',reagent,K.CONDITION);inputs.append((y,R.CONDITION))
        y=g.product(key+'_out',label+'; preparation includes unremoved companion species',K.MOLECULAR,state=state,chemical_state=chemical_state)
        g.react(key,label,S.GENERIC,inputs,(y,));return y
    def join(self,key,inputs,label,state=None,chemical_state=None):
        g=self.g;y=g.product(key+'_out',label+'; reaction components remain until an explicit separation',K.MOLECULAR,state=state,chemical_state=chemical_state)
        g.react(key,label,S.GENERIC,tuple((x,R.TEMPLATE if i==0 else R.STRUCTURAL) for i,x in enumerate(inputs)),(y,));return y
    def clean(self,key,x,label):
        g=self.g;z=g.input(key+'_selection',label,K.SELECTION)
        y=g.product(key+'_kept','Source physically retained fraction: '+label,K.MOLECULAR)
        w=g.product(key+'_removed','Source physically removed fractions; exact composition/purity unassigned',K.WASTE)
        g.react(key,label,S.CLEANUP,((x,R.TEMPLATE),(z,R.SELECTION)),(y,w));g.discard(w,'Only the source physical separation; no perfect selectivity asserted.');return y
    def split(self,key,x,label,other_label,other=Disposition.UNRESOLVED):
        g=self.g;y=g.product(key+'_production',label,K.MOLECULAR)
        z=g.product(key+'_other',other_label,K.MOLECULAR)
        g.react(key,'Source physical allocation: '+label,S.SPLIT,((x,R.TEMPLATE),),(y,z));g.terminals[z.id]=Terminal(z.id,other,other_label);return y,z
    def qc(self,key,x,label):
        g=self.g;y=g.product(key+'_production','Production preparation retained independently from analytical aliquot',K.MOLECULAR)
        z=g.product(key+'_assay',label,K.QC)
        g.react(key,label,S.SPLIT,((x,R.TEMPLATE),),(y,z));g.terminals[z.id]=Terminal(z.id,Disposition.QC,'Analytical aliquot does not become the sequencing substrate.');return y
    def opaque(self,key,x,label):
        g=self.g;y=g.product(key+'_out',label,K.MOLECULAR)
        g.react(key,label,S.MAGIC,((x,R.TEMPLATE),),(y,),reason='Selected originals establish this operational continuation; full molecular inventory/geometry remains unassembled.');return y

def finish(g,x,readout):
    g.sequence(x,SequencingPlan(acquisition='unspecified',limitation=readout))
    g.declare_reviewed_scope(tuple(g.materials)+tuple(g.steps),SOURCE)
    g.add_source_item(SourceItem('unassembled','Unassembled source inputs, molecular geometry and directly linked targets',LIMIT));return g

def render_page():
    return head(TITLE)+'<main class="wrap"><h1>'+TITLE+'</h1>'+info(LIMIT)+'<p><a href="01_chem.html">Research notes</a></p>'+workflow_panel(workflow_from_protocol_graph(protocol_graph(),ENDPOINT))+'</main>'


from checks import Source,have
from storage import sources_root
import illumina as il
TITLE='NOIR — platform-specific ligated molecular tags and linear copying'
SOURCE='PublishedPMC4535617 Methods2.3/2.4; published dsv010supp.pdf Methodspp1–2,TablesS2/S3/footnotes'
LIMIT='Private incomplete source graph. NOIRconsensussoftwareisnotDNAchemistry. SourceIonN12 andIllumina15-baseBDHVBDHVBDHVBDH moleculebarcodes aredifferentpreparations,notoneassumed12bpwindow. Moleculartags versus sampleindexes areseparatestructuredinventories. Doublestrand-labelreplacementexperimentbypassesthe10-cyclelinearroute. Fullcohesiveends/junctions/modifiedprimerpools andhistoricalplatformgeometryremainunassembled;no perfecteliminationofearlyPCRerrors.'
ENDPOINT='ion_load_out'
PATTERN='BDHVBDHVBDHVBDH'
@dataclass(frozen=True)
class BarcodeClass:
    platform:str
    pattern:str
    def __post_init__(self):
        if (self.platform,self.pattern) not in {('Ion','N'*12),('Illumina',PATTERN)}:raise ValueError('Keepexactsourceplatform-specificbarcodepattern')
    @property
    def semantic(self):return feature('noir.'+self.platform.lower()+'.molecule','umi','random' if self.platform=='Ion' else 'degenerate',note='Source-declaredmolecularbarcode;notallPCRdaughtererrorscorrectedbychemistry')
    def copied(self):return (self.platform,self.pattern,'newDNAcopy;tagidentityretained,consensuslatercomputed')

def ill_adapter():
    return Construct([Segment('sourceRead1adapter',il.TRUSEQ_READ1),Segment('sourcefixedAGCTA','AGCTA'),Segment('Illumina15baseconstrainedmolecularbarcode',PATTERN,placeholder=True,feature=BarcodeClass('Illumina',PATTERN).semantic),Segment('sourcecommonend','GTACATATTGTCGTTAGAACGCG')])

def sample_inventory(platform):
    sc=Scene();sc.strand('index',[Segment('source5ntindividualindex' if platform=='Ion' else 'source6ntPCRsampleindex','N'*(5 if platform=='Ion' else 6),placeholder=True,feature=feature('noir.'+platform.lower()+'.sample','sample_index','unknown',note='Memberandwholeprimerplacementunassembled'))],label='Identifierinventoryonly,notaseparatefreeindexreagent');return sc

def prepare(g,o,prefix,dna,platform):
    x=o.op(prefix+'_digest',dna,'Sourceset-specificrestrictiondigestion12µL;30°C10min/37°C20min/45°C10min/65°C20min','Set1AlwNI+Alw26I,Set2EarI+NcoIorSetKCEarI+NmuCI;actualtargetallocationunassembled')
    a=g.input(prefix+'_adapter','Sourcemixedtarget-adaptorstock;exactstickyendallocationunassembled',K.OLIGO,state=ill_adapter() if platform=='Illumina' else None,chemical_state=BarcodeClass(platform,PATTERN if platform=='Illumina' else 'N'*12))
    x=o.join(prefix+'_ligation',(x,a),'E.coliDNAligase/NADsource24µLreaction;4×(16°C1h/13°C1h/10°C1h),molecularbarcodeligatedbeforecopying')
    x=o.clean(prefix+'_clean1',x,'Source1.2xAMPureligationcleanup1')
    x=o.clean(prefix+'_clean2',x,'Source1.2xAMPureligationcleanup2')
    return x

def protocol_graph():
    g=ProtocolGraph(TITLE,authority=SOURCE);o=Operations(g)
    for platform in ('Ion','Illumina'):
        k=platform.lower();dna=g.input(k+'_dna','Independentplatform-specificsourcegDNA/cfDNApreparation;notinventedsplitofonephysicalsample',K.MOLECULAR,primary=True)
        x=prepare(g,o,k,dna,platform)
        p=g.input(k+'_region_primer','Sourceregion-specificlinearprimerpool;exacttargetallocationandPSlinkagesunassembled',K.OLIGO)
        x=o.join(k+'_linear_mix',(x,p),'DissolvesourcepurifiedDNAin20µLQ5linear-amplificationmixture;originalDNA/unusedadaptersremain')
        x=o.clean(k+'_bead_remove',x,'SourceexplicitremovalofAMPurebeadsbeforethermalcycling;retainsolution,notnewselectiveDNAcleanup')
        x=o.op(k+'_linear',x,'Q5HotStart10linearcycles98°C10s/68°C10s/72°C30s,newtaggedDNAcopiesretainligatedmolecularidentity','Q5HotStart')
        x=o.clean(k+'_linear_clean',x,'Source1.2xAMPurepurificationoflinearproducts')
        f=g.input(k+'_pcr_f','SourcePGM/ProtonorindexedIlluminaPCRprimer,wholememberunassembled',K.OLIGO)
        r=g.input(k+'_pcr_r','SourceplatformPCRprimer;sampleindexmeaningsinstructuredinventory,wholebodyunassembled',K.OLIGO,chemical_state=sample_inventory(platform))
        x=o.join(k+'_pcr_mix',(x,f,r),'DissolvepurifiedDNAin20µLPlatinumTaqHighFidelityPCRmixture;Q5comparisonisaseparatesourcevariant')
        x=o.clean(k+'_pcr_bead_remove',x,'SourceexplicitremovalofAMPurebeadsbeforePCR;retainsolution')
        x=o.op(k+'_pcr',x,'Source30cyclesPCR95°C15s/60°C1minTP53or63°CKRAS/CTNNB1;targetallocationunreviewed','PlatinumTaqHighFidelity')
        x=o.opaque(k+'_workup_choice',x,'SourcepostPCRrouteallocationunassembled:Ionbeadorgel;Illuminaprintedtwice1xAMPureandalsoIllumina-gelrecommendation;sourceconditionselectionnotforced')
        x=o.clean(k+'_purify',x,'ActualsourcephysicalpostPCRpurification,platform/target-specificbeadorMinElutegelrouteunassigned')
        if platform=='Ion':x=o.opaque('ion_load',x,'SourceIonPIOT2-200v2orv3emulsionPCR/beadenrichmentOneTouch→IonPISequencing200v2orv3Proton;exacthistoricaltemplategeometryunassembled')
        else:x=o.opaque('illumina_load',x,'SourceIlluminaMiSeqsingle-endlibraryloading;exactindex/readprimerlandingunassembled')
        g.sequence(x,SequencingPlan(acquisition='unspecified',limitation='Published'+platform+'acquisitionroute;fullgeometry/cyclesorvariantallocationunassembled. Computationalbarcodefilteringisnotphysicaldiscard.'))
    dna=g.input('ds_dna','Independentsourceheteroduplexdouble-strand-labelvalidationexperiment;notallclinicalinput',K.MOLECULAR,primary=True)
    x=prepare(g,o,'ds',dna,'Ion')
    a=g.input('ds_t_pcr_a','SourceT_PCR_A',K.OLIGO,state=Construct([Segment('sourceprimer','CCATCTCATCCCTGCGTGTC')]))
    b=g.input('ds_region','SourceDSvalidationregion-specificprimerpool;wholebodiesunassembled',K.OLIGO)
    x=o.join('ds_replace',(x,a,b),'PlatinumTaqHighFidelitykit/PyrococcusGB-D72°C20or30minreplacement-synthesisexperiment;newstranddiffersfromnativemismatchedpartner')
    x=o.op('ds_pcr',x,'Source30cyclesPCRwithout10-cyclelinear-amplificationstep;nativeheteroduplexbasesnotcomputationallycorrected','PlatinumTaqHighFidelity')
    x=o.opaque('ds_finish',x,'Sourcevalidationlibraryphysicalworkup/loadingdetailsnotfullyassignedinselectedDSsection')
    g.sequence(x,SequencingPlan(acquisition='unspecified',limitation='Publishedheteroduplexvalidationsequenced/analyzedusingmoleculartags;completeworkup/rungeometryunassembled.'))
    g.declare_reviewed_scope(tuple(g.materials)+tuple(g.steps),SOURCE);g.add_source_item(SourceItem('unassembled','Remainingwholeoligos,endgeometry,platform/cleanupallocationandlinkedsourcekitrevisions',LIMIT));return g

def selftest():
    from checks import Check
    c=Check();g=protocol_graph();ion=BarcodeClass('Ion','N'*12);ill=BarcodeClass('Illumina',PATTERN)
    c('Platformmolecularbarcodewidthsremain12and15',(len(ion.pattern),len(ill.pattern)),(12,15))
    c.raises('CannotreinterpretIlluminaasfreeN12',lambda:BarcodeClass('Illumina','N'*12))
    c('Copiesretaintheoriginaltagidentity',ion.copied()[:2],('Ion','N'*12))
    c('Sourcephysicallyremovesbeadsbeforelinearthermalcycling',g.steps['ion_linear'].inputs[0].material,'ion_bead_remove_kept')
    c('PCRonlyafterpostlinearphysicalpurification',g.steps['ion_pcr_mix'].inputs[0].material,'ion_linear_clean_kept')
    c('DSreplacementdoesnotrunthroughlinear10cycles','ds_linear' in g.steps,False)
    c('IonandIlluminaoriginaltemplatesindependent',g.steps['ion_digest'].inputs[0].material!=g.steps['illumina_digest'].inputs[0].material,True)
    src=Source(sources_root()/'noir-sequencing__10.1093+dnares+dsv010'/'dnares_dsv010_supp_supp_dsv010_dsv010supp.pdf.txt','PublishedSI TableS3ILMadS00andmolecularbarcodefootnote')
    if have(c,src,label='Sourceplatformbarcodetranscription'):
        t=src.path.read_text();c('WholeILMadS00sourcebodyretained',''.join(s.top for s in ill_adapter().segments) in t,True)
        c('Sourceexplicitlydeclares15basedegeneratemoleculartag','"'+PATTERN+'" sequence is molecular barcode.' in t,True)
    c('Incompletegeometrynopromotion',g.publishable,False);c.report()
