"""Published 2OMe-seq Method1 high/low-dNTP pair; private partial chemistry."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Scene,Segment,revcomp
from illumina import TRUSEQ_SMALL_RNA_RA3,TRUSEQ_SMALL_RNA_RA5
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,SourceItem,SourceCoverage,SequencingPlan,Terminal,Disposition
from adapter_ligation import AdapterLigationGeometry
from pcr import PCRGeometry
TITLE='2OMe-seq Method1 without fragmentation'
SOURCE='PublishedPMC5388417 Methods70/Fig1;publishedSI File006 Method2 lines215–277 inspectedonlyasseparatealternative. Twooriginals; nofullsource/absenceaudit.'
RT_BODY=revcomp(TRUSEQ_SMALL_RNA_RA3)
ADAPTER=revcomp(TRUSEQ_SMALL_RNA_RA5)

def inventory(name):
    sc=Scene()
    if name=='rt':sc.strand('rt',[Segment('source_handle',RT_BODY),Segment('random_hexamer','N'*6,placeholder=True,note='RTprimingdegeneracy,notUMI')])
    else:sc.strand('adapter',[Segment('source_adapter',ADAPTER)],mod5='P',mod3='C3 spacer')
    return sc

def known_ends():
    sc=inventory('adapter');sc.strand('cdna',[Segment('source_cdna','[cDNA sequence unassembled]',placeholder=True)],mod3='OH')
    sc.ligation_end('adapter','source_adapter',kind='single-strand',polarity='5-prime',owner='adapter',text='Sourceadapter5P;3C3blocksoppositeend')
    sc.ligation_end('cdna','source_cdna',kind='single-strand',polarity='3-prime',owner='substrate',text='SourcecDNA3OH;wholejunctionunassembled');return sc

def unknown():
    sc=Scene();sc.strand('unknown',[Segment('unassembled','[unassembled]',placeholder=True)]);return sc

def protocol_graph():
    g=ProtocolGraph(TITLE,authority=SOURCE,representation='published_instance');stock=g.input('rna_stock','SourcetotalRNAstock;massunassignedbeyondtwo2ugallocations',K.MOLECULAR,primary=True)
    arms={a:g.product(a+'_rna','Source2ugRNAallocatedfor'+a+'dNTPMethod1library',K.MOLECULAR)for a in ('high','low')};rem=g.product('rna_remainder','UntransferredRNAstock;amountunassigned',K.MOLECULAR);g.react('allocate','Sourcehigh/lowlibraryallocations',S.SPLIT,((stock,R.TEMPLATE),),(*arms.values(),rem));g.mark_terminal(Terminal(rem.id,Disposition.UNRESOLVED,'Untransferredstock; source storage/fate unassigned'))
    for arm,x in arms.items():
        def key(k):return arm+'_'+k
        def clean(k,label):
            nonlocal x
            b=g.input(key(k+'_condition'),label,K.CONDITION);y=g.product(key(k+'_kept'),'SourcephysicallyretainedDNApreparation;yield/carryoverunassigned',K.MOLECULAR);w=g.product(key(k+'_waste'),'Actualsourcewash/gel/cleanupunretainedfraction',K.WASTE);g.react(key(k),label,S.CLEANUP,((x,R.TEMPLATE),(b,R.CONDITION)),(y,w));g.discard(w,'Physicalfractionnotperfectmolecularselectivity');x=y
        spike=g.input(key('spike'),'Source100nginvitroRNAspikepoolwhereused;selectedspike-containingcase,wholemixturecompositionunassembled',K.MOLECULAR);y=g.product(key('spiked'),'SourcebiologicalRNAplusstudysyntheticRNAcontribution',K.MOLECULAR);g.react(key('spike_mix'),'Source100ngspikeallocation',S.GENERIC,((x,R.TEMPLATE),(spike,R.SPIKE_IN)),(y,),changes_molecule=False);x=y
        primer=g.input(key('rt_primer'),'SourceRA3reversecomplement/randomhexamerRTprimer',K.OLIGO,state=inventory('rt'));enz=g.input(key('amv'),'SourceAMVRT.5Uul/20ulRT',K.CATALYST);nts=g.input(key('dntps'),'Source'+('1mM'if arm=='high'else'.004mM')+'finaldNTPasstated;per-base-versus-totalconventionunassigned',K.MOLECULAR);b=g.input(key('rt_buffer'),'Source42C30min;RTbuffer/primerconcentrationunassigned',K.CONDITION);y=g.product(key('rt_product'),'SourcecDNApopulation/RNA/freeprimerreaction;lowdNTPstopsnotperfectsite-specificconversion',K.MOLECULAR);g.react(key('rt'),'Sourcehigh/lowdNTPRT',S.REVERSE_TRANSCRIPTION,((x,R.TEMPLATE),(primer,R.STRUCTURAL),(nts,R.STRUCTURAL),(enz,R.CATALYST),(b,R.CONDITION)),(y,));x=y
        heat=g.input(key('rt_heat'),'Source95C10minAMVinactivation',K.CONDITION);y=g.product(key('rt_stopped'),'SourceheatedRTmixture',K.MOLECULAR);g.react(key('rt_stop'),'SourceRTenzymeheatstop',S.THERMAL_TREATMENT,((x,R.TEMPLATE),(heat,R.CONDITION)),(y,));x=y
        rnases=g.input(key('rnases'),'Source1ulRNaseH+1ulRNaseCocktail37C30min',K.CATALYST);y=g.product(key('rna_digested'),'SourcecDNAandRNAcleavage/freeprimercompanionsbeforephysicalcleanup',K.MOLECULAR);g.react(key('digest_rna'),'SourceRNAtemplateenzymaticdegradation',S.DIGESTION,((x,R.TEMPLATE),(rnases,R.CATALYST)),(y,));x=y
        clean('ampure','SourceAMPureXPrecoverscDNA/removesexcessprimer;notperfectselection')
        clean('gel','Source10percentTBEureaPAGE70–300ntslice')
        elute=g.input(key('diffusion'),'SourceTE10mMTris/1mMEDTA37C16hmoderateshake',K.CONDITION);eluate=g.product(key('gel_eluate'),'SourcesolublecDNArecoveredfromgelslice;recoveryunassigned',K.MOLECULAR);res=g.product(key('gel_residue'),'Remaininggelslice/unrecoveredDNA;laterfateunassigned',K.MOLECULAR);g.react(key('diffuse'),'SourcepassivegelDNArecovery',S.EXTRACTION,((x,R.TEMPLATE),(elute,R.CONDITION)),(eluate,res));g.mark_terminal(Terminal(res.id,Disposition.UNRESOLVED,'Sourcegelrecoverydoesnotestablishcompleterecovery'));x=eluate;clean('precipitate','Sourceethanolprecipitation')
        adapter=g.input(key('adapter'),'Source5P/3C3reverseRA5DNAadapter',K.ADAPTER,state=inventory('adapter'),chemical_state={'five_prime':'P','three_prime':'C3 spacer','polymer':'DNA'})
        circ=g.input(key('circligase'),'Source200UCircLigaseII65C6h;othercofactoramountsunassembled',K.CATALYST);y=g.product(key('adapted'),'Sourceadapter-ligatedandunligatedcDNAclasses;sealedproductunassembled',K.MOLECULAR);g.react(key('ligate'),'Sourceintermolecularadapter5PtocDNA3OHligation',S.LIGATION,((x,R.TEMPLATE),(adapter,R.STRUCTURAL),(circ,R.CATALYST)),(y,),geometry=AdapterLigationGeometry.known('single-strand',known_ends(),statement='PublishedMethods70whole26nt5P/3C3adapterandcDNA3OH;single-strandendinterfaceknown,sealedproduct/RTstoppositionsunassembled.'));x=y
        mix=g.input(key('pcr_mix'),'Aggregate actual source PCR supplies with standardTruSeq primerpair; packaging/stockcount/wholeorders/enzyme identity unassigned',K.MOLECULAR);f=g.product(key('pcr_f'),'Forward primer component of aggregate PCR supplies; packaging/ordering unassigned',K.OLIGO);r=g.product(key('pcr_r'),'Reverse primer component of aggregate PCR supplies; packaging/ordering unassigned',K.OLIGO);e=g.product(key('pcr_e'),'Catalyst component of aggregate PCR supplies; identity/packaging unassigned',K.CATALYST);b=g.product(key('pcr_b'),'Buffer component of aggregate PCR supplies; packaging unassigned',K.CONDITION);pn=g.product(key('pcr_n'),'Required nucleotide contribution from aggregate PCR supplies; formulation/packaging unassigned',K.MOLECULAR);g.react(key('pcr_components'),'Expose aggregate actual PCR supplies; packaging and stockcount unassigned, not an asserted bulk premix',S.GENERIC,((mix,R.TEMPLATE),),(f,r,e,b,pn),changes_molecule=False)
        y=g.product(key('pcr_product'),'Source15cycleindexedampliconmixture;fullprimer/strandgeometryunassembled',K.MOLECULAR);g.react(key('pcr'),'SourceMethod1exact15cycles,notMethod2eighteen',S.PCR,((x,R.TEMPLATE),(f,R.STRUCTURAL),(r,R.STRUCTURAL),(e,R.CATALYST),(b,R.CONDITION),(pn,R.STRUCTURAL)),(y,),geometry=PCRGeometry.unavailable(unknown(),reason='SourceTruSeqPCRpairnamedbutexactindexprimerorderings/landingandampliconsunassembled.'));x=y;clean('pcr_clean','SourceAMPureXPexcessprimerremoval')
        selected=g.product(key('pool_aliquot'),'Sourcecontributionnormalizedforequimolarpool;volumeunassigned',K.MOLECULAR);remaining=g.product(key('pool_remainder'),'UntransferredpurifiedPCRstockifany;volume/lateruseunassigned',K.MOLECULAR);g.react(key('pool_allocate'),'Sourceequimolarpoolstockallocation',S.SPLIT,((x,R.TEMPLATE),),(selected,remaining));g.mark_terminal(Terminal(remaining.id,Disposition.UNRESOLVED,'Normalizationvolumeanduntransferredstockfateunassembled'));arms[arm]=selected
    pool=g.product('equimolar_pool','Sourceequimolarhigh/lowlibrarypool;normalization/volumesunassigned',K.MOLECULAR);g.react('pool','Sourceequimolarlibrarypool',S.POOL,tuple((x,R.TEMPLATE)for x in arms.values()),(pool,));selected=g.product('load_selected','Selected source pooled-library loading contribution; amountunassigned',K.MOLECULAR);remainder=g.product('load_remainder','Untransferred prepared pool; amountunassigned',K.MOLECULAR);g.react('load_allocate','Select loading aliquot from actualequimolarpool',S.SPLIT,((pool,R.TEMPLATE),),(selected,remainder),changes_molecule=False);g.mark_terminal(Terminal(remainder.id,Disposition.UNRESOLVED,'Untransferred pool, notsilentlydiscarded'));lib=g.product('library','SourceMethod1pooledNGSlibrary;fullreadgeometryunassembled',K.LIBRARY);g.react('identity','PhysicalNGSlibraryidentity',S.GENERIC,((selected,R.TEMPLATE),),(lib,),changes_molecule=False);g.sequence(lib,SequencingPlan(limitation='SourceNextSeq500. Exactrun/indexmasksunassigned. RandomRT6NisnotUMI. SourceMethod2fragmentation/50–150ntcDNA/18PCRcyclesremainsaseparatealternative,notmixed.'));g.add_source_item(SourceItem('selected',SOURCE,'Knownsingle-strandligationendsnotwholeproduct/completeRTvalidation.'));return _publication_scope(g)

def selftest():
    from checks import Check,Source,have
    from storage import sources_root
    c=Check();g=protocol_graph()
    c('WholeRTpattern',RT_BODY+'N'*6,'CCTTGGCACCCGAGAATTCCANNNNNN')
    c('WholemodifiedDNAadapter','/5P/'+ADAPTER+'/3C3/','/5P/GATCGTCGGACTGTAGAACTCTGAAC/3C3/')
    c('Knownligationendroles',[(e.owner,e.polarity)for e in known_ends().ligation_ends],[('adapter','5-prime'),('substrate','3-prime')])
    c('NoMethod1fragmentation',any(x.kind==S.FRAGMENTATION for x in g.steps.values()),False)
    c('BothRNaseworkupsbeforeligation',g.steps['low_ligate'].inputs[0].material,'low_precipitate_kept')
    c('RTdNTPsmolecularsubstrates',[(i.material,i.role)for i in g.steps['low_rt'].inputs if i.material=='low_dntps'],[('low_dntps',R.STRUCTURAL)])
    c('Poolconsumesnormalizedaliquots',[i.material for i in g.steps['pool'].inputs],['high_pool_aliquot','low_pool_aliquot'])
    c('UntransferredPCRstockpreserved',g.terminals['low_pool_remainder'].disposition,Disposition.UNRESOLVED)
    c('Physicalpoolcontainsbothlibraries',len(g.steps['pool'].inputs),2)
    src=Source(sources_root()/'2ome-seq__10.1093+nar+gkw810'/'nar_gkw810_PMC5388417.xml.txt','https://pmc-oa-opendata.s3.amazonaws.com/PMC5388417.1/PMC5388417.1.xml;Methods70')
    if have(c,src):
        t=src.path.read_text();c('WholeRTorderingagainstoriginal',RT_BODY+'N'*6 in t,True);c('Wholeadapteragainstoriginal',ADAPTER in t,True)
    c('Untransferred RNA stock has no invented storage disposition',g.terminals['rna_remainder'].disposition,Disposition.UNRESOLVED);c('Loading consumes selected pool contribution',g.steps['identity'].inputs[0].material,'load_selected');c('PCR nucleotide class derives from same aggregate supplies',g.producer_of('high_pcr_n').id,'high_pcr_components');c('Selected source-qualified publication gate',g.publishable,True);c.report()



PAGE_FILE='2ome-seq.html'
# Source-defined high/low experiment has two distinct preparations. These literal
# suffixes expand that source pair, never discover/auto-account graph nodes.
REVIEWED_ARM_TARGETS=('spike','spike_mix','rt_primer','amv','dntps','rt_buffer','rt','rt_heat','rt_stop','rnases','digest_rna','ampure_condition','ampure','gel_condition','gel','diffusion','diffuse','precipitate_condition','precipitate','adapter','circligase','ligate','pcr_mix','pcr_components','pcr','pcr_clean_condition','pcr_clean','pool_allocate')
def _publication_scope(g):
    g.cover_source_item(SourceCoverage('selected','modelled','library'))
    g.declare_reviewed_scope(('rna_stock','allocate','pool','load_allocate','identity','library'),SOURCE)
    for arm in ('high','low'):
        g.declare_reviewed_scope(tuple(arm+'_'+target for target in REVIEWED_ARM_TARGETS),SOURCE)
    return g.require_publishable()

def render_page():
    from page import head,caveat,table
    from chemdraw import panel
    from protocol_views import assemble_page,WorkflowSection,PlanLimitationSection
    g=protocol_graph()
    return assemble_page(g,[head('2OMe-seq'),'<main class="wrap" data-protocol-complete="true"><h1>2OMe-seq</h1><p>Compare reverse-transcription stops at high and low nucleotide concentrations to map candidate RNA 2′-O-methylation sites.</p><p><a href="https://doi.org/10.1093/nar/gkw810">Defining published article</a> · <a href="01_chem.html">Research note</a></p>',
      caveat('Selected published Method1, without fragmentation: two separate 2-µg RNA preparations, shown with the optional source spike-in mixture. Native RNA/cDNA populations, exact RT stops, complete ligation/PCR products, index assignments and sequencing-primer landing remain unassembled. Low-dNTP stopping is not a claim of perfect site specificity or quantitative conversion.'),
      '<h2>Preparation</h2>',table(('Stage','Source operation'),(('Reverse transcribe','High 1 mM or low 0.004 mM dNTP as stated; AMV at42°C for30min. The concentration convention is unassigned.'),('Recover cDNA','Heat stop, RNaseH/cocktail digestion, AMPure,70–300nt PAGE and passive diffusion/ethanol recovery.'),('Ligate','CircLigaseII joins the source 5′-phosphorylated,3′-C3-blocked DNA adapter to cDNA3′OH. Full sealed products remain unknown.'),('Amplify and load','15PCRcycles, AMPure, normalized contributions pooled equimolarly; select a loading contribution for NextSeq500.'))),
      '<h2>Known molecular parts</h2>',panel(inventory('rt').rows(),caption='Whole source RT primer: constant handle and random priming hexamer, not a UMI'),panel(inventory('adapter').rows(),caption='Whole source DNA adapter;5′ phosphate and3′ C3 block'),
      '<h2>Read setup and limits</h2><p>Source acquisition uses NextSeq500. Complete run/index masks and primer landing remain unassembled. The source states that reads begin one nucleotide downstream of the RT stop; stop enrichment and site inference are computational comparisons, not physical removal of every unrelated molecule.</p><p>Fragmented Method2 uses different RNA/cDNA selections and18PCRcycles. It is a separate source alternative and is not merged into this unfragmented15-cycle preparation.</p>',
      '<details><summary>Complete preparation, stocks and unresolved chemistry</summary>',WorkflowSection('High/low physical preparations'),PlanLimitationSection('library'),'</details></main>'])

if __name__=='__main__':selftest()
