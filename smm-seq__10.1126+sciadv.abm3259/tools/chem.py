"""Published selected preparation; private partial molecular model, not a complete audit."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,SourceItem,SequencingPlan,Terminal,Disposition
from chemdraw import Scene,Segment
from adapter_ligation import AdapterLigationGeometry
from pcr import PCRGeometry
from circular import CircularizationGeometry

TITLE='SMM-seq selected pulse-RCA route'
AUTHORITY='PublishedPMC8993124 main76/82/134–140;originalPDFp4verified;publishedSI tablesS1–S3/FigS1–S4captions.2originals;historicalkits/predecessorsdeferred.'
INITIAL='Source QuickgDNA-extracted human genomicDNA; sourcecell/treatmentidentityunassigned, not a synthetic calibration molecule'
READOUT='SourceNovaSeqPE150. SourceUID/strand-family consensus and>=7-readfiltersarecomputational, not physical DNA removal or guaranteed independent error-free copies. SeparateundilutedUSER-conventionalcontrolis not thisselectedpulseRCAroute.'
DEFERRED='HistoricalNEB/Bioroncomponentoriginals;fullHpyCH4IIIadaptercut/ligation/dumbbell/RCA/PCRgeometry; sourcepulsegoto3vsnarrativedenaturationdifference preserved, notsilentlyfixed.'

def inventory_state(label):return None
def inventory_chemistry(label):return None

from illumina import NEBNEXT_ARM_READ2,TRUSEQ_READ1
from chemdraw import feature
RCA_PRIMER='GTAGGGAAAGAGTGTAGACTGGAGTTC'
def hairpin_sequence():return 'TCTTCTACAGT'+'N'*6+'A'+NEBNEXT_ARM_READ2+'U'+TRUSEQ_READ1
def inventory_chemistry(label):
    if label=='P5_HP6N':return {'polymer':'DNA','sequence':hairpin_sequence(),'deoxy_u_positions':(hairpin_sequence().index('U'),),'uid_start':11,'uid_length':6,'ends':'source end modifications not assigned'}
    return None
def inventory_state(label):
    sc=Scene()
    if label=='P5_HP6N':sc.strand('P5_HP6N',[Segment('source prefix','TCTTCTACAGT'),Segment('source6UID','N'*6,placeholder=True,feature=feature('smm_uid','umi','random')),Segment('source adapter DNA','A'+NEBNEXT_ARM_READ2+'U'+TRUSEQ_READ1,note='Internal deoxyU retained for sourceRCA; ordering only, notfold/cleavagegeometry')]);return sc
    if label.startswith('P5-RCA'):sc.strand('P5_RCA',[Segment('source RCA primer',RCA_PRIMER)]);return sc
    return None

SOURCE_ORDERINGS={'hairpin_oligo':{'source':'sciadv.abm3259_PMC8993124.xml.txt:136; sciadv.abm3259_PMC8993124.pdf p4','raw':'TCTTC TACAGT NNNNNN AGATCG GAAGAG CACACG TCTGAA CTCCAG TC /ideoxyU/ ACACTC TTTCCC TACACG ACGCTC TTCCGA TCT','modifications':'Internal DNAdeoxyU only; terminal phosphorylation/PS not invented.'}}
SOURCE_OLIGO_EXCEPTIONS={'hairpin_oligo':'Published main136/originalPDFp4 verified whole hairpin; text matcher splits spaced TC /ideoxyU/ into separate hits. Local source-backed whole-ordering check retains DNAU.'}
def unknown_scene():
    sc=Scene();sc.strand('unassembled',[Segment('geometry','[source geometry unassembled]',placeholder=True)]);return sc

def protocol_graph():
    g=ProtocolGraph(TITLE,authority=AUTHORITY,representation='published_instance')
    x=g.input('input',INITIAL,K.MOLECULAR,primary=True)
    def op(key,label,kind,reagent,product,oligos=()):
        nonlocal x
        enzyme=kind in {S.DIGESTION,S.END_REPAIR,S.EXTENSION,S.TAILING,S.CIRCULARIZATION,S.LIGATION,S.PCR,S.DEPROTEINIZATION,S.TAGMENTATION,S.TRANSCRIPTION,S.REVERSE_TRANSCRIPTION,S.CHEMICAL_CONVERSION}
        c=g.input(key+'_stock',reagent,K.CATALYST if enzyme else K.CONDITION)
        inputs=[(x,R.TEMPLATE),(c,R.CATALYST if enzyme else R.CONDITION)]
        for i,a in enumerate(oligos):
            o=g.input(key+'_oligo_'+str(i),a,K.OLIGO,state=inventory_state(a),chemical_state=inventory_chemistry(a));inputs.append((o,R.STRUCTURAL))
        kw={'changes_molecule':False} if kind==S.GENERIC else {}
        if kind==S.MAGIC:kw['reason']='Actual published component; full source-specific molecular transition unassembled.'
        if kind==S.CIRCULARIZATION:kw['geometry']=CircularizationGeometry.unavailable('Source intramolecular closure known; full substrate and closure geometry unassembled, not a claim of source absence.')
        if kind==S.LIGATION:kw['geometry']=AdapterLigationGeometry('unavailable','unavailable',(unknown_scene(),),'Source ligation operation known; full ends and sealed junction unassembled.')
        if kind==S.PCR:kw['geometry']=PCRGeometry.unavailable(unknown_scene(),reason='Source primer inventories/whole amplicon geometry not assembled in this selected model.')
        y=g.product(key+'_product',product,K.MOLECULAR);g.react(key,label,kind,inputs,(y,),**kw);x=y
    def clean(key,label):
        nonlocal x
        c=g.input(key+'_condition',label,K.CONDITION);y=g.product(key+'_kept','Physical source recovered preparation; carryover and recovery unquantified',K.MOLECULAR);w=g.product(key+'_waste','Actual source physically separated fractions',K.WASTE)
        g.react(key,label,S.CLEANUP,((x,R.TEMPLATE),(c,R.CONDITION)),(y,w));g.discard(w,'Physical workup, not a perfect molecular specificity assertion.');x=y
    def aliquot(key,label):
        nonlocal x
        y=g.product(key+'_used',label,K.MOLECULAR);r=g.product(key+'_rest','Physicalsource remainder; mass/fateunassigned',K.MOLECULAR);g.react(key,label,S.SPLIT,((x,R.TEMPLATE),),(y,r));g.mark_terminal(Terminal(r.id,Disposition.UNRESOLVED,'Not implicitly discarded or counted as sequenced.'));x=y
    def qc(key,label):
        nonlocal x
        y=g.product(key+'_retained','Production preparation retained after independent analytical allocation',K.MOLECULAR);q=g.product(key+'_assay',label,K.QC);g.react(key,label,S.SPLIT,((x,R.TEMPLATE),),(y,q));g.mark_terminal(Terminal(q.id,Disposition.QC,'Analytical aliquot; not returned to production.'));x=y
    qc('input_quant','SourceQubit extractedDNAQC')
    op('digest','AluI/MluCIdouble digest37Covernight',S.DIGESTION,'SourceAluI/MluCI','Source restriction fragments and unprocessed DNA; genomic cleavage coordinates unassembled')
    clean('digest_clean','Source1.5xAMPure')
    op('endprep','SourceFig1A/main82 endrepair/dAtail; historicalkitinternalsunassembled',S.END_REPAIR,'SourceNEBNextUltraII;internalconditionsunreviewed','Sourceend-conditionedDNA plusunprocessedmolecules/reagents')
    hp=g.input('hairpin_oligo','P5_HP6N source ordering inventory',K.OLIGO,state=inventory_state('P5_HP6N'),chemical_state=inventory_chemistry('P5_HP6N'));cond=g.input('fold_conditions','Source86ul10uMhairpin0.1xTE;95C5min37C5min',K.CONDITION);fold=g.product('folded_hp','Source self-annealedhairpin preparation; exactfoldgeometryunassembled',K.MOLECULAR);g.react('hairpin_fold','Sourceintramolecularannealing',S.FOLDING,((hp,R.TEMPLATE),(cond,R.CONDITION)),(fold,))
    k=g.input('hp_klenow','SourceKlenowexo-minus10U/10ulCutSmart/2ul10mMdNTP37C30min',K.CATALYST);filled=g.product('filled_hp','Source5primeextension-filledhairpinmixture; random6copiedstem/fullgeometryunassembled',K.MOLECULAR);g.react('hairpin_fill','Sourcefillinhairpin',S.EXTENSION,((fold,R.TEMPLATE),(k,R.CATALYST)),(filled,))
    c=g.input('hp_workup','SourceQIAquickNucleotideRemoval workup; finalelute100ulEB',K.CONDITION);pre=g.product('hp_preclean','Sourcephysicallyrecoveredfilledhairpin',K.MOLECULAR);firstw=g.product('hp_first_waste','FirstNucleotideRemovalphysicalfractions',K.WASTE);g.react('hairpin_first_cleanup','SourcecleanupbeforeHpyCH4III',S.CLEANUP,((filled,R.TEMPLATE),(c,R.CONDITION)),(pre,firstw));g.discard(firstw,'Actualfirstadapterworkup.');e=g.input('hp_hpy','SourceHpyCH4III10U',K.CATALYST);cut=g.product('cut_hp_mix','SourcecutadapterandDNAcleavageproducts; no inventedprecisecutjunction',K.MOLECULAR);g.react('hairpin_cut','SourceHpyCH4IIIafteractualfirstcleanup',S.DIGESTION,((pre,R.TEMPLATE),(e,R.CATALYST)),(cut,));ready=g.product('ready_hp','Sourcephysicallyrecoveredadapterpreparation',K.ADAPTER);w=g.product('hp_waste','SourceNucleotideRemovalphysicalfractions',K.WASTE);g.react('hairpin_workup','SourceNucleotideRemovalfinalcleanup',S.CLEANUP,((cut,R.TEMPLATE),(c,R.CONDITION)),(ready,w));g.discard(w,'Actualphysicalworkup; no perfectadapterpurityasserted.')
    enzyme=g.input('ligation_enzyme','HistoricalNEBNextUltraIIligase',K.CATALYST);y=g.product('ligated','Sourcedumbbell-likeconstructswithintactU plus partial/unligatedproducts',K.MOLECULAR);g.react('ligation','SourceP5_HP6Nreplaceskitadapter; fullendgeometryunassembled',S.LIGATION,((x,R.TEMPLATE),(ready,R.STRUCTURAL),(enzyme,R.CATALYST)),(y,),geometry=AdapterLigationGeometry('unavailable','unavailable',(unknown_scene(),),'Sourcehairpinligationknown;fulladaptercleavage/genomicends/dumbbellgeometryunassembled.'));x=y
    clean('size','SourceAMPuredouble-sidedsize selection; ratiosunassigned')
    qc('qubit','SourceQubitdumbbellquantification')
    qc('bioanalyzer','SourceBioanalyzerHSQC')
    op('dilute','SourceD=M*N_A*2/10**25; intended~8amolDNAperul',S.GENERIC,'Diluentidentityandabsolutevolumenotassigned','Source dilutedligation preparation; no exactmoleculeyield guarantee')
    aliquot('rca_input','Source1ul dilutedtemplate into20ul pulseRCA; restretainedoutsidechosenreaction')
    op('pulse_rca','Source92C2min(step1),92C30s(step2),60C30s(step3),65C150s(step4),goto3ninetimes;4Chold',S.EXTENSION,'BIORONstranddisplacementpolymeraseHS25U/2ulbuffer/1ul10mMdNTP/0.6ul100mMMgCl2','Source pulse-RCAproduct mixture containing concatemerizedrepresentations ofbothoriginalstrands;fullcopycounts/ends/fidelityunassembled',('P5-RCA1ul concentrationunassigned',))
    clean('rca_clean','Source1.5xAMPure/elute23ulTE')
    op('pcr','Whole23ulRCAproduct;98C30s,[98C10s65C75s]goto2eighttimes;65C5min',S.PCR,'NEBNextUltraIIQ5MasterMix25ul;listed50ultotal','SourcePCR libraryofindividualRCAcopies plusoriginaltemplates/primers/remnants; no USER imposed',('SourceP5dualindexoligo listed1ul; exacthistoricalinventoryunassembled','SourceP7dualindexoligo listed1ul; allocationwordingretained'))
    clean('pcr_clean','Source0.7xAMPure/elute30ulTE')
    qc('final_quant','SourceQubitpostPCRquantification')
    others=g.input('other_libraries','Sourceotherindependentlypreparedsamplelibraries; poolmembers/massesunassigned',K.LIBRARY,primary=True);pool=g.product('pool','Physicalsourcesequencingpool',K.MOLECULAR);g.react('pool_samples','Sourcepooling',S.POOL,((x,R.TEMPLATE),(others,R.TEMPLATE)),(pool,));x=pool
    library=g.product('library','Source sequencing preparation; exact full-library geometry remains unassembled',K.LIBRARY)
    g.react('loading','Physical source loading preparation',S.GENERIC,((x,R.TEMPLATE),),(library,),changes_molecule=False)
    g.sequence(library,SequencingPlan(limitation=READOUT))
    g.add_source_item(SourceItem('pending',DEFERRED,'Selected original sections read; linked-source/whole molecular geometry review incomplete. No accepted absence claim.'))
    return g

def selftest():
    from checks import Check
    c=Check();g=protocol_graph()
    c('PulseRCAreceivesselecteddilutionaliquot',g.steps['pulse_rca'].inputs[0].material,'rca_input_used')
    c('USERnotinsertedintoSMMroute',any('user'in k.lower() for k in g.steps),False)
    c('PCRusescleanedRCAproducts',g.steps['pcr'].inputs[0].material,'rca_clean_kept')
    c('Actualtemplate remainder retained',g.terminals['rca_input_rest'].disposition,Disposition.UNRESOLVED)
    c('Sourcehairpinwholeorderingtranscription',hairpin_sequence(),'TCTTCTACAGTNNNNNNAGATCGGAAGAGCACACGTCTGAACTCCAGTCUACACTCTTTCCCTACACGACGCTCTTCCGATCT')
    c('SourceRCAprimertranscription',RCA_PRIMER,'GTAGGGAAAGAGTGTAGACTGGAGTTC')
    from checks import Source,have
    from storage import sources_root
    import re
    src=Source(sources_root()/'smm-seq__10.1126+sciadv.abm3259'/'sciadv.abm3259_PMC8993124.xml.txt','Publishedmain136/fullhairpin')
    if have(c,src,label='Publishedwholehairpin'):
        t=src.path.read_text();raw=t.split('5′-TCTTC',1)[1].split('-3′',1)[0];raw='TCTTC'+raw;raw=re.sub(r'\s+','',raw).replace('/ideoxyU/','U');c('Wholemodifiedhairpin source transcription',hairpin_sequence(),raw)
    c('Private partial model cannot be publicly promoted',bool(g.publication_problems()),True)
    c.report()
if __name__=='__main__':selftest()
