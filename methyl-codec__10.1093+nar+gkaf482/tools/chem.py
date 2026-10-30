"""Published selected preparation; private partial molecular model, not a complete audit."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,SourceItem,SequencingPlan,Terminal,Disposition
from chemdraw import Scene,Segment
from adapter_ligation import AdapterLigationGeometry
from pcr import PCRGeometry
from circular import CircularizationGeometry

TITLE='Methyl-CODEC selected hydroxy-dCTP route'
AUTHORITY='Published2025 PMC12135180 main75–99/111/141;publishedSI Methodsfirstpage andoriginalTableS1p9 inspected; originalp8FigS7 viewed while locating table, no chemistry adopted;2originals. No genericCODEC constructor inherited.'
INITIAL='SourceNA12878genomicDNA;20ngfragmentedDNAusedforlibrary; selected5hydroxymethyl-dCTP-copy route, otheranalogues/HCT116controls independent'
READOUT='SelectedNA12878sourceNovaSeqS4PE169,noindexcycles; HCT116NovaSeqXPE161/eightindexcycles is a separate source cohort. U in chemical DNA is reported as T afterPCR. Pairedcopy/genetic/methylation inference doesnotguarantee two-product recovery orperfectconversion.'
DEFERRED='Full12TableS1 oligo chemical orderings andquadruplex/USER/endgeometry unassembled although originaltablevisuallyread; source16ulpostconversioneluate vs20ulPCRinput unresolved; historicalkit andCODECprecursor originals notreviewed.'

def inventory_state(label):return None

def unknown_scene():
    sc=Scene();sc.strand('unassembled',[Segment('geometry','[source geometry unassembled]',placeholder=True)]);return sc

def protocol_graph():
    g=ProtocolGraph(TITLE,authority=AUTHORITY,representation='published_instance')
    x=g.input('input',INITIAL,K.MOLECULAR,primary=True)
    def op(key,label,kind,reagent,product,oligos=()):
        nonlocal x
        enzyme=kind in {S.DIGESTION,S.END_REPAIR,S.EXTENSION,S.TAILING,S.CIRCULARIZATION,S.LIGATION,S.PCR,S.DEPROTEINIZATION,S.TAGMENTATION,S.CHEMICAL_CONVERSION}
        c=g.input(key+'_stock',reagent,K.CATALYST if enzyme else K.CONDITION)
        inputs=[(x,R.TEMPLATE),(c,R.CATALYST if enzyme else R.CONDITION)]
        for i,a in enumerate(oligos):
            o=g.input(key+'_oligo_'+str(i),a,K.OLIGO,state=inventory_state(a));inputs.append((o,R.STRUCTURAL))
        kw={}
        if kind==S.MAGIC:kw['reason']='Actual published component; full source-specific molecular transition unassembled.'
        if kind==S.CIRCULARIZATION:kw['geometry']=CircularizationGeometry.unavailable('Source intramolecular closure known; full substrate and closure geometry unassembled, not a claim of source absence.')
        if kind==S.LIGATION:kw['geometry']=AdapterLigationGeometry('unavailable','unavailable',(unknown_scene(),),'Source ligation operation known; full ends and sealed junction unassembled.')
        if kind==S.PCR:kw['geometry']=PCRGeometry.unavailable(unknown_scene(),reason='Source primer inventories/whole amplicon geometry not assembled in this selected model.')
        y=g.product(key+'_product',product,K.MOLECULAR);g.react(key,label,kind,inputs,(y,),**kw);x=y
    def clean(key,label):
        nonlocal x
        c=g.input(key+'_condition',label,K.CONDITION);y=g.product(key+'_kept','Physical source recovered preparation; carryover and recovery unquantified',K.MOLECULAR);w=g.product(key+'_waste','Actual source physically separated fractions',K.WASTE)
        g.react(key,label,S.CLEANUP,((x,R.TEMPLATE),(c,R.CONDITION)),(y,w));g.discard(w,'Physical workup, not a perfect molecular specificity assertion.');x=y
    def qc(key,label):
        nonlocal x
        y=g.product(key+'_retained','Production preparation retained after independent analytical allocation',K.MOLECULAR);q=g.product(key+'_assay',label,K.QC);g.react(key,label,S.SPLIT,((x,R.TEMPLATE),),(y,q));g.mark_terminal(Terminal(q.id,Disposition.QC,'Analytical aliquot; not returned to production.'));x=y
    op('shear','Covarismean150bp',S.FRAGMENTATION,'Source ultrasonication','DNAfragments with unfragmentedDNA/endsunassigned')
    clean('shear_clean','SourceAMPurecleanup/elutelowTE')
    op('endprep','50ulDNA+7ulbuffer+3ulenzyme20C30min65C30min',S.END_REPAIR,'NEBNextUltraIIE7645 EndPrep','Source end-repaired/dA-tailed mixture')
    clean('end_clean','Source1.8xAMPure/elute60ul')
    op('ligate','20C1h',S.LIGATION,'30ulLigationMasterMix/1ulEnhancer/3ul5primeDeadenylase','Sourcequadruplexadapter-ligatedandunligatedDNA; fullgeometryunassembled',('2.5ul500nMMethylCODECquadruplex;eacholigo5uM40ullowTE100mMNaCl85C3min→20C→RT12h; sourceTableS1mC/5P/terminalPS/DDD preserved as knownunassembled inventory',))
    clean('ligation_clean','Source0.9xAMPure/elute40ul')
    op('user','USER1.5ul37C30min',S.DIGESTION,'USERenzyme','SourceUSER-cleaved adapter substrates plusunprocessedDNA/fragments/reagents; exactcleavagegeometryunassembled')
    op('copy','Source41.5ulUSERmixture+10ulbuffer+1ulPhi29+47ulwater+2ul10mMdNTP;30C20min',S.EXTENSION,'Phi29M0269L/dATP/dGTP/dTTP/5hydroxymethyl-dCTP equally','Source strand-displaced protected copies linked to originalDNA and otherproducts; no universal perfectcopy/duplex yield guarantee')
    clean('copy_clean','Source0.75xAMPure/elute28ul')
    op('protect','TET2/BGToxidation37C1h',S.CHEMICAL_CONVERSION,'10ulTET2buffer/1uloxidationsupplement/1ulBGT/printed1DTT(unitunassigned)/4ulTET2/5ul1to1250FeIIsulfate','SourceTET2/BGTmodifiednative/protectedcytosineclasses plus unprocessedDNA/reagents')
    op('stop','StopReagent1ul37C30min',S.QUENCHING,'SourceStopReagent','Stoppedoxidationmixture')
    clean('oxidation_clean','Source0.75xAMPure/elute16ul')
    op('denature','Add4ulformamide;85C10min/ice',S.THERMAL_TREATMENT,'Formamide/thermaldenaturation','Source20uldenaturedDNAmixture; bothstrandclassesremain')
    op('deaminate','Add68ulwater/10ulbuffer/1ulAPOBEC/1ulBSA;37C3h',S.CHEMICAL_CONVERSION,'SourceAPOBEC/E7125S','SourceconvertednativeDNAcontaininguracil and protectedcopyclasses; conversionfractions/offproducts unassigned')
    clean('converted_clean','Source0.75xAMPure/elute16ul; mainPCRstates20ulinput,volume bridgeunassigned')
    op('pcr','Source30sinitialdenaturation(temperatureunassigned);8–10cycles98C20s60C30s65C2min;65C5min',S.PCR,'Source25ulQ5UHotStartM0515 inlisted50ulreaction','SourcePCRDNAreportingchemicalUasT; nativeprotectedmethylationinformation encodedbyconversion,notretainedasidenticalPCRDNAmarks',('KAPAKK2623forwardprimer component ofone5ulsourcepremix','KAPAKK2623reverseprimer component ofsame5ulsourcepremix;notsecondaddition'))
    clean('pcr_clean','Sourceconsecutive0.65xAMPuretwice/final20uleluate')
    qc('qubit','SourceQubitdsDNAHSQC')
    qc('agilent','SourceAgilentHighSensitivityDNAQC')
    library=g.product('library','Source sequencing preparation; exact full-library geometry remains unassembled',K.LIBRARY)
    g.react('loading','Physical source loading preparation',S.GENERIC,((x,R.TEMPLATE),),(library,),changes_molecule=False)
    g.sequence(library,SequencingPlan(limitation=READOUT))
    g.add_source_item(SourceItem('pending',DEFERRED,'Selected original sections read; linked-source/whole molecular geometry review incomplete. No accepted absence claim.'))
    return g

def selftest():
    from checks import Check
    c=Check();g=protocol_graph()
    c('USERafterligationcleanup',g.steps['user'].inputs[0].material,'ligation_clean_kept')
    c('Phi29usesUSERproduct',g.steps['copy'].inputs[0].material,'user_product')
    c('Oxidationbeforedeamination',g.steps['denature'].inputs[0].material,'oxidation_clean_kept')
    c('PCRusesconvertedDNA',g.steps['pcr'].inputs[0].material,'converted_clean_kept')
    c('QCisindependent',g.terminals['agilent_assay'].disposition,Disposition.QC)
    c('Private partial model cannot be publicly promoted',bool(g.publication_problems()),True)
    c.report()
if __name__=='__main__':selftest()
