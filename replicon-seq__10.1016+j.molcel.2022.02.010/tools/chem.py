"""Published selected preparation; private partial molecular model, not a complete audit."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,SourceItem,SequencingPlan,Terminal,Disposition
from chemdraw import Scene,Segment
from adapter_ligation import AdapterLigationGeometry
from pcr import PCRGeometry
from circular import CircularizationGeometry

TITLE='Replicon-seq selected Mcm4-MNase/BrdU route'
AUTHORITY='AssociatedauthorMSPMC8995386 Methods128–144/limitations102;publishedfinalandSI/vendororiginals unreviewed. Oneoriginal;preprint excluded.'
INITIAL='SourceW303-derivedyeastwithMcm4-MNasefusion; selectedreplicationDNAroute, not4thiouracilnascentRNAcontrol'
READOUT='SourceSQKLSK109orLSK110libraries,750ng–1ugloadedMinIONFLOMIN106R9>48h; exactkitrunmapping/adaptergeometryunassembled. DNAscentBrdUscore>=0.5 filtering is computational, notphysicalenrichment; unlabelled parentalDNAco-sequences. ExactMNaseendposition maybiasreadout perauthoredlimitation.'
DEFERRED='Published final version andsourceSI1docx (gated), historicalLSK109/110 originals, fusionconstruct/end/adapter geometry. Saleh2021citedChECcomponent notblindly inherited; noSNV/markcallingperformance guarantee.'

def inventory_state(label):return None

def unknown_scene():
    sc=Scene();sc.strand('unassembled',[Segment('geometry','[source geometry unassembled]',placeholder=True)]);return sc

def protocol_graph():
    g=ProtocolGraph(TITLE,authority=AUTHORITY,representation='published_instance')
    x=g.input('input',INITIAL,K.BIOLOGICAL,primary=True)
    def op(key,label,kind,reagent,product,oligos=()):
        nonlocal x
        enzyme=kind in {S.DIGESTION,S.END_REPAIR,S.EXTENSION,S.TAILING,S.CIRCULARIZATION,S.LIGATION,S.PCR,S.DEPROTEINIZATION,S.TAGMENTATION,S.CHEMICAL_CONVERSION} or key=='wall'
        c=g.input(key+'_stock',reagent,K.CATALYST if enzyme else K.CONDITION)
        inputs=[(x,R.TEMPLATE),(c,R.CATALYST if enzyme else R.CONDITION)]
        if key=='wall':
            buf=g.input('wall_buffer','Source250ulZymolyasebuffer:1Msorbitol/1mMEGTA/10mMbetaME/50mMTrispH7.5',K.CONDITION);inputs.append((buf,R.CONDITION))
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
    op('label','SourceYPD25C/OD6000.35;alpha-factor2h45min;BrdU400ugml30minbeforereleaseand30minrelease',S.MAGIC,'Alpha-factorhourly5/5/1ugml;BrdUinYPD','SourceS-phasecellswithnascentBrdU-DNAandunlabelledparentalDNA; no perfectlabelspecificity')
    clean('harvest','100mlculture4000rpm1min/iceEDTA20mM-EGTA1mMwash/re-pellet/liquidnitrogenfreeze;thawice10minEDTA-EGTA/10000g30sec')
    clean('wash','Three1mlSolutionAwashes;15mMTrispH7.5/80mMKCl/0.1mMEGTA/PIC/spermidine/spermine')
    op('permeabilize','570ulA30C1min+30ul2percentdigitonin30C5min/ice2min',S.PERMEABILIZATION,'Source digitonin','Permeabilizedsourcecellpreparation; notpurifiednascentDNA')
    op('cut','AddCaCl2to2mM;ice10seconds',S.DIGESTION,'Source alreadycell-associatedMcm4MNase (notanaddedfreeMNasestock);addedCaCl2','Sourcefork-proximalcutDNAincludingdaughterandparentalendsplusuncutDNA; exactstrandendsnotpredicted')
    op('quench','Add610ul2xstopbuffer',S.QUENCHING,'400mMNaCl/20mMEDTA/8mMEGTA','Quenchedsourcecell/DNApreparation')
    clean('pellet','10000g20secsourcecellpellet')
    op('wall','250ul30C3min',S.PERMEABILIZATION,'Zymolyase100T5mg','Sourcespheroplastmixturepluswallfragments; no implicitDNAextraction')
    clean('spheroplast_wash','TwoZymolyasebufferwashes/oneRNasebufferwash')
    op('rna','500ulbuffer+10ulRNasecocktail37C2h',S.DIGESTION,'RNAsecocktail','SourceDNA/RNAcleavageproductsandmatrix')
    op('protein','SDS2percent/50mgChelex/10ul20mgmlPK55C1h',S.DEPROTEINIZATION,'ProteinaseK/SDS/Chelex','SourceDNA/proteinfragments/RNAremnants/resin/reagents')
    clean('extract','Twophenolchloroformphaselockextractions/gentleinversion/13000rpm5min/aqueouspourtransfer;ethanolprecipitate/resuspendovernight40ulwater')
    op('kit','Selected1.2–1.8ugDNAhistoricalLSK109orLSK110repair/endprep/A-tail/adaptercomponent',S.MAGIC,'Source10min20Crepair/endprep;30minRTA-tail;AMPureelution37C30min;completekitinternalsunreviewed','SourcePCR-freeNanoporeDNApreparation; sourceforkcleavageisnotposthoccomputationalfragmentation')
    library=g.product('library','Source sequencing preparation; exact full-library geometry remains unassembled',K.LIBRARY)
    g.react('loading','Physical source loading preparation',S.GENERIC,((x,R.TEMPLATE),),(library,),changes_molecule=False)
    g.sequence(library,SequencingPlan(limitation=READOUT))
    g.add_source_item(SourceItem('pending',DEFERRED,'Selected original sections read; linked-source/whole molecular geometry review incomplete. No accepted absence claim.'))
    return g

def selftest():
    from checks import Check
    c=Check();g=protocol_graph()
    c('MNaseonlyafterpermeabilization',g.steps['cut'].inputs[0].material,'permeabilize_product')
    c('Stopbeforecellpelleting',g.steps['pellet'].inputs[0].material,'quench_product')
    c('RNAcleavagebeforeproteinworkup',g.steps['protein'].inputs[0].material,'rna_product')
    c('Physicalextractionafterproteolysis',g.steps['extract'].inputs[0].material,'protein_product')
    c('SourceKitafterDNAphysicalrecovery',g.steps['kit'].inputs[0].material,'extract_kept')
    c('Wallhydrolysis is not DNA/RNA nuclease digestion',g.steps['wall'].kind,S.PERMEABILIZATION)
    c('Zymolyase catalyst role explicit',g.steps['wall'].inputs[1].role,R.CATALYST)
    c('Private partial model cannot be publicly promoted',bool(g.publication_problems()),True)
    c.report()
if __name__=='__main__':selftest()
