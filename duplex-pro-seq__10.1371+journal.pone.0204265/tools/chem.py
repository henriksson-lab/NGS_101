"""Published selected preparation; private partial molecular model, not a complete audit."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,SourceItem,SequencingPlan,Terminal,Disposition
from chemdraw import Scene,Segment
from adapter_ligation import AdapterLigationGeometry
from pcr import PCRGeometry
from circular import CircularizationGeometry

TITLE='Duplex Proximity Sequencing — targeted droplet route'
AUTHORITY='Published PLoSONE2018 PMC6168144 main Results70–99/Methods194–226 and S2Fig caption284–286. DNA Pro-Seq, not RNA run-on PRO-seq. S1/S2/S3workbook headers administratively inspected; four-original cap overrun recorded.'
INITIAL='Source non-denatured plasma cfDNA preparation, desired genomic template copies selected; never predenatured to discard duplex pairing'
READOUT='Source MiSeq overlapping PE reads, customR1/index1/index2/R2; index2 before turnaround and omit dark cycles. Intensity files/non-pass-filter reads retained. DS index checks and mixed-base masking are computational, not physical product purification. Whole Genome PEG-loop single-extension route is distinct and unmodeled.'
DEFERRED='Published S2Fig original TIFF gated; actual gene-specific primer and custom read-primer identities/whole PEG complex geometry unassembled. No additional source acquisition after administrative four-original overrun.'

from illumina import P7,P5
# The first printed arm is 3prime-to-5prime; reversal, not reverse complementation, gives ordering orientation.
LINKED_LEFT_3TO5='AGTCAACGTCGTCTTCTTCCACATCTAGAGCCACCAGCGGCATAGTAA'
LINKED_RIGHT_5TO3=P5+'TCCCTCCTATCATGGACAC'
UNIVERSAL_REVERSE=P7+'GCGACGGTTAGACGAACGGTACG'
def linked_arms_5to3():return (LINKED_LEFT_3TO5[::-1],LINKED_RIGHT_5TO3)
def inventory_chemistry(label):
    return {'linker':'Sp9 PEG','connection':'5prime-to-5prime','arms_5to3':linked_arms_5to3()} if 'two-arm5primePEG' in label else None
def inventory_state(label):
    sc=Scene()
    if 'two-arm5primePEG' in label:
        for name,seq in zip(('left arm','right arm'),linked_arms_5to3()):
            sc.strand(name,[Segment(name,seq,note='Source two-arm inventory; both5prime ends joined by nonnucleotide Sp9 PEG, not a DNA concatemer')],mod5='PEG Sp9')
        return sc
    if 'universalreverse400nM' in label:
        sc.strand('reverse',[Segment('universal reverse',UNIVERSAL_REVERSE)]);return sc
    return None

def unknown_scene():
    sc=Scene();sc.strand('unassembled',[Segment('geometry','[source geometry unassembled]',placeholder=True)]);return sc

def protocol_graph():
    g=ProtocolGraph(TITLE,authority=AUTHORITY,representation='published_instance')
    x=g.input('input',INITIAL,K.MOLECULAR,primary=True)
    def op(key,label,kind,reagent,product,oligos=()):
        nonlocal x
        enzyme=kind in {S.DIGESTION,S.END_REPAIR,S.EXTENSION,S.TAILING,S.CIRCULARIZATION,S.LIGATION,S.PCR,S.DEPROTEINIZATION}
        c=g.input(key+'_stock',reagent,K.CATALYST if enzyme else K.CONDITION)
        inputs=[(x,R.TEMPLATE),(c,R.CATALYST if enzyme else R.CONDITION)]
        for i,a in enumerate(oligos):
            o=g.input(key+'_oligo_'+str(i),a,K.OLIGO,state=inventory_state(a),chemical_state=inventory_chemistry(a));inputs.append((o,R.STRUCTURAL))
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
    clean('extract','Source plasma2000g10min4C/QIAampCirculatingNucleicAcidKit;2x100ul0.1xTE10min/20000g3min;Monarchcleanup/elute15ul0.1xTDTE')
    qc('copy_quant','SourceCOG5/ALB qPCR input-copy assay with serial dilutions and reference standard')
    op('partition','Source40ul droplet preparation; zero-or-one template occupancy expectation, not exact occupancy',S.MAGIC,'RainDropSource/ThunderBolts; printed approximately8,000,0000droplets retained as source typography','Droplet compartments with source templates, background DNA and full reaction components')
    op('pcr','Q5dropletPCR98C30s;38x98C10s/60C30s/72C30s;1C/s;4Chold',S.PCR,'Q5HotStart0.02Uul/dNTP0.2mM/Q5buffer/RDTstabilizer','PEG-linked two-primer complexes plus unlinked strands, unamplified background and reaction components; no guarantee each complex contains both original senses',('Gene-specific index1forward25nM and index2forward25nM; sequences unassigned','Gene-specificreverse50nM, exact sequence unassigned','Source universalreverse400nM, literal in note','Source200nM two-arm5primePEGSp9-linked universal primer; literal orientation in note'))
    op('break','Source62.5ul destabilizer; source emulsions broken',S.MAGIC,'RainDance destabilizationreagents','Released linked complexes/unlinked DNA and emulsion remnants')
    clean('post_break','AMPure0.8:1/elute20ul0.1xTE')
    op('unlinked_digest','50ul37C1h/70C20min',S.DIGESTION,'T7exonuclease6.7U/RecJf41.7U/NEBuffer4','Source linked-complex-enriched mixture and digestion fragments; no perfect specificity or duplex certainty imposed')
    clean('digest_clean','AMPure1.6:1/elute25ul0.1xTE')
    qc('kapa','SourceKAPALibraryQuantKitQC')
    op('denature','18ullibrary+2ul1NNaOH5minRT thenice',S.THERMAL_TREATMENT,'SourceNaOH','Denatured sequencing preparation retaining PEG-linked structure and reaction components')
    phix=g.input('phix','Source separately denatured PhiX,5percent libraryconcentration',K.MOLECULAR,primary=True);pool=g.product('phix_pool','Physical loading pool with PhiX',K.MOLECULAR);g.react('phix_add','SourcePhiX addition',S.POOL,((x,R.TEMPLATE),(phix,R.TEMPLATE)),(pool,));x=pool
    op('neutralize','Source2ul1NHCl/HT1to600ul,5.5pMfinal',S.QUENCHING,'SourceHCl/HT1buffer','Neutralized loading preparation')
    library=g.product('library','Source sequencing preparation; exact full-library geometry remains unassembled',K.LIBRARY)
    g.react('loading','Physical source loading preparation',S.GENERIC,((x,R.TEMPLATE),),(library,),changes_molecule=False)
    g.sequence(library,SequencingPlan(limitation=READOUT))
    g.add_source_item(SourceItem('pending',DEFERRED,'Selected original sections read; linked-source/whole molecular geometry review incomplete. No accepted absence claim.'))
    return g

def selftest():
    from checks import Check
    c=Check();g=protocol_graph()
    c('Duplex retained before droplet PCR',g.steps['pcr'].inputs[0].material,'partition_product')
    c('Unlinked digestion follows emulsion removal',g.steps['unlinked_digest'].inputs[0].material,'post_break_kept')
    c('Cleavage products physically removed after digest',g.steps['digest_clean'].inputs[0].material,'unlinked_digest_product')
    c('Loading denaturation follows separate quantification',g.steps['denature'].inputs[0].material,'kapa_retained')
    c('Printed first arm is reversed not complemented',linked_arms_5to3()[0],'AATGATACGGCGACCACCGAGATCTACACCTTCTTCTGCTGCAACTGA')
    c('Printed second arm source transcription',linked_arms_5to3()[1],'AATGATACGGCGACCACCGAGATCTACACTCCCTCCTATCATGGACAC')
    c('Universal reverse source transcription',UNIVERSAL_REVERSE,'CAAGCAGAAGACGGCATACGAGATGCGACGGTTAGACGAACGGTACG')
    c('PhiX mixed before neutralization',g.steps['neutralize'].inputs[0].material,'phix_pool')
    c('Private partial model cannot be publicly promoted',bool(g.publication_problems()),True)
    c.report()
if __name__=='__main__':selftest()
