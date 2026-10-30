"""Published selected preparation; private partial molecular model, not a complete audit."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,SourceItem,SequencingPlan,Terminal,Disposition
from chemdraw import Scene,Segment
from adapter_ligation import AdapterLigationGeometry
from pcr import PCRGeometry
from circular import CircularizationGeometry

TITLE='Circle-damage-seq — selected genomic CPD route'
AUTHORITY='Published Science Advances2021 PMC8324051 main Methods139–185; original PDF p11 adapters inspected. Deaminated branch documented in note, not silently merged.'
INITIAL='Selected 1ug genomic DNA from UVB1000J/m2 fibroblasts, isolated immediately; genomic CPD route'
READOUT='Source HiSeq2500 PE150. Separate computational 2x80 hard clipping and divergent alignment analysis do not remove physical molecules. Source historical i5/i7 exact inventories and whole read geometry unassembled.'
DEFERRED='Published supplementary materials not in fetched figure-only zip; historical NEB adapter/primer originals; deaminated CPD and E.coli DpnI pilot variants remain distinct.'

from illumina import NEBNEXT_ARM_READ2,TRUSEQ_READ1

def adapter_sequences():
    # Published right arm deletes one C in the canonical C run; do not silently canonicalize.
    right=TRUSEQ_READ1[:11]+TRUSEQ_READ1[12:]
    t=NEBNEXT_ARM_READ2+'U'+right
    return (t,t[:-1]+'C')

def adapter_orderings():
    return tuple({'sequence':seq,'five_prime':'P','deoxy_u_positions':(seq.index('U'),),'phosphorothioate_after':(len(seq)-2,)} for seq in adapter_sequences())

def adapter_inventory():
    sc=Scene()
    for name,seq in zip(('printed T variant','printed C variant'),adapter_sequences()):
        sc.strand(name,[Segment(name,seq,note='Published internal deoxyU; phosphorothioate bond before final nucleotide; full folding unassembled')],mod5='P')
    return sc

def unknown_scene():
    sc=Scene();sc.strand('unassembled',[Segment('geometry','[source geometry unassembled]',placeholder=True)]);return sc

def protocol_graph():
    g=ProtocolGraph(TITLE,authority=AUTHORITY,representation='published_instance')
    x=g.input('input',INITIAL,K.MOLECULAR,primary=True)
    def op(key,label,kind,reagent,product,oligos=()):
        nonlocal x
        enzyme=kind in {S.DIGESTION,S.END_REPAIR,S.EXTENSION,S.TAILING,S.CIRCULARIZATION,S.LIGATION,S.PCR,S.DEPROTEINIZATION,S.CHEMICAL_CONVERSION}
        c=g.input(key+'_stock',reagent,K.CATALYST if enzyme else K.CONDITION)
        inputs=[(x,R.TEMPLATE),(c,R.CATALYST if enzyme else R.CONDITION)]
        for i,a in enumerate(oligos):
            o=g.input(key+'_oligo_'+str(i),a,K.OLIGO,state=adapter_inventory() if key=='adapter' else None,chemical_state=adapter_orderings() if key=='adapter' else None);inputs.append((o,R.STRUCTURAL))
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
    op('shear','CovarisE220 4C to300–400bp',S.FRAGMENTATION,'Source sonication','Genomic fragments, lesions and unfragmented DNA')
    clean('shear_clean','DCC5 cleanup/elute42ul')
    op('repair','50ul repair24C30min/70C20min',S.END_REPAIR,'RNaseH2U/T4DNApolymerase3U/T4PNK10U/dNTPs/T4ligasebuffer','End-repair mixture including repaired/unprocessed DNA and RNA fragments')
    clean('repair_clean','DCC5 elute100ul')
    qc('dna_quant','Source NanoDrop concentration assay; allocation unspecified')
    op('circle','Selected1ug blunt DNA;200ul16Covernight',S.CIRCULARIZATION,'T4ligase1600U','Circularized dsDNA and residual linear/unclosed DNA; exact distribution unassigned')
    clean('circle_clean','DCC5 elute42ul')
    op('linear_digest','50ul37C30min/70C20min',S.DIGESTION,'PlasmidSafeATPDependentDNase20U','Source circle-enriched mixture plus linear cleavage products; no perfect digestion asserted')
    clean('linear_clean','1.8xAMPure90ul/elute35ul')
    op('lesion_open','Selected CPD circles40ul37C20min',S.DIGESTION,'T4PDG10U/APE110U/NEBuffer4/BSA','CPD-opened circular DNA and other processed/unprocessed species; cleavage coordinates unassembled')
    clean('open_clean','1.8xAMPure72ul/elute25ul')
    op('photorepair','50ul;365nmUVA15cm90minRT',S.CHEMICAL_CONVERSION,'E.coliphrBphotolyase0.75ug/UVA','Source photolyase repair mixture; no quantitative lesion-removal guarantee')
    clean('photo_clean','1.8xAMPure90ul/elute32ul')
    op('s1','40ul4minRT',S.DIGESTION,'S1nuclease5U','Source gap-cleaved DNA plus fragments and remaining substrates')
    op('quench','EDTA/Tris then70C10min',S.QUENCHING,'2ul0.5MEDTA/1ul1MTrispH8','Quenched cleavage mixture')
    clean('s1_clean','Source1.8xAMPure72ul/elute48ul; stated volume retained without reinterpretation')
    op('atail','55ul37C30min/70C30min',S.TAILING,'Klenowexo-minus10U/dATP0.2mM/NEBNextbuffer','Source dA-tailed DNA mixture')
    op('adapter','20C15min; source states93ul although listed55+30+1+3=89ul',S.LIGATION,'UltraIILigationMasterMix30ul/Enhancer1ul', 'Adapter-ligated and unligated DNA mixture; source T/C adapter alternatives not resolved into one perfect junction',('Source1.5uM adapter3ul: both printed5P/deoxyU/terminalPS T-or-C inventories in note',))
    op('user','USER3U37C20min',S.DIGESTION,'USERenzyme','Source deoxyU-cleaved adapter mixture including cleavage fragments and unprocessed substrates')
    clean('size','1xAMPure96ul then0.5x300–1000bp selection/elute25ul')
    op('pcr','12cycles98C10s/65C75s;98C30s initial/65C5min final; stated50ul but listed55ul',S.PCR,'SourceQ5MasterMix25ul','Source amplified library mixture',('Sourcei5primer2.5ul, exact historical identity unassembled','Sourcei7primer2.5ul, exact historical identity unassembled'))
    clean('pcr_clean','1xAMPure50ul/elute20ul')
    qc('qubit','SourceQubit3QC')
    qc('bioanalyzer','SourceBioanalyzerQC')
    qc('kapa','SourceKAPAlibraryquantificationQC')
    library=g.product('library','Source sequencing preparation; exact full-library geometry remains unassembled',K.LIBRARY)
    g.react('loading','Physical source loading preparation',S.GENERIC,((x,R.TEMPLATE),),(library,),changes_molecule=False)
    g.sequence(library,SequencingPlan(limitation=READOUT))
    g.add_source_item(SourceItem('pending',DEFERRED,'Selected original sections read; linked-source/whole molecular geometry review incomplete. No accepted absence claim.'))
    return g

def selftest():
    from checks import Check
    c=Check();g=protocol_graph()
    c('CPD opening before photolyase',g.steps['photorepair'].inputs[0].material,'open_clean_kept')
    c('S1 follows photorepair cleanup',g.steps['s1'].inputs[0].material,'photo_clean_kept')
    c('Chemical USER follows adapter ligation',g.steps['user'].inputs[0].material,'adapter_product')
    c('PCR follows size selection',g.steps['pcr'].inputs[0].material,'size_kept')
    c('Diagnostic QC remains separate',g.terminals['kapa_assay'].disposition,Disposition.QC)
    c('Printed T adapter transcription',adapter_sequences()[0],'GATCGGAAGAGCACACGTCTGAACTCCAGTCUACACTCTTTCCTACACGACGCTCTTCCGATCT')
    c('Printed C adapter alternative retained',adapter_sequences()[1][-1],'C')
    c('Internal DNA uracil retained',adapter_sequences()[0].count('U'),1)
    c('Private partial model cannot be publicly promoted',bool(g.publication_problems()),True)
    c.report()
if __name__=='__main__':selftest()
