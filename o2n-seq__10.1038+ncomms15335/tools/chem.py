"""Published selected preparation; private partial molecular model, not a complete audit."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,SourceItem,SequencingPlan,Terminal,Disposition
from chemdraw import Scene,Segment
from adapter_ligation import AdapterLigationGeometry
from pcr import PCRGeometry
from circular import CircularizationGeometry

TITLE='o2n-seq — nick-controlled single-strand circle copying'
AUTHORITY='Published NatureCommunications2017 PMC5458117 main Methods121–151/Results65–67/analysis165–172; published SI Fig1 caption read. Full SI graphical geometry deferred rather than inferred from caption.'
INITIAL='Source1–4ug genomic DNA; selected standard uncaptured o2n route, distinct from comparison Cir-seq/AmpliSeq'
READOUT='Source o2n HiSeq2500 PE125, tandem-copy PE consensus after computational initial5base/internaladapter trimming. Single original strand circle copied, not proof of independent original duplex agreement. ComparatorCir-seqPE250 and AmpliSeqPE150 are not adopted.'
DEFERRED='Published SI Fig1 original graphical layout; exact second-strand primer/o2n-blocker and NEXTflex514102 historical primer/adapter geometry; third-party vendor originals; optional0.42MbRoche capture is separate from selected uncaptured route.'

FIRST_ADAPTER=('GATCAGTCGTACGTGCTTACTCTCAATAGCAGCTT','GTGGGCAGTCGGTGAACGACTGAUCT')
def inventory_chemistry(label):
    for seq in FIRST_ADAPTER:
        if seq in label:return {'sequence':seq,'five_prime':'P','deoxy_u_positions':tuple(i for i,b in enumerate(seq) if b=='U')}
    return None
def inventory_state(label):
    for seq in FIRST_ADAPTER:
        if seq in label:
            sc=Scene();sc.strand('source adapter',[Segment('DNA adapter',seq,note='Source DNA with internal deoxyU where printed; Y/end geometry unassembled')],mod5='P');return sc
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
    op('shear','CovarisS220130ulAE;600s4C; source selected gel route',S.FRAGMENTATION,'Duty10percent/intensity5/100cyclesburst','Source approximately100bpDNA fragments and other sizes')
    clean('size','OligoCleanConcentrator then4percentagarose80V70min/60–120bp bands/MinEluteGelExtraction6xQG;900sno-gel alternative not merged')
    op('endprep','65ul20C30min/65C30min',S.END_REPAIR,'NEBNextUltraE7442S EnzymeMix3ul/buffer6.5ul/DNA55.5ul','Source repaired/dA-tailed DNA mixture')
    op('adapter','Source add15ulE7445Smix/2uladapter/1ulenhancer20C30min/65C5min/ice',S.LIGATION,'NEBNextUltraLigationE7445S','Source first-adapter ligated and unligated DNA; exact sealed geometry unassembled',('5P-GATCAGTCGTACGTGCTTACTCTCAATAGCAGCTT','5P-GTGGGCAGTCGGTGAACGACTGAUCT; internalU is source DNA deoxyU nicking substrate'))
    clean('lig_clean','MinEluteReactionCleanup;65C5min/ice;1.8xAMPure')
    qc('dsquant','Qubit2dsDNAHS assay')
    op('denature','Selected30ulDNA<300ng95C3min/ice3min',S.THERMAL_TREATMENT,'Thermaldenaturation','Source separated strands and residual duplexes; complete separation not guaranteed')
    op('circle','CircLigaseII60C14h/25C2h/95C2min/ice3min',S.CIRCULARIZATION,'3ul10xbuffer/1.5ul50mMMnCl2/1.5ulCircLigaseII stated100U','Source ssDNA circles and unclosed substrates; final volume arithmetic retained as printed')
    op('linear_digest','37C1h',S.DIGESTION,'ExonucleaseI1ul/ExonucleaseIII1ul/Fpg1ul','Source circles and processed/unprocessed linear/damaged strands plus fragments')
    clean('circle_clean','MinEluteReactionCleanup3xERC')
    qc('ssquant','SourceQubitssDNA assay')
    op('primer_anneal','DNA15.8ul/primer1ul/NEBuffer4 2ul;95C3min/45C5min/ice3min',S.HYBRIDIZATION,'Thermalannealing', 'Source primer-associated circle mixture; exact source primer identity not assigned',('Source second-strand primer1ul, exact sequence unassembled',))
    op('second_strand','25C30min/75C20min',S.EXTENSION,'T4DNApolymerase0.5ul/dNTP0.5ul2.5mM/BSA0.2ul100x','Source circle-directed second-strand products and unprocessed substrates')
    op('nick','USER1ul37C30min/50C4min/ice',S.DIGESTION,'USERenzyme','Source nicked circle preparation plus unprocessed species; no imported generic RCA product')
    clean('nick_clean','1.5xAMPure')
    op('displace','20ul60C1h',S.EXTENSION,'DNA16ul/Isothermalbuffer2ul/dNTP1ul2.5mM/Bst2WarmStart1ul','Source nick-controlled strand-displacement products; full end geometry/copy count distribution unassembled')
    clean('displace_clean','1xAMPure')
    op('endprep2','20C30min/65C30min',S.END_REPAIR,'RecoveredDNA18.5ul/buffer2.17ul/EnzymeMix1ul','Source repaired/dA-tailed displacement product mixture')
    op('adapter2','20C30min',S.LIGATION,'BluntTAligase5ul/enhancer0.33ul','Source NEXTflex-ligated/unligated DNA',('BiooNEXTflex514102 barcodeadapter0.17ul; exact index unassigned',))
    clean('adapter2_clean','0.8xAMPure')
    op('pcr','98C45s;5–11cycles98C15s/65C30s/72C1min;72C4min/4Chold',S.PCR,'KAPAHiFiHotStartReadyMix25ul in50ul', 'Source amplified sequencing library mixture',('NEXTflex514102 PrimerMix1ul forward component; historical exactsequence unassembled','NEXTflex514102 PrimerMix1ul reverse component; same physical mixture, not a second1ul addition'))
    clean('final_clean','Source twice0.8xAMPure or2percentgel size selection; chosen branch unassigned')
    library=g.product('library','Source sequencing preparation; exact full-library geometry remains unassembled',K.LIBRARY)
    g.react('loading','Physical source loading preparation',S.GENERIC,((x,R.TEMPLATE),),(library,),changes_molecule=False)
    g.sequence(library,SequencingPlan(limitation=READOUT))
    g.add_source_item(SourceItem('pending',DEFERRED,'Selected original sections read; linked-source/whole molecular geometry review incomplete. No accepted absence claim.'))
    return g

def selftest():
    from checks import Check
    c=Check();g=protocol_graph()
    c('Linear digest immediately follows ss circularization',g.steps['linear_digest'].inputs[0].material,'circle_product')
    c('Second strand follows annealing',g.steps['second_strand'].inputs[0].material,'primer_anneal_product')
    c('USER nicking precedes strand displacement',g.steps['displace'].inputs[0].material,'nick_clean_kept')
    c('Historical second adapter separate from first',g.steps['adapter2'].inputs[0].material,'endprep2_product')
    c('First adapter source DNA sequence',FIRST_ADAPTER[0],'GATCAGTCGTACGTGCTTACTCTCAATAGCAGCTT')
    c('Nicking arm retains DNA uracil',FIRST_ADAPTER[1],'GTGGGCAGTCGGTGAACGACTGAUCT')
    c('Private partial model cannot be publicly promoted',bool(g.publication_problems()),True)
    c.report()
if __name__=='__main__':selftest()
