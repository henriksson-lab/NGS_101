"""Published ID-seq two-stage reporter barcoding; incomplete source-local model."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
import illumina as il
from chemdraw import Construct,Scene,Segment,feature,revcomp,panel,workflow_from_protocol_graph,workflow_panel
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,Terminal,Disposition,SourceItem,SequencingPlan
from pcr import PCRGeometry
from page import head,info
PAGE_FILE='chem.html'
TITLE='ID-seq — source producer and two-stage reporter barcoding (partial)'
SOURCE='Published DOI10.1038/s41467-018-04761-0, PMC6008431 Methods sample barcoding; published SI Note1/Fig1 and SupplementaryData3'
LIMIT='Private partial model: antibody conjugation, modified released-dsDNA ends, complete PCR product strands and sequencing-primer/run geometry remain unassembled. A production template oligo is not silently equated with its manufactured antibody-linked dsDNA product.'
WORKBOOK='s41467-018-04761-0_supp_41467_2018_4761_MOESM6_ESM.xlsx.txt'
FORWARD=il.P5[:19]
ANCHOR='ATCAGTCAACAGATAAGCG'
PRODUCTION_REVERSE='TCGCTTATCTGTTGACTGAT'
WELL_BARCODE='AGAGTCGCGT'
WELL_HANDLE='GATCGGAAGAGCACACGTC'
PLATE_PRIMER_INDEX='ATCAGT'

def production_template():
    return Construct([Segment('source P5 prefix',il.P5[:25]),Segment('Read1 region',il.TRUSEQ_READ1),
        Segment('FLD0001 diversity spacer','A'),Segment('reporter UMI','N'*15,placeholder=True,feature=feature('idseq.umi','umi','random')),
        Segment('FLD0001 antibody barcode','GTATCGTCGT',feature=feature('idseq.antibody','feature_barcode','fixed')),
        Segment('template anchor',ANCHOR)],'Source production template FLD0001; not released dsDNA')

def well_primer():
    return Construct([Segment('second reverse handle',revcomp(WELL_HANDLE)),
        Segment('FLD0151 well barcode reverse strand',revcomp(WELL_BARCODE),feature=feature('idseq.well','sample_index','fixed')),
        Segment('production reverse segment',PRODUCTION_REVERSE)],'Source well reverse primer FLD0151')

def plate_primer():
    return Construct([Segment('P7',il.P7),Segment('source plate25 index in primer orientation',PLATE_PRIMER_INDEX,'idx',feature=feature('idseq.plate','sample_index','fixed')),
        Segment('source truncated Read2 region',il.TRUSEQ_READ2[:-1])],'Source plate reverse primer25')

def unknown_scene():
    sc=Scene();sc.strand('unknown',[Segment('released dsDNA ends and complete product','[unassembled molecular geometry]',placeholder=True)])
    sc.note('unknown','Source-scope assembly boundary; not a claim of missing disclosure.')
    return sc

def cleanup(g,key,x,label):
    sel=g.input(key+'_selection',label,K.SELECTION)
    kept=g.product(key+'_kept','Retained DNA after '+label,K.MOLECULAR)
    waste=g.product(key+'_waste','Physically removed '+key+' fractions',K.WASTE)
    g.react(key,label,S.CLEANUP,((x,R.TEMPLATE),(sel,R.SELECTION)),(kept,waste));g.discard(waste,'Removed only at this source-described separation')
    return kept

def protocol_graph():
    g=ProtocolGraph(TITLE,authority=SOURCE,representation='published_instance')
    # Actual source physical producer interfaces; full modified-dsDNA/read assembly remains a required source item below.
    templates=g.input('production_templates','Actual source production-oligo family; source reference components shown separately, no study member or N3-end position assigned',K.OLIGO,primary=True,chemical_state={'identity':'Source UMI15 and antibodyBarcode10; exact manufactured termini and selected family members unassembled'})
    manufactured=g.product('manufactured_dna','Actual N3-dsDNA production output; full modified end/body relationship unassembled',K.MOLECULAR)
    g.react('produce_dna','Source-reported N3-dsDNA production interface',S.MAGIC,((templates,R.TEMPLATE),),(manufactured,),reason='Defining main positively reports N3-dsDNA production and cites detailed2016 original. Exact source production chemistry, N3 position and modified duplex ends remain unassembled; no predecessor recipe inherited.')
    def reaction(key,x,label,kind,extras=()):
        y=g.product(key+'_product','Actual '+label+' preparation; complete chemical structures/side products unassembled',K.MOLECULAR)
        g.react(key,label,kind,((x,R.TEMPLATE),)+tuple(extras),(y,),changes_molecule=kind is not S.GENERIC);return y
    def mol(key,label):return(g.input(key,label,K.MOLECULAR),R.STRUCTURAL)
    def cond(key,label):return(g.input(key,label,K.CONDITION),R.CONDITION)
    dna_tco=reaction('dna_functionalize',manufactured,'DBCO-PEG12-TCO functionalization of N3-dsDNA',S.CHEMICAL_CONVERSION,(mol('dbco_tco','Actual DBCO-PEG12-TCO molecular reagent'),cond('dna_environment','Actual source functionalization environment; exact formulation unassembled')))
    antibodies=g.input('antibody_family','Actual antibody family for source reporter conjugates',K.MOLECULAR)
    antibody_tz=reaction('antibody_functionalize',antibodies,'Antibody NHS-s-s-PEG4-tetrazine functionalization at1:10',S.CHEMICAL_CONVERSION,(mol('tetrazine_linker','Actual NHS-s-s-PEG4-tetrazine disulfide-linker reagent'),cond('bbs','Actual50mM BBS pH8.4/150mMNaCl')))
    conjugated=g.product('conjugated','Actual antibody–dsDNA conjugate family plus reaction material; protein-DNA complex, not pure oligo',K.MOLECULAR)
    g.react('conjugate','Purified functionalized Ab:DNA4:1,4hroomtemperature',S.CHEMICAL_CONVERSION,((dna_tco,R.TEMPLATE),(antibody_tz,R.STRUCTURAL),cond('conjugation_buffer','Actual BBS pH8.4/roomtemperature4h')),(conjugated,))
    quenched=reaction('conjugate_quench',conjugated,'Quench with excess3,6-diphenyltetrazine',S.QUENCHING,(mol('quench_reagent','Actual excess3,6-diphenyltetrazine'),cond('quench_condition','Source quenching operation, no extra reagent dose')))
    # Source reports pooling the conjugate family before this affinity purification; count/individual members unassigned.
    pooled=reaction('conjugate_pool',quenched,'Source family-array conjugate pooling before affinity purification; full member allocation unassembled',S.GENERIC)
    affinity=cleanup(g,'conjugate_affinity',pooled,'Actual pooled conjugate ProtA/G capture overnight/BBS and PBS washes')
    eluate=reaction('conjugate_elute',affinity,'Source0.1McitratepH2.3 elution then immediate Tris8.8 neutralization',S.CHEMICAL_CONVERSION,(mol('citrate','Actual citrate elution contribution'),mol('neutralization_tris','Actual Tris-HCl neutralization contribution'),cond('elution_condition','Source pH2.3 then8.8 neutralization sequence')))
    purified=cleanup(g,'zeba_1',eluate,'Actual first Zeba spin exchange into PBS7.4')
    purified=cleanup(g,'zeba_2',purified,'Actual second Zeba spin exchange into PBS7.4')
    conj_qc=g.product('conjugate_qc','Actual source conjugate-DNA agarose analytical check; aliquot amount unassigned',K.QC)
    conjugate_stock=g.product('conjugate_stock','Actual conjugate stock after analytical allocation',K.MOLECULAR)
    g.react('conjugate_qc_allocate','Actual analytical check apart from staining stock',S.SPLIT,((purified,R.TEMPLATE),),(conj_qc,conjugate_stock));g.mark_terminal(Terminal(conj_qc.id,Disposition.QC,'Source actual gel check, not a purity calculation'))
    cells=g.input('cells','Actual keratinocyte well collection; individual wells/count/genotype assignments unassembled',K.BIOLOGICAL)
    fixed=reaction('fix',cells,'Source4%PFA/PBS15minRT and PBS washes',S.CHEMICAL_CONVERSION,(mol('pfa','Actual paraformaldehyde molecular contribution'),cond('fix_condition','Actual source15minRT fixation')))
    blocked=reaction('cell_block',fixed,'Permeabilize/block30min with source blocking mixture',S.BINDING,(mol('cell_block_mix','Actual0.5xprotein-freeblockingbuffer/PBS/0.1%Triton/200ngmlsssalmonDNA carrier'),))
    stained=g.product('stained','Actual well-array antibody–DNA staining material; bound/unbound fractions unassigned',K.MOLECULAR)
    g.react('stain','Actual0.1µg/mlconjugate staining overnight4°C',S.BINDING,((conjugate_stock,R.TEMPLATE),(blocked,R.ANALYTE)),(stained,))
    washed=cleanup(g,'stain_wash',stained,'Source PBS3short/3x15min/3short washes')
    first=g.product('first_well','Actual one washed well contribution before release',K.MOLECULAR);other=g.product('other_wells','Actual remaining washed well collection; count/identities unassigned',K.MOLECULAR)
    g.react('well_array','Actual distinct well contributions BEFORE release and PCR; not copies of a completed product',S.SPLIT,((washed,R.TEMPLATE),),(first,other))
    releases={}
    for key,well in (('a',first),('b',other)):
        released=g.product(key+'_released','Actual released reporter-DNA stock:20–50µl per source well; entire end/linker chemistry unassembled',K.MOLECULAR)
        residual=g.product(key+'_cell_residual','Actual residual cell/protein material after collected reporter release',K.MOLECULAR)
        g.react(key+'_release','Actual10mMDTT/BBS8.4,90minRT; mix every30min and collect DNA',S.CHEMICAL_CONVERSION,((well,R.TEMPLATE),mol(key+'_dtt','Actual source DTT molecular reductant'),cond(key+'_release_condition','Source20–50µlperwell release-buffer allocation')),(released,residual));g.mark_terminal(Terminal(residual.id,Disposition.UNRESOLVED,'Residual physical cell/protein fraction; source disposal unassigned'))
        selected=g.product(key+'_pcr_input','Actual8–15µlperwell release-stock allocation into PCR1',K.MOLECULAR);stock=g.product(key+'_release_stock','Untransferred released-DNA stock class; perwell volume unassigned',K.MOLECULAR)
        g.react(key+'_release_allocate','Allocate actualPCR1 template before independent well reactions',S.SPLIT,((released,R.TEMPLATE),),(selected,stock));g.mark_terminal(Terminal(stock.id,Disposition.UNRESOLVED,'Source untransferred release stock'));releases[key]=selected
    fw=g.input('forward','Source19nt forward primer',K.OLIGO,state=Construct([Segment('forward',FORWARD)]))
    enzyme=g.input('pfu','Actual PCR1 well-array PFU polymerase allocations; stock packaging/count unassigned',K.CATALYST)
    nts=g.input('pcr1_nts','Actual PCR1 dNTP allocations at0.2mM',K.MOLECULAR)
    buffer_stock=g.input('pcr1_buffer_stock','Actual PCR1 PFU buffer allocations; same-stock salts/Tris/Triton/BSA components',K.MOLECULAR)
    buffer=g.product('pcr1_buffer','Copresent actual PCR1 PFU buffer condition components',K.CONDITION)
    bsa=g.product('pcr1_bsa','Copresent actual PFU-buffer BSA0.1mg/ml; no fresh added dose',K.MOLECULAR)
    g.react('pcr1_buffer_inventory','Components of the same actual PFU-buffer allocations; no physical separation',S.GENERIC,((buffer_stock,R.TEMPLATE),),(buffer,bsa),changes_molecule=False)
    branches=[]
    for key in ('a','b'):
        x=releases[key]
        rv=g.input(key+'_well_primer','Distinct well-specific reverse primer; FLD0151 inventory is an example, not assigned to this symbolic well',K.OLIGO)
        spike=g.input(key+'_spike','Source spike-in DNA barcodes; identities not assigned',K.MOLECULAR)
        y=g.product(key+'_barcoded','Well-barcoded reporters and source spike-ins in '+key,K.MOLECULAR)
        g.react(key+'_well_pcr','Actual25µlPCR1 per independent well class; source repeat instruction retained, exact cycle interpretation unassigned',S.PCR,
            ((x,R.TEMPLATE),(spike,R.SPIKE_IN),(fw,R.STRUCTURAL),(rv,R.STRUCTURAL),(enzyme,R.CATALYST),(nts,R.STRUCTURAL),(buffer,R.CONDITION),(bsa,R.STRUCTURAL)),(y,),
            geometry=PCRGeometry.unavailable(unknown_scene(),reason='Full manufactured reporter and well-specific product strands unassembled; printed primer inventory is separate.'))
        branches.append(y)
    pool=g.product('pool','Pooled well-barcoded products from one plate',K.MOLECULAR)
    g.react('plate_pool','Pool well-specific PCR products per plate',S.POOL,tuple((x,R.TEMPLATE) for x in branches),(pool,))
    clean=cleanup(g,'first_column',pool,'Qiagen PCR purification column')
    exo=g.input('exoi','ExonucleaseI',K.CATALYST)
    dig=g.product('digested','Pooled DNA with primer digestion products still present',K.MOLECULAR)
    g.react('primer_digest','ExonucleaseI37C30min',S.DIGESTION,((clean,R.TEMPLATE),(exo,R.CATALYST)),(dig,))
    heat=g.input('heat','80C20min inactivation',K.CONDITION)
    inactive=g.product('inactive','Heat-inactivated reaction; digestion products not physically discarded',K.MOLECULAR)
    g.react('exo_inactivation','Heat-inactivate ExonucleaseI',S.THERMAL_TREATMENT,((dig,R.TEMPLATE),(heat,R.CONDITION)),(inactive,))
    allocated=g.product('plate_pcr_input','Actual15–17µlpostExoIstock allocation into PCR2',K.MOLECULAR)
    remainder=g.product('plate_pcr_stock','UntransferredpostExoIstock; initialcolumn30µl elution plus lateradditions, actualvolumeunassigned',K.MOLECULAR)
    g.react('plate_pcr_allocate','ActualPCR2template allocation; do not normalize unknown postdigestvolume',S.SPLIT,((inactive,R.TEMPLATE),),(allocated,remainder));g.mark_terminal(Terminal(remainder.id,Disposition.UNRESOLVED,'Actual untransferred plate stock'))
    rv=g.input('plate_reverse','Actual plate reverse-primer class; plate25 remains reference inventory, experimental index member unassigned',K.OLIGO)
    pol2=g.input('pcr2_pfu','Actual PCR2 PFU1µl source allocation; no repeated first-stage dose',K.CATALYST)
    nt2=g.input('pcr2_nts','Actual PCR2 dNTP0.2mM contribution',K.MOLECULAR)
    buffer2_stock=g.input('pcr2_buffer_stock','Actual PCR2 PFU buffer allocation',K.MOLECULAR)
    buf2=g.product('pcr2_buffer','Copresent actual PCR2 PFU buffer condition components',K.CONDITION)
    bsa2=g.product('pcr2_bsa','Copresent molecular BSA component of same actual PCR2 buffer, no extra dose',K.MOLECULAR)
    g.react('pcr2_buffer_inventory','Components of the same actual PCR2 buffer allocation; no physical separation',S.GENERIC,((buffer2_stock,R.TEMPLATE),),(buf2,bsa2),changes_molecule=False)
    indexed=g.product('indexed','Plate-index PCR reaction; complete duplex unassembled',K.MOLECULAR)
    g.react('plate_pcr','PCR2 adds plate index and second adapter',S.PCR,((allocated,R.TEMPLATE),(fw,R.STRUCTURAL),(rv,R.STRUCTURAL),(pol2,R.CATALYST),(nt2,R.STRUCTURAL),(buf2,R.CONDITION),(bsa2,R.STRUCTURAL)),(indexed,),
        geometry=PCRGeometry.unavailable(unknown_scene(),reason='Plate25 is separate source inventory only; actual historical plate-index member, full pooled product/end assembly and selected experimental well remain unassigned.'))
    second=cleanup(g,'second_column',indexed,'Qiagen PCR purification after PCR2')
    finaldig=g.product('final_digest','Second ExonucleaseI reaction including digestion products',K.MOLECULAR)
    g.react('final_primer_digest','ExonucleaseI45min',S.DIGESTION,((second,R.TEMPLATE),(exo,R.CATALYST)),(finaldig,))
    selected=cleanup(g,'pkis_size_selection',finaldig,'PKIS branch only: Zymo selection for fragments greater than150bp')
    out=cleanup(g,'minelute',selected,'Final Qiagen MinElute purification')
    qa=g.product('final_qa','Actual source BioAnalyzer analytical allocation from10µl finalMinElute; amountunassigned',K.QC)
    stock=g.product('final_stock','Actual library stock after analytical allocation; whole molecules unassembled',K.LIBRARY)
    g.react('qa_allocate','Actual selected BioAnalyzer QA branch; sourcealternativeagarose branch not merged',S.SPLIT,((out,R.TEMPLATE),),(qa,stock));g.mark_terminal(Terminal(qa.id,Disposition.QC,'Actual analyticalallocation, no imaginedproductionreturn'))
    loading=g.product('loading','Actual NextSeq500 loading contribution; runendnumber/cyclesunassigned',K.LIBRARY)
    rest=g.product('library_stock','Untransferredpreparedlibrary stock; amount/fateunassigned',K.LIBRARY)
    g.react('loading_allocate','Actual source sequencing allocation',S.SPLIT,((stock,R.TEMPLATE),),(loading,rest));g.mark_terminal(Terminal(rest.id,Disposition.UNRESOLVED,'Prepared library remainder'))
    g.sequence(loading,SequencingPlan(acquisition='unspecified',limitation='SourceNextSeq500 and approximately185bp DNAproducts. Numberofacquiredreadends, cyclecounts/indexmask/fullprimerlanding unassembled; no scID-run inheritance.'))
    g.declare_reviewed_scope(REVIEWED_TARGETS,SOURCE)
    g.add_source_item(SourceItem('unassembled_source','Remaining source and full reporter/read assembly',LIMIT))
    return g

REVIEWED_TARGETS=tuple("""production_templates produce_dna manufactured_dna dbco_tco dna_environment dna_functionalize dna_functionalize_product antibody_family tetrazine_linker bbs antibody_functionalize antibody_functionalize_product conjugation_buffer conjugate conjugated quench_reagent quench_condition conjugate_quench conjugate_quench_product conjugate_pool conjugate_pool_product conjugate_affinity_selection conjugate_affinity conjugate_affinity_kept conjugate_affinity_waste citrate neutralization_tris elution_condition conjugate_elute conjugate_elute_product zeba_1_selection zeba_1 zeba_1_kept zeba_1_waste zeba_2_selection zeba_2 zeba_2_kept zeba_2_waste conjugate_qc_allocate conjugate_qc conjugate_stock cells pfa fix_condition fix fix_product cell_block_mix cell_block cell_block_product stain stained stain_wash_selection stain_wash stain_wash_kept stain_wash_waste well_array first_well other_wells
a_dtt a_release_condition a_release a_released a_cell_residual a_release_allocate a_pcr_input a_release_stock b_dtt b_release_condition b_release b_released b_cell_residual b_release_allocate b_pcr_input b_release_stock
forward pfu pcr1_nts pcr1_buffer_stock pcr1_buffer_inventory pcr1_buffer pcr1_bsa a_well_primer a_spike a_well_pcr a_barcoded b_well_primer b_spike b_well_pcr b_barcoded plate_pool pool first_column_selection first_column first_column_kept first_column_waste exoi primer_digest digested heat exo_inactivation inactive plate_pcr_allocate plate_pcr_input plate_pcr_stock plate_reverse pcr2_pfu pcr2_nts pcr2_buffer_stock pcr2_buffer_inventory pcr2_buffer pcr2_bsa plate_pcr indexed second_column_selection second_column second_column_kept second_column_waste final_primer_digest final_digest pkis_size_selection_selection pkis_size_selection pkis_size_selection_kept pkis_size_selection_waste minelute_selection minelute minelute_kept minelute_waste qa_allocate final_qa final_stock loading_allocate loading library_stock""".split())

def production_binding_reference():
    # Source template's19ntanchor pairs19of20production-reverse bases; no manufactured-end identity is asserted.
    template=production_template();sc=Scene();sc.strand('template',list(template),label='Printed FLD0001 template reference, not manufactured modified dsDNA')
    sc.anneal('production reverse',[Segment('source reverse5prime extra T',PRODUCTION_REVERSE[:1]),Segment('source reverse3prime anchor',PRODUCTION_REVERSE[1:])],to='template',pair=('source reverse3prime anchor','template anchor'),unpaired=('source reverse5prime extra T',))
    return sc

def render_page():
    g=protocol_graph()
    def inventory(x):
        sc=Scene();sc.strand('oligo',x.segments,label=x.name);return sc.rows()
    return '\n'.join([head(TITLE),'<main class="wrap"><h1>'+TITLE+'</h1>',info(LIMIT),'<p><a href="01_id_seq.html">Research notes</a></p>',
        info('Actualphysicalproducerextended, but required fullreporter/readassembly sourceitem remains unmet. No complete-protocol marker or publication promotion.'),panel(production_binding_reference().rows(),caption='Known source reference relation:19bp anchor match and one unpaired5′T on production reverse primer; actual manufactured N3-end/duplex and full read assembly remain unresolved.'),workflow_panel(workflow_from_protocol_graph(g,'loading')),panel(inventory(production_template()),caption='Verbatim production-template component example FLD0001; not manufactured dsDNA geometry'),
        panel(inventory(well_primer()),caption='Source well-primer FLD0151 inventory; source-field orientation preserved'),'</main>'])

def selftest():
    from checks import Check
    from storage import sources_root
    c=Check();g=protocol_graph()
    c('Private incomplete model',g.publishable,False)
    c('No false completion marker','data-protocol-complete=' in render_page(),False)
    c('Well barcoding uses actualrelease-stock allocation',g.steps['a_well_pcr'].inputs[0].material,'a_pcr_input')
    c('Perplate pooling follows both independent wellPCRs',[i.material for i in g.steps['plate_pool'].inputs],['a_barcoded','b_barcoded'])
    c('Digest debris remains through heat inactivation',g.steps['exo_inactivation'].inputs[0].material,'digested')
    p=sources_root()/Path(__file__).resolve().parents[1].name/WORKBOOK
    if p.exists():
        source=''.join(p.read_text().split())
        for x in (production_template(),well_primer(),plate_primer()):c('Published Data3 transcription '+x.name,x.top() in source,True)
        c('Published forward transcription',FORWARD in source,True)
    else:c.skip('Published Data3 transcriptions','External original unavailable')
    production_binding_reference()
    c('Required fullassembly source item remains unmet',g.publication_problems(),["authoritative source item 'unassembled_source' is not accounted for"])
    c('Actualplateindex remains unassigned',g.materials['plate_reverse'].state,None)
    c('ReleaseDNA is physically produced',g.producer_of('a_released').id,'a_release')
    c('Collected DNA and residual cells stay distinct',g.terminals['a_cell_residual'].disposition,Disposition.UNRESOLVED)
    c.report(exit_on_fail=True)
