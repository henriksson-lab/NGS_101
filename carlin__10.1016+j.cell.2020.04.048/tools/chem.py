"""CARLIN accepted-manuscript recovery variants;gated final primer original not substituted."""
from dataclasses import dataclass
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Scene,Segment,feature,workflow_from_protocol_graph,workflow_panel
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,Terminal,Disposition,SourceItem,SequencingPlan
from pcr import PCRGeometry
from page import head,info
PAGE_FILE='chem.html'
TITLE='CARLIN — genomic,bulk RNA and single-cell array recovery'
SOURCE='Conditional accepted author manuscript PMC7529102 Methods design/assembly193–203,bulk232–236,singlecell238–244;publishedarticleHTTP403,S1primeroriginalgated;no supersededpreprintSI primers imported'
LIMIT='Private accepted-manuscript source-version model. Exact CARLIN array/primers,source-specific modified strand/end geometry and historical Chromium kit originals remain unassembled;published main/supplement applicability unresolved. Recorder engineering precedes recovery;allele creation is not inferred from the library PCR.'
LINEAGE=feature('carlin.array','lineage_barcode','unknown',note='Inherited edited array shared by descendants;not unique current-cell or molecule identity.')
CELL=feature('carlin.cell','cell_barcode','unknown')
SCUMI=feature('carlin.sc_umi','umi','unknown')
BULKUMI=feature('carlin.bulk_umi','umi','random')

@dataclass(frozen=True)
class RecorderContext:
    target_sites:int=10
    preserved_gestalt_guide_designs:int=6
    newly_designed_guides:int=4
    def rows(self):
        sc=Scene();sc.strand('array',[Segment('CARLIN edited recorder','[source/sample-dependent array allele]',placeholder=True,feature=LINEAGE)])
        sc.note('array','Source array has10sense-oriented20bp target sites with3bpPAM/4bp linkers;fullsequence/editedlength unassembled. GFP3primeUTR recorder creation/induction is upstream of this assay.');return sc.rows()

@dataclass(frozen=True)
class BiotinPrimer:
    source_name:str='10X-CARLIN_1-bio'
    modification:str='biotin'
    position:int|None=None
    def rows(self):
        sc=Scene();sc.strand('primer',[Segment(self.source_name,'[full primer and attachment position unassembled]',placeholder=True)])
        sc.note('primer','Source biotin-bearing primer;exact attachment position/strand geometry awaits actual S1Table. Biotin is not automatically copied into all new strands.');return sc.rows()

@dataclass(frozen=True)
class BiotinProduct:
    beads:bool=False
    biotin_origin:str='source biotin-primer-derived strand;full duplex and product population unassembled'
    def rows(self):
        sc=Scene();sc.strand('mixture',[Segment('source target preparation','[heterogeneous CARLIN DNA and reaction species]',placeholder=True)],label='bead-associated preparation' if self.beads else 'source biotin-PCR preparation')
        sc.note('mixture','Biotin origin: '+self.biotin_origin+';no claim of biotin on every complementary copy or perfect capture.');return sc.rows()

def feature_inventory(single_cell=False,umi_only=False,index_width=None):
    sc=Scene()
    if not umi_only and not index_width:sc.strand('lineage',[Segment('inherited recorder','[edited allele]',placeholder=True,feature=LINEAGE)],label='Component class only;not a joined full molecule')
    if single_cell:
        sc.strand('cell',[Segment('Chromium cellID','[width/placement unassembled]',placeholder=True,feature=CELL)])
        sc.strand('umi',[Segment('Chromium UMI','[width/placement unassembled]',placeholder=True,feature=SCUMI)])
    elif umi_only:sc.strand('umi',[Segment('source bulk recovery UMI','N'*10,placeholder=True,feature=BULKUMI)],label='Known component only;full primer and placement unassembled')
    if index_width:
        sc.strand('index',[Segment('source library index','N'*index_width,'idx',placeholder=True,feature=feature('carlin.library_index','sample_index','unknown'))],label='Run index width known,actual primer/member unassembled')
    return sc

def variant_graph(variant):
    if variant not in ('single_cell_v2','single_cell_v3','bulk_dna','bulk_rna_in_vitro','bulk_rna_in_vivo'):raise ValueError('Choose a source-defined CARLIN recovery variant')
    g=ProtocolGraph(TITLE+' / '+variant,authority=SOURCE)
    def inp(key,label,kind=K.CATALYST,state=None):return g.input(key,label,kind,state=state)
    def op(key,x,label,kind,reagents=(),geometry=None,state=None):
        y=g.product(key+'_product',label+' preparation;residual species retained until physical workup',K.MOLECULAR,state=state);kw={} if geometry is None else {'geometry':geometry}
        g.react(key,label,kind,((x,R.TEMPLATE),)+tuple(reagents),(y,),**kw);return y
    def clean(key,x,label,state=None):
        sel=inp(key+'_selection',label,K.SELECTION);y=g.product(key+'_kept','Source physically retained DNA preparation',K.MOLECULAR,state=state);w=g.product(key+'_waste','Removed source workup fractions',K.WASTE)
        g.react(key,label,S.CLEANUP,((x,R.TEMPLATE),(sel,R.SELECTION)),(y,w));g.discard(w,'Source physical separation;purity/yield/carryover not quantified');return y
    def pcr(key,x,forward,reverse,label,polymerase,state=None):
        f=inp(key+'_forward',forward,K.OLIGO);r=inp(key+'_reverse',reverse,K.OLIGO)
        enzyme=inp(key+'_enzyme',polymerase)
        return op(key,x,label,S.PCR,((f,R.STRUCTURAL),(r,R.STRUCTURAL),(enzyme,R.CATALYST)),PCRGeometry.unavailable(Scene(),reason='Actual sourceS1primer original is gated;no whole polymer/primer geometry filled from preprint or another assay.'),state)
    def take(key,x,label,remainder):
        chosen=g.product(key+'_chosen',label,K.MOLECULAR);unused=g.product(key+'_unused',remainder,K.MOLECULAR)
        g.react(key,label,S.SPLIT,((x,R.TEMPLATE),),(chosen,unused));g.mark_terminal(Terminal(unused.id,Disposition.UNRESOLVED,'Untransferred source fraction not presumed disposed.'));return chosen
    if variant.startswith('single_cell'):
        cdna=g.input('cdna','Source Chromium3prime'+('v2' if variant.endswith('v2') else 'v3')+' amplified whole cDNA at sourceStep2.3;upstream capture/RT historicalkit original unreviewed',K.MOLECULAR,primary=True,state=feature_inventory(True))
        target=g.product('target_aliquot','3–5ng amplified whole-cDNA mixture for CARLIN capture',K.MOLECULAR)
        gene=g.product('gene_aliquot','Distinct source amplifiedcDNA aliquot for whole-transcriptome library',K.MOLECULAR)
        unused=g.product('unused_cdna','Unallocated source cDNA preparation',K.MOLECULAR)
        g.react('library_allocation','Allocate independent aliquots;not physical CARLIN purification',S.SPLIT,((cdna,R.TEMPLATE),),(target,gene,unused));g.mark_terminal(Terminal(unused.id,Disposition.UNRESOLVED,'Amounts outside source3–5ng target allocation not assigned disposal.'))
        kit=inp('gene_kit','Actual source Chromiumv2/v3 remaining library kit steps unreviewed',K.CONDITION)
        gene_product=g.product('gene_library','Source whole-transcriptome library;exact strands/end geometry unassembled',K.MOLECULAR)
        g.react('transcriptome_kit','Opaque historical source transcriptome library continuation',S.MAGIC,((gene,R.TEMPLATE),(kit,R.CONDITION)),(gene_product,),reason='Main explicitly names historicalv2/v3 protocol and amplifiedcDNA;actual linked kit revisions not read.No current source-selecting shared recipe substituted.')
        g.sequence(gene_product,SequencingPlan(acquired_roles=('Read 1','Index 1','Read 2'),read_lengths={'Read 1':28,'Index 1':8,'Read 2':91},limitation='Source NextSeq500 transcriptome run;complete product/primer-binding geometry and cell/UMI subdivisions unassembled.'))
        bio=inp('biotin_forward','Source10X-CARLIN_1-bio;exactmodifiedoligo not located',K.OLIGO,BiotinPrimer())
        p5=inp('p5_pr1','SourceP5-PR1;full sequence unassembled',K.OLIGO)
        enzyme=inp('sc_polymerase','Polymerase supplied with10Xv3kit used forALL source10XPCRs,includingv2capture inputs')
        x=op('target_pcr1',target,'15cycle initial biotin-targeted CARLINPCR',S.PCR,((bio,R.STRUCTURAL),(p5,R.STRUCTURAL),(enzyme,R.CATALYST)),PCRGeometry.unavailable(Scene(),reason='GatedsourceTable1 prevents full primer/binding/end reconstruction;source biotin modification class known.'),BiotinProduct())
        x=clean('target_spri',x,'0.8×SPRIselect cleanup',BiotinProduct())
        beads=inp('dynabeads','SourceDynabeads kilobaseBINDER3h roomtemperature roller',K.SELECTION)
        bound=g.product('bound','Source biotin-DNA/Dynabead retained fraction',K.MOLECULAR,state=BiotinProduct(True));w=g.product('unbound','Source supernatant/wash fractions',K.WASTE)
        g.react('biotin_capture','Capture biotinPCR products;magnetremove supernatant/wash',S.AFFINITY_CAPTURE,((x,R.TEMPLATE),(beads,R.SELECTION)),(bound,w));g.discard(w,'Actual source magnetic separation/washes;capture purity not assumed')
        half=take('half_beads',bound,'Half sourceCARLIN-tagged Dynabeads for nestedPCR','Remaining bound beads not transferred')
        nested=inp('nested_forward','Source10X-CARLIN_2;sequence/modifications unassembled',K.OLIGO)
        x=op('target_pcr2',half,'15cycle nestedPCR from bound template',S.PCR,((nested,R.STRUCTURAL),(p5,R.STRUCTURAL),(enzyme,R.CATALYST)),PCRGeometry.unavailable(Scene(),reason='Exactsource nestedprimer/duplex geometry unassembled.No copied biotin on every new strand asserted.'),BiotinProduct(True))
        x=clean('remove_beads',x,'Magnetic Dynabead removal;retain library-containing supernatant')
        x=clean('nested_spri',x,'0.8×SPRIselect purification')
        x=take('index_aliquot',x,'One tenth nestedPCR product for indexing','Untransferred nestedPCR products')
        si=inp('si_primer','SourceSIprimer;sequence unassembled',K.OLIGO)
        i7=inp('chromium_i7','Source Chromiumi7 sampleprimer;member/full sequence unassembled',K.OLIGO,state=feature_inventory(umi_only=False,index_width=8))
        x=op('index_pcr',x,'8cycle source indexingPCR',S.PCR,((si,R.STRUCTURAL),(i7,R.STRUCTURAL),(enzyme,R.CATALYST)),PCRGeometry.unavailable(Scene(),reason='Historicalsource sampleprimer/member and final CARLINproduct/priming geometry unassembled.'))
        g.sequence(x,SequencingPlan(acquired_roles=('Read 1','Index 1','Read 2'),read_lengths={'Read 1':28,'Index 1':8,'Read 2':350},limitation='Source MiSeq500cyclev2kit with10%PhiX;named acquisition preserved.Fullcell/UMI/array read spans and priming positions unassembled.'))
    elif variant=='bulk_dna':
        dna=g.input('dna','Source genomicDNA of pooled cells,up to250ng',K.MOLECULAR,primary=True,state=RecorderContext())
        primer=inp('umi_primer','SourceNGS_UIM_D_F (source spelling),10ntUMI withCARLIN_fwd2/adaptor;full oligo unassembled',K.OLIGO,state=feature_inventory(umi_only=True))
        enzyme=inp('umi_extension_reagents','Source extension polymerase/conditions unassembled;not copied from bulkRNA')
        x=op('umi_extension',dna,'Singleextension UMI tagging;five anneal/extend ramps WITHOUT repeated denaturation',S.EXTENSION,((primer,R.STRUCTURAL),(enzyme,R.CATALYST)))
        x=clean('extension_ampure',x,'AMPureXP purification;source ratio unassigned')
        x=pcr('bulk_pcr',x,'NGS2_F','NGS1_R','35cycle genomic targetPCR','Source gDNAPCR enzyme identity not assigned in inspectedparagraph')
        x=pcr('bulk_index',x,'P5','NGS2R_I#','10cycle genomic indexingPCR','Source bulkindexPCR reagent identity unassembled')
        g.sequence(x,SequencingPlan(acquired_roles=('Read 1','Index 1','Read 2'),read_lengths={'Read 1':250,'Index 1':6,'Read 2':250},limitation='Bulksection sourceMiSeq250/6/250;exact gDNAfinalworkup/primer/indexmember geometry unassembled.SourceRNA-only cleanup ratios are not silently assigned here.'))
    else:
        rna=g.input('rna','Source pooled-cell totalRNA,up to1ug;upstream isolation outsidegraph',K.MOLECULAR,primary=True,state=RecorderContext())
        primer=inp('bulk_rt_primer','RT_Bulk_CARLIN gene-specific source10ntUMI primer;fulloligo unassembled',K.OLIGO,state=feature_inventory(umi_only=True))
        ssiii=inp('ssiii','SourceSuperScriptIII')
        x=op('bulk_rt',rna,'Gene-specific UMI reverse transcription',S.REVERSE_TRANSCRIPTION,((primer,R.STRUCTURAL),(ssiii,R.CATALYST)))
        if variant.endswith('in_vitro'):
            x=clean('rt_cleanup',x,'Sourceinvitro1×AMPure cDNA purification')
            x=pcr('rna_pcr',x,'NGS2_F','NGS1_Bulk_R','35cycle invitroRNA targetPCR','SourcePlatinumTaq')
            x=clean('rna_pcr_cleanup',x,'0.8×AMPure cleanup')
            x=pcr('rna_index',x,'P5','NGS2R_I#','10cycle invitroRNA indexingPCR','Source indexingPCR mixture unassembled')
        else:
            x=clean('rt_cleanup',x,'Sourceinvivo purifiedcDNA;exact initialcleanup ratio not assigned')
            x=pcr('rna_pcr1',x,'NGS1_F','NGS1_Bulk_R','15cycle invivoRNA firstPCR','SourceQ5HighFidelity')
            x=clean('rna_pcr1_cleanup',x,'0.8×AMPure cleanup')
            x=take('rna_half',x,'Half firstPCR product for secondPCR','Untransferred firstPCR material')
            x=pcr('rna_pcr2',x,'NGS2_F','NGS1_Bulk_R','15cycle invivoRNA nestedPCR','SourceQ5HighFidelity')
            x=clean('rna_pcr2_cleanup',x,'0.8×AMPure cleanup')
            x=pcr('rna_index',x,'P5','NGS2R_I#','9cycle invivoRNA indexingPCR','Source indexingPCR mixture unassembled')
        x=clean('final_rna_cleanup',x,'Source bulkRNAfinal0.8×AMPure purification')
        g.sequence(x,SequencingPlan(acquired_roles=('Read 1','Index 1','Read 2'),read_lengths={'Read 1':250,'Index 1':6,'Read 2':250},limitation='SourcebulkRNA MiSeq500cyclev2kit,10%PhiX;full molecular priming/index geometry unassembled.10ntbulkUMI is independent of Chromiumcell/UMI identity.'))
    g.declare_reviewed_scope(tuple(g.materials)+tuple(g.steps),SOURCE);g.add_source_item(SourceItem('actual_source_table','Actual final-version CARLINTableS1/array/primer original',LIMIT));return g

def protocol_graph():return variant_graph('single_cell_v2')

def render_page():
    variants=('single_cell_v2','single_cell_v3','bulk_dna','bulk_rna_in_vitro','bulk_rna_in_vivo')
    def panels(v):
        g=variant_graph(v);ends=('index_pcr_product','gene_library') if v.startswith('single_cell') else ('bulk_index_product',) if v=='bulk_dna' else ('final_rna_cleanup_kept',)
        return ''.join(workflow_panel(workflow_from_protocol_graph(g,e)) for e in ends)
    return '\n'.join([head(TITLE),'<main class="wrap"><h1>'+TITLE+'</h1>',info(LIMIT),'<p><a href="01_carlin.html">Research notes</a></p>',
        ''.join('<h2>'+v+'</h2>'+panels(v) for v in variants),'</main>'])

def selftest():
    from checks import Check
    from storage import sources_root
    c=Check();g=protocol_graph();dna=variant_graph('bulk_dna');vitro=variant_graph('bulk_rna_in_vitro');vivo=variant_graph('bulk_rna_in_vivo');v3=variant_graph('single_cell_v3')
    c('Allvariants private',any(x.publishable for x in (g,dna,vitro,vivo,v3)),False);c('No finished marker','data-protocol-complete=' in render_page(),False)
    c('TargetedPCR uses wholecDNA aliquot',g.steps['target_pcr1'].inputs[0].material,'target_aliquot')
    c('Biotin belongs to actual source primer state',g.materials['biotin_forward'].state.modification,'biotin')
    c('NestedPCR takes only sourcehalfbeads',g.steps['target_pcr2'].inputs[0].material,'half_beads_chosen')
    c('Magnet removal precedes postnestedSPRI',g.steps['nested_spri'].inputs[0].material,'remove_beads_kept')
    c('IndexPCR takes onlyonetenth',g.steps['index_pcr'].inputs[0].material,'index_aliquot_chosen')
    c('gDNAUMItagging is extension notPCR',dna.steps['umi_extension'].kind,S.EXTENSION)
    c('InvivoRNA nestedPCR uses allocatedhalf',vivo.steps['rna_pcr2'].inputs[0].material,'rna_half_chosen')
    c('InvitroRNA doesnot silentlyuseQ5nestedroute','rna_pcr2' in vitro.steps,False)
    p=sources_root()/Path(__file__).resolve().parents[1].name/'j.cell.2020.04.048_PMC7529102.xml.txt'
    if p.exists():
        text=' '.join(p.read_text().split())
        c('Selectedsource explicitlyspecifiesno-repeatdenaturation','five cycles without a denaturalization step' in text,True)
        c('Selectedsource biotincapture andhalfbeads','Half of the CARLIN-tagged Dynabeads' in text and '3h at room temperature' in text,True)
    else:c.skip('Source operation transcriptions','External original unavailable')
    c.report()
