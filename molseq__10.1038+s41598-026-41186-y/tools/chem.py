"""Published MOLseq: source-owned photocleavable letters and partial RNA library.

The local linker descriptor preserves a non-nucleotide boundary. It is not an
atomic photolysis model or a claim of quantitative cleavage/ligation in cells.
"""
from dataclasses import dataclass
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Scene,Segment,revcomp,feature
from nextera import ADAPTOR_S5,ADAPTOR_S7
from protocol_graph import (ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,
                            SourceItem,SourceCoverage,SequencingPlan,Terminal,Disposition)
from pcr import PCRGeometry
HANDLE='AGAGAATGGATAGGTTGTGT'
LETTERS={'L1':HANDLE,'L2':'ATGGAAGAGATGTGTTAGGT','L3':'TAGGTTGTGTAGAGAATGGA','L4':'TGTGTTAGGTATGGAAGAGA'}
SPLINTS={'12':'TCTCTTCCATACACA','21':'TCCATTCTCTACCTA','31':'TCCATTCTCTTCCAT','42':'TCTCTTCCATTCTCT'}
SMART3='ACGAGCATCAGCAGCATACGA'
LD_REVERSE='AAGCAGTGGTATCAACGCAGAGT'
TSO='AAGCAGTGGTATCAACGCAGAGTACATrGrG+G'
ENRICH_FORWARD=ADAPTOR_S5+SMART3
ENRICH_REVERSE=ADAPTOR_S7
RT_PRIMER='/5Phos/'+HANDLE+'N'*20+'T'*30
AUTHORITY='Published PMC13136315 main46–70/100–110/124–168, actual MethodsPDFp9; published SI Figs1/2/4 and Table1 actualp6–8. Two originals; preprint and AIP excluded.'

@dataclass(frozen=True)
class PhotoLetter:
    """Ordered DNA cap / non-nucleotide iSpPC / retained letter, from SI Table1."""
    name: str
    cap: str
    body: str
    linker: str='/iSpPC/'
    def __post_init__(self):
        if not self.cap or not self.body or set(self.cap+self.body)-set('ACGT'):
            raise ValueError('Explicit DNA parts required on either side of PC')
        if self.linker!='/iSpPC/':
            raise ValueError('This source-owned descriptor only covers the cited iSpPC reagent')
    @property
    def ordering(self):
        return self.cap+self.linker+self.body
    @property
    def nucleotide_count(self):
        return len(self.cap)+len(self.body)
    @property
    def identifier_feature(self):
        return feature('molseq.'+self.name.lower(),'spatial_barcode','fixed',group='molseq.spatial',part=self.name.lower())
    def cleaved_class(self):
        return {'identifier_feature':self.identifier_feature,'retained_DNA':self.body,'retained_5prime':'P',
                'separated_cap_DNA':self.cap,'linker_fragment_chemistry':'unassembled',
                'scope':'Source SI Fig2 cleaved class only; intact/off-target/damaged classes and yields unassigned'}

def letter(name):
    return PhotoLetter(name,'T'*20,LETTERS[name])

def literal_scene(seq,mod5='',mod3=''):
    s=Scene();s.strand('oligo',[Segment('source ordering',seq)],mod5=mod5,mod3=mod3);return s

def primer_scene():
    s=Scene();s.strand('primer',[Segment('source handle',HANDLE),
        Segment('source random window','N'*20,placeholder=True,note='Purpose unassigned, not automatically UMI'),
        Segment('polyT','T'*30)],mod5='P');return s

def local_interface(incoming,prior,splint):
    """Exact declared local nick geometry; upstream PC/full cDNA remain outside view."""
    if revcomp(splint)!=incoming[-5:]+prior[:10]:
        raise ValueError('Source splint does not bridge supplied incoming3end and prior5end')
    s=Scene();s.strand('incoming',[Segment('incoming head',incoming[:-5],feature=feature('molseq.incoming_head','spatial_barcode','fixed',group='molseq.incoming',part='head')),
        Segment('incoming tail',incoming[-5:],feature=feature('molseq.incoming_tail','spatial_barcode','fixed',group='molseq.incoming',part='tail'))],
        mod5='upstream T20–iSpPC cap retained in chemical_state; partial DNA interface only',mod3='OH')
    s.anneal('splint',[Segment('prior complement',splint[:10]),Segment('incoming complement',splint[10:])],
             to='incoming',pair=('incoming complement','incoming tail'),unpaired=('prior complement',))
    s.anneal('prior',[Segment('prior head',prior[:10]),Segment('prior remainder',prior[10:])],
             to='splint',pair=('prior head','prior complement'),unpaired=('prior remainder',),mod5='P')
    s.ligation_end('incoming','incoming tail',kind='splinted',polarity='3-prime',owner='adapter')
    s.ligation_end('prior','prior head',kind='splinted',polarity='5-prime',owner='substrate')
    return s

def unknown_scene():
    s=Scene();s.strand('unknown',[Segment('unassembled','[whole molecular product unassembled]',placeholder=True)]);return s

def protocol_graph():
    g=ProtocolGraph('MOLseq published spatial transcriptome preparation',authority=AUTHORITY,representation='published_instance')
    x=g.input('cells','Source plated adjacent U2OS/MEF regions; plating/precoating outside selected graph',K.BIOLOGICAL,primary=True)
    def inp(key,label,kind=K.CONDITION,role=R.CONDITION,**kw):
        return(g.input(key,label,kind,**kw),role)
    def op(key,label,kind,inputs=(),product='Source reaction mixture retaining incomplete substrates and free reagents until workup',**kw):
        nonlocal x
        if kind==S.GENERIC:kw['changes_molecule']=False
        y=g.product(key+'_product',product,K.MOLECULAR);g.react(key,label,kind,((x,R.TEMPLATE),)+tuple(inputs),(y,),**kw);x=y
    def clean(key,label):
        nonlocal x
        y=g.product(key+'_kept','Physically retained source preparation; recovery/carryover unassigned',K.MOLECULAR)
        w=g.product(key+'_waste','Source removed wash/cleanup fractions',K.WASTE)
        g.react(key,label,S.CLEANUP,((x,R.TEMPLATE),inp(key+'_condition',label)),(y,w));g.discard(w,'Source physical separation, not perfect molecular specificity.');x=y
    def pcr(key,label,forward=None,reverse=None,enzyme='Source catalyst component; historical identity unassembled'):
        bulk=g.input(key+'_stock','Aggregate actual source PCR supplies; packaging/stockcount/quantities/internal historical formulation unassigned',K.MOLECULAR)
        e=g.product(key+'_enzyme',enzyme,K.CATALYST);n=g.product(key+'_nts','Nucleotide component of aggregate actual PCR supplies; packaging unassigned; pool/dose unassigned',K.MOLECULAR)
        f=g.product(key+'_forward','Forward-primer component of aggregate actual PCR supplies; packaging unassigned',K.OLIGO,state=literal_scene(forward)if forward else None)
        r=g.product(key+'_reverse','Reverse-primer component of aggregate actual PCR supplies; packaging unassigned',K.OLIGO,state=literal_scene(reverse)if reverse else None)
        b=g.product(key+'_buffer','Condition component of aggregate actual PCR supplies; packaging unassigned',K.CONDITION)
        g.react(key+'_components','Expose aggregate actual PCR supplies; packaging/stockcount unassigned, not an asserted prepared premix',S.GENERIC,((bulk,R.TEMPLATE),),(e,n,f,r,b),changes_molecule=False)
        op(key,label,S.PCR,((e,R.CATALYST),(n,R.STRUCTURAL),(f,R.STRUCTURAL),(r,R.STRUCTURAL),(b,R.CONDITION)),
           geometry=PCRGeometry.unavailable(unknown_scene(),reason='Full indexed source products/landing and cycle programmes unassembled; no absence claim.'))
    op('fix','Source4percentPFAin1xDPBS; duration unassigned',S.CROSS_LINKING,(inp('pfa','SourceDPBS/fixation environment; duration unassigned'),inp('pfa_reagent','Actual source4percentPFA molecular contribution',K.MOLECULAR,R.STRUCTURAL)))
    primer=g.input('rt_primer','Source1uM '+RT_PRIMER,K.OLIGO,state=primer_scene(),chemical_state={'source_ordering':RT_PRIMER,'N20_purpose':'unassigned'})
    op('hybridize','SourcepolyTprimer1uM/2xSSC37Covernight',S.HYBRIDIZATION,((primer,R.STRUCTURAL),inp('hybridization','SourceSSC/temperature')))
    clean('primer_wash','Source2x30percentformamide/2xSSC47C30min,2xDPBS,2xwater,2xSSC; remove unbound/non-specific primer')
    template=g.product('rt_template','Same retained cellularRNA component',K.MOLECULAR);p=g.product('retained_primer','Co-retained sourcepolyTprimer component; occupancy unassigned, no fresh addition',K.OLIGO,state=primer_scene())
    g.react('retained_components','Same-tube components, no physical separation',S.GENERIC,((x,R.TEMPLATE),),(template,p),changes_molecule=False);x=template
    tso=g.input('tso','Source2.5uM '+TSO,K.OLIGO,chemical_state={'source_ordering':TSO,'terminal_modifications':'rG/rG/LNA-G'})
    op('rt','SourceMaximaHMinus10Uul37Covernight withTSO, two distinct RNase inhibitors600Uml each',S.REVERSE_TRANSCRIPTION,
       ((p,R.STRUCTURAL),(tso,R.STRUCTURAL),inp('rt_nts','Source1.875mMdNTPmix(25mMeachstock); convention unassigned',K.MOLECULAR,R.STRUCTURAL),
        inp('maxima','SourceMaximaHMinus',K.CATALYST,R.CATALYST),inp('rt_conditions','Source1xMaximaRTbuffer'),inp('murine_inhibitor','Actual source murineRNaseinhibitor600Uml contribution',K.MOLECULAR,R.STRUCTURAL),inp('rnasin_plus','Actual source RNasinPlus600Uml contribution, distinct from murine inhibitor',K.MOLECULAR,R.STRUCTURAL)),product='SourcecDNA/RNA/template-switch classes and residual reagents; complete products unassembled')
    clean('rt_wash','Source2xSSCtwice/PBStwice afterRT')
    # Local descriptor checks before building the opaque source-specific schedule.
    local_interface(LETTERS['L2'],HANDLE,SPLINTS['21'])
    local_interface(LETTERS['L4'],LETTERS['L2'],SPLINTS['42'])
    letters=g.input('letter_aliquots','Source per-round letter aliquots0.4uM; separate staged supplies, not an invented one-pot mixture',K.ADAPTER,
        chemical_state={'source_photo_letters':tuple(letter(k)for k in LETTERS),'whole_spatial_schedule':'unassembled','source_orderings':tuple(letter(k).ordering for k in LETTERS)})
    splints=g.input('splint_aliquots','Source per-round splint supplies0.4uM; actual occupancy/pool membership/schedule unassembled',K.OLIGO,chemical_state={'source_table_orderings':SPLINTS})
    op('spatial_barcoding','Source iterative light-controlled splinted ligation/365nm cleavage/washes; whole round schedule unassembled',S.MAGIC,
       ((letters,R.STRUCTURAL),(splints,R.STRUCTURAL),inp('barcoding_atp','Actual source ATP2mM supplementation; stock-versus-final convention unassigned',K.MOLECULAR,R.STRUCTURAL),inp('barcoding_dtt','Actual source DTT2mM supplementation; stock-versus-final convention unassigned',K.MOLECULAR,R.STRUCTURAL),inp('quickligase','SourceQuickLigase40Uul',K.CATALYST,R.CATALYST),
        inp('barcoding_conditions','Source1xligationbuffer/.1percentTriton,RT1h dark;2xstock LITERAL132nMTris/20mMMg/15percentPEG,ATP/DTT contributions exposed separately;DMD3mWmm2/30s forspatialexperiments')),
       product='Source spatially barcoded cDNA plus failed/off-target/intact/cleaved populations; cap/linker wash fates and whole molecules unassembled',
       reason='Source light exposes5P on retained letter; exact chronological letter/splint/ROI schedule, full PC-bearing products and byproduct chemistry not assembled. Local21/42 interfaces do not establish a full co-culture barcode path.')
    tag=g.input('smart3_tag','Source finalSmartSeq3sequence '+SMART3,K.OLIGO,state=literal_scene(SMART3))
    op('final_tag','Source final ligation ofSmartSeq3sequence tobarcodedcDNA',S.MAGIC,((tag,R.STRUCTURAL),inp('final_tag_preparation','Source final ligation preparation; exact historical stock/splint unassembled',K.MOLECULAR,R.STRUCTURAL)),reason='Main168/SIFig1 source finaltag positive; whole target end/splint and sealed construct unassembled, not filled from another SmartSeq recipe.')
    solution=g.product('extracted_solution','Source50ul aspiratedDNA-containingextraction solution; recovery/carryover unassigned',K.MOLECULAR)
    residual=g.product('cell_residual','Unaspiratedfixedcell/surfacefraction; laterfateunassigned',K.MOLECULAR)
    g.react('extract','Source50ulwater/.1percentTriton93C4min,aspirate and coolRT',S.EXTRACTION,
            ((x,R.TEMPLATE),inp('extraction_condition','Sourceheatedextraction')),(solution,residual))
    g.mark_terminal(Terminal(residual.id,Disposition.UNRESOLVED,'Sourceaspiratessolution; cell/surfaceresidualnotimplicitlysequenced.'));x=solution
    pcr('ld_pcr','SourceLD-PCR perTakaraSmartseqV4 instructions; exact cycle count unassembled',SMART3,LD_REVERSE,'SourceSeqAmpDNApolymerase from same preparedPCRstock')
    clean('ld_clean','Source0.8xAMPureXP')
    op('tagment','SourceNexteraXTFC-131-1024tagmentation; exact loadedoligos/amounts/wholeproducts unassembled',S.MAGIC,
       (inp('xt_preparation','Source selectedhistoricalNexteraXTtagmentationpreparation',K.MOLECULAR,R.STRUCTURAL),),reason='Main168/SIFig1 establishes tagmentation followed by5prime selection; no modern constructor inherited.')
    pcr('enrich','Source5primeenrichmentPCR, not unbiased wholetranscriptome molecule recovery',ENRICH_FORWARD,ENRICH_REVERSE)
    clean('enrich_clean','Source0.8xAMPureXP')
    pcr('index_pcr','Source furtherPCR addsP5/i5/P7/i7; actual index identities/wholeprimers unassembled')
    clean('index_clean','Source0.8xAMPureXP')
    used=g.product('load_selected','Actual selected libraryloading contribution; amountunassigned',K.MOLECULAR)
    rest=g.product('load_remainder','Untransferred preparedlibrary stock; amountunassigned',K.MOLECULAR)
    g.react('loading_allocation','Select physical sequencing contribution without whole-stock assumption',S.SPLIT,((x,R.TEMPLATE),),(used,rest),changes_molecule=False)
    g.mark_terminal(Terminal(rest.id,Disposition.UNRESOLVED,'Untransferred productionstock'));x=used
    lib=g.product('library','SourceMOLseqsequencingpreparation; full molecules/read cycle mapping unassembled',K.LIBRARY)
    g.react('load','SourceMiSeqloading',S.GENERIC,((x,R.TEMPLATE),),(lib,),changes_molecule=False)
    g.sequence(lib,SequencingPlan(limitation='SourceMS-102-3001MiSeqv3(150-cycle) and literal150paired-endreads both retained; precise runmask unassigned. Read1 used forbarcodefiltering,Read2 separately aligned; computation does not physically remove molecules. N20purpose not automaticallyUMI. Full spatial schedules/wholePCRproducts/readlanding unresolved.'))
    g.add_source_item(SourceItem('selected',AUTHORITY,'Literal132nMTris confirmed in mainXML+originalPDFp9; intended concentration unresolved. Source run-description statements not silently repaired.'))
    return _publication_scope(g)

def selftest():
    from checks import Check,Source,have
    from storage import sources_root
    c=Check();g=protocol_graph()
    c('PC is non-nucleotide, not added to DNA length',letter('L2').nucleotide_count,40)
    c('Cleaved class retains downstream DNA with5P',letter('L2').cleaved_class()['retained_DNA'],LETTERS['L2'])
    c('Separated cap persists as a distinct class',letter('L2').cleaved_class()['separated_cap_DNA'],'T'*20)
    c('RT uses retained primer',g.steps['rt'].inputs[0].material,'rt_template')
    c('FinalTag follows actualbarcodingboundary',g.steps['final_tag'].inputs[0].material,'spatial_barcoding_product')
    root=sources_root()/Path(__file__).resolve().parents[1].name
    source=Source(root/'s41598-026-41186-y_supp_41598_2026_41186_MOESM1_ESM.pdf.txt','https://doi.org/10.1038/s41598-026-41186-y')
    if have(c,source):
        text=''.join(source.path.read_text().split())
        for index,name in enumerate(LETTERS,1):
            following=('Letter'+str(index+1)+'(L'+str(index+1)+')')if index<4 else'Letter1(L1)forfluorescentdetection'
            row=text.split('Letter'+str(index)+'('+name+')',1)[1].split(following,1)[0]
            c('Whole modified sourceletter row '+name,letter(name).ordering,row)
        rtrow=text.split('poly-Tprimerforreverse',1)[1].split('poly-Tprimerforbarcodesize',1)[0].replace('transcription','')
        c('Whole sourceRT row with interleaved column label removed',RT_PRIMER,rtrow)
        c('Whole sourceTSO',TSO in text,True)
        for name,seq in SPLINTS.items():c('Whole source splint '+name,seq in text,True)
        for name,seq in [('SmartSeq3',SMART3),('LDreverse',LD_REVERSE),('enrichForward',ENRICH_FORWARD),('enrichReverse',ENRICH_REVERSE)]:c('Whole sourceprimer '+name,seq in text,True)
    c('Loading uses selected contribution',g.steps['load'].inputs[0].material,'load_selected')
    c('Distinct actual RT inhibitor proteins contribute to RT',all(k in {i.material for i in g.steps['rt'].inputs} for k in ('murine_inhibitor','rnasin_plus')),True)
    c('Source-qualified publication gate',g.publishable,True);c.report()



PAGE_FILE='molseq.html'
REVIEWED_TARGETS=(
 'cells','pfa','pfa_reagent','fix','rt_primer','hybridization','hybridize','primer_wash_condition','primer_wash','retained_components','tso','rt_nts','maxima','rt_conditions','murine_inhibitor','rnasin_plus','rt','rt_wash_condition','rt_wash',
 'letter_aliquots','splint_aliquots','barcoding_atp','barcoding_dtt','quickligase','barcoding_conditions','spatial_barcoding','smart3_tag','final_tag_preparation','final_tag','extraction_condition','extract',
 'ld_pcr_stock','ld_pcr_components','ld_pcr','ld_clean_condition','ld_clean','xt_preparation','tagment',
 'enrich_stock','enrich_components','enrich','enrich_clean_condition','enrich_clean','index_pcr_stock','index_pcr_components','index_pcr','index_clean_condition','index_clean','loading_allocation','load','library')
def _publication_scope(g):
    g.cover_source_item(SourceCoverage('selected','modelled','spatial_barcoding'))
    g.declare_reviewed_scope(REVIEWED_TARGETS,AUTHORITY)
    return g.require_publishable()

def render_page():
    from page import head,caveat,table
    from chemdraw import panel
    from protocol_views import assemble_page,WorkflowSection,PlanLimitationSection
    g=protocol_graph()
    return assemble_page(g,[head('MOLseq'),'<main class="wrap" data-protocol-complete="true"><h1>MOLseq</h1><p>Light-controlled DNA-letter ligation records spatial selections on cDNA before transcriptome library preparation.</p><p><a href="https://doi.org/10.1038/s41598-026-41186-y">Defining published article</a> · <a href="01_molseq.html">Research note</a></p>',
      caveat('Selected published U2OS/MEF spatial-transcriptome preparation. Whole occupied cellular RNA/cDNA, chronological letter/splint/ROI schedule, complete photocleavage byproducts, final PCR molecules and sequencing-primer landings remain unassembled. Local nick views below describe source-supported reagent interfaces; they are not a reconstructed full spatial barcode or a quantitative ligation/photolysis claim.'),
      '<h2>Preparation</h2>',table(('Stage','Source operation'),(('Capture and copy','Fix cells, hybridize the source phosphorylated poly(T) primer, wash, then Maxima RT with TSO and two distinct RNase-inhibitor contributions.'),('Record spatial selections','Staged photocleavable DNA letters and splints undergo QuickLigase reactions,365nm illumination and washes; source whole ROI schedule remains compound.'),('Recover and prepare','Ligate the final SmartSeq3 tag; heat/extract the DNA-containing solution. LD-PCR,0.8xAMPure, historical NexteraXT tagmentation,5′ enrichment PCR, indexPCR and cleanups follow.'),('Read','MiSeq acquisition; Read1 supplies barcode filters and corresponding Read2 maps to the transcriptome. Those filters do not physically remove molecules.'))),
      '<h2>Known molecular parts</h2>',panel(primer_scene().rows(),caption='Whole source poly(T) primer; N20 purpose unassigned, not automatically a UMI'),
      table(('Source part','Chemical scope'),tuple((k,letter(k).ordering)for k in ('L1','L2','L3','L4'))),
      '<p>The iSpPC linker is a non-nucleotide boundary. The source cleaved class retains the downstream letter with5′ phosphate; separated cap and incompletely processed classes are not silently erased. Letter identities are spatial barcode components.</p>',
      '<h2>Published values retained</h2><p>Methods prints132nM Tris-HCl for the2x ligation buffer. This literal value is confirmed in the original PDF; its intended concentration remains unresolved and is not silently changed to mM. The source also names MiSeqv3(150-cycle) and “150 paired-end reads”; the precise run mask remains unassigned. The model does not normalize either statement into an inferred run.</p>',
      '<details><summary>Full preparation and unresolved molecular transitions</summary>',WorkflowSection('Source physical workflow'),PlanLimitationSection('library'),'</details></main>'])

if __name__=='__main__':selftest()
