"""MPX source-qualified selected preparation, with explicit molecular limits."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from protocol_graph import ProtocolGraph, MaterialKind as K, InputRole as R, StepKind as S, SourceItem,SourceCoverage, SequencingPlan, Terminal, Disposition
from chemdraw import Scene,Segment
from adapter_ligation import AdapterLigationGeometry
from circular import CircularizationGeometry
from pcr import PCRGeometry
TITLE='Molecular Pixelation selected main-article route'
AUTHORITY='PublishedPMC11166577 main85/183–239 and SI MOESM1 Steps1–4/Appendix1. Third original MOESM3 administrative header only; chemistry unused. Capped pass, not a complete source audit.'
READOUT='SourceNextSeq2000PE44+78,0.65nM/15percentPhiX. Eight-nucleotide sample indexes named but exact index-read masks/primer landing unassembled. AOC UMI, protein identity and two DNA-pixel UPIs have different meanings; UPIs identify molecular neighborhoods, not predetermined cell barcodes. Pixelator spatial inference and filtering are computational, not physical removal.'
SOURCE_DIFFERENCES='Publishedmain209/213/217 uses30min AOC staining,15000cells and30min USER; publishedSI1 Steps2/3 uses45min AOC staining,20000cells and15min USER. This graph selects main values; SI-only PCR cleanup allocation is separately identified. No claim of a corrected recipe or experimental equivalence.'

def unknown_state():
    sc=Scene();sc.strand('unassembled',[Segment('unassembled','[source-known operation; molecular geometry unassembled]',placeholder=True)]);return sc

# Literal curator-reviewed selected published preparation, not graph enumeration.
REVIEWED_INPUTS = '''cells fix_stock fix_wash_buffer block_stock block_wash_buffer aocs stain_buffer aoc_wash_buffer secondary_stock secondary_wash_buffer
 a_template a_padlock a_ligase a_ex a_rca_mix a_sap b_template b_padlock b_ligase b_ex b_rca_mix b_sap
 a_hyb_buffer a_gfp a_hyb_wash_buffer a_gap_fill_ligate_stock a_gap_wash_buffer user_stock user_wash_buffer
 b_hyb_buffer b_gfp b_hyb_wash_buffer b_gap_fill_ligate_stock b_gap_wash_buffer lambda_stock pcr_stock pcr_oligo_0 pcr_oligo_1 ampure1_buffer ampure2_buffer phix loading_buffer'''.split()
REVIEWED_OPERATIONS = '''fix fix_wash block block_wash stain aoc_wash secondary secondary_wash a_circle a_exonuclease a_circle_aliquot a_rca_components a_rca a_rsap a_pixel_aliquot b_circle b_exonuclease b_circle_aliquot b_rca_components b_rca b_rsap b_pixel_aliquot pixelation_input a_hybridize a_hyb_wash a_gap_fill_ligate_components a_gap_fill_ligate a_gap_wash user user_wash b_hybridize b_hyb_wash b_gap_fill_ligate_components b_gap_fill_ligate b_gap_wash cell_count ex_input lambda pcr_components pcr cleanup_input ampure1 ampure2 qubit dilute interface library'''.split()

def protocol_graph():
    g=ProtocolGraph(TITLE,authority=AUTHORITY,representation='published_instance')
    x=g.input('cells','Source cell suspension in PBS, selected fixed-cell assay; no capping/rituximab/uropod perturbation silently imposed',K.BIOLOGICAL,primary=True)
    def op(key,label,kind,stock,product,oligos=()):
        nonlocal x
        enzymatic=kind in {S.DIGESTION,S.CIRCULARIZATION,S.LIGATION,S.EXTENSION,S.PCR,S.MAGIC}
        supplied_components = key in {'a_gap_fill_ligate','b_gap_fill_ligate','pcr'}
        reagent=g.input(key+'_stock',stock,K.MOLECULAR if supplied_components else K.CATALYST if enzymatic else (K.MOLECULAR if kind==S.BINDING else K.CONDITION))
        inputs=[(x,R.TEMPLATE)]
        if supplied_components:
            enzyme=g.product(key+'_enzyme',('Q5 HotStart HiFi polymerase component of the actual master mix' if key=='pcr' else 'Source Taq DNA ligase40U and T4 DNA polymerase3U components of the same gap-fill supplies'),K.CATALYST)
            nts=g.product(key+'_nts',('Co-present Q5 master mix dNTP substrates; source concentration unassigned' if key=='pcr' else 'Actual source100uM dNTP substrates of the same gap-fill supplies; per-base convention unassigned'),K.MOLECULAR)
            buffer=g.product(key+'_buffer','Co-present source reaction buffer component',K.CONDITION)
            components=[enzyme,nts,buffer]
            inputs.extend(((enzyme,R.CATALYST),(nts,R.STRUCTURAL),(buffer,R.CONDITION)))
            if key!='pcr':
                nad=g.product(key+'_nad','Actual source 0.5mM NAD+ cofactor from the same gap-fill supplies',K.MOLECULAR)
                components.append(nad);inputs.append((nad,R.STRUCTURAL))
            g.react(key+'_components','Co-present source supply functions; packaging and stock count unassigned, no physical separation or extra doses',S.GENERIC,((reagent,R.TEMPLATE),),components,changes_molecule=False)
        else:
            inputs.append((reagent,R.CATALYST if enzymatic else (R.STRUCTURAL if kind==S.BINDING else R.CONDITION)))
        for i,name in enumerate(oligos):inputs.append((g.input(key+'_oligo_'+str(i),name,K.OLIGO),R.STRUCTURAL))
        kw={}
        if kind==S.MAGIC:kw['reason']='Source chemistry positive; precise combined substrate/product relation unassembled.'
        if kind==S.GENERIC:kw['changes_molecule']=False
        if kind==S.PCR:kw['geometry']=PCRGeometry.unavailable(unknown_state(),reason='Actual main ILM_p5_PCR/ILM_p7_PCR and eight-base sample indexes; whole sequences/amplicon unassembled in this pass.')
        out=g.product(key+'_product',product,K.MOLECULAR)
        g.react(key,label,kind,inputs,(out,),**kw);x=out
    def wash(key,label):
        nonlocal x
        c=g.input(key+'_buffer',label,K.CONDITION);out=g.product(key+'_kept','Physically retained cell-associated preparation, recovery/carryover unquantified',K.MOLECULAR);w=g.product(key+'_waste','Actual removed supernatants/washes, content not perfectly assigned',K.WASTE)
        g.react(key,label,S.CLEANUP,((x,R.TEMPLATE),(c,R.CONDITION)),(out,w));g.discard(w,'Physical source workup; not complete molecular selectivity.');x=out
    def aliquot(key,label,remainder):
        nonlocal x
        out=g.product(key+'_used',label,K.MOLECULAR);rest=g.product(key+'_rest',remainder,K.MOLECULAR)
        g.react(key,label,S.SPLIT,((x,R.TEMPLATE),),(out,rest));g.mark_terminal(Terminal(rest.id,Disposition.UNRESOLVED,'Source stock remaining outside selected reaction.'));x=out
    def qc(key,label):
        nonlocal x
        out=g.product(key+'_kept','Production stock after independent analytical allocation',K.MOLECULAR);q=g.product(key+'_qc',label,K.QC)
        g.react(key,label,S.SPLIT,((x,R.TEMPLATE),),(out,q));g.mark_terminal(Terminal(q.id,Disposition.QC,'Actual source analytical allocation not returned to production.'));x=out
    op('fix','Main1percentPFA/PBS15minRT',S.CROSS_LINKING,'SourcePFA','Sourcefixedcells; exactcrosslinksunassembled')
    wash('fix_wash','MainonePBSwash')
    op('block','Mainblocking/quenching15min4C',S.QUENCHING,'1percentFBS/.1percentBSA/1mgmlssDNA/50ugmlhumanIgG/125mMglycine/4mMEDTA/.04percentProClin300/PBS','Sourceblockedfixedcell preparation')
    wash('block_wash','MainPBSwash')
    aocs=g.input('aocs','Source preconjugated validated antibody-oligo cocktail,80AOCs each5ugml; complete chemical orderings and DBCO-azide coupling unassembled',K.MOLECULAR)
    stain=g.input('stain_buffer','Main.2percentBSA/2mMEDTA/PBS;50ul30min4C',K.CONDITION);stained=g.product('stained','Sourcecell-boundAOCs plusunboundAOCs; no perfect binding/specificity assertion',K.MOLECULAR)
    g.react('stain','Main30min AOC binding; SI45min alternative not substituted',S.BINDING,((x,R.TEMPLATE),(aocs,R.STRUCTURAL),(stain,R.CONDITION)),(stained,));x=stained
    wash('aoc_wash','Mainthreewashes50mMNaCl/1mMEDTA/20mMTrispH8')
    op('secondary','Main20ugmlsecondaryantibody/.2percentBSA/2mMEDTA/PBS30min37C',S.BINDING,'Source secondary antibody','Stabilizedcell-boundAOCs; upstream conjugate geometry unassembled')
    wash('secondary_wash','Main two wash-buffer washes')
    # Pixel preparation is an independent reagent branch; intact cells are not circularized.
    pixels=[]
    for letter,ntps in [('a','Source0.75mM dAUGC or dNTPs; exact per-pixel formulation unassigned'),('b','Source0.75mM dAUGC or dNTPs; exact per-pixel formulation unassigned')]:
        template=g.input(letter+'_template','Source300nM template oligo; exact sequence and UPI variable window unassembled',K.OLIGO)
        pad=g.input(letter+'_padlock','Source200nM padlock; exact source sequence/end modification unassembled',K.OLIGO)
        lig=g.input(letter+'_ligase','Source400UT4ligase/1mMATP/33mMTrisAc7.9/10mMMgAc/66mMKAc/.1percentTween20/1mMDTT50ul37C30min→75C10min',K.CATALYST)
        circle=g.product(letter+'_circle_mix','Source circularized templates plus unclosed/linear substrates; exact closure geometry unassembled',K.MOLECULAR)
        g.react(letter+'_circle','Source independent pixel-template circularization',S.CIRCULARIZATION,((template,R.TEMPLATE),(pad,R.STRUCTURAL),(lig,R.CATALYST)),(circle,),geometry=CircularizationGeometry.unavailable('Source circle operation; selected oligo source tables not read within cap.'))
        ex=g.input(letter+'_ex','SourceExoI10U/ExoIII20U37C30min→85C20min',K.CATALYST);dig=g.product(letter+'_ex_product','Source circle-enriched mixture including digested linear DNA and reaction companions',K.MOLECULAR)
        g.react(letter+'_exonuclease','Sourcetemplate exonuclease treatment, not physical cleanup',S.DIGESTION,((circle,R.TEMPLATE),(ex,R.CATALYST)),(dig,))
        used_circle=g.product(letter+'_circle_used','Sourcecircle-preparationaliquotgivingfinal5nMtemplatein75ulRCA; actualvolume/yieldunassigned',K.MOLECULAR)
        circle_rest=g.product(letter+'_circle_rest','Untransferred50ulcircle-preparationstock; exactremainingvolume/fateunassigned',K.MOLECULAR)
        g.react(letter+'_circle_aliquot','SourcecirclepreparationstockallocatedtoRCA; noassumedquantitativeyield',S.SPLIT,((dig,R.TEMPLATE),),(used_circle,circle_rest));g.mark_terminal(Terminal(circle_rest.id,Disposition.UNRESOLVED,'Circlepreparationnotimplicitlyconsumedinfull.'))
        rcac=g.input(letter+'_rca_mix','Source aggregate RCA supplies: phi29 7.5U/'+ntps+'/75ul; packaging and stock count unassigned',K.MOLECULAR)
        enzyme=g.product(letter+'_rca_enzyme','Phi29 component of the same source RCA supplies',K.CATALYST)
        nts=g.product(letter+'_rca_nts',ntps+'; nucleotide substrates of the same source RCA supplies',K.MOLECULAR)
        buffer=g.product(letter+'_rca_buffer','Same source RCA buffer as template preparation; no fresh inferred formulation',K.CONDITION)
        g.react(letter+'_rca_components','Co-present RCA supply functions; no extra phi29/nucleotide doses or physical separation',S.GENERIC,((rcac,R.TEMPLATE),),(enzyme,nts,buffer),changes_molecule=False)
        rca=g.product(letter+'_rca_product','Source RCA DNA-pixel preparation; repeats/neighborhood UPIs known conceptually, exact copies/ends unassembled',K.MOLECULAR)
        g.react(letter+'_rca','SourceRCA pixel production37C4min→65C10min',S.EXTENSION,((used_circle,R.TEMPLATE),(enzyme,R.CATALYST),(nts,R.STRUCTURAL),(buffer,R.CONDITION)),(rca,))
        sap=g.input(letter+'_sap','SourcerSAP1U37C20min→65C5min',K.CATALYST);ready=g.product(letter+'_ready','Sourcepixelreagentmixture afterfree-nucleotide inactivation; not purified DNA alone',K.MOLECULAR)
        g.react(letter+'_rsap','SourcerSAPfree-nucleotide inactivation',S.DEPHOSPHORYLATION,((rca,R.TEMPLATE),(sap,R.CATALYST)),(ready,))
        used_pixel=g.product(letter+'_pixel_used','Sourcepixelstockaliquotgiving1nMfinalin50ulhybridization; exacttransfer/yieldunassigned',K.MOLECULAR)
        pixel_rest=g.product(letter+'_pixel_rest','Untransferred75ulRCA-derivedpixelstockafterrSAP; exactremainingvolume/fateunassigned',K.MOLECULAR)
        g.react(letter+'_pixel_aliquot','Sourcepixelstockallocatedtohybridization; noentire-stockconsumption',S.SPLIT,((ready,R.TEMPLATE),),(used_pixel,pixel_rest));g.mark_terminal(Terminal(pixel_rest.id,Disposition.UNRESOLVED,'Remainingpixelstockoutsidechosenassay.'));pixels.append(used_pixel)
    aliquot('pixelation_input','Main15000stainedcells selectedfor50ulpixelation','Other stained cells/stock outside selected pixelation; exact count unassigned')
    for letter,pixel in zip(('a','b'),pixels):
        hyb=g.input(letter+'_hyb_buffer','Main50ul:1nMpixels/2uM'+letter.upper()+'-GFP/1ugulshearedsalmonDNA/300mMNaCl/15mMMgCl2/20mMTris8/.05percentTween20;55C15min',K.CONDITION)
        primer=g.input(letter+'_gfp','Source'+letter.upper()+'-GFP oligo; exact modified ordering unassembled',K.OLIGO)
        joined=g.product(letter+'_hyb_product','Sourcepixel/AOC neighborhoodhybridizationmixture; precise binding geometry unassembled',K.MOLECULAR)
        g.react(letter+'_hybridize','Sourcepixel'+letter.upper()+' hybridization',S.HYBRIDIZATION,((x,R.TEMPLATE),(pixel,R.STRUCTURAL),(primer,R.STRUCTURAL),(hyb,R.CONDITION)),(joined,));x=joined
        wash(letter+'_hyb_wash','Main two wash-buffer washes')
        # Simultaneous fill plus ligation remains opaque rather than falsely claiming one exact sealed junction.
        op(letter+'_gap_fill_ligate','Main simultaneous gap-fill ligation37C20min',S.MAGIC,'SourceTaqDNAligase40U/T4DNApolymerase3U/100uMdNTP/.5mMNAD+/1xrCutSmart50ul','Source'+letter.upper()+'-UPI-bearingAOCreporters plus incomplete/unmodifiedclasses; copy/ligationgeometryunassembled')
        wash(letter+'_gap_wash','Main one wash-buffer wash')
        if letter=='a':
            op('user','MainA-pixel USER1U/50ulwashbuffer37C30min; SI15min not substituted',S.DIGESTION,'SourceUSER','SourceA-pixel cleavage products and A-UPI-bearingAOCreporters, not physically purified')
            wash('user_wash','Mainonewash afterUSER; actualphysicalremoval followscleavage')
    qc('cell_count','Sourcehemocytometer analytical cell allocation')
    aliquot('ex_input','Mainfewerthan1000cells into15ulexonucleasereaction','Remaining pixelated cells outside selected PCR; not counted as sequenced')
    op('lambda','MainlambdaExo10U/100uMdNTP/1mMNAD+/1xrCutSmart37C30min→75C10min',S.DIGESTION,'Sourcelambdaexonuclease','Source exonuclease-treated reporter mixture including cleavage products; selection not perfectphysicalpurification')
    op('pcr','Mainwhole15ulExoproduct into40ulQ5HotStartHiFiPCR/.4uMprimers',S.PCR,'SourceQ5HotStartHiFi master mix; actualPCRcycles unassigned inmainselectedscope','Sourceindexedampliconmixture with AOCUMI/proteinID/twoUPI meanings, exact wholelibrary unassembled',('ILM_p5_PCR with source8ntsampleindex; wholemodifiedordering unassembled','ILM_p7_PCR with source8ntsampleindex; wholemodifiedordering unassembled'))
    aliquot('cleanup_input','SI1Step4 uses30ul of40ulPCR;source-qualified cleanupallocation','SI1remaining10ulPCR may be stored−20C; stocknotimplicitlydiscarded')
    wash('ampure1','SI1Step4first1.5xAMPure/2x180ul70percentEtOH/elute30ulTE')
    wash('ampure2','SI1Step4second1.5xAMPure/2xEtOH/elute20ulTE; mainstatesdoubleAMPure')
    qc('qubit','MainQubitHS analytical aliquot')
    phix=g.input('phix','SourcePhiXDNAsequencinglibrary15percentphysicalspike; actualmass/volumeunassigned',K.MOLECULAR)
    diluent=g.input('loading_buffer','Main0.65nMloading; sourcebufferandexacttransferunassigned',K.CONDITION)
    mixture=g.product('loading_mixture','SourcephysicalMPX+PhiXsequencingmixture; PhiXDNAretainedinloadedmaterial',K.MOLECULAR)
    g.react('dilute','Mainphysical15percentPhiXspikeandloadingdilution',S.GENERIC,((x,R.TEMPLATE),(phix,R.SPIKE_IN),(diluent,R.CONDITION)),(mixture,),changes_molecule=False);x=mixture
    lib=g.product('library','SourceMPXNGSlibrary; complete molecular and read-primer geometry unresolved',K.LIBRARY)
    g.react('interface','Source library interface',S.GENERIC,((x,R.TEMPLATE),),(lib,),changes_molecule=False);g.sequence(lib,SequencingPlan(limitation=READOUT))
    g.add_source_item(SourceItem('selected',AUTHORITY,'Sourcesreadwithin3originalcap; no complete source/absence review.'))
    g.add_source_item(SourceItem('differences',SOURCE_DIFFERENCES,'Publishedmain and SI differing instructions preserved explicitly.'))
    g.cover_source_item(SourceCoverage('selected','modelled','library'))
    g.cover_source_item(SourceCoverage('differences','modelled','stain'))
    g.declare_reviewed_scope(REVIEWED_INPUTS, AUTHORITY+' Selected actual molecular supplies and prepared-stock boundaries; unknown whole oligos remain unassembled.')
    g.declare_reviewed_scope(REVIEWED_OPERATIONS, AUTHORITY+' Actual selected physical preparation, allocations, QC and sequencing contribution; no assertion of complete occupied products or yield.')
    return g

def selftest():
    from checks import Check
    c=Check();g=protocol_graph()
    c('A-pixelUSERbeforephysicalwash',g.steps['user_wash'].inputs[0].material,'user_product')
    c('Bpixel followsactualAcleavageandwash',g.steps['b_hybridize'].inputs[0].material,'user_wash_kept')
    c('Exonuclease followsactualchosenaliquot',g.steps['lambda'].inputs[0].material,'ex_input_used')
    c('PCRwholeexoproduct',g.steps['pcr'].inputs[0].material,'lambda_product')
    c('PCRcleanupstockretained',g.terminals['cleanup_input_rest'].disposition,Disposition.UNRESOLVED)
    c('Bothpixelbranches independentofcells',all(g.steps[k+'_circle'].inputs[0].material==k+'_template' for k in ('a','b')),True)
    c('Fouractualpixelstockremainders',sorted(k for k in g.terminals if k.endswith('_circle_rest')or k.endswith('_pixel_rest')),['a_circle_rest','a_pixel_rest','b_circle_rest','b_pixel_rest'])
    c('PhiXisphysicalDNAspike',[(i.material,i.role)for i in g.steps['dilute'].inputs if i.material=='phix'],[('phix',R.SPIKE_IN)])
    c('RCAconsumesselectedcirclestock',g.steps['a_rca'].inputs[0].material,'a_circle_used')
    c('Hybridizationconsumesselectedpixelstock',g.steps['a_hybridize'].inputs[1].material,'a_pixel_used')
    c('NoBpixelUSERinvented',any(k=='b_user' for k in g.steps),False)
    for letter in ('a','b'):
        c('RCA consumes nucleotide component of same actual supplies '+letter,
          (g.steps[letter+'_rca_components'].inputs[0].material,
           [(i.material,i.role) for i in g.steps[letter+'_rca'].inputs if i.material==letter+'_rca_nts']),
          (letter+'_rca_mix',[(letter+'_rca_nts',R.STRUCTURAL)]))
        key=letter+'_gap_fill_ligate'
        c('Gap-fill consumes same-source NAD and nucleotide components '+letter,
          {i.material for i in g.steps[key].inputs if i.role==R.STRUCTURAL},
          {key+'_nts',key+'_nad'})
    c('Q5 master polymerase and NT functions derive from one source addition',
      (g.steps['pcr_components'].inputs[0].material,
       [(i.material,i.role)for i in g.steps['pcr'].inputs if i.material in {'pcr_enzyme','pcr_nts'}]),
      ('pcr_stock',[('pcr_enzyme',R.CATALYST),('pcr_nts',R.STRUCTURAL)]))
    c('Selected reviewed preparation passes normal publication gate',g.publication_problems(),[])
    g.input('unreviewed_addition','Regression-only unreviewed stock',K.CONDITION)
    c('New stock cannot auto-publish',any('unreviewed_addition' in p and 'source-inventory' in p for p in g.publication_problems()),True)
    c.report()


PAGE_FILE = 'molecular_pixelation.html'
def render_page():
    from page import head,info
    from protocol_views import assemble_page,WorkflowSection,PlanLimitationSection
    g=protocol_graph().require_publishable()
    return assemble_page(g,[head('Molecular pixelation'),'<div class="wrap" data-protocol-complete="true"><h1>Molecular pixelation</h1>',
      '<p>Antibody–oligo conjugates label cell-surface proteins. Two successive DNA-pixel reactions copy neighborhood identifiers onto the reporters, producing sequencing reads used to reconstruct protein organization.</p>',
      info('Source: <a href="https://doi.org/10.1038/s41592-024-02268-9">Molecular pixelation: spatial proteomics of single cells by sequencing</a>, published main and Supplementary Protocol. This page selects the main-article fixed-cell route. Antibody conjugates are prepared stocks; CD3 capping, rituximab and uropod perturbations are separate experiments.'),
      '<h2>Published instruction differences</h2>',info(SOURCE_DIFFERENCES),
      '<h2>Molecular components</h2><p>The antibody cocktail contains 80 conjugates. Each independent A/B pixel factory circularizes source template/padlock oligos, treats linear substrates with exonucleases, allocates circle preparation to phi29 RCA and treats remaining free nucleotides with rSAP. The assay then uses selected portions of these pixel stocks; it does not circularize the cells or consume every prepared stock in full.</p>',
      '<p>Protein identifiers, antibody-oligo UMIs, pixel UPIs and eight-base sample indexes have distinct meanings. UPIs describe molecular neighborhoods, not predetermined cell barcodes. Exact source oligo bodies, variable-window coordinates, occupied complexes, sealed junctions and primer landing are unassembled here. No illustrative whole construct is substituted for them.</p>',
      '<h2>Preparation and reads</h2><p>Fix, block and stain cells, then stabilize bound conjugates with secondary antibody. Hybridize A pixels, perform gap filling and ligation, cleave A pixels with USER and wash. Repeat hybridization and gap filling with B pixels. Count cells and allocate fewer than 1,000 to the lambda-exonuclease/PCR branch. Purify the selected PCR portion twice and quantify before physical PhiX addition.</p>',
      info('The source sequence run is NextSeq 2000 paired-end 44 + 78 cycles with 15% PhiX. Full read/index masks remain unresolved. SI Step 4 supplies the separately qualified 30-of-40-µl PCR cleanup allocation. Pixelator graph reconstruction and filtering are computational, not additional physical purification.'),PlanLimitationSection('library'),
      '<details><summary>Full preparation, pixel factories and chemistry warnings</summary>',WorkflowSection('Source-qualified main-article preparation'),'</details>',
      '<p><a href="01_mpx.html">Source notes and model limits</a></p></div>'])
if __name__=='__main__':selftest()
