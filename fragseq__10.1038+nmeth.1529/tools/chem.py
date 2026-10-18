"""FragSeq owns its three RNA-end branches and historical SOLiD boundary."""
from adapter_ligation import AdapterLigationGeometry
from chemdraw import Construct, Scene, feature
from pcr import PCRGeometry
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind, Disposition, Terminal
from protocol_model import molecule, role
TITLE='FragSeq — nuclease accessibility profiling of RNA'
NOTES='01_fragseq.html'
SOURCE='Underwood et al., Nature Methods (2010), <a href="https://doi.org/10.1038/nmeth.1529">doi:10.1038/nmeth.1529</a>.'
SUMMARY='Nuclease P1 generates5′monophosphate/3′OH RNA ends; parallel untreated and PNK controls distinguish accessible cuts from endogenous breaks.'
CAVEAT=('Purified refolded nuclear RNA reports in-vitro nuclease accessibility. Historical SOLiD adapter/RT/PCR species and exact run cycles are not reconstructed. '
        'The published software filter FASTA contains computational masking patterns, including padded Ns; those entries are not physical ordered adapters.')
SAMPLE=feature('fragseq_sample_index','sample_index','unknown')
FINAL_LIBRARY=Construct([role('SOLiD P1-derived library end; physical geometry unassigned'),role('RNA-fragment-derived cDNA'),
    role('SOLiD P2/PCR barcode-derived end; selected barcode geometry unassigned',molecular_feature=SAMPLE)],name='FragSeq source-limited SOLiD endpoint')
SEQ_PRIMERS=()
READ_LENGTHS={'Read 1':50}
FINAL_CAPTION='Mapped RNA fragment termini define cutting scores; size selection and end-compatible ligation select which termini are represented.'
SEQUENCING_INTRO='Published study uses ABI SOLiD3 and individually PCR-barcoded libraries.'
SEQUENCING_ENDING='Only reads reaching the3′adapter establish both fragment ends; other reads establish the sampled5′end.'
SEQUENCING_UNAVAILABLE='Historical SOLiD sequencing-primer species, color-space ligation-cycle recipe and selected barcodes remain unassigned; the published biological read length is50bp.'
def oligos(): return ()

def rna_ends(label, *, ready, selected=True):
    sc=Scene();sc.strand('RNA',[role('20–100nt RNA fragment; exact sequence variable' if selected else 'RNA cleavage fragment; size not yet selected')],label=label,
        mod5='P' if ready else '',mod3='OH' if ready else '')
    if not ready: sc.note('RNA','End classes remain mixed; ligation requires5′P/3′OH and does not repair incompatible ends.')
    return sc

def protocol_graph():
    g=ProtocolGraph(TITLE,authority=SOURCE)
    rna=g.input('input','purified refolded mouse nuclear RNA',MaterialKind.BIOLOGICAL,state=molecule('nuclear RNA',(role('folded nuclear RNA; heterogeneous ends'),)),primary=True)
    treated=g.product('treated_branch','P1 digestion aliquot',MaterialKind.MOLECULAR,state=rna.state)
    control=g.product('control_branch','undigested aliquot',MaterialKind.MOLECULAR,state=rna.state)
    pnk=g.product('pnk_branch','PNK control aliquot',MaterialKind.MOLECULAR,state=rna.state)
    g.react('split','Prepare separate P1,undigested and PNK control samples',StepKind.SPLIT,((rna,InputRole.TEMPLATE),),(treated,control,pnk),changes_molecule=False)
    enzyme=g.input('p1_enzyme','nuclease P1: limited cleavage of accessible RNA',MaterialKind.CATALYST)
    cut=g.product('p1_fragments','P1 fragments with5′monophosphate/3′OH ends',MaterialKind.MOLECULAR,state=rna_ends('representative P1-generated compatible end class',ready=True,selected=False))
    g.react('p1_digest','Partially digest accessible RNA with P1',StepKind.DIGESTION,((treated,InputRole.TEMPLATE),(enzyme,InputRole.CATALYST)),(cut,))
    reviewed=['input','split','p1_enzyme','p1_digest']
    library_species=[]
    for prefix,template,label in (('treated',cut,'P1-treated'),('control',control,'undigested'),('pnk',pnk,'PNK control')):
        if prefix=='pnk':
            kinase=g.input('pnk_enzyme','polynucleotide kinase and ATP',MaterialKind.CATALYST)
            repaired=g.product('pnk_repaired','PNK-processed endogenous-break control RNA',MaterialKind.MOLECULAR,state=rna_ends('PNK-compatible end class',ready=True,selected=False))
            g.react('pnk_control','Phosphorylate5′ends and remove3′cyclic phosphate in control',StepKind.END_REPAIR,((template,InputRole.TEMPLATE),(kinase,InputRole.CATALYST)),(repaired,))
            template=repaired;reviewed.extend(('pnk_enzyme','pnk_control'))
        gel=g.input(prefix+'_gel','FlashPAGE isolate20–100nt RNA',MaterialKind.SELECTION)
        selected=g.product(prefix+'_selected',label+'20–100nt RNA',MaterialKind.MOLECULAR,state=rna_ends(label+' size-selected RNA',ready=prefix in ('treated','pnk')))
        offsize=g.product(prefix+'_offsize','RNA outside excised gel window',MaterialKind.WASTE)
        g.react(prefix+'_size_select','Gel isolate20–100nt '+label+' RNA',StepKind.SELECTION,((template,InputRole.TEMPLATE),(gel,InputRole.SELECTION)),(selected,offsize),changes_molecule=False)
        g.discard(offsize,'not included in the physically excised20–100nt gel fraction')
        reviewed.extend((prefix+'_gel',prefix+'_size_select'))
        adapter=g.input(prefix+'_adapters','historical SOLiD small-RNA adapters; physical species unassigned',MaterialKind.ADAPTER)
        ligase=g.input(prefix+'_ligase','small-RNA library ligation reagents; identities unassigned',MaterialKind.CATALYST)
        ligated=g.product(prefix+'_ligated',label+' adapter-ligated RNA',MaterialKind.MOLECULAR,state=molecule('adapter-ligated RNA',(role('5′RNA adapter; native strand/end unassigned'),role('compatible5′P/3′OH RNA fragment'),role('3′RNA adapter; native strand/end unassigned'))))
        not_ligated=g.product(prefix+'_unligated',label+' molecules not entering adapter-positive library',MaterialKind.MOLECULAR)
        g.react(prefix+'_ligate','Ligate adapters to end-compatible '+label+' RNA',StepKind.LIGATION,((selected,InputRole.TEMPLATE),(adapter,InputRole.STRUCTURAL),(ligase,InputRole.CATALYST)),(ligated,not_ligated),geometry=AdapterLigationGeometry.unavailable('Article requires5′P and3′OH on RNA. Physical adapter strands/modifications and native ligating geometry remain undeclared; software filter motifs are not ordered reagents.'))
        g.mark_terminal(Terminal(not_ligated.id,Disposition.UNRESOLVED,'End incompatibility/non-ligation prevents library entry; it does not establish physical removal. Exact molecular recovery is unassigned.'))
        rtprimer=g.input(prefix+'_rt_primer','SOLiD small-RNA RT primer; physical species unassigned',MaterialKind.OLIGO)
        rt=g.input(prefix+'_rt_enzyme','source-delegated reverse-transcription reagents',MaterialKind.CATALYST)
        cdna=g.product(prefix+'_cdna',label+' first-strand cDNA',MaterialKind.MOLECULAR,state=molecule('actual first-strand cDNA',(role('3′RNA-adapter-derived RT-primer end'),role('RNA-fragment cDNA'),role('copied5′RNA-adapter end'))))
        residual_rna=g.product(prefix+'_rt_template_rna',label+' RNA templates after RT; physical fate unresolved',MaterialKind.MOLECULAR,state=ligated.state)
        g.react(prefix+'_rt','Reverse-transcribe '+label+' adapter-ligated RNA',StepKind.REVERSE_TRANSCRIPTION,((ligated,InputRole.TEMPLATE),(rtprimer,InputRole.STRUCTURAL),(rt,InputRole.CATALYST)),(cdna,residual_rna))
        g.mark_terminal(Terminal(residual_rna.id,Disposition.UNRESOLVED,'RT produces cDNA but does not physically remove template RNA. This graph follows the cDNA target; RNA recovery or degradation through later preparation is unassigned.'))
        f=g.input(prefix+'_forward','SOLiD library forward PCR primer; species unassigned',MaterialKind.OLIGO)
        r=g.input(prefix+'_reverse','sample-barcoding SOLiD PCR primer; selected member unassigned',MaterialKind.OLIGO)
        pol=g.input(prefix+'_polymerase','source-delegated SOLiD PCR reagents',MaterialKind.CATALYST)
        library=g.product(prefix+'_pcr_product',label+' barcoded PCR product',MaterialKind.MOLECULAR,state=FINAL_LIBRARY)
        g.react(prefix+'_pcr','PCR-barcode '+label+' library',StepKind.PCR,((cdna,InputRole.TEMPLATE),(f,InputRole.STRUCTURAL),(r,InputRole.STRUCTURAL),(pol,InputRole.CATALYST)),(library,),geometry=PCRGeometry.unavailable(cdna.state,reason='Tracked first-strand cDNA is the PCR target; remaining RNA material is separately unresolved. Historical primers, first-copy geometry and completed SOLiD endpoint remain undeclared.'))
        gel_after=g.input(prefix+'_pcr_gel','6%nativeTBEgel;25–100nt insert window',MaterialKind.SELECTION)
        gel_library=g.product(prefix+'_library',label+' barcoded gel-selected library',MaterialKind.MOLECULAR,state=library.state)
        gel_waste=g.product(prefix+'_pcr_gel_waste','material outside the excised nativegel window',MaterialKind.WASTE)
        g.react(prefix+'_pcr_gel_select','Select barcoded products on6%nativeTBEgel',StepKind.SELECTION,((library,InputRole.TEMPLATE),(gel_after,InputRole.SELECTION)),(gel_library,gel_waste),changes_molecule=False)
        g.discard(gel_waste,'published SupplementaryMethods selects25–100nt inserts afterPCR; excluded gel material is not sequenced')
        library_species.append((prefix,gel_library))
        reviewed.extend((prefix+'_pcr_gel',prefix+'_pcr_gel_select'))
        reviewed.extend((prefix+'_adapters',prefix+'_ligase',prefix+'_ligate',prefix+'_rt_primer',prefix+'_rt_enzyme',prefix+'_rt',prefix+'_forward',prefix+'_reverse',prefix+'_polymerase',prefix+'_pcr',prefix+'_library'))
    pooling=g.input('equimolar_pooling','equimolar pooling of separately barcoded libraries',MaterialKind.CONDITION)
    pooled_scene=Scene()
    for prefix,library in library_species:
        pooled_scene.strand(prefix,[role(prefix+' separately barcoded library species; native outer geometry unassigned')],label=prefix+' species in one physical pool')
    pool=g.product('pooled_libraries','equimolar pool of independent barcoded libraries',MaterialKind.MOLECULAR,state=pooled_scene)
    g.react('pool_libraries','Pool individually barcoded gel-selected libraries equimolarly',StepKind.POOL,tuple((library,InputRole.TEMPLATE) for _,library in library_species)+((pooling,InputRole.CONDITION),),(pool,),changes_molecule=False)
    emulsion=g.input('solid_emulsion_pcr','source SOLiD emulsionPCR; native bead/primer geometry unassigned',MaterialKind.CONDITION)
    sequenced=g.product('solid_library','SOLiD emulsionPCR-prepared pooled library species',MaterialKind.LIBRARY,state=pool.state)
    g.react('solid_emulsion','Prepare the source SOLiD emulsionPCR pool',StepKind.MAGIC,((pool,InputRole.TEMPLATE),(emulsion,InputRole.CONDITION)),(sequenced,),reason='Published SupplementaryMethods names emulsionPCR before SOLiD3 sequencing; native physical bead/primer/product geometry is not reconstructed here.')
    g.sequence(sequenced,SequencingPlan(acquired_roles=('Read 1',),read_lengths={'Read 1':50},limitation='Published SupplementaryMethods reports ABI SOLiD3 50bp reads. This states biological read length, not a reconstructed color-space ligation-cycle recipe or index acquisition. '+SEQUENCING_ENDING))
    reviewed.extend(('equimolar_pooling','pool_libraries','solid_emulsion_pcr','solid_emulsion','solid_library'))
    g.declare_reviewed_scope(reviewed,SOURCE)
    return g.require_publishable()
