"""ARM-seq pretreatment graph; historical kit library remains an explicit boundary."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'lib'))
from chemdraw import Construct, Segment
from protocol_graph import (ProtocolGraph, MaterialKind, InputRole, StepKind,
                            Terminal, Disposition, SequencingPlan)
from protocol_views import assemble_page, WorkflowSection
from page import head, info
PAGE_FILE="arm-seq.html"
SOURCE = 'Cozen et al. 2015 doi:10.1038/nmeth.3508, Methods Treatment of RNA with AlkB; library preparation'
BOUNDARY = ('The source names NEBNext Small RNA Library Prep Kit for Illumina; applicable '
            'historical kit identity/revision and oligos have not been reviewed. '
            'MiSeq paired-end output is reported, but no exact library or run lengths are derived.')


def rna_state(name):
    # A variable population, not a fabricated example tRNA or residue-specific modification map.
    return Construct([Segment(name, '[variable RNA population]', placeholder=True)], name=name)


def protocol_graph():
    g=ProtocolGraph('ARM-seq — matched treatment and control library preparations',authority=SOURCE+'; associated author manuscript, whole VOR equivalence not established',representation='published_instance')
    source=g.input('source_rna','Actual selected bulk RNA sample stock; transcript/modified-site identities and total stock amount unassigned',MaterialKind.MOLECULAR,primary=True)
    treated=g.product('treated_input','Actual 50 µg RNA allocation for AlkB treatment',MaterialKind.MOLECULAR)
    control=g.product('control_input','Actual separate 50 µg RNA allocation for matched storage-buffer control',MaterialKind.MOLECULAR)
    remaining=g.product('source_stock','Unallocated source RNA stock class; amount/presence unassigned, not claimed discarded',MaterialKind.MOLECULAR)
    g.react('source_allocate','Allocate actual matched 50 µg treatment/control inputs before reactions',StepKind.SPLIT,((source,InputRole.TEMPLATE),),(treated,control,remaining))
    g.mark_terminal(Terminal(remaining.id,Disposition.UNRESOLVED,'Source stock outside selected matched preparations; amount/fate unassigned'))
    for prefix,template in (('treated_',treated),('control_',control)):
        def reagent(key,label,kind=MaterialKind.MOLECULAR,role=InputRole.STRUCTURAL):
            return(g.input(prefix+key,label,kind),role)
        def op(key,label,kind,extras=(),chemical_state=None):
            nonlocal template
            product=g.product(prefix+key+'_product','Actual '+label+' material; whole occupied molecules, recovery and residual modification fractions unassigned',MaterialKind.MOLECULAR,chemical_state=chemical_state)
            g.react(prefix+key,label,kind,((template,InputRole.TEMPLATE),)+tuple(extras),(product,));template=product
        def cleanup(key,label,extras=()):
            nonlocal template
            retained=g.product(prefix+key+'_kept','Actual retained '+label+' material; purity/recovery unassigned',MaterialKind.MOLECULAR)
            waste=g.product(prefix+key+'_waste','Actual removed '+label+' fractions',MaterialKind.WASTE)
            g.react(prefix+key,label,StepKind.CLEANUP,((template,InputRole.TEMPLATE),)+tuple(extras),(retained,waste));g.discard(waste,'Physical source workup reject, not a quantitative recovery assertion');template=retained
        common=(reagent('iron','Actual 75 µM ferrous ammonium sulfate; molecular cofactor contribution'),reagent('akg','Actual 1 mM alpha-ketoglutarate'),reagent('ascorbate','Actual 2 mM sodium ascorbate'),reagent('bsa','Actual 50 µg/ml BSA carrier'),reagent('reaction_buffer','Actual fresh HEPES-KOH pH8 reaction buffer; 200 µl reaction,37°C100min',MaterialKind.CONDITION,InputRole.CONDITION))
        if prefix=='treated_':
            enzyme=reagent('alkb','Actual 50 µg wild-type E.coli AlkB',MaterialKind.CATALYST,InputRole.CATALYST)
            op('demethylate','AlkB treatment37°C100min; residual modification fractions unassigned',StepKind.CHEMICAL_CONVERSION,common+(enzyme,),{'scope':'Source AlkB-sensitive methylation removal, not universal repair or site-by-site DNA mutation prediction'})
        else:
            op('incubate','Matched control incubation37°C100min with storage buffer instead of AlkB',StepKind.THERMAL_TREATMENT,common+(reagent('storage_buffer','Actual source AlkB storage-buffer contribution, no AlkB enzyme'),),{'scope':'Untreated source control; no demethylation activity assigned'})
        op('stop','Add actual 200 µl EDTA/ammonium-acetate stop mixture',StepKind.QUENCHING,(reagent('stop_reagent','Actual source 11 mM EDTA/200 mM ammonium acetate stop-buffer components'),reagent('stop_condition','Source stopping operation with the same actual buffer addition, not another dose',MaterialKind.CONDITION,InputRole.CONDITION)))
        cleanup('recover','Phenol extraction/ethanol precipitation and water resuspension',(reagent('phenol','Actual phenol extraction contribution'),reagent('ethanol','Actual ethanol precipitation contribution'),reagent('water','Actual pellet resuspension water')))
        cleanup('mirvana','MirVana selection of RNA below200nt',(reagent('mirvana_reagent','Actual historical MirVana miRNA Isolation Kit contribution; internal recipe unassembled',MaterialKind.SELECTION,InputRole.SELECTION),))
        cleanup('concentrate','RNA Clean and Concentrate-25 recovery; source describes RNA concentrated to25µg',(reagent('concentration_column','Actual source concentration/column materials',MaterialKind.SELECTION,InputRole.SELECTION),))
        selected=g.product(prefix+'dnase_input','Actual 10 µg allocation from source-described25µg RNA stock',MaterialKind.MOLECULAR)
        stock=g.product(prefix+'small_rna_stock','Nominal remaining15µg source-described small-RNA stock; actual stock origin/yield relationship not normalized',MaterialKind.MOLECULAR)
        g.react(prefix+'dnase_allocate','Actual 25→10µg allocation; nominal15µg untransferred stock',StepKind.SPLIT,((template,InputRole.TEMPLATE),),(selected,stock));g.mark_terminal(Terminal(stock.id,Disposition.UNRESOLVED,'Untransferred source RNA stock, no disposal or recovery yield inferred'));template=selected
        op('dnase','Actual DNase I treatment of10µg RNA allocation',StepKind.DIGESTION,(reagent('dnase_enzyme','Actual source DNase I',MaterialKind.CATALYST,InputRole.CATALYST),))
        cleanup('dnase_cleanup','Actual post-DNase column cleanup',(reagent('dnase_column','Actual source RNA column cleanup materials',MaterialKind.SELECTION,InputRole.SELECTION),))
        one=g.product(prefix+'kit_input','Actual 1 µg RNA allocation for historical NEBNext Small RNA Library Prep Kit for Illumina',MaterialKind.MOLECULAR)
        rest=g.product(prefix+'post_dnase_stock','Untransferred post-DNase RNA stock; nominal9µg before unassigned cleanup losses, actual recovered mass unassigned',MaterialKind.MOLECULAR)
        g.react(prefix+'kit_allocate','Allocate actual1µg kit input; preserve remaining RNA stock without inventing yield',StepKind.SPLIT,((template,InputRole.TEMPLATE),),(one,rest));g.mark_terminal(Terminal(rest.id,Disposition.UNRESOLVED,'Source stock remainder, no false exact recovery mass'));template=one
        kit=g.input(prefix+'kit_reagents','Actual named historical NEBNext Small RNA Library Prep Kit for Illumina contributions; edition, components and packaging unassembled',MaterialKind.MOLECULAR)
        dna=g.product(prefix+'kit_library','Actual source DNA-library preparation from1µgRNA using named historical kit; whole molecules, kit reactions and oligos unassembled',MaterialKind.LIBRARY,chemical_state={'source_positive_output':'Source reports prepared DNA libraries sized/QCed/sequenced; this is not RNA directly sequenced','unassembled':'Historical kit adapter/end/RT/PCR/primer/index and reagent internals not supplied from memory or a modern kit','validation_exclusion':'Source separate radiolabeled primer extension is not this library recipe'})
        g.react(prefix+'kit_prepare','Actual historical-kit RNA→DNA-library interface',StepKind.MAGIC,((template,InputRole.TEMPLATE),(kit,InputRole.STRUCTURAL)),(dna,),reason='Defining accepted manuscript positively reports 1µg kit input and resulting DNA libraries. Historical kit internals remain unreviewed; no ligase/adapter/RT/PCR default or vendor revision inferred.');template=dna
        cleanup('egel','Actual2%SizeSelect agarose E-Gel retention corresponding to18–120ntRNA inserts',(reagent('gel_material','Actual source E-Gel and50bp ladder contributions; not a claim that whole library DNA is18–120nt',MaterialKind.SELECTION,InputRole.SELECTION),))
        cleanup('library_column','Actual column cleanup and concentration before source analytical dilutions',(reagent('library_column_material','Actual source library column/concentration materials',MaterialKind.SELECTION,InputRole.SELECTION),))
        qa=g.product(prefix+'bioanalyzer_input','Actual source library dilution/analytical allocation for BioAnalyzer; amount unassigned',MaterialKind.QC)
        library=g.product(prefix+'qa_stock','Actual prepared library stock after analytical allocation; amount unassigned',MaterialKind.LIBRARY)
        g.react(prefix+'qa_allocate','Allocate actual BioAnalyzer dilution apart from production stock',StepKind.SPLIT,((template,InputRole.TEMPLATE),),(qa,library))
        qa_reagent=g.input(prefix+'qa_reagents','Actual Agilent High Sensitivity DNA kit analytical contributions; composition unassembled',MaterialKind.MOLECULAR)
        qa_result=g.product(prefix+'qa_readout','Actual BioAnalyzer analytical readout; not production library input',MaterialKind.QC)
        g.react(prefix+'qa','Actual source BioAnalyzer traces',StepKind.MAGIC,((qa,InputRole.TEMPLATE),(qa_reagent,InputRole.STRUCTURAL)),(qa_result,),reason='Source-positive analytical assay; detailed analytical reagent chemistry unassembled.');g.mark_terminal(Terminal(qa_result.id,Disposition.QC,'Source analytical allocation'))
        load=g.product(prefix+'loading','Actual selected MiSeq library loading contribution; amount unassigned',MaterialKind.LIBRARY)
        stored=g.product(prefix+'library_stock','Untransferred prepared library stock; amount/fate unassigned',MaterialKind.LIBRARY)
        g.react(prefix+'load_allocate','Actual source sequencing allocation',StepKind.SPLIT,((library,InputRole.TEMPLATE),),(load,stored));g.mark_terminal(Terminal(stored.id,Disposition.UNRESOLVED,'Prepared library remainder, not discarded'))
        g.sequence(load,SequencingPlan(acquired_roles=('Read 1','Read 2'),limitation='Source MiSeq paired-end. Acquired read lengths, index/read masks, barcode meaning/width and full primer landing remain unassembled. Historical kit DNA production is source-positive; inner chemistry is not reconstructed.'))
    g.declare_reviewed_scope(('source_rna','source_allocate','treated_input','control_input','source_stock'),SOURCE)
    for prefix in ('treated_','control_'):g.declare_reviewed_scope(tuple(prefix+k for k in SHARED_TARGETS),SOURCE)
    g.declare_reviewed_scope(('treated_alkb','treated_demethylate','treated_demethylate_product','control_storage_buffer','control_incubate','control_incubate_product'),SOURCE)
    return g.require_publishable()

# Curator-written target inventory for the inspected source route. No graph/node enumeration.
SHARED_TARGETS=tuple("""iron akg ascorbate bsa reaction_buffer stop_reagent stop_condition stop stop_product phenol ethanol water recover recover_kept recover_waste mirvana_reagent mirvana mirvana_kept mirvana_waste concentration_column concentrate concentrate_kept concentrate_waste dnase_allocate dnase_input small_rna_stock dnase_enzyme dnase dnase_product dnase_column dnase_cleanup dnase_cleanup_kept dnase_cleanup_waste kit_allocate kit_input post_dnase_stock kit_reagents kit_prepare kit_library gel_material egel egel_kept egel_waste library_column_material library_column library_column_kept library_column_waste qa_allocate bioanalyzer_input qa_stock qa_reagents qa qa_readout load_allocate loading library_stock""".split())

def render_page():
    g=protocol_graph()
    return assemble_page(g,(head('ARM-seq — modification-sensitive RNA recovery'),'<main class="wrap" data-protocol-complete="true"><h1>ARM-seq</h1><p>Wild-type AlkB treatment improves recovery of RNA carrying susceptible RT-blocking methylations. Matched treated and storage-buffer controls yield distinct small-RNA libraries; differential abundance does not establish a universal per-base mutation signal.</p><p><a href="https://doi.org/10.1038/nmeth.3508">Cozen et al. (2015)</a> · <a href="01_arm_seq.html">Research note</a></p>',info('Source-qualified associated author-manuscript preparation; whole version-of-record equivalence is not established. The named historical NEBNext kit is an explicit RNA-to-DNA-library handoff. Its edition, adapters, ends, RT/PCR recipe and primer/index identities remain unassembled, with no modern kit defaults. Exact occupied RNA and DNA molecules are not fabricated.'),'<h2>Preparation and source allocations</h2><p>Separate 50 µg treatment/control allocations receive the same fresh cofactors and molecular BSA; the control receives storage buffer instead of AlkB. Extraction, MirVana selection and concentration lead to the source-described 25 µg stock. Ten micrograms enter DNase treatment, and one microgram enters the named kit. Untransferred stocks remain separate; cleanup yield is not inferred from these labels.</p>','<h2>Library selection, quality assessment and reads</h2><p>The actual source DNA libraries undergo 2% SizeSelect E-Gel selection for products corresponding to 18–120 nt RNA inserts. BioAnalyzer dilutions are analytical material, separate from loading and library stock. The source reports MiSeq paired-end sequencing, with acquired cycle counts, index masks and exact read landing unassigned. The radiolabeled primer-extension validation assay is not reused as the library recipe.</p>','<details><summary>Full source preparation, kit boundaries and physical fates</summary>',WorkflowSection('Matched source preparations and separate DNA-library outputs'),'</details></main>'))

def selftest():
    from checks import Check
    c=Check();g=protocol_graph()
    c('Matched source preparations allocated before treatment',tuple(g.steps['source_allocate'].outputs),('treated_input','control_input','source_stock'))
    c('Control does not receive AlkB','control_alkb' in g.materials,False)
    c('BSA remains molecular material',g.materials['treated_bsa'].kind,MaterialKind.MOLECULAR)
    c('Actual source10µg allocation retains stock',g.terminals['treated_small_rna_stock'].disposition,Disposition.UNRESOLVED)
    c('Kit input derives from postDNase cleanup',g.steps['treated_kit_allocate'].inputs[0].material,'treated_dnase_cleanup_kept')
    c('Historical kit is explicit opaque source handoff',g.steps['treated_kit_prepare'].kind,StepKind.MAGIC)
    c('Actual kit output is DNA library',g.materials['treated_kit_library'].kind,MaterialKind.LIBRARY)
    c('Actual analytical fractions terminate as QC',tuple(g.terminals[k+'qa_readout'].disposition for k in ('treated_','control_')),(Disposition.QC,Disposition.QC))
    c('Two actual DNA library sequencing outputs',len(g.sequencing),2)
    c('Run lengths remain unassigned',tuple(bool(p.read_lengths) for p in g.sequencing.values()),(False,False))
    c('Normal publication gate',g.publication_problems(),[])
    from checks import Source,have
    from storage import sources_root
    source=Source(sources_root()/'arm-seq__10.1038+nmeth.3508'/'nmeth.3508_PMC4553111.pdf.txt','Accepted-manuscript sourceMethods351–402 only')
    if have(c,source,label='Source stock/control/kit assertions'):
        text=' '.join(source.path.read_text().split())
        c('Source control uses storage buffer','using AlkB storage buffer instead of AlkB' in text,True)
        c('Source historical kit and actual stock quantities','1 µg was used as input for NEBNext Small RNA Library Prep Kit' in text,True)
    c.report(exit_on_fail=True)

if __name__ == '__main__':
    selftest()
