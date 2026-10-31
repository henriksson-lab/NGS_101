"""Source-owned ChIP-BS-seq operational graph; unresolved molecular geometry is explicit."""
from pathlib import Path
import sys
HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE), str(HERE.parents[1] / 'lib')]
from chemdraw import Scene, workflow_from_protocol_graph, workflow_panel
from protocol_graph import ProtocolGraph, MaterialKind as K, InputRole as R, StepKind as S, SequencingPlan, SourceItem, Terminal, Disposition
from protocol_model import molecule, role
from adapter_ligation import AdapterLigationGeometry
from pcr import PCRGeometry
from page import head, info

PAGE_FILE = 'chip-bs-seq.html'
SOURCE = 'Brinkman et al.2012 DOI10.1101/gr.133728.111; published PMC3371717 Methods140–150 and Fig1A.'


def protocol_graph():
    g = ProtocolGraph('ChIP-BS-seq — captured DNA to bisulfite library', authority=SOURCE, representation='published_instance')
    def operation(key, label, kind, material, reagent_label, output_label=None, additional=(), molecular_reagent=False, **kw):
        enzyme = kind in {S.END_REPAIR, S.TAILING}
        reagent = g.input(key + '_reagent', reagent_label, K.CATALYST if enzyme else (K.MOLECULAR if molecular_reagent else K.CONDITION))
        product = g.product(key + '_material', output_label or label + '; retained mixture and carryover unquantified', K.MOLECULAR)
        g.react(key, label, kind, ((material, R.TEMPLATE), (reagent, R.CATALYST if enzyme else (R.STRUCTURAL if molecular_reagent else R.CONDITION)))+tuple(additional), (product,), **kw)
        return product
    def recover(key, material, label, *, molecular_reagent=False):
        reagent = g.input(key + '_procedure', label, K.MOLECULAR if molecular_reagent else K.SELECTION)
        retained = g.product(key, 'Physically recovered source preparation; exact recovery and purity unassigned', K.MOLECULAR)
        waste = g.product(key + '_waste', 'Source-separated fractions; exact molecular composition unassigned', K.WASTE)
        g.react(key + '_step', label, S.CLEANUP, ((material, R.TEMPLATE), (reagent, R.STRUCTURAL if molecular_reagent else R.SELECTION)), (retained, waste))
        g.discard(waste, 'Source physical separation; no ideal depletion or recovery claim.')
        return retained
    dna = g.input('chromatin', 'Source H3K27me3-associated chromatin; harvesting/ChIP predecessor unreviewed', K.BIOLOGICAL, primary=True)
    dna = operation('chip', 'Source anti-H3K27me3 ChIP and captured DNA isolation', S.MAGIC, dna, 'Actual anti-H3K27me3 molecular supply; cited Denissov2007 internals unreviewed', molecular_reagent=True, reason='Actual source identifies ChIP; exact upstream crosslinking, fragment and isolation geometry remains unreviewed. MethylCap is a separate precursor branch, not an antibody substitute.')
    selected = g.product('chip_100ng','Actual100ng captured ChIP DNA allocated for source library preparation',K.MOLECULAR)
    unused = g.product('chip_stock','Untransferred captured DNA if any; amount/lateruse unassigned',K.MOLECULAR)
    g.react('chip_allocation','Source100ng capturedDNA input allocation',S.SPLIT,((dna,R.TEMPLATE),),(selected,unused))
    g.mark_terminal(Terminal(unused.id,Disposition.UNRESOLVED,'No positive remaining DNA amount inferred.'))
    dna=selected
    def atp_buffer(key):
        stock=g.input(key+'_stock','Actual1xT4LigaseBuffer-with-ATP working contribution',K.MOLECULAR)
        buffer=g.product(key+'_buffer','Buffer component of same actual source working reagent',K.CONDITION)
        atp=g.product(key+'_atp','ATP co-present in same actual buffer; dose unassigned',K.MOLECULAR)
        g.react(key+'_components','Expose co-present ATP/buffer, not fresh extra ATP',S.GENERIC,((stock,R.STRUCTURAL),),(buffer,atp),changes_molecule=False)
        return ((buffer,R.CONDITION),(atp,R.STRUCTURAL))
    repair_nt=g.input('repair_nt','Actual0.125mMeach dATP,dCTP,dGTP,dTTP contribution; stock packaging unassigned',K.MOLECULAR)
    repair_buffer=atp_buffer('repair_buffer')
    dna = operation('repair', 'Source T4 DNA polymerase, Klenow and PNK end repair', S.END_REPAIR, dna, 'Actual6UT4DNApolymerase/2.5UKlenow/20UPNK catalyst contributions; stockcount unassigned',additional=((repair_nt,R.STRUCTURAL),)+repair_buffer)
    dna = recover('repair_cleanup', dna, 'Phenol:chloroform:isoamyl alcohol extraction and ethanol precipitation',molecular_reagent=True)
    tail_nt=g.input('tail_datp','Actual0.5mM dATP molecular contribution',K.MOLECULAR)
    tail_buffer=g.input('tail_buffer','Actual1xNEBBuffer2',K.CONDITION)
    dna = operation('a_tail', 'Source exo-minus Klenow adenylation', S.TAILING, dna, 'Actual10UKlenowexo-minus catalyst',additional=((tail_nt,R.STRUCTURAL),(tail_buffer,R.CONDITION)))
    dna = recover('tail_cleanup', dna, 'Phenol extraction and ethanol precipitation',molecular_reagent=True)
    adapter = g.input('methylated_adapter', 'Preannealed ATDBio Illumina genomic adapters: all source cytosines supplied as5mC; exact bodies/end geometry unassembled', K.ADAPTER)
    ligation_buffer=atp_buffer('ligation_buffer')
    ligase = g.input('ligase', 'Actual2000UT4DNA ligase;16°C16–20h', K.CATALYST)
    ligated = g.product('ligated', 'Adapter-ligated DNA plus free methylated adapters, ligase and possible unligated DNA', K.MOLECULAR)
    boundary = AdapterLigationGeometry('unavailable', 'unavailable', (molecule("Source adapter boundary", (role("Methylated adapter and substrate end geometry unassembled"),)),), 'Historical methylated adapter bodies, substrate ends and complete junction geometry not yet assembled; not a source-absence claim.')
    g.react('ligate', 'Ligate source methylated adapters before bisulfite treatment', S.LIGATION, ((dna, R.TEMPLATE), (adapter, R.STRUCTURAL), (ligase, R.CATALYST))+ligation_buffer, (ligated,), geometry=boundary)
    dna = recover('ligation_cleanup', ligated, 'Phenol extraction and ethanol precipitation',molecular_reagent=True)
    carrier = g.input('carrier', 'Source50ng sheared/dephosphorylated E.coli K12 DNA carrier; no assumed adapter or experimental barcode', K.MOLECULAR)
    mixed = g.product('carrier_mix', 'Adapter-ligated ChIP DNA and physically distinct carrier DNA', K.MOLECULAR)
    g.react('add_carrier', 'Add source carrier for gel isolation and bisulfite conversion', S.GENERIC, ((dna, R.TEMPLATE), (carrier, R.SPIKE_IN)), (mixed,), changes_molecule=False)
    low = g.product('low', 'Low-size source fraction:80–280bp insert plus possible carrier/carryover', K.MOLECULAR)
    high = g.product('high', 'High-size source fraction:280–430bp insert plus possible carrier/carryover', K.MOLECULAR)
    rejected = g.product('gel_reject', 'Physically excluded gel fractions; composition unassigned', K.WASTE)
    gel = g.input('gel', 'Source2.5%Nusieve gel; cut apparent200–400bp and400–550bp regions separately, MinElute extraction', K.SELECTION)
    g.react('size_split', 'Physically isolate low and high libraries; do not pool molecules before conversion', S.CLEANUP, ((mixed, R.TEMPLATE), (gel, R.SELECTION)), (low, high, rejected))
    g.discard(rejected, 'Actual unexcised gel fractions; migration ranges are not exact final molecular lengths.')
    for branch, sample in (('low', low), ('high', high)):
        # Diagnostic PCR consumes an aliquot; it is never the preparative bisulfite input.
        def diagnostic(key, material, stage):
            main=g.product(key+'_retained','Retained source preparation after analytical consumption; actual volume unassigned',K.MOLECULAR)
            aliquot=g.product(key+'_aliquot','Actual three0.3ul template aliquots for source diagnostic cycle-comparison reactions',K.MOLECULAR)
            g.react(key,'Allocate actual diagnosticPCR material before separate retained preparation'+stage,S.SPLIT,((material,R.TEMPLATE),),(main,aliquot))
            supplies=g.input(key+'_supplies','Actual three10ul diagnosticPCR reagent inventories; packaging/stockcount unassigned',K.MOLECULAR)
            f=g.product(key+'_f','LPX1.1 diagnostic primer0.3uM; body unreviewed',K.OLIGO)
            r=g.product(key+'_r','LPX2.1 diagnostic primer0.3uM; body unreviewed',K.OLIGO)
            cat=g.product(key+'_cat','Actual0.5UPfuTurboCx per diagnostic10ul reaction',K.CATALYST)
            nt=g.product(key+'_nt','Actual0.25mMeach diagnostic dNTP contribution',K.MOLECULAR)
            buffer=g.product(key+'_buffer','Actual1xTurboCx buffer contribution',K.CONDITION)
            g.react(key+'_components','Expose actual diagnostic supply classes; not a documented premix',S.GENERIC,((supplies,R.STRUCTURAL),),(f,r,cat,nt,buffer),changes_molecule=False)
            qc=g.product(key+'_qc','Actual diagnostic PCR collection for analytical gel, not preparative material',K.QC)
            g.react(key+'_amplify','Source three independent diagnostic10ul reactions; cycle comparison'+('15/19/22' if 'after' in stage else '10/15/20'),S.PCR,((aliquot,R.TEMPLATE),(f,R.STRUCTURAL),(r,R.STRUCTURAL),(cat,R.CATALYST),(nt,R.STRUCTURAL),(buffer,R.CONDITION)),(qc,),geometry=PCRGeometry.unavailable(Scene(),reason='Source diagnostic reaction recipe known; exact LPX/whole-template geometry unreviewed.'))
            g.mark_terminal(Terminal(qc.id,Disposition.QC,'Source analytical PCR/gel only; not carried into bisulfite or preparative PCR.'))
            return main
        sample = diagnostic(branch + '_pre_qc', sample, ' before bisulfite')
        sample = operation(branch + '_bisulfite', 'Source two successive5h EpiTect treatments, FFPE program', S.MAGIC, sample, 'Historical EpiTect Bisulfite kit; exact edition/intermediate strand and degradation chemistry unassembled', output_label='Bisulfite-processed insert and carrier plus retained methylated-adapter regions; possible unconverted/degraded molecules unquantified', reason='Source conversion is known; exact complete polymer/adduct/recovery geometry and historical kit operations not assembled. PhysicalU is distinct from later copiedT.')
        sample = diagnostic(branch + '_post_qc', sample, ' after bisulfite')
        pcr_input = g.product(branch + '_pcr_input', 'Source24µl across eight3µl preparative PCR aliquots', K.MOLECULAR)
        backup = g.product(branch + '_backup', 'Source untransferred bisulfite DNA after diagnostic consumption and24ul preparative allocation; remaining volume unassigned, stored−80C', K.MOLECULAR)
        g.react(branch + '_pcr_split', 'Retain backup separately from preparative PCR aliquots', S.SPLIT, ((sample, R.TEMPLATE),), (pcr_input, backup))
        g.mark_terminal(Terminal(backup.id, Disposition.STORED, 'Source backup stored at−80°C; not chemically destroyed or added to the library.'))
        f = g.input(branch + '_lpx1', 'Source Illumina LPX1.1; full body unreviewed', K.OLIGO)
        r = g.input(branch + '_lpx2', 'Source Illumina LPX2.1; full body unreviewed', K.OLIGO)
        nts=g.input(branch+'_pcr_nt','Actual0.25mMeachdNTP for each25ul preparativePCR; stockcount unassigned',K.MOLECULAR)
        buffer=g.input(branch+'_pcr_buffer','Actual1xTurboCx buffer for preparativePCR',K.CONDITION)
        enzyme = g.input(branch + '_polymerase', 'Source1.25UPfuTurboCx Hotstart per25ul preparative reaction; cycle count chosen from analytical gels, not universally22', K.CATALYST)
        product = g.product(branch + '_pcr', 'Actual eight-PCR reaction collection with new DNA, residual original templates/carrier and reagents; no one-product cloning', K.MOLECULAR)
        g.react(branch + '_amplify', 'Source allocate3ul per eight independent25ul preparativePCR reactions BEFORE amplification; modified templates distinct from copies', S.PCR, ((pcr_input, R.TEMPLATE), (f, R.STRUCTURAL), (r, R.STRUCTURAL), (enzyme, R.CATALYST),(nts,R.STRUCTURAL),(buffer,R.CONDITION)), (product,), geometry=PCRGeometry.unavailable(molecule("Source PCR boundary", (role("LPX primer/full-template geometry unassembled"),)), reason='Exact historical LPX primers, methylated adapter topology and full product strands unreviewed.'))
        product = recover(branch + '_column', product, 'MinElute PCR purification')
        product = recover(branch + '_final_gel', product, 'Source final gel purification without added carrier; no perfect carrier depletion asserted')
        lib = g.product(branch + '_library', 'Source'+branch+' library after gel/qc/Quant-iT; exact duplex/end geometry unassembled', K.LIBRARY)
        gel_qa=g.product(branch+'_final_gel_qa','Source final analytical-gel material; amount unassigned',K.QC)
        quant_qa=g.product(branch+'_quant_qa','Source Quant-iT dsDNA analytical material; amount unassigned',K.QC)
        g.react(branch+'_assessment','Separate actual final-gel/Quant-iT allocations from retained production library',S.SPLIT,((product,R.TEMPLATE),),(gel_qa,quant_qa,lib))
        g.mark_terminal(Terminal(gel_qa.id,Disposition.QC,'Source analytical gel material, not loaded productionDNA.'))
        g.mark_terminal(Terminal(quant_qa.id,Disposition.QC,'Source analytical Quant-iT material, not productionDNA.'))
        g.sequence(lib, SequencingPlan(limitation='Source GenomeAnalyzerIIx or HiSeq2000 standard36-base protocol; exact acquired index/read-primer geometry unassembled. Low/high reads combined computationally only after sequencing.'))
    g.add_source_item(SourceItem('predecessors', 'Unreviewed directly cited originals and historical kit geometry', 'Denissov2007 ChIP;Smith2009/Gu2010/Gu2011 library preparation; exact source versions/components not inherited.'))
    return g


def render_page():
    g = protocol_graph()
    return '\n'.join([head('ChIP-BS-seq'), '<main class="wrap"><h1>ChIP-BS-seq</h1>', '<p><a href="01_chip_bs.html">Research notes</a></p>', info(SOURCE), info('Private incomplete operational model. Original cytosine modification, bisulfite-treated template and copied DNA are distinct. Historical primers and full substrate/end geometry remain unassembled; no source-absence or complete validation claim.'), workflow_panel(workflow_from_protocol_graph(g, 'low_library')), workflow_panel(workflow_from_protocol_graph(g, 'high_library')), '</main>'])


def selftest():
    from checks import Check
    c = Check(); g = protocol_graph()
    c('Low and high preparations have separate sequencing endpoints', set(g.sequencing), {'low_library', 'high_library'})
    c('Carrier is a spike-in, not an adapter or barcode', g.steps['add_carrier'].inputs[1].role, R.SPIKE_IN)
    c('Methylated adapter is installed before carrier and conversion', g.steps['ligation_cleanup_step'].inputs[0].material, 'ligated')
    c('Analytical pre-PCR is not preparative bisulfite template', g.steps['low_bisulfite'].inputs[0].material, 'low_pre_qc_retained')
    c('Postconversion QC remains distinct from preparative PCR', g.steps['high_pcr_split'].inputs[0].material, 'high_post_qc_retained')
    c('Stored backup has an explicit nonsequencing fate', g.terminals['low_backup'].disposition, Disposition.STORED)
    c('End repair receives actual selected100ng ChIPDNA',g.steps['repair'].inputs[0].material,'chip_100ng')
    c('Source ATP derives from actual working buffer',g.producer_of('repair_buffer_atp').id,'repair_buffer_components')
    c('Diagnostic PCR material does not feed preparative template',g.steps['low_bisulfite'].inputs[0].material,'low_pre_qc_retained')
    c('Final analytical fractions are QC, not production library',tuple(g.terminals[k].disposition for k in ('low_final_gel_qa','low_quant_qa')),(Disposition.QC,Disposition.QC))
    c('Unassembled geometry prevents public completion', bool(g.publication_problems()))
    c('Private preview has no finished marker', 'data-protocol-complete=' not in render_page())
    c.report()
