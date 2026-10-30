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
    def operation(key, label, kind, material, reagent_label, output_label=None, **kw):
        enzyme = kind in {S.END_REPAIR, S.TAILING}
        reagent = g.input(key + '_reagent', reagent_label, K.CATALYST if enzyme else K.CONDITION)
        product = g.product(key + '_material', output_label or label + '; retained mixture and carryover unquantified', K.MOLECULAR)
        g.react(key, label, kind, ((material, R.TEMPLATE), (reagent, R.CATALYST if enzyme else R.CONDITION)), (product,), **kw)
        return product
    def recover(key, material, label):
        reagent = g.input(key + '_procedure', label, K.SELECTION)
        retained = g.product(key, 'Physically recovered source preparation; exact recovery and purity unassigned', K.MOLECULAR)
        waste = g.product(key + '_waste', 'Source-separated fractions; exact molecular composition unassigned', K.WASTE)
        g.react(key + '_step', label, S.CLEANUP, ((material, R.TEMPLATE), (reagent, R.SELECTION)), (retained, waste))
        g.discard(waste, 'Source physical separation; no ideal depletion or recovery claim.')
        return retained
    dna = g.input('chromatin', 'Source H3K27me3-associated chromatin; harvesting/ChIP predecessor unreviewed', K.BIOLOGICAL, primary=True)
    dna = operation('chip', 'Source anti-H3K27me3 ChIP and captured DNA isolation', S.MAGIC, dna, 'Source anti-H3K27me3 and cited Denissov2007 procedure', reason='Actual source identifies ChIP; exact upstream crosslinking, fragment and isolation geometry remains unreviewed. MethylCap is a separate precursor branch, not an antibody substitute.')
    dna = operation('repair', 'Source T4 DNA polymerase, Klenow and PNK end repair', S.END_REPAIR, dna, 'Source end-repair enzymes, four dNTPs and ATP;20°C30min')
    dna = recover('repair_cleanup', dna, 'Phenol:chloroform:isoamyl alcohol extraction and ethanol precipitation')
    dna = operation('a_tail', 'Source exo-minus Klenow adenylation', S.TAILING, dna, 'Klenow exo-minus/dATP/NEB buffer2;37°C30min')
    dna = recover('tail_cleanup', dna, 'Phenol extraction and ethanol precipitation')
    adapter = g.input('methylated_adapter', 'Preannealed ATDBio Illumina genomic adapters: all source cytosines supplied as5mC; exact bodies/end geometry unassembled', K.ADAPTER)
    ligase = g.input('ligase', 'T4 DNA ligase/ATP;16°C16–20h', K.CATALYST)
    ligated = g.product('ligated', 'Adapter-ligated DNA plus free methylated adapters, ligase and possible unligated DNA', K.MOLECULAR)
    boundary = AdapterLigationGeometry('unavailable', 'unavailable', (molecule("Source adapter boundary", (role("Methylated adapter and substrate end geometry unassembled"),)),), 'Historical methylated adapter bodies, substrate ends and complete junction geometry not yet assembled; not a source-absence claim.')
    g.react('ligate', 'Ligate source methylated adapters before bisulfite treatment', S.LIGATION, ((dna, R.TEMPLATE), (adapter, R.STRUCTURAL), (ligase, R.CATALYST)), (ligated,), geometry=boundary)
    dna = recover('ligation_cleanup', ligated, 'Phenol extraction and ethanol precipitation')
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
            main = g.product(key + '_retained', 'Retained source'+stage+' preparation', K.MOLECULAR)
            aliquot = g.product(key + '_aliquot', 'Source analytical PCR aliquot'+stage, K.QC)
            g.react(key, 'Aliquot analytical PCR separately from retained material'+stage, S.SPLIT, ((material, R.TEMPLATE),), (main, aliquot))
            g.mark_terminal(Terminal(aliquot.id, Disposition.QC, 'Analytical PCR/gel assessment only; not carried forward as the preparative template.'))
            return main
        sample = diagnostic(branch + '_pre_qc', sample, ' before bisulfite')
        sample = operation(branch + '_bisulfite', 'Source two successive5h EpiTect treatments, FFPE program', S.MAGIC, sample, 'Historical EpiTect Bisulfite kit; exact edition/intermediate strand and degradation chemistry unassembled', output_label='Bisulfite-processed insert and carrier plus retained methylated-adapter regions; possible unconverted/degraded molecules unquantified', reason='Source conversion is known; exact complete polymer/adduct/recovery geometry and historical kit operations not assembled. PhysicalU is distinct from later copiedT.')
        sample = diagnostic(branch + '_post_qc', sample, ' after bisulfite')
        pcr_input = g.product(branch + '_pcr_input', 'Source24µl across eight3µl preparative PCR aliquots', K.MOLECULAR)
        backup = g.product(branch + '_backup', 'Source remainder of40µl bisulfite DNA stored at−80°C', K.MOLECULAR)
        g.react(branch + '_pcr_split', 'Retain backup separately from preparative PCR aliquots', S.SPLIT, ((sample, R.TEMPLATE),), (pcr_input, backup))
        g.mark_terminal(Terminal(backup.id, Disposition.STORED, 'Source backup stored at−80°C; not chemically destroyed or added to the library.'))
        f = g.input(branch + '_lpx1', 'Source Illumina LPX1.1; full body unreviewed', K.OLIGO)
        r = g.input(branch + '_lpx2', 'Source Illumina LPX2.1; full body unreviewed', K.OLIGO)
        enzyme = g.input(branch + '_polymerase', 'Source PfuTurboCx Hotstart; cycle count chosen from analytical gels, not universally22', K.CATALYST)
        product = g.product(branch + '_pcr', 'New PCR DNA plus residual original bisulfite templates, carrier and reaction material', K.MOLECULAR)
        g.react(branch + '_amplify', 'Source eight preparative25µl reactions; original modified templates remain distinct from copied DNA', S.PCR, ((pcr_input, R.TEMPLATE), (f, R.STRUCTURAL), (r, R.STRUCTURAL), (enzyme, R.CATALYST)), (product,), geometry=PCRGeometry.unavailable(molecule("Source PCR boundary", (role("LPX primer/full-template geometry unassembled"),)), reason='Exact historical LPX primers, methylated adapter topology and full product strands unreviewed.'))
        product = recover(branch + '_column', product, 'MinElute PCR purification')
        product = recover(branch + '_final_gel', product, 'Source final gel purification without added carrier; no perfect carrier depletion asserted')
        lib = g.product(branch + '_library', 'Source'+branch+' library after gel/qc/Quant-iT; exact duplex/end geometry unassembled', K.LIBRARY)
        g.react(branch + '_assessment', 'Source analytical gel and Quant-iT dsDNA assessment', S.GENERIC, ((product, R.TEMPLATE),), (lib,), changes_molecule=False)
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
    c('Unassembled geometry prevents public completion', bool(g.publication_problems()))
    c('Private preview has no finished marker', 'data-protocol-complete=' not in render_page())
    c.report()
