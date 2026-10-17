"""Source-owned original TUC-seq model (Riml et al. 2017, Angewandte SI)."""
from chemdraw import Construct, Scene, Segment, oligo
from pcr import PCRGeometry
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind
from protocol_model import molecule, role
import seqprimers as sp

TITLE='TUC-seq — original thiouridine-to-cytidine conversion'
NOTES='01_tuc_seq.html'
SOURCE='<a href="https://doi.org/10.1002/anie.201707465">Riml et al., Angewandte Chemie International Edition (2017)</a> and its Supporting Information (anie201707465-sup-0001). Later TUC-seq DUAL and protocol chapters are version-scoped comparisons.'
SUMMARY='HEK293 cells are pulse-labelled with 4sU; total RNA is treated with OsO4 in ammonium chloride, which converts 4sU into native cytidine; random-hexamer cDNA is amplified with gene-specific barcoded primers for seven transcripts, and the gel-purified amplicon pool is sequenced to count U-to-C changes. This native-C product differs from the substituted cytidine analogue in TimeLapse chemistry.'
SI='Riml et al. 2017 Supporting Information (anie201707465-sup-0001)'
CAVEAT=('Sources: the 2017 Supporting Information (metabolic labelling, RNA isolation, tRNA and amplicon methods, Supporting Fig. 1 and Table 4) and the PubMed abstract. The main article, whose Experimental Section the SI cites for the 4sU-to-C conversion, has not been read. '
        'Conversion conditions are therefore partial: Supporting Fig. 1 gives OsO4 450 µM in 180 mM ammonium chloride for synthetic oligoribonucleotides, and the amplicon section gives 3 h at room temperature for cellular RNA; pH, the cellular-RNA OsO4 concentration and the post-conversion purification are unassigned. '
        'The barcode overhangs on the gene-specific PCR primers are not printed, and library preparation and Illumina HiSeq sequencing were a GATC Biotech service with no adapter, index or read-length details; that step is a warned boundary. '
        'Published cross-reference discrepancy: the SI amplicon section (page 15) says the primer sequences "are listed in Supporting Table S1", but Supporting Table 1 lists T-to-C mutation frequencies and the amplicon primers are printed in Supporting Table 4 ("Amplicon Sequencing"). The model uses the Table 4 primers. '
        'The tRNA-Val stem-loop RT-PCR branch ends in cloning and Sanger sequencing of 40 clones and is not drawn. '
        'INFERRED — RNA classes distinguished by pulse-label incorporation and a possible unconverted labelled class express population possibilities, not physical fractions or quantitative conversion yields.')
RECODING_SOURCE=(SI + ': Supporting Fig. 1 legend (OsO4 450 µM in 180 mM ammonium chloride, synthetic 17- and 42-nt 4sU RNAs) and the amplicon-sequencing section ("treated with OsO4 as described above for 3 h at room temperature"); the SI defers the conversion protocol to the unread main-text Experimental Section.')
AMPLICON_PRIMERS=(
    ('CcnA2', 'CCAGAAGTAGCAGAGTTTGT', 'TTGAGGAGAGAAACACCATG'),
    ('CcnB1', 'ACATCGAAGCATGCTAAGAT', 'CTATGCAGCAGATTCTCCAT'),
    ('CcnD1', 'GAGGGCAGTTTTCTAATGGA', 'ATCAAGGGGAGATTGCATTT'),
    ('CcnE1', 'CTGATGAAGATGCACACAAC', 'CTTTTGTTGTTGTGGGAGTC'),
    ('p21/Cdkn1A', 'CTTGAGTGGGGTTATCTCTG', 'ATATTCAGCATTGTGGGAGG'),
    ('eGFP', 'CCATCTTCTTCAAGGACGAC', 'TACTTGTACAGCTCGTCCAT'),
    ('PCNA', 'ACCAAACCAGGAGAAAGTTT', 'TCCTTCTTCATCCTCGATCT'),
)
FINAL_LIBRARY=Construct((role('GATC Biotech left adapter end; unassigned'),role('forward sample barcode overhang (not printed)'),
                         role('gene-specific amplicon with U-to-C readout'),role('reverse sample barcode overhang (not printed)'),
                         role('GATC Biotech right adapter end; unassigned')),
                        name='original TUC-seq barcoded amplicon library (service library preparation unassigned)')
FINAL_CAPTION='Barcoded gene-specific amplicons as far as the SI describes them. The barcode overhangs, the service-provider library adapters, indices and loaded run primers are not printed.'
SEQ_PRIMERS=()
SEQUENCING_ENDING='SI: reads were demultiplexed by sample-specific forward and reverse barcodes (flexbar), aligned to amplicon references (bowtie2) and U-to-C positions called within 146 bases of the amplicon ends. Illumina HiSeq (GATC Biotech); read length not stated.'
SEQUENCING_UNAVAILABLE='The service-provider library architecture and the acquired read layout are not printed in the SI; the barcode overhang sequences are also not printed.'
READ_SPECS=(sp.read_spec('Read 1','unassigned','sample barcode then gene-specific amplicon with U-to-C conversions; exact run/read primer unassigned','undetermined'),)


def oligos():
    for gene, fw, rv in AMPLICON_PRIMERS:
        yield oligo(f'{gene}_fw (Supporting Table 4; barcode overhang not printed)', [Segment('gene-specific', fw)])
        yield oligo(f'{gene}_rv (Supporting Table 4; barcode overhang not printed)', [Segment('gene-specific', rv)])


def rna_pool(converted=False):
    sc=Scene()
    sc.strand('old',[role('RNA without a pulse-installed label')],label='RNA without a pulse-installed label')
    mark=Segment('native cytidine formed from 4sU','C','r1') if converted else Segment('pulse-labelled 4sU','4sU',placeholder=True,length_bp=1)
    sc.strand('new',[role('RNA flank before the marked site; synthesis time unassigned'),mark,role('RNA continuation after the marked site; synthesis time unassigned')],
        label='RNA with a representative pulse-derived converted site' if converted else 'RNA with a representative pulse-incorporated 4sU site')
    sc.note('new','Representative site only; conversion completeness and fraction sizes are not asserted.')
    if converted:
        sc.strand('possible_unconverted',[role('possible pulse-labelled RNA with unconverted 4sU',inferred=True)],label='INFERRED — possible unconverted labelled RNA')
    return sc

def _primer_set(direction):
    col = 1 if direction == 'fw' else 2
    return [Construct([role(f'{"forward" if direction == "fw" else "reverse"} sample barcode overhang (not printed)'),
                       Segment(f'{row[0]}_{direction}', row[col])], name=f'{row[0]}_{direction}')
            for row in AMPLICON_PRIMERS]


def protocol_graph():
    graph=ProtocolGraph(TITLE,authority=SI)
    cells=graph.input('cells','HEK293 (293-HA-StrepIII-eGFP) cells, 5 × 10^5 per well, Dox-inducible eGFP',MaterialKind.BIOLOGICAL,state=molecule('cells',(role('cellular RNA populations'),)),primary=True)
    label=graph.input('four_su','4-thiouridine 0.05 or 0.1 mM in DMEM, 1 h; doxycycline 1 µg/ml added after 30 min',MaterialKind.MOLECULAR)
    polymerase=graph.input('cellular_rna_polymerase','cellular RNA polymerases',MaterialKind.CATALYST)
    labelled=graph.product('labelled_rna','RNA classes distinguished by pulse-label incorporation',MaterialKind.MOLECULAR,state=rna_pool())
    graph.react('pulse_label','Metabolically label newly synthesized RNA',StepKind.TRANSCRIPTION,
        ((cells,InputRole.TEMPLATE),(label,InputRole.STRUCTURAL),(polymerase,InputRole.CATALYST)),(labelled,))
    trizol=graph.input('rna_isolation','TRIzol isolation, DNase I (NEB), phenol–chloroform extraction, 0.3 M NaAc precipitation',MaterialKind.CONDITION)
    ready=graph.product('rna_for_conversion','DNase-treated total RNA',MaterialKind.MOLECULAR,state=labelled.state)
    graph.react('rna_isolation_step','Isolate and DNase-treat total RNA',StepKind.EXTRACTION,
        ((labelled,InputRole.TEMPLATE),(trizol,InputRole.CONDITION)),(ready,))
    chemistry=graph.input('recoding_chemistry','OsO4 in ammonium chloride buffer, 3 h room temperature (Supporting Fig. 1: 450 µM OsO4, 180 mM NH4Cl for synthetic RNAs); pH and cellular-RNA concentration unassigned',MaterialKind.CONDITION,source=RECODING_SOURCE)
    recoded=graph.product('recoded_rna','RNA with representative native cytidine conversion',MaterialKind.MOLECULAR,state=rna_pool(True))
    graph.react('recode','Convert 4sU to native cytidine with OsO4/NH4Cl',StepKind.CHEMICAL_CONVERSION,
        ((ready,InputRole.TEMPLATE),(chemistry,InputRole.CONDITION)),(recoded,))
    hexamers=graph.input('random_hexamers','random hexamer primers',MaterialKind.OLIGO)
    goscript=graph.input('goscript','GoScript Reverse Transcriptase (Promega), after unspecified post-conversion purification',MaterialKind.CATALYST)
    cdna_state=Construct((role('random-hexamer-primed cDNA; converted sites read as G opposite C'),),name='TUC-seq cDNA (role level)')
    cdna=graph.product('cdna','random-primed first-strand cDNA',MaterialKind.MOLECULAR,state=cdna_state)
    graph.react('reverse_transcribe','Reverse-transcribe with random hexamers',StepKind.REVERSE_TRANSCRIPTION,
        ((recoded,InputRole.TEMPLATE),(hexamers,InputRole.STRUCTURAL),(goscript,InputRole.CATALYST)),(cdna,))
    fwd=graph.input('amplicon_forward','forward gene-specific primers (_fw) with sample barcode overhangs for CcnA2, CcnB1, CcnD1, CcnE1, p21/Cdkn1A, eGFP and PCNA (cores in Supporting Table 4; overhangs not printed)',MaterialKind.OLIGO,
        state=_primer_set('fw'),source=SI+' Supporting Table 4 (Amplicon Sequencing)')
    rev=graph.input('amplicon_reverse','reverse gene-specific primers (_rv) with sample barcode overhangs (cores in Supporting Table 4; overhangs not printed)',MaterialKind.OLIGO,
        state=_primer_set('rv'),source=SI+' Supporting Table 4 (Amplicon Sequencing)')
    pcr_pol=graph.input('pcr_polymerase','"standard PCR conditions" (polymerase and cycles not stated)',MaterialKind.CATALYST)
    amp_state=Construct((role('forward sample barcode overhang (not printed)'),role('gene-specific amplicon'),role('reverse sample barcode overhang (not printed)')),name='barcoded TUC-seq amplicon (role level)')
    amplicons=graph.product('amplicons','barcoded gene-specific amplicons',MaterialKind.MOLECULAR,state=amp_state)
    graph.react('amplicon_pcr','Amplify selected targets with barcoded gene-specific primers',StepKind.PCR,
        ((cdna,InputRole.TEMPLATE),(fwd,InputRole.STRUCTURAL),(rev,InputRole.STRUCTURAL),(pcr_pol,InputRole.CATALYST)),(amplicons,),
        geometry=PCRGeometry.unavailable(amp_state,reason='Supporting Table 4 prints only the gene-specific primer cores; barcode overhangs, target sequences and amplicon lengths are not printed, so primer landing cannot be drawn.'))
    gel=graph.input('agarose_gel','2% agarose gel separation, gel purification and equimolar pooling',MaterialKind.SELECTION)
    pool=graph.product('amplicon_pool','gel-purified, equimolar amplicon pool',MaterialKind.MOLECULAR,state=amp_state)
    gel_waste=graph.product('gel_excluded','gel material outside the excised amplicon bands',MaterialKind.WASTE)
    graph.react('gel_purify_pool','Gel-purify amplicons and pool equimolarly',StepKind.SELECTION,
        ((amplicons,InputRole.TEMPLATE),(gel,InputRole.SELECTION)),(pool,gel_waste))
    graph.discard(gel_waste,'left in the agarose gel')
    library=graph.product('library',FINAL_LIBRARY.name,MaterialKind.LIBRARY,state=FINAL_LIBRARY,inferred=True)
    graph.react('library_preparation','INFERRED — GATC Biotech library preparation from the amplicon pool (service; adapters unassigned)',StepKind.MAGIC,
        ((pool,InputRole.TEMPLATE),),(library,),reason='SI: "Library preparation from the amplicon pool and sequencing using the Illumina HighSeq platform was performed by GATC Biotech"; no adapter, index, end-processing or cycle information is printed.')
    graph.sequence(library,SequencingPlan(library=FINAL_LIBRARY,limitation=SEQUENCING_UNAVAILABLE))
    graph.declare_reviewed_scope(('cells','four_su','cellular_rna_polymerase','pulse_label','rna_isolation','rna_isolation_step',
        'recoding_chemistry','recode','random_hexamers','goscript','reverse_transcribe','amplicon_forward','amplicon_reverse','pcr_polymerase',
        'amplicon_pcr','agarose_gel','gel_purify_pool','library_preparation','library'),graph.authority)
    return graph.require_publishable()
