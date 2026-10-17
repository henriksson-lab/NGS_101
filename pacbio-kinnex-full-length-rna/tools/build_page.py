#!/usr/bin/env python3
from pathlib import Path
import sys
H=Path(__file__).resolve().parent; sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import kinnex as K
import seqprimers as sp
from chemdraw import annotation_rows, panel, protocol_graph_panels, strand_row
from dumbbell import dumbbell_rows
from page import head, info
OUT=H.parent/"kinnex-full-length-rna.html"

def render():
    cycle=K.polymerase_cycle()
    graph=K.protocol_graph()
    parts=[head(K.TITLE),'<div class="wrap" data-protocol-complete="true">',f'<h1>{K.TITLE}</h1>',
           '<p class="research-notes"><a href="01_kinnex-full-length-rna.html">Research notes</a></p>',
           info(K.SOURCE),info(K.SUMMARY),f'<div class="caveat">{K.CAVEAT}</div>']
    parts.extend(['<h2>Source-known analysis references</h2>',
                  '<p>PacBio’s <a href="https://www.pacb.com/wp-content/uploads/Kinnex_barcoded_adapter_indexes.fasta_.zip">Kinnex index FASTA</a> prints the four terminal-adapter barcode cores: '
                  + ', '.join(f'{name} <code>{bases}</code>' for name,bases in K.KINNEX_TERMINAL_BARCODES.items()) + '. '
                  'Its <a href="https://downloads.pacbcloud.com/public/dataset/Kinnex-full-length-RNA/REF-MAS_adapters/MAS-Seq_Adapter_v3/mas8_primers.fasta">MAS8 FASTA</a> '
                  'provides the segmentation recognition motifs below; these short references do not establish full physical adapter or PCR-primer sequences.</p>',
                  '<p>' + ', '.join(f'{name} <code>{bases}</code>' for name,bases in K.MAS8_REFERENCE_MOTIFS.items()) + '</p>'])
    parts.append('<p>The <a href="https://www.pacb.com/wp-content/uploads/SMRT-Link-User-Guide-v13.1.pdf">SMRT Link v13.1 guide</a>, p. 101, gives the Iso-Seq v2 primer family. '
                 f'Example bc01 5′ trimming reference: <code>{K.ISOSEQ_V2_BC01_5P_REFERENCE}</code>; 3′ trimming reference: <code>{K.ISOSEQ_V2_3P_REFERENCE}</code>. '
                 'The chosen bc01–bc12 depends on the sample. Page 100 specifies that the 3′ reference excludes poly(T) and is reverse-complemented; these records are for primer removal, not full physical oligo ordering.</p>')
    parts.extend(['<h2>Reaction workflow</h2>',
                  protocol_graph_panels(graph, cls="long")])
    parts.extend(['<h2>Sequencing-primer and polymerase binding</h2>',
        panel(dumbbell_rows(K.SMRTBELL),cls="long",caption="INFERRED — the Kinnex sequencing primer (103-238-700 Rev09 section 7) anneals in a terminal adapter whose bases have not been verified from the checked guide. Its exact binding site and extension direction remain unresolved."),
        '<h2>HiFi passes and segmentation</h2>',
        panel([strand_row(cycle),*annotation_rows(cycle)],cls="long",caption="INFERRED — one schematic polymerase circuit traverses the eight-insert array in both orientations. Consensus is formed first; the seven ordered segmentation junctions then delimit the eight original cDNAs."),
        sp.read_spec_table([sp.read_spec('Sequencing output',None,
            'HiFi consensus of an eight-cDNA segmented SMRTbell array','include')],
            'PacBio records repeated polymerase passes; cycle count is molecule- and movie-dependent.'),'</div>'])
    return '\n'.join(parts)

if __name__=='__main__': OUT.write_text(render(),encoding='utf-8'); print(f'wrote {OUT}')
