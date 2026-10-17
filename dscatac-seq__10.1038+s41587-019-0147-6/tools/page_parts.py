"""Sparse individual dscATAC-seq and dsciATAC-seq renderers."""
from __future__ import annotations
import sys
from pathlib import Path
HERE=Path(__file__).resolve().parent; sys.path[:0]=[str(HERE),str(HERE.parents[1]/"lib")]
import dscatac as D
import nextera as nx
import seqprimers as sp
from chemdraw import (Construct,Scene,Segment,annotation_rows,oligo,panel,
                      protocol_graph_panels,strand_row)
from page import head,info,table

def seg(n,t,g=None,**kw): return Segment(n,t,g,**kw)
def render_page(protocol, *, graph_factory=None):
    graph=graph_factory() if graph_factory else D.protocol_graph(protocol)
    sci=protocol=="dsciATAC-seq"
    notes=("01_dscatac-seq.html" if not sci else
           "../dscatac-seq__10.1038+s41587-019-0147-6/01_dscatac-seq.html")
    bead=oligo("published 3' portion of bead oligo",D.bead_published_segments())
    tn5=oligo(("barcoded Tn5 Read-1 adaptor" if sci else "Tn5 Read-1 adaptor"),D.tn5_r1(protocol))
    frag=D.filled_fragment(protocol); final=D.final_library(protocol)
    r1=[seg("BC1","N"*7,"cbc",placeholder=True),seg("phase","NN",placeholder=True),
        seg("constant 1",D.C1),seg("BC2","N"*7,"cbc",placeholder=True),seg("constant 2",D.C2),
        seg("BC3","N"*7,"cbc",placeholder=True),seg("s5",nx.S5,"s5")]
    if sci:r1.append(seg("Tn5 barcode","N"*6,"cbc",placeholder=True))
    r1 += [seg("ME",nx.ME,"me"),seg("genome","XXXXXXXX...",placeholder=True)]
    rows=[("Read 1 (118)","bead barcode"+(" + Tn5 barcode" if sci else "")+" + genomic DNA"),
          ("Index 1 (8)","sample index"),("Read 2 (40)","genomic DNA from the s7 end")]
    return "\n".join([head(f"{protocol} library chemistry"),'<div class="wrap" data-protocol-complete="true">',
      f'<h1>{protocol} &mdash; droplet-barcoded single-cell ATAC-seq</h1>',
      f'<p class="research-notes"><a href="{notes}">Research notes</a></p>',
      info('Defining source: <a href="https://doi.org/10.1038/s41587-019-0147-6">Lareau et al., <i>Nature Biotechnology</i> (2019)</a>.'),
      '<div class="caveat"><b>Bead and primer sequence boundary.</b> The paper publishes the bead barcode/constant/s5 portion but leaves its linker cell empty. Dotted 5\' sequence marks bases unresolved from the checked originals. Bio-Rad guides identify the custom Read 1 primer as catalog 16005986 but do not print its bases; the publication-era 2018 guide and reagent insert remain unavailable.</div>',
      '<div class="caveat"><b>Droplet programme boundary.</b> The Methods give the in-droplet programme as 37 &deg;C 30 min, 85 &deg;C 10 min, 72 &deg;C 5 min, then 8 PCR cycles; the Supplementary Note says bead oligos are cleaved from the bead. The cleavage chemistry and the purpose of the 85 &deg;C step are not stated (warned boundary), and the 72 &deg;C step is modelled as Tn5 gap fill by inference. Beads are super-loaded, so a droplet may hold more than one bead; such barcodes are merged computationally (BAP), not removed. ' + ('Cells are also loaded above Poisson and distinguished by their Tn5 barcode. ' if sci else '') + 'The droplet index primer is drawn from the Supplementary Table 6 N70X layout; the Methods name only Index Kit 12009360, and the second-PCR primers are not named.</div>',
      f'<h2>Key oligos</h2><seq>{tn5}{bead}</seq>',
      '<h2>Reaction workflow</h2>',protocol_graph_panels(graph,cls="long"),
      '<h2>Molecular details</h2>',
      ('<h3>(1) Barcode Tn5 in wells, tagment, then pool intact cells</h3>' if sci else '<h3>(1) Tagment intact cells or nuclei in bulk</h3>'),
      panel([*Scene.duplex(list(frag)).rows(),*annotation_rows(frag)],cls="long",caption=("Gap-filled representative fragment; the s5 end carries the first, Tn5 barcode." if sci else "Gap-filled representative s5/s7 fragment.")),
      '<h3>(2) Co-encapsulate tagmented cells with super-loaded barcode beads</h3>',
      panel(D.bead_priming_scene().rows(),cls="long",caption="The released bead oligo primes through its published 3' s5 sequence."),
      '<h3>(3) Droplet PCR copies the bead barcode onto s5/s7 fragments</h3>',
      panel([strand_row(final),*annotation_rows(final)],cls="long",caption="INFERRED — complete supported inner structure; dotted bead/P5/Read-1 arm stops where the checked article, Supplementary Table 6 and Bio-Rad guide stop specifying bases."),
      '<h2>Read layout</h2>',sp.unavailable_diagram(final,"the custom Read-1 and library-primer bases are unresolved from the checked article, Supplementary Table 6 and Bio-Rad guide"),panel([strand_row(Construct(r1),prefix="Read 1  ",suffix="")],cls="long",caption="Custom Read 1 begins at BC1; its bases remain unresolved from the checked article, Supplementary Table 6 and Bio-Rad guide."),
      table(("Read","Reports"),rows),
      sp.read_spec_table([
          sp.read_spec("Read 1",118,
              "bead barcode" + (", Tn5 barcode" if sci else "") + ", then genomic DNA",
              "extract_then_include"),
          sp.read_spec("Index 1 (i7)",8,"sample index","demultiplex"),
          sp.read_spec("Read 2",40,"genomic DNA from the s7 end","include")],
          "Read roles and cycle counts are source-declared; custom primer bases remain unresolved from the checked originals, and the publication-era guide is unavailable."),
      '</div>'])
