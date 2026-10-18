from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent
sys.path[:0]=[str(HERE),str(HERE.parents[1]/"lib")]
import chem
from checks import Check
c=Check();g=chem.protocol_graph()
for p in ('s_phase','g1'):
    c(p+' WGA geometry uses native pre-amplification',g.steps[p+'_wga_step'].geometry.template is g.materials[p+'_preamp'].state,True)
    c(p+' library PCR uses allocated bulk template',g.steps[p+'_library_pcr'].geometry.template is g.materials[p+'_bulk_input'].state,True)
    c(p+' diagnostic PCR stays terminal QC',g.terminals[p+'_qpcr_product'].disposition.value,'QC')
    c(p+' no inferred index-cycle length','Index 1 (i7)' in g.sequencing[p+'_library'].read_lengths,False)
c.report()
