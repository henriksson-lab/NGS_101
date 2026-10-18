from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent
sys.path[:0]=[str(HERE),str(HERE.parents[1]/"lib")]
import scrrbs as S
from checks import Check
c=Check();g=S.protocol_graph()
c("firstPCR geometry uses recovered U-bearing input",g.steps["library_pcr"].geometry.template is g.materials["purified"].state,True)
c("secondPCR geometry uses actual source transfer",g.steps["second_pcr"].geometry.template is g.materials["second_template"].state,True)
c("no conventional index acquisition supplied",g.sequencing["library"].read_lengths,{})
c("tRNA is independent of DNA in recovered scene",len(g.materials["purified"].state.strands),5)
c.report()
