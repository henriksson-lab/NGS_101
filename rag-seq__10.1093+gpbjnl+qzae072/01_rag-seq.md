# RAG-seq

[Published defining article](https://doi.org/10.1093/gpbjnl/qzae072).

🟢 Selected branch: the published detailed Methods describe RAG-seq3.0. Versions1.0/2.0 remain components/earlier branches of this entry. We do not combine an invented1.0 RT with3.0 downstream choices. A selected single-cell sample is isolated from source PBS/1%BSA into2µllysisbuffer (0.5%Triton,4URiboLock,DNasebuffer,0.4UDNaseI),20°C30min. Add1µl0.12U/µlthermolabilePK,25°C30min then55°C10min. Add50pmolNSR and5fmolP5-oligo-dT,0.8µl25mMdNTP andwater;72°C3min/ice. Add5µlRTmixture containing100UMaximaHMinus,10URiboLock,2µl5Mbetaine;25°C5min/37°C15min/40°C90min/70°C15min.

🟢 Tn5 adaptor factory: equimolar100µM T7-Tn5-MEB and5P ME-rev,94°C5min then−1°C/min to10°C; actual Tn5 loading delegates cited predecessor53 and remains unassembled. Whole upper oligo `GATCACTAATACGACTCACTATAGGGTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG`; lower `/Phos/CTGTCTCTTATACACATCT`. Compute exact19bp ME duplex locally, without inheriting SHERRY/Tn5 loading formulation. Select1µl produced loaded-stock aliquot; stock remainder retained as unresolved.

🟢 Add10µltagmentationmixture:1µltransposome,4µl5×TD(50mMTris/25mMMgCl2/50%DMF/4.25mMATP),4.5µl40%PEG8000/0.5µlwater;55°C30min. Add1µl0.4%SDS,55°C7min. Bst3.0(1µl8U/µl) and1×Q5mix fill/extend72°C15min/80°C5min. Fresh3µl2×Q5,4.5pmolT7customB andcustomP5 each;94°C2min then2cycles low40°C/2cycles60°C plus n=11 forsinglecell (8for10cells/1ng,5for100cells/10ng), sourceextensionadds10s percycle. Whole occupied RNA/hybrid junctions, products and strand-transfer distribution remain unassembled.

🟢 DASH:27published18S/28Stargets (not whole sgRNA-scaffold orderings). Source T7IVT guide production, purification/equalpool/storage are opaque precursor boundary. Cas9:sgRNA1:2,1×buffer3.1,37°C15min;5µlcomplex plus50µlPCRproduct37°C2h. PK1µl37°C15min/55°C10min;1×AMPure, elute5µl **do not discard beads**. DNA/beads continue into5µlIVTmix (1µl10×buffer,2µl10mMNTP,0.5µlRNaseinhibitor,0.5µl0.1U/µlYIPP,1µl50U/µlT7polymerase),37°C10–16h. Add45µlTrispH8→55µl then55µlHXP(20%PEG/2.5MNaCl); eluteRNA11µl.

🟢 SecondRT: whole11µlRNA+1µl10mMdNTP+1µl10µMP5primer65°C5min/ice≥2min;4µl5×buffer/1µlinhibitor/1µlSSIII200U/µl/1µlwater25°C5min/50°C1h/70°C15min. Select10µlcDNA and retain remainder. FinalPCR50µl:0.5µl10µMP7index,0.5µl10µMP5,25µl2×Q5,14µlwater;3–5cycles98°C30s/60°C20s/72°C2min.0.8×AMPure/elute20µl/Qubit/NovaSeq6000. Exact actual study index and whole library/read landing remain unassembled.

🟢 Whole TableS2: NSR `TCCGATCTGANNNNNNN` (last6N underlined in actual workbook; the precedingN is not called UMI). P5-oligodT `CCTACACGACGCTCTTCCGATCTCT`+T30+VN. CustomP5 `CCTACACGACGCTCTTCCGATCTCT`; T7customB `GATCACTAATACGACTCACTATAGGGTCTCG`. FinalP5 `AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCTCT`; finalP7 `CAAGCAGAAGACGGCATACGAGAT`+[8ntindex]+`GTCTCGTGGGCTCGGAGATGT`. NSR binding sequences are source-degenerate classes rather than a claimed enumerated pool. Smart-seq3 comparators in the same table are not inherited.

🟡 Deferred: directly cited Tn5loading predecessor53, NSR pool design/application, whole sgRNA scaffolds, exact occupied hybrids/junctions/physical read masks. Three originals closed; no source absence or full validation claim.
