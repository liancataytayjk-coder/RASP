#!/usr/bin/env python3
"""Stage 3: write BEAST 1.10.4 XML (GTR+G, UCLD relaxed clock, birth-death, 3 Normal secondary calibrations,
monophyly constraints, constraint-compatible UPGMA starting tree). Usage: 03_make_xml.py out.xml chain_length log_every prefix"""
import sys, csv, itertools
out, chain, every, prefix = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
al = {}; n = None
for l in open("results/1_alignment/COI_aligned.fasta"):
    l = l.strip()
    if l.startswith(">"): n = l[1:]; al[n] = ""
    else: al[n] += l
ing = [k for k in al if k not in ("P_memnon", "P_machaon")]
assert len(ing) == 40 and len(al) == 42
# ---- constraint-compatible starting tree (UPGMA on p-distances, ingroup), heights in Ma
def pd(a, b):
    d = t = 0
    for x, y in zip(al[a], al[b]):
        if x in "ACGT" and y in "ACGT": t += 1; d += x != y
    return d / t
class Nd:
    def __init__(s, name=None, kids=(), h=0.0): s.name, s.kids, s.h = name, list(kids), h
cl = {frozenset([k]): Nd(k) for k in ing}
D = {}
keys = list(cl)
def dist(a, b):
    return sum(pd(x, y) for x in a for y in b) / (len(a) * len(b))
while len(cl) > 1:
    ks = list(cl); best = None
    for a, b in itertools.combinations(ks, 2):
        d = dist(a, b)
        if best is None or d < best[0]: best = (d, a, b)
    d, a, b = best
    cl[a | b] = Nd(kids=[cl.pop(a), cl.pop(b)], h=d / 2)
    # (cl.pop above removes a,b)
root = next(iter(cl.values()))
sc = 17.0 / root.h
def fix(nd):
    if not nd.kids: nd.h = 0.0; return 0.0
    nd.h *= sc
    m = max(fix(k) for k in nd.kids)
    if nd.h < m + 0.05: nd.h = m + 0.05
    return nd.h
fix(root)
def nw(nd, ph):
    if not nd.kids: return f"{nd.name}:{ph - nd.h:.5f}"
    return "(" + ",".join(nw(k, nd.h) for k in nd.kids) + f"):{ph - nd.h:.5f}"
ingroup_h = root.h
mem_h = max(19.41, ingroup_h + 0.5); mach_h = max(21.25, mem_h + 0.5)
s_ing = nw(root, mem_h).rstrip(";")
s_mem = f"({s_ing},P_memnon:{mem_h:.5f}):{mach_h - mem_h:.5f}"
start = f"({s_mem},P_machaon:{mach_h:.5f});"
# ---- XML
N = len(al); X = []
w = X.append
w('<?xml version="1.0" standalone="yes"?>\n<beast>')
w('  <taxa id="taxa">'); [w(f'    <taxon id="{k}"/>') for k in al]; w('  </taxa>')
w('  <taxa id="Achillides">'); [w(f'    <taxon idref="{k}"/>') for k in ing]; w('  </taxa>')
w('  <taxa id="Achillides_memnon">'); [w(f'    <taxon idref="{k}"/>') for k in ing + ["P_memnon"]]; w('  </taxa>')
w('  <alignment id="alignment" dataType="nucleotide">')
for k, v in al.items(): w(f'    <sequence><taxon idref="{k}"/>{v.replace("-","?")}</sequence>')
w('  </alignment>')
w('  <patterns id="patterns" from="1" strip="false"><alignment idref="alignment"/></patterns>')
w(f'  <newick id="startingTree" usingDates="false" units="years">{start}</newick>')
w('''  <treeModel id="treeModel">
    <newick idref="startingTree"/>
    <rootHeight><parameter id="treeModel.rootHeight"/></rootHeight>
    <nodeHeights internalNodes="true"><parameter id="treeModel.internalNodeHeights"/></nodeHeights>
    <nodeHeights internalNodes="true" rootNode="true"><parameter id="treeModel.allInternalNodeHeights"/></nodeHeights>
  </treeModel>
  <birthDeathModel id="birthDeath" units="years">
    <birthMinusDeathRate><parameter id="birthDeath.meanGrowthRate" value="0.2" lower="0.0" upper="Infinity"/></birthMinusDeathRate>
    <relativeDeathRate><parameter id="birthDeath.relativeDeathRate" value="0.5" lower="0.0" upper="1.0"/></relativeDeathRate>
  </birthDeathModel>
  <speciationLikelihood id="speciation">
    <model><birthDeathModel idref="birthDeath"/></model>
    <speciesTree><treeModel idref="treeModel"/></speciesTree>
  </speciationLikelihood>''')
w('  <monophylyStatistic id="monophyly(Achillides)"><mrca><taxa idref="Achillides"/></mrca><treeModel idref="treeModel"/></monophylyStatistic>')
w('  <monophylyStatistic id="monophyly(Achillides_memnon)"><mrca><taxa idref="Achillides_memnon"/></mrca><treeModel idref="treeModel"/></monophylyStatistic>')
w('  <tmrcaStatistic id="tmrca(Achillides)" includeStem="false"><mrca><taxa idref="Achillides"/></mrca><treeModel idref="treeModel"/></tmrcaStatistic>')
w('  <tmrcaStatistic id="tmrca(Achillides_memnon)" includeStem="false"><mrca><taxa idref="Achillides_memnon"/></mrca><treeModel idref="treeModel"/></tmrcaStatistic>')
w(f'''  <discretizedBranchRates id="branchRates">
    <treeModel idref="treeModel"/>
    <distribution><logNormalDistributionModel meanInRealSpace="true">
      <mean><parameter id="ucld.mean" value="0.01" lower="0.0" upper="Infinity"/></mean>
      <stdev><parameter id="ucld.stdev" value="0.3" lower="0.0" upper="Infinity"/></stdev>
    </logNormalDistributionModel></distribution>
    <rateCategories><parameter id="branchRates.categories" dimension="{2*N-2}"/></rateCategories>
  </discretizedBranchRates>
  <rateStatistic id="meanRate" name="meanRate" mode="mean" internal="true" external="true"><treeModel idref="treeModel"/><discretizedBranchRates idref="branchRates"/></rateStatistic>
  <rateStatistic id="coefficientOfVariation" name="coefficientOfVariation" mode="coefficientOfVariation" internal="true" external="true"><treeModel idref="treeModel"/><discretizedBranchRates idref="branchRates"/></rateStatistic>
  <gtrModel id="gtr">
    <frequencies><frequencyModel dataType="nucleotide"><frequencies><parameter id="frequencies" value="0.25 0.25 0.25 0.25"/></frequencies></frequencyModel></frequencies>
    <rateAC><parameter id="ac" value="1.0" lower="0.0"/></rateAC>
    <rateAG><parameter id="ag" value="1.0" lower="0.0"/></rateAG>
    <rateAT><parameter id="at" value="1.0" lower="0.0"/></rateAT>
    <rateCG><parameter id="cg" value="1.0" lower="0.0"/></rateCG>
    <rateGT><parameter id="gt" value="1.0" lower="0.0"/></rateGT>
  </gtrModel>
  <siteModel id="siteModel">
    <substitutionModel><gtrModel idref="gtr"/></substitutionModel>
    <gammaShape gammaCategories="4"><parameter id="alpha" value="0.5" lower="0.0"/></gammaShape>
  </siteModel>
  <treeDataLikelihood id="treeLikelihood" useAmbiguities="false">
    <partition><patterns idref="patterns"/><siteModel idref="siteModel"/></partition>
    <treeModel idref="treeModel"/>
    <discretizedBranchRates idref="branchRates"/>
  </treeDataLikelihood>
  <operators id="operators" optimizationSchedule="default">
    <scaleOperator scaleFactor="0.5" weight="4"><parameter idref="ac"/></scaleOperator>
    <scaleOperator scaleFactor="0.5" weight="4"><parameter idref="ag"/></scaleOperator>
    <scaleOperator scaleFactor="0.5" weight="4"><parameter idref="at"/></scaleOperator>
    <scaleOperator scaleFactor="0.5" weight="4"><parameter idref="cg"/></scaleOperator>
    <scaleOperator scaleFactor="0.5" weight="4"><parameter idref="gt"/></scaleOperator>
    <deltaExchange delta="0.02" weight="3"><parameter idref="frequencies"/></deltaExchange>
    <scaleOperator scaleFactor="0.75" weight="3"><parameter idref="alpha"/></scaleOperator>
    <scaleOperator scaleFactor="0.75" weight="3"><parameter idref="ucld.mean"/></scaleOperator>
    <scaleOperator scaleFactor="0.75" weight="3"><parameter idref="ucld.stdev"/></scaleOperator>
    <scaleOperator scaleFactor="0.75" weight="3"><parameter idref="birthDeath.meanGrowthRate"/></scaleOperator>
    <randomWalkOperator windowSize="0.1" weight="3" autoOptimize="false"><parameter idref="birthDeath.relativeDeathRate"/></randomWalkOperator>
    <subtreeSlide size="2.0" gaussian="true" weight="30"><treeModel idref="treeModel"/></subtreeSlide>
    <narrowExchange weight="30"><treeModel idref="treeModel"/></narrowExchange>
    <wideExchange weight="10"><treeModel idref="treeModel"/></wideExchange>
    <wilsonBalding weight="10"><treeModel idref="treeModel"/></wilsonBalding>
    <scaleOperator scaleFactor="0.75" weight="3"><parameter idref="treeModel.rootHeight"/></scaleOperator>
    <uniformOperator weight="30"><parameter idref="treeModel.internalNodeHeights"/></uniformOperator>
    <upDownOperator scaleFactor="0.75" weight="3">
      <up><parameter idref="ucld.mean"/></up>
      <down><parameter idref="treeModel.allInternalNodeHeights"/></down>
    </upDownOperator>
    <swapOperator size="1" weight="10" autoOptimize="false"><parameter idref="branchRates.categories"/></swapOperator>
    <uniformIntegerOperator weight="10"><parameter idref="branchRates.categories"/></uniformIntegerOperator>
  </operators>
  <mcmc id="mcmc" chainLength="{chain}" autoOptimize="true" operatorAnalysis="{prefix}.ops">
    <posterior id="posterior">
      <prior id="prior">
        <booleanLikelihood><monophylyStatistic idref="monophyly(Achillides)"/></booleanLikelihood>
        <booleanLikelihood><monophylyStatistic idref="monophyly(Achillides_memnon)"/></booleanLikelihood>
        <normalPrior mean="17.0" stdev="2.1"><statistic idref="tmrca(Achillides)"/></normalPrior>
        <normalPrior mean="19.41" stdev="2.25"><statistic idref="tmrca(Achillides_memnon)"/></normalPrior>
        <normalPrior mean="21.25" stdev="2.45"><parameter idref="treeModel.rootHeight"/></normalPrior>
        <uniformPrior lower="0.0" upper="1.0"><parameter idref="ucld.mean"/></uniformPrior>
        <exponentialPrior mean="0.3333333333333333" offset="0.0"><parameter idref="ucld.stdev"/></exponentialPrior>
        <uniformPrior lower="0.0" upper="10.0"><parameter idref="birthDeath.meanGrowthRate"/></uniformPrior>
        <uniformPrior lower="0.0" upper="1.0"><parameter idref="birthDeath.relativeDeathRate"/></uniformPrior>
        <exponentialPrior mean="1.0" offset="0.0"><parameter idref="ac"/></exponentialPrior>
        <exponentialPrior mean="1.0" offset="0.0"><parameter idref="ag"/></exponentialPrior>
        <exponentialPrior mean="1.0" offset="0.0"><parameter idref="at"/></exponentialPrior>
        <exponentialPrior mean="1.0" offset="0.0"><parameter idref="cg"/></exponentialPrior>
        <exponentialPrior mean="1.0" offset="0.0"><parameter idref="gt"/></exponentialPrior>
        <exponentialPrior mean="0.5" offset="0.0"><parameter idref="alpha"/></exponentialPrior>
        <uniformPrior lower="0.0" upper="1.0"><parameter idref="frequencies"/></uniformPrior>
        <speciationLikelihood idref="speciation"/>
      </prior>
      <likelihood id="likelihood"><treeDataLikelihood idref="treeLikelihood"/></likelihood>
    </posterior>
    <operators idref="operators"/>
    <log id="screenLog" logEvery="{max(every*10,1000)}">
      <column label="Posterior" dp="4" width="12"><posterior idref="posterior"/></column>
      <column label="Prior" dp="4" width="12"><prior idref="prior"/></column>
      <column label="Likelihood" dp="4" width="12"><likelihood idref="likelihood"/></column>
      <column label="rootHeight" sf="6" width="12"><parameter idref="treeModel.rootHeight"/></column>
      <column label="ucld.mean" sf="6" width="12"><parameter idref="ucld.mean"/></column>
    </log>
    <log id="fileLog" logEvery="{every}" fileName="{prefix}.log" overwrite="false">
      <posterior idref="posterior"/><prior idref="prior"/><likelihood idref="likelihood"/>
      <parameter idref="treeModel.rootHeight"/>
      <tmrcaStatistic idref="tmrca(Achillides)"/><tmrcaStatistic idref="tmrca(Achillides_memnon)"/>
      <parameter idref="birthDeath.meanGrowthRate"/><parameter idref="birthDeath.relativeDeathRate"/>
      <parameter idref="ac"/><parameter idref="ag"/><parameter idref="at"/><parameter idref="cg"/><parameter idref="gt"/>
      <parameter idref="frequencies"/><parameter idref="alpha"/>
      <parameter idref="ucld.mean"/><parameter idref="ucld.stdev"/>
      <rateStatistic idref="meanRate"/><rateStatistic idref="coefficientOfVariation"/>
      <speciationLikelihood idref="speciation"/><treeDataLikelihood idref="treeLikelihood"/>
    </log>
    <logTree id="treeFileLog" logEvery="{every}" nexusFormat="true" fileName="{prefix}.trees" sortTranslationTable="true">
      <treeModel idref="treeModel"/>
      <trait name="rate" tag="rate"><discretizedBranchRates idref="branchRates"/></trait>
      <posterior idref="posterior"/>
    </logTree>
  </mcmc>
  <report><property name="timer"><mcmc idref="mcmc"/></property></report>
</beast>''')
open(out, "w").write("\n".join(X) + "\n")
print("wrote", out, "| start tree crown", round(ingroup_h, 2))
