# coding: utf-8
import sys
sys.path.append("../")
import calvados

# 1 synuclein
model1 = calvados.CalvadosModel(file_res='../residues.csv')
model1.add_proteins(["aS_WT", "Gly10"], file_fasta='../sequences.fasta')
system1 = calvados.OMMsystem(model1, n_chains=[100, 1000], box=50)
ommrunner1 = calvados.OMMrunner(system1)

#
### 20 synuclein
#system1_1 = calvados.OMMsystem(model1, n_chains=100, box=[200, 200, 200], name='aS_WT_100_L200')
#ommrunner1_1 = calvados.OMMrunner(system1_1)
#ommrunner1_1.run(0.1)
#
### 1 synuclein +  1 pLK
#model2 = calvados.CalvadosModel()
#model2.add_proteins(["aS_WT", "pLK"])
#system2 = calvados.OMMsystem(model2)
##ommrunner2 = calvados.OMMrunner(system2)
#
## 30 synuclein +  2 pLK
#system2_2 = calvados.OMMsystem(model2, n_chains=[30,2])
##ommrunner2_2 = calvados.OMMrunner(system2_2)
