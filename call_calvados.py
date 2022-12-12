# coding: utf-8
import sys
import calvados

# 1 synuclein
model1 = calvados.CalvadosModel()
model1.add_proteins("aS_WT")
system1 = calvados.OMMsystem(model1)
#ommrunner1 = calvados.OMMrunner(system1)

## 20 synuclein
system1_1 = calvados.OMMsystem(model1, n_chains=20)
#ommrunner1_1 = calvados.OMMrunner(system1_1)

## 1 synuclein +  1 pLK
model2 = calvados.CalvadosModel()
model2.add_proteins(["aS_WT", "pLK"])
system2 = calvados.OMMsystem(model2)
#ommrunner2 = calvados.OMMrunner(system2)

# 30 synuclein +  2 pLK
system2_2 = calvados.OMMsystem(model2, n_chains=[30,2])
#ommrunner2_2 = calvados.OMMrunner(system2_2)
