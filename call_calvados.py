# coding: utf-8
import sys
import calvados

# 1 synuclein
model1 = calvados.CalvadosModel()
model1.add_proteins("aS_WT", file_fasta="sequences.fasta")
system1 = calvados.OMMsystem(model1, box=100)
ommrunner1 = calvados.OMMrunner(system1)
#ommrunner1.run(time=0.01)
ommrunner1.run(steps=10000)

# 20 synuclein in cubic box
system_multi = calvados.OMMsystem(model1, n_chains=20, box=200, name='aS_WT_20_L200')
ommrunner_multi = calvados.OMMrunner(system_multi)
#ommrunner_multi.run(steps=10000)

# 20 synuclein in slab configuration
system_slab = calvados.OMMsystem(model1, n_chains=20, box=[10, 10, 200], name='aS_WT_20_slab')
ommrunner_slab = calvados.OMMrunner(system_slab)
ommrunner_slab.run(steps=100000)

## 1 synuclein +  1 pLK
#model2 = calvados.CalvadosModel()
#model2.add_proteins(["aS_WT", "pLK"])
#system2 = calvados.OMMsystem(model2)
#ommrunner2 = calvados.OMMrunner(system2)

# 30 synuclein +  2 pLK
#system2_2 = calvados.OMMsystem(model2, n_chains=[30,2])
#ommrunner2_2 = calvados.OMMrunner(system2_2)
