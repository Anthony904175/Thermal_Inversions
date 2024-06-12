#Import modules
#import sys
#sys.path.insert(1, "/Users/aosborne3/Desktop/research_tutorials/n2v")
import n2v
import psi4
import matplotlib.pyplot as plt
import numpy as np

import matplotlib as mpl
mpl.rcParams["font.size"] = 11
mpl.rcParams["font.family"] = "sans-serif"
mpl.rcParams["axes.edgecolor"] = "#eae8e9" 
#Set options within Psi4.
#Saving the jk object prevents n2v to build it again. 
psi4.set_options({"save_jk" : True})

#n2v is driven by psi4's reference option. Make sure you set it accordingly. 
psi4.set_options({"reference" : "rhf"})

#Set memory and avoid inheriting options from a different calculation. 
psi4.set_memory(int(2.50e9))
psi4.core.clean()
# Let us focus on the Neon atom.

# Define Psi4 geometries. Symmetries need to be set to C1!
Ne = psi4.geometry(
"""
0 1
Ne 0.0 0.0 0.0
noreorient
nocom
units bohr
symmetry c1
""")

# Perform a DFT calculation.
e, wfn = psi4.energy("svwn/cc-pvtz", return_wfn=True, molecule=Ne)


# Define inverter objects for each molecule. Simply use the wnf object from psi4 as an argument.
ine = n2v.Inverter.from_wfn(wfn)

# Define a grid to visualize all of the generated potentials.
# Since we are using radially symmetric atoms, we focus on the space from the origin all the way to 10 a.u.
npoints=10002
x = np.linspace(0,10,npoints)[:,None]
y = np.zeros_like(x)
z = y
grid = np.concatenate((x,y,z), axis=1).T
#let us invert the kohn-sham equations:
vxc_inverted = ine.invert('direct', grid=grid, correction=False)
#Compare the performance of direct and reversed calculation of vxc. 
vxc0 = ine.eng.grid.vxc(Da=wfn.Da().np, Db=wfn.Db().np, grid=grid)
print(vxc0)
#---------------> Plot
#fig, ax = plt.subplots(1,1, dpi=100)


#---------------> Data
#ax.plot(x, vxc0, label="Forward DFT")
#ax.plot(x, vxc_inverted, label="Inverted DFT")

#---------------> Formatting
#ax.set_xlim(0.001,10)
#ax.set_ylim(-10,1)
#ax.set_ylabel('Density (e/$b_0^3$)')
#ax.set_xscale('log')
#ax.set_xlabel('z ($b_0$)')
#ax.legend()

#fig.suptitle("Neon Atom Vxc")
