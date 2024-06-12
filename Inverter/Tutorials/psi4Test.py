import psi4

# Set options within Psi4.
psi4.set_options({"save_jk" : True})
psi4.set_options({"reference" : "rhf"})
psi4.set_memory(int(2.50e9))
psi4.core.clean()

# Define the molecule (Neon atom in this case).
neon = psi4.geometry("""
0 1
Ne 0.0 0.0 0.0
noreorient
nocom
units bohr
symmetry c1
""")

# Perform a Hartree-Fock (HF) calculation.
e_hf = psi4.energy("hf/cc-pvtz", molecule=neon)

print("Hartree-Fock Energy:", e_hf)
