#
from ase.io.espresso import read_espresso_ph
from ase.io.extxyz import read_extxyz
from ase.io import read, write
from ase import Atoms

from pathlib import Path
import os
import argparse
import pandas as pd


import numpy as np
from ase.thermochemistry import HarmonicThermo
from ase.vibrations import Vibrations
from ase import units


#  units = units.create_units("2014")
#  print(units)
THZ2EV = 0.004135665538



parser = argparse.ArgumentParser(description="Give something ...")
parser.add_argument("-extxyz_path", type=str, required=True, help="..")
args = parser.parse_args()
extxyz_path = args.extxyz_path


if "isolated" in extxyz_path:
    keyword = "isolated"
elif "polymeric" in extxyz_path:
    keyword = "polymeric"
else:
    keyword = "test"


keyword += "_24atoms_3x3x3"

atoms_list = read(extxyz_path, index=":")

df_helmholtz = pd.DataFrame()
temperature_list = range(0, 1001, 10)
df_helmholtz["Temperature"] = temperature_list
for i, atoms in enumerate(atoms_list):
    #  name = atoms.info["label"]
    name = "structure_0"
    potentialenergy = atoms.get_potential_energy()
    BASE_DIR = Path(f"{keyword}")
    print(name)

    #  try:
    ph_path = f"{BASE_DIR}/{name}/ph.{name}.out"
    phout = read_espresso_ph(open(ph_path))
    vib_energies = []
    for i in phout.keys():
        vib_energies += [phout[i]["freqs"] * THZ2EV]
    #  for i in range(1, 15):
    #      if i == 1:
    #          vib_energies = np.array(phout[i]["freqs"])
    #      else:
    #          vib_energies += np.array(phout[i]["freqs"])

    vib_energies = np.concatenate(vib_energies)
        #  break
    #  vib_energies *= THZ2EV
    #  vib_energies = phout[3]["freqs"]
    #  print(vib_energies * THZ2EV)
    #  quit()
        #  print(vib_energies)
    #  except:
        #  continue

    #  quit()

    #  vib_energies = [mode for mode in np.real(vib_energies) if mode > 1e-3]
    #  vib_energies = [complex(1.0e-8, 0) if energy < 1.0e-4 else energy for energy in vib_energies]

    #  print(vib_energies)
    # get free energy from ase
    #  free_energy_class = HarmonicThermo(vib_energies, potentialenergy=0.)#potentialenergy)
    free_energy_class = HarmonicThermo(vib_energies, potentialenergy=potentialenergy)
    #  df_gibbs["Temperature"] = temperature_list
    helm_holtz_energies = []
    for temperature in temperature_list:
        #  print(temperature)
        if temperature == 0:
            temperature = 1e-5
            #  temperature = 1
        helm_holtz_energy = free_energy_class.get_helmholtz_energy(temperature, verbose=True)
        helm_holtz_energies += [helm_holtz_energy / len(atoms)]
    df_helmholtz[f"structure_{i}"] = helm_holtz_energies
    #  df_gibbs[f"structure_{i}"] = gibbs

#  df_helmholtz.to_csv(f"{BASE_DIR}/helmholtz_{keyword}.csv")
df_helmholtz.to_csv(f"helmholtz_{keyword}.csv")
#  df_gibbs.to_csv(f"gibbs_{keyword}.csv")

