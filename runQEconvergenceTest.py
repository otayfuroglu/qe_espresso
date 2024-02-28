from ase.calculators.espresso import Espresso, EspressoProfile

from ase.io.extxyz import read_extxyz
from ase.io import read, write
from ase import Atoms

from pathlib import Path
import os, sys
import argparse



def getQECalculator(name, calculation, tmpdir, func, ecutwfc, ecutrho):
    pseudopotentials = {"H": "H.pbe-rrkjus_psl.1.0.0.UPF",
                        "Li": "Li.pbe-s-kjpaw_psl.1.0.0.UPF",
                        "Al": "Al.pbe-n-kjpaw_psl.1.0.0.UPF"
                       }
    # set env
    espresso_profile = EspressoProfile(binary=f'mpirun -n {nprocs} pw.x'.split(),
                                       pseudo_dir=Path.home()/Path("qe_pseudopotentials/PBE_efficiency/"))


    input_data = {"control":   {'prefix': name,
                                'calculation': calculation,
                                "etot_conv_thr": 1.0e-7,
                                "forc_conv_thr": 5.0e-4,
                                "outdir": tmpdir,
                                #  "outdir": "/tmp/" + name,
                                #  "verbosity": 'high',
                               },
                  "system":    {'ibrav': 0,
                                'nosym': True,
                                "input_dft": func,
                                'ecutwfc': ecutwfc,
                                'ecutrho': ecutrho,
                               },
                  "electrons": {'conv_thr': 1.0e-8,
                                'electron_maxstep': 200,
                                'mixing_beta': 0.7,
                               }
                 }

    calc = Espresso(
        profile=espresso_profile,
        #  ecutwfc=25,
        kpts = (6, 6, 6),
        pseudopotentials=pseudopotentials,
        input_data = input_data,
        tstress=True, tprnfor=True, nosym=True,)
    return calc


parser = argparse.ArgumentParser(description="Give something ...")
parser.add_argument("-geoms_path", type=str, required=True)
parser.add_argument("-func", type=str, required=True)
parser.add_argument("-nprocs", type=int, required=True)
parser.add_argument("-tmpdir", type=str, required=True)
args = parser.parse_args()


geoms_path = args.geoms_path
func = args.func.upper()
nprocs = args.nprocs
tmpdir = args.tmpdir
atoms = read(geoms_path)


if "isolated" in geoms_path:
    keyword = "isolated"
elif "polymeric" in geoms_path:
    keyword = "polymeric"

try:
    name = atoms.info["label"]
except:
    name = f"structure_{idx}"

calculation = 'scf'

OUT_DIR = f"{keyword}_{func}/{name}"
Path(OUT_DIR).mkdir(parents=True, exist_ok=True)
os.chdir(OUT_DIR)
#  cwd = os.getcwd()

ecutwfc = 120
testpar_list = range(400, 801, 100)
testpar_name = "ecutrho"
with open(f"{testpar_name}_convergence_test.csv", "w") as fl:
    print(f"{testpar_name},Energy", file=fl)
    for testpar in testpar_list:
        calc = getQECalculator(name=name,
                               calculation=calculation,
                               tmpdir=tmpdir,
                               func=func,
                               ecutwfc=120,
                               ecutrho=testpar)
        atoms.calc = calc
        pot_e = atoms.get_potential_energy()
        print(f"{testpar},{pot_e}", file=fl)
        fl.flush()
