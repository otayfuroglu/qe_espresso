from ase.calculators.espresso import Espresso, EspressoProfile

from ase.io.extxyz import read_extxyz
from ase.io import read, write
from ase import Atoms

from pathlib import Path
import os, sys
import argparse



def getQECalculator():
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
                                'ecutwfc': 100,
                                'ecutrho': 400,
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


def runPh(name):
    input_ph = f"""
    &INPUTPH
      outdir = tmpdir
      prefix = '{name}'
      tr2_ph = 1d-14
      ldisp = .true.
      recover = .true.
      nq1 = 3
      nq2 = 3
      nq3 = 3
      fildyn = '{name}.dyn'
    /
    """

    fl = open(f"ph.{name}.in", "w")
    print(input_ph, file=fl)
    fl.close()
    os.system(f"mpirun -n {nprocs} ph.x -i ph.{name}.in > ph.{name}.out")

def runQ2r(name):
    input_q2r = f"""
    &INPUT
      fildyn = '{name}.dyn'
      zasr = 'crystal'
      flfrc = '{name}.fc'
    /
    """

    fl = open(f"q2r.{name}.in", "w")
    print(input_q2r, file=fl)
    fl.close()
    os.system(f"mpirun -n {nprocs}  q2r.x -i q2r.{name}.in > q2r.{name}.out")

def runMetadyn(name):
    input_metadyn = f"""
    &INPUT
      asr = 'crystal'
      flfrc = '{name}.fc'
      flfrq = '{name}.freq'
      flvec = '{name}.modes'
    !  loto_2d = .true.
      q_in_band_form = .true.
    /
    5
    0.500 0.500 0.500   20 ! L
    0.000 0.000 0.000   20 ! G
    0.500 0.000 0.500   20 ! X
    0.375 0.375 0.750   20 ! K
    0.000 0.000 0.000    1 ! G
    """

    fl = open(f"metadyn.{name}.in", "w")
    print(input_metadyn, file=fl)
    fl.close()
    os.system(f"mpirun -n {nprocs} metadyn.x -i metadyn.{name}.in > metadyn.{name}.out")


parser = argparse.ArgumentParser(description="Give something ...")
parser.add_argument("-calc_type", type=str, required=True)
parser.add_argument("-geoms_path", type=str, required=True)
parser.add_argument("-func", type=str, required=True)
parser.add_argument("-idx", type=int, required=True)
parser.add_argument("-nprocs", type=int, required=True)
parser.add_argument("-tmpdir", type=str, required=True)
args = parser.parse_args()


calc_type = args.calc_type
geoms_path = args.geoms_path
func = args.func.upper()
idx = args.idx
nprocs = args.nprocs
tmpdir = args.tmpdir
atoms = read(geoms_path, index=idx)


if "isolated" in geoms_path:
    keyword = "isolated"
elif "polymeric" in geoms_path:
    keyword = "polymeric"

try:
    name = atoms.info["label"]
except:
    name = f"structure_{idx}"

if calc_type == "sp":
    calculation = 'scf'
elif calc_type == "opt":
    calculation = 'vc-relax'
elif calc_type == "freq":
    calculation = 'scf'
elif calc_type == "opt_freq":
    calculation = 'vc-relax'

calc = getQECalculator()
OUT_DIR = f"{keyword}_{calc_type}_{func}/{name}"

if not os.path.exists(OUT_DIR):

    # create outpur folders
    # easy create nested forder wether it is exists
    Path(OUT_DIR).mkdir(parents=True, exist_ok=True)
    # change to local scratch directory
    #  os.chdir(TMP_DIR)

    cwd = os.getcwd()
    os.chdir(OUT_DIR)

    atoms.calc = calc
    atoms.get_potential_energy()
    if calc_type == "opt" or calc_type == "opt_freq":
        print("Readiding optimized structure to start freq calculations")
        atoms = read("espresso.pwo")
        atoms.info["label"] = name

    if calc_type == "freq" or calc_type == "opt_freq":
        runPh(name)
        #  runQ2r(name)
        #  runMetadyn(name)

    os.chdir(cwd)
    if calc_type == "opt" or calc_type == "opt_freq":
        write(f"qe_{args.calc_type}_{keyword}.extxyz", atoms, append=True)
else:
    print(f"{idx} calcultaion is alrady done")

