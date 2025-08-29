import os
import yaml
import shutil
import argparse
import logging
import subprocess
import mdtraj as md
import re


def rename_atoms_mdtraj(pdb_in, pdb_out):
    traj = md.load(pdb_in)
    topology = traj.topology
    element_counts = {}
    atom_renames_list = []

    for atom in topology.atoms:
        element = atom.element.symbol if atom.element else atom.name[0]
        element_counts[element] = element_counts.get(element, 0) + 1
        new_name = f"{element}{element_counts[element]}"
        atom.name = new_name
        atom_renames_list.append(new_name)

    traj.save_pdb(pdb_out)
    return atom_renames_list


def rename_itp(itp_in, itp_out, atom_renames_list):
    new_lines = []
    in_atoms_section = False
    atom_index = 0

    for line in open(itp_in):
        stripped = line.strip()

        if stripped.startswith("[ atoms ]"):
            in_atoms_section = True
            new_lines.append(line)
            continue
        elif stripped.startswith("[") and "]" in stripped and not stripped.startswith("[ atoms ]"):
            in_atoms_section = False
            new_lines.append(line)
            continue

        if in_atoms_section and stripped and not stripped.startswith(";"):
            parts = line.split()
            if atom_index < len(atom_renames_list):
                parts[4] = atom_renames_list[atom_index]
            else:
                raise ValueError(f"Too few renamed atoms provided for {itp_in}")
            atom_index += 1

            # Ensure 8 parts for formatting
            padded = parts + [''] * (8 - len(parts))
            new_lines.append("{:<6} {:<9} {:>2} {:<4} {:<4} {:>5} {:>15} {:>15}\n".format(*padded))
        else:
            new_lines.append(line)

    with open(itp_out, 'w') as f:
        f.writelines(new_lines)


def rename_molecule_files(pdb_path, itp_path):
    pdb_out = pdb_path.replace(".pdb", "_renamed.pdb")
    itp_out = itp_path.replace(".itp", "_renamed.itp") if itp_path else None

    #print(f"Renaming atoms for PDB: {pdb_path} -> {pdb_out}")
    atom_renames_list = rename_atoms_mdtraj(pdb_path, pdb_out)

    if itp_out and os.path.exists(itp_path):
        #print(f"Renaming atoms for ITP: {itp_path} -> {itp_out}")
        rename_itp(itp_path, itp_out, atom_renames_list)
    else:
        itp_out = itp_path  # fallback: no rename

    return pdb_out, itp_out


# The rest of the script remains unchanged, but included here for completeness:

def load_config(path: str) -> dict:
    with open(path, 'r') as f:
        config = yaml.safe_load(f)

    config_dir = os.path.abspath(os.path.dirname(path))

    for key, value in config.items():
        if isinstance(value, str) and not os.path.isabs(value):
            if any(value.endswith(ext) for ext in [".pdb", ".itp", ".mdp", ".top", ".gro", ".yaml", ".yml", ".txt", ".log"]):
                config[key] = os.path.abspath(os.path.join(config_dir, value))

    return config


def run_command(command):
    result = subprocess.run(command, shell=True, capture_output=True, text=True)
    logging.info(result.stdout)
    if result.stderr:
        logging.error(result.stderr)
    if result.returncode != 0:
        raise RuntimeError(f"Command failed: {command}")


def write_packmol_input(path, sim_name, n_solute, n_solvent, config):
    num_solute = n_solute - 1
    box_size = config["box_size"]

    with open(path, "w") as f:
        f.write("tolerance 2.0\nfiletype pdb\n")
        f.write(f"output {sim_name}_mixture.pdb\n\n")

        # Alchemical ligand
        f.write(f"""structure {os.path.basename(config["alchemical_ligand_pdb"])}
  number 1
  inside box 0. 0. 0. {box_size}. {box_size}. {box_size}.
end structure\n\n""")

        if num_solute > 0:
            f.write(f"""structure {os.path.basename(config["solute_pdb"])}
  number {num_solute}
  inside box 0. 0. 0. {box_size}. {box_size}. {box_size}.
end structure\n\n""")

        if n_solvent > 0:
            f.write(f"""structure {os.path.basename(config["solvent_pdb"])}
  number {n_solvent}
  inside box 0. 0. 0. {box_size}. {box_size}. {box_size}.
end structure\n""")


def write_topology(path, config, sim_name, alchemical_ligand_name, solute_name, n_solute, solvent_name, n_solvent):
    num_solute = n_solute - 1
    with open(path, "w") as f:
        f.write("[ defaults ]\n1 2 yes 0.5 0.833333\n\n")
        f.write(f'#include "{config["all_atomtypes_itp"]}"\n')
        f.write(f'#include "{config["alchemical_ligand_itp"]}"\n')
        f.write(f'#include "{config["solute_itp"]}"\n')
        f.write(f'#include "{config["solvent_itp"]}"\n')

        f.write("\n[ system ]\nPackmol system\n\n[ molecules ]\n")
        f.write(f"{alchemical_ligand_name} 1\n")
        if num_solute > 0:
            f.write(f"{solute_name} {num_solute}\n")
        if n_solvent > 0:
            f.write(f"{solvent_name} {n_solvent}\n")


def parse_arguments():
    parser = argparse.ArgumentParser(description="CFMC simulation setup")
    parser.add_argument("--config", "-y", default="config.yaml", help="YAML config file")
    parser.add_argument("--parent_dir", "-p", default="sim_data")
    parser.add_argument("--output_cmd_directory", "-c", default="cmds")
    parser.add_argument("--log_file", "-l", default="setup.log")
    parser.add_argument("--skip", "-s", default="",
                        help="Comma-separated roles to skip renaming: e.g. 'solvent,solute'. Default: none")
    return parser.parse_args()


def main():
    args = parse_arguments()
    config = load_config(args.config)

    skip_list = [x.strip().lower() for x in args.skip.split(",") if x.strip()]
    print(f"Skipping renaming for roles: {skip_list}")

    logging.basicConfig(filename=args.log_file, level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

    base_dir = os.getcwd()
    os.makedirs(args.parent_dir, exist_ok=True)
    os.makedirs(args.output_cmd_directory, exist_ok=True)

    mdp_files = [
        config["em_mdp"],
        config["nvt_mdp"],
        config["npt_mdp"],
        config["prod_mdp"],
    ]
    itp_files = [
        config["all_atomtypes_itp"],
        config["alchemical_ligand_itp"],
        config["solute_itp"],
        config["solvent_itp"],
    ]

    for role in ["alchemical_ligand", "solute", "solvent"]:
        if role in skip_list:
            print(f"Skipping renaming for role '{role}' as per user request")
            continue

        pdb_path = config.get(f"{role}_pdb")
        itp_path = config.get(f"{role}_itp", None)

        if pdb_path and os.path.exists(pdb_path):
            pdb_renamed, itp_renamed = rename_molecule_files(pdb_path, itp_path)
            config[f"{role}_pdb"] = pdb_renamed
            if itp_renamed:
                config[f"{role}_itp"] = itp_renamed
        else:
            print(f"Warning: {role} pdb file {pdb_path} not found or not specified.")

    for num_solute in range(1, config["n_solute"] + 1):
        condition_dir = os.path.join(args.parent_dir, f"{config['solute_name']}_{num_solute}")
        os.makedirs(condition_dir, exist_ok=True)

        for rep in range(1, config["n_replicates"] + 1):
            sim_name = f"rep{rep}"
            sim_dir = os.path.abspath(os.path.join(condition_dir, sim_name))
            os.makedirs(sim_dir, exist_ok=True)

            logging.info(f"Setting up {sim_name} in {condition_dir}: {num_solute} {config['solute_name']}, 1 transforming, {config['n_solvent']} {config['solvent_name']}")

            inp_file = os.path.join(sim_dir, "packmol.inp")
            write_packmol_input(inp_file, sim_name, num_solute, config['n_solvent'], config)

            shutil.copy(config["alchemical_ligand_pdb"], sim_dir)
            shutil.copy(config["solute_pdb"], sim_dir)
            shutil.copy(config["solvent_pdb"], sim_dir)

            for mdp in mdp_files:
                shutil.copy(mdp, sim_dir)

            os.chdir(sim_dir)
            run_command("packmol < packmol.inp")
            run_command(f"gmx editconf -f {sim_name}_mixture.pdb -o {sim_name}_final.gro -box {config['box_size']} {config['box_size']} {config['box_size']}")
            write_topology(os.path.join(sim_dir, f"{sim_name}.top"),
                           config,
                           sim_name,
                           config["alchemical_ligand_name"],
                           config["solute_name"],
                           num_solute,
                           config["solvent_name"],
                           config["n_solvent"])
            os.chdir(base_dir)

            script_filename = f"{sim_name}_{config['solute_name']}{num_solute}.sh"
            script_path = os.path.join(base_dir, args.output_cmd_directory, script_filename)
            with open(script_path, "w") as script:
                script.write(f"""#!/bin/bash
#SBATCH --job-name=cfmc_rep{rep}_{num_solute}
#SBATCH --output=cfmc_rep{rep}_{num_solute}_%j.out
#SBATCH --error=cfmc_rep{rep}_{num_solute}_%j.err
#SBATCH --account=ucb-general
#SBATCH --partition=amilan
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks={config["nt"]}
#SBATCH --constraint=ib
#SBATCH --qos=normal

export OMPI_MCA_btl="self,openib,vader,tcp"
export OMPI_MCA_pml="ob1"
module purge
ml gcc
ml openmpi/5.0.6
ml anaconda
conda activate {config["conda_env"]}

cd "{sim_dir}"
gmx grompp -f em.mdp -c {sim_name}_final.gro -p {sim_name}.top -o em.tpr -maxwarn 1
gmx mdrun -deffnm em -nt {config["nt"]}

gmx grompp -f nvt.mdp -c em.gro -p {sim_name}.top -o nvt.tpr -maxwarn 1
gmx mdrun -deffnm nvt -nt {config["nt"]}

gmx grompp -f npt.mdp -c nvt.gro -p {sim_name}.top -o npt.tpr -maxwarn 1
gmx mdrun -deffnm npt -nt {config["nt"]}

gmx grompp -f prod.mdp -c npt.gro -p {sim_name}.top -o prod.tpr -maxwarn 1
gmx mdrun -deffnm prod -nt {config["nt"]}
""")
            os.chmod(script_path, 0o755)


if __name__ == "__main__":
    main()

