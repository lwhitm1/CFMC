#!/usr/bin/env python3
import os
import argparse
from datetime import datetime

stages = ["em", "nvt", "npt", "prod"]
LOG_FILE = "restart_script_generation.log"

def log(msg):
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    with open(LOG_FILE, "a") as f:
        f.write(f"{timestamp} | {msg}\n")

def parse_args():
    parser = argparse.ArgumentParser(description="Generate restart scripts from existing job scripts.")
    parser.add_argument("--input_dir", "-i", default="cmds", help="Directory containing original job scripts.")
    parser.add_argument("--output_dir", "-o", default="restart_cmds", help="Directory to write restart scripts.")
    parser.add_argument("--nt", type=int, default=8, help="Number of threads to use in restart scripts (overrides parsed value).")
    return parser.parse_args()

def extract_common_parts(lines):
    header = []
    cd_line = None
    for line in lines:
        if line.startswith("cd "):
            cd_line = line.strip()
            break
        header.append(line.strip())
    return header, cd_line

def get_prev_stage(stage):
    idx = stages.index(stage)
    return stages[idx - 1] if idx > 0 else None

def extract_top_name(cd_line):
    last_dir = os.path.basename(cd_line.strip().split()[-1].strip('"'))
    return f"{last_dir}.top"

def is_stage_completed(sim_dir, stage):
    log_path = os.path.join(sim_dir, f"{stage}.log")
    if not os.path.exists(log_path):
        return False
    try:
        with open(log_path, "rb") as f:
            f.seek(-2048, os.SEEK_END)
            tail = f.read().decode(errors='ignore')
            return "Finished mdrun on" in tail
    except Exception:
        return False

def determine_restart_stage(sim_dir):
    for i in reversed(range(len(stages))):
        stage = stages[i]
        gro_path = os.path.join(sim_dir, f"{stage}.gro")

        if os.path.exists(gro_path):
            if stage == "prod":
                if not is_stage_completed(sim_dir, "prod"):
                    return "prod"
                else:
                    return None
            else:
                return stages[i + 1] if i + 1 < len(stages) else None
    return None

def write_restart_script(script_path, header, cd_line, top_name, nt, start_stage):
    lines = []
    lines.extend(header)
    lines.append("")
    lines.append(cd_line)
    lines.append("")

    lines.append(f"gmx mdrun -s {start_stage}.tpr -deffnm {start_stage} -cpi {start_stage}.cpt -nt {nt}")
    lines.append("")

    remaining = stages[stages.index(start_stage)+1:]
    for stage in remaining:
        prev = get_prev_stage(stage)
        lines.append(f"gmx grompp -f {stage}.mdp -c {prev}.gro -p {top_name} -o {stage}.tpr -maxwarn 1")
        lines.append(f"gmx mdrun -deffnm {stage} -nt {nt}")
        lines.append("")

    with open(script_path, "w") as f:
        f.write("\n".join(lines))
    os.chmod(script_path, 0o755)

def main():
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    log("=== New run: generating restart scripts ===")

    for filename in os.listdir(args.input_dir):
        if not filename.endswith(".sh"):
            continue

        in_path = os.path.join(args.input_dir, filename)
        with open(in_path, "r") as f:
            lines = f.readlines()

        header, cd_line = extract_common_parts(lines)
        if cd_line is None:
            msg = f"{filename} | No cd line found | Skipping"
            print(f"Skipping {filename}: no cd line found.")
            log(msg)
            continue

        sim_dir = cd_line.split("cd")[-1].strip().strip('"')
        sim_name = os.path.basename(sim_dir.rstrip("/"))

        if not os.path.isdir(sim_dir):
            msg = f"{filename} | {sim_name} | Directory '{sim_dir}' not found | Skipping"
            print(f"Skipping {filename}: simulation directory '{sim_dir}' not found.")
            log(msg)
            continue

        start_stage = determine_restart_stage(sim_dir)
        if start_stage is None:
            msg = f"{filename} | {sim_name} | Already complete | Skipping"
            print(f"Skipping {filename}: simulation already completed.")
            log(msg)
            continue

        top_name = extract_top_name(cd_line)

        # Priority: command-line arg > script header > fallback
        if args.nt is not None:
            nt = args.nt
        else:
            nt = 4  # fallback
            for line in header:
                if "--ntasks=" in line:
                    try:
                        nt = int(line.strip().split("=")[-1])
                    except ValueError:
                        pass

        out_path = os.path.join(args.output_dir, filename)
        write_restart_script(out_path, header, cd_line, top_name, nt, start_stage)

        msg = f"{filename} | {sim_name} | Restart from: {start_stage} | Threads: {nt} | Script generated"
        print(f"Generated restart script for {filename}, starting from '{start_stage}' with {nt} threads.")
        log(msg)

    print(f"\nAll restart scripts written to: {args.output_dir}/")
    log("=== Finished script generation ===\n")

if __name__ == "__main__":
    main()

