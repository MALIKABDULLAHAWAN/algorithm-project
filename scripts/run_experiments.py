import argparse
import csv
import os
import subprocess
from pathlib import Path


def run(cmd):
    print("RUN:", " ".join(cmd))
    p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    print(p.stdout)
    return p.returncode


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--bench', default=str(Path('build') / 'bench'), help='Path to bench executable')
    ap.add_argument('--out', default='results.csv', help='CSV output file (append)')
    ap.add_argument('--threads', type=int, nargs='+', default=[1,2,4,8])
    ap.add_argument('--gen', choices=['er','ba','rmat'], default='er')
    ap.add_argument('--n', type=int, default=20000)
    ap.add_argument('--p', type=float, default=0.0005)
    ap.add_argument('--m0', type=int, default=4)
    ap.add_argument('--m_attach', type=int, default=2)
    ap.add_argument('--scale', type=int, default=16)
    ap.add_argument('--edgefactor', type=int, default=16)
    ap.add_argument('--delta', type=float, nargs='+', default=[2.0,4.0,8.0])
    ap.add_argument('--seed', type=int, default=42)
    ap.add_argument('--input', default='', help='Edge-list file (overrides generator)')
    ap.add_argument('--skip-baseline', action='store_true', help='Skip Dijkstra baseline run')
    args = ap.parse_args()

    if args.out != 'results.csv':
        parent = Path(args.out).parent
        if parent and str(parent) != '.':
            os.makedirs(parent, exist_ok=True)

    def build_gen_args():
        gen_args = ['--gen', args.gen]
        if args.gen == 'er':
            gen_args += ['--n', str(args.n), '--p', str(args.p)]
        elif args.gen == 'ba':
            gen_args += ['--n', str(args.n), '--m0', str(args.m0), '--m_attach', str(args.m_attach)]
        elif args.gen == 'rmat':
            gen_args += ['--scale', str(args.scale), '--edgefactor', str(args.edgefactor)]
        return gen_args

    def build_input_args():
        return ['--input', args.input, '--weights','1','--one_based','1','--undirected','1']

    # Run Dijkstra baseline first (for speedup comparison)
    if not args.skip_baseline:
        cmd = [args.bench, '--algo', 'dijkstra', '--threads', '1', '--csv', args.out, '--seed', str(args.seed)]
        if args.input:
            cmd += build_input_args()
        else:
            cmd += build_gen_args()
        rc = run(cmd)
        if rc != 0:
            raise SystemExit(rc)

    # Run delta-stepping with varying threads and delta
    for t in args.threads:
        for d in args.delta:
            cmd = [args.bench,
                   '--algo','delta',
                   '--threads', str(t),
                   '--delta', str(d),
                   '--csv', args.out,
                   '--seed', str(args.seed)]
            if args.input:
                cmd += build_input_args()
            else:
                cmd += build_gen_args()
            rc = run(cmd)
            if rc != 0:
                raise SystemExit(rc)


if __name__ == '__main__':
    main()
