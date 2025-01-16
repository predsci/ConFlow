#!/usr/bin/env python3

import os
from pathlib import Path
import pandas as pd
import argparse

# Version 1.0.0

def argParsing():
    parser = argparse.ArgumentParser(
        description='conflow_realizations: Generates a HipFT flow input file spanning multiple realizations of ConFlow.'
    )
    parser.add_argument(
        'rootdir',
        type=str,
        help='Path of directory containing the ConFlow realization runs.'
    )
    parser.add_argument(
        '-o',
        dest='oFile',
        help='Name of output csv file.',
        default='conflow_realizations.csv'
    )
    parser.add_argument(
        '-outdir',
        dest='outdir',
        type=str,
        help='Directory where the soft links will be created.',
        required=False
    )
    parser.add_argument(
        '-nr',
        dest='nr',
        type=int,
        help='Number of realizations to be considered.',
        required=False
    )
    parser.add_argument(
        '-reps',
        dest='replicates',
        type=int,
        help='Number of repeats for each realization.',
        default=3
    )

    return parser.parse_args()

def run(args):
    if args.rootdir is None:
        print("Please provide the path of the directory containing the realizations.")
        return

    rootdir = Path(args.rootdir).resolve()
    args.outdir = args.outdir or os.path.join(os.getcwd(), 'conflow_realizations')
    outdir = Path(args.outdir).resolve()
    
    if args.outdir:
        os.makedirs(outdir, exist_ok=True)

    output_file = outdir / args.oFile
    with open(output_file, 'w') as ofile:
        ofile.write("   TIME(JD), VTFILENAME, VPFILENAME\n")

        idx = 1
        fidx = 0

        for dir_entry in rootdir.iterdir():
            if not dir_entry.is_dir():
                continue

            if args.nr and fidx >= args.nr:
                break

            fidx += 1
            for rep in range(args.replicates):
                print(f"Processing Realization: {dir_entry} - Replicate {rep+1}")

                for filevt, filevp in zip(sorted(dir_entry.glob('vt*.h5')),sorted(dir_entry.glob('vp*.h5'))):
                    vtpath = outdir / f'vt{idx:0>6}.h5'
                    vppath = outdir / f'vp{idx:0>6}.h5'

                    for file in [vtpath, vppath]:
                        if file.exists():
                            file.unlink()

                    vtpath.symlink_to(filevt)
                    vppath.symlink_to(filevp)

                    time = (1.0 / 96.0) * (idx - 1)
                    ofile.write('%11.5f,vt%06d.h5,vp%06d.h5\n' % (time, idx, idx))
                    idx += 1

    print(f"List file is saved in: {output_file}")

def main():
    args = argParsing()
    print(args)
    run(args)

if __name__ == '__main__':
    main()
