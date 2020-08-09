#!/usr/bin/env python

import argparse
from rdkit import Chem


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('iname', type=str)
    parser.add_argument('--margin', '-m', type=float, default=10)
    args = parser.parse_args()

    iname = args.iname
    margin = args.margin
    coords = []
    mol = Chem.MolFromMolFile(iname, sanitize=False)
    conf = mol.GetConformer(0)
    coords = conf.GetPositions()
    dims = [round(x, 3) for x in coords.max(axis=0) - coords.min(axis=0) + margin]
    center = [round(x, 3) for x in coords.mean(axis=0)]
    x, y, z = center
    print(f'--center_x {x} --center_y {y} --center_z {z} --size_x {dims[0]} --size_y {dims[1]} --size_z {dims[2]}')
