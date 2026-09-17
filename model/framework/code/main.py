#!/usr/bin/env python3
import os
import sys 

from molecule_generation import load_model_from_directory
from molecule_generation.utils.cli_utils import (
    setup_logging,
    supress_tensorflow_warnings,
)
import csv
import random
from tqdm import tqdm
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs
from rdkit.Chem.Scaffolds import MurckoScaffold

os.environ["CUDA_VISIBLE_DEVICES"] = "-1"

ROOT = os.path.dirname(os.path.abspath(__file__))
BLOCKS_LIST = os.path.join(ROOT, "..", "..", "checkpoints", "fragments_from_enamine.smi")

N_SAMPLES = 1000


def get_murcko_scaffold(smiles):
    molecule = Chem.MolFromSmiles(smiles)
    scaffold = MurckoScaffold.GetScaffoldForMol(molecule)
    scaffold_smiles = Chem.MolToSmiles(scaffold)
    return scaffold_smiles


def read_blocks():
    blocks_list = []
    with open(BLOCKS_LIST, "r") as f:
        reader = csv.reader(f, delimiter="\t")
        for r in reader:
            blocks_list += [r[0]]
    return blocks_list

def read_smiles(input_file):
    smiles = []
    with open(input_file, "r") as f:
        reader = csv.reader(f)
        next(reader)
        for r in reader:
            smiles += [r[0]]
    print("These are the SMILES: ", smiles)
    return smiles

def tanimoto_calc(smi1, smi2):
    mol1 = Chem.MolFromSmiles(smi1)
    mol2 = Chem.MolFromSmiles(smi2)
    fp1 = AllChem.GetMorganFingerprintAsBitVect(mol1, 3, nBits=2048)
    fp2 = AllChem.GetMorganFingerprintAsBitVect(mol2, 3, nBits=2048)
    s = round(DataStructs.TanimotoSimilarity(fp1, fp2), 3)
    return s

def scaffold_based_sampling_unique(query_scaffold, blocks_list, model, target=N_SAMPLES, max_rounds=5, rng=None):
    # A single encode/decode call against `target` random fragments was found to return
    # heavily duplicated molecules (up to ~82% duplicates on some scaffolds, confirmed
    # empirically) with no dedup anywhere in this pipeline. Draw fresh, non-overlapping
    # fragment batches across multiple rounds (reusing the already-loaded model) until
    # `target` unique molecules are collected or the round budget / fragment pool runs out.
    rng = rng or random.Random()
    pool = list(blocks_list)
    rng.shuffle(pool)

    seen = set()
    results = []
    idx = 0
    for _ in range(max_rounds):
        if len(results) >= target or idx >= len(pool):
            break
        batch = pool[idx: idx + target]
        idx += len(batch)

        embeddings = model.encode(batch)
        decoded = model.decode(embeddings, scaffolds=[query_scaffold] * len(batch))
        for o in decoded:
            if not o:
                continue
            mol = Chem.MolFromSmiles(o)
            key = Chem.MolToSmiles(mol) if mol is not None else o
            if key in seen:
                continue
            seen.add(key)
            results.append(o)
            if len(results) >= target:
                break

    results = results[:target]
    if len(results) < target:
        results = results + [""] * (target - len(results))
    return results

def main() -> None:
    supress_tensorflow_warnings()
    setup_logging()

    blocks_list = read_blocks()

    input_file = sys.argv[1]
    output_file = sys.argv[2]

    smiles_list = read_smiles(input_file=input_file)
    scaff_list = [get_murcko_scaffold(smi) for smi in smiles_list]

    model_directory = os.path.abspath(os.path.join(ROOT, "..", "..", "checkpoints", "MODEL_DIR"))
    rng = random.Random(42)

    R = []
    with load_model_from_directory(model_directory) as model:
        for smi in tqdm(scaff_list):
            decoded = scaffold_based_sampling_unique(smi, blocks_list, model, rng=rng)
            R += [decoded]

    with open(output_file, "w") as f:
        writer = csv.writer(f)
        header = ["smi_{:03d}".format(i) for i in range(N_SAMPLES)]
        writer.writerow(header)
        for r in R:
            writer.writerow(r)


if __name__ == "__main__":
    main()
