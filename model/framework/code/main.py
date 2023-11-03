#!/usr/bin/env python3
import os
import sys 

from molecule_generation import load_model_from_directory
from molecule_generation.utils.cli_utils import (
    setup_logging,
    supress_tensorflow_warnings,
)
import csv
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs

os.environ["CUDA_VISIBLE_DEVICES"] = "-1"


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

def scaffold_based_sampling(smiles_list):
    ROOT = os.path.dirname(os.path.abspath(__file__))
    model_directory = os.path.abspath(os.path.join(ROOT,"..","..","checkpoints","MODEL_DIR"))
    with load_model_from_directory(model_directory) as model:
        print(type(model))
        print(smiles_list)
        model.num_workers = 1
        embeddings = model.encode(smiles_list)
        print(embeddings)
        decoded = model.decode(embeddings, scaffolds=["CN", "CCC"])
    return decoded


def main() -> None:
    supress_tensorflow_warnings()
    setup_logging()

    # Replace 'path/to/model_directory' with the actual absolute path to your model directory.

	
    input_file = sys.argv[1]
    output_file = sys.argv[2]

    smiles_list = read_smiles(input_file=input_file)
    decoded = scaffold_based_sampling(smiles_list)

    print(decoded)


if __name__ == "__main__":
    main()
