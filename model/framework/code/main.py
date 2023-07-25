#!/usr/bin/env python3
import os
import sys 

from molecule_generation import load_model_from_directory
from molecule_generation.utils.cli_utils import (
    setup_logging,
    supress_tensorflow_warnings,
)
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs

def tanimoto_calc(smi1, smi2):
    mol1 = Chem.MolFromSmiles(smi1)
    mol2 = Chem.MolFromSmiles(smi2)
    fp1 = AllChem.GetMorganFingerprintAsBitVect(mol1, 3, nBits=2048)
    fp2 = AllChem.GetMorganFingerprintAsBitVect(mol2, 3, nBits=2048)
    s = round(DataStructs.TanimotoSimilarity(fp1, fp2), 3)
    return s

def smiles_samples_to_csv(model_dir: str, num_samples: int, input_sample: str, output_file: str, **model_kwargs) -> None:
    with load_model_from_directory(model_dir) as model:
        samples = model.sample(num_samples)
        similarity_scores = []
        for sample in samples:
            similarity_score = tanimoto_calc(sample, input_sample)
            similarity_scores.append(similarity_score)
    dict_smiles = {'smiles': samples, 'similarity_score': similarity_scores}
    df = pd.DataFrame(dict_smiles)
    df = df.sort_values(by=['similarity_score'], ascending=False)
    with open(output_file, 'w') as f:
        df['smiles'].iloc[0:100].to_csv(output_file, sep='\t', index=False)

def main() -> None:
    supress_tensorflow_warnings()
    setup_logging()

    # Replace 'path/to/model_directory' with the actual absolute path to your model directory.
    ROOT = os.path.dirname(os.path.abspath(__file__))
    model_directory = os.path.abspath(os.path.join(ROOT,"..","..","checkpoints","MODEL_DIR"))
	
    input_file = sys.argv[1]
    output_file = sys.argv[2]

    num_sample = 1000

    input_sample_df = pd.read_csv(input_file)
    input_sample = input_sample_df['smiles'][0]

    smiles_samples_to_csv(model_dir=model_directory,
                          num_samples=num_sample,
                          input_sample=input_sample,
                          output_file=output_file,
                          beam_size=None,
                          seed=None,
                          num_workers=None)


if __name__ == "__main__":
    main()
