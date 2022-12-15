#!/usr/bin/env python3
import argparse

from molecule_generation import load_model_from_directory
from molecule_generation.utils.cli_utils import (
    get_model_loading_parser,
    setup_logging,
    supress_tensorflow_warnings,
)
import csv
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs

def tanimoto_calc(smi1, smi2):
    mol1 = Chem.MolFromSmiles(smi1)
    mol2 = Chem.MolFromSmiles(smi2)
    fp1 = AllChem.GetMorganFingerprintAsBitVect(mol1, 3, nBits=2048)
    fp2 = AllChem.GetMorganFingerprintAsBitVect(mol2, 3, nBits=2048)
    s = round(DataStructs.TanimotoSimilarity(fp1,fp2),3)
    return s


def smiles_samples_to_csv(model_dir: str, num_samples: int, input_sample: str, output_file: str, **model_kwargs) -> None:
    with load_model_from_directory(model_dir, **model_kwargs) as model:
        samples = model.sample(num_samples)
        similarity_scores = []
        for sample in samples:
            similarity_score = tanimoto_calc(sample, input_sample)
            similarity_scores.append(similarity_score)
    dict_smiles = {'smiles': samples, 'similarity_score': similarity_scores}
    df = pd.DataFrame(dict_smiles )
    df= df.sort_values(by=['similarity_score'],  ascending=False)
    with open(output_file, 'w') as f:
        df.iloc[0:99].to_csv(output_file, sep='\t', index=False)
      
def get_argparser() -> argparse.ArgumentParser:
    parser = get_model_loading_parser(description="Sample SMILES strings from a trained model.")
    #parser.add_argument("NUM_SAMPLES", type=str, help="Number of samples to generate.")
    parser.add_argument("INPUT_SAMPLE_csv", type=str, help="Input sample to calculate similarity")
    parser.add_argument("OUT_FILE", type=str, help="Output File directory.")
    parser.add_argument(
        "--beam-size", dest="beam_size", type=int, help="Beam size to use during decoding."
    )
    return parser


def run_from_args(args: argparse.Namespace) -> None:
    model_kwargs = {key: getattr(args, key) for key in ["beam_size", "seed", "num_workers"]}
    #num_samples_df = pd.read_csv(args.NUM_SAMPLES)
    num_sample = 1000
    
    input_sample_df = pd.read_csv(args.INPUT_SAMPLE_csv)
    input_sample = input_sample_df['smiles'][0]

    smiles_samples_to_csv(model_dir=args.MODEL_DIR,
        num_samples= num_sample,
        input_sample = input_sample,
        output_file = args.OUT_FILE, 
        **{key: value for (key, value) in model_kwargs.items() if value is not None},)

def main() -> None:
    supress_tensorflow_warnings()
    setup_logging()

    run_from_args(get_argparser().parse_args())


if __name__ == "__main__":
    main()
