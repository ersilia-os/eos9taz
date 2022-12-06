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


def smiles_samples_to_csv(model_dir: str, num_samples: int, output_file: str, **model_kwargs) -> None:
    with load_model_from_directory(model_dir, **model_kwargs) as model:
        samples = model.sample(num_samples)
        df = pd.DataFrame(samples, columns =['smiles'])
    with open(output_file, 'w') as f:
        df.to_csv(output_file)
      
def get_argparser() -> argparse.ArgumentParser:
    parser = get_model_loading_parser(description="Sample SMILES strings from a trained model.")
    parser.add_argument("NUM_SAMPLES", type=str, help="Number of samples to generate.")
    parser.add_argument("OUT_FILE", type=str, help="Output File directory.")
    parser.add_argument(
        "--beam-size", dest="beam_size", type=int, help="Beam size to use during decoding."
    )
    return parser


def run_from_args(args: argparse.Namespace) -> None:
    model_kwargs = {key: getattr(args, key) for key in ["beam_size", "seed", "num_workers"]}
    num_samples_df = pd.read_csv(args.NUM_SAMPLES)
    num_sample = num_samples_df['num_samples'][0]

    smiles_samples_to_csv(model_dir=args.MODEL_DIR,
        num_samples= num_sample,
        output_file = args.OUT_FILE,
        **{key: value for (key, value) in model_kwargs.items() if value is not None},)

def main() -> None:
    supress_tensorflow_warnings()
    setup_logging()

    run_from_args(get_argparser().parse_args())


if __name__ == "__main__":
    main()
