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
from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator
from rdkit.Chem.Scaffolds import MurckoScaffold

os.environ["CUDA_VISIBLE_DEVICES"] = "-1"

ROOT = os.path.dirname(os.path.abspath(__file__))
BLOCKS_LIST = os.path.join(ROOT, "..", "..", "checkpoints", "fragments_from_enamine.smi")

N_SAMPLES = 100


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

def scaffold_based_sampling_unique(query_scaffold, blocks_list, model, target=N_SAMPLES, max_rounds=5, rng=None, exclude_smiles=None):
    # A single encode/decode call against `target` random fragments was found to return
    # heavily duplicated molecules (up to ~82% duplicates on some scaffolds, confirmed
    # empirically) with no dedup anywhere in this pipeline. Draw fresh, non-overlapping
    # fragment batches across multiple rounds (reusing the already-loaded model) until
    # `target` unique molecules are collected or the round budget / fragment pool runs out.
    #
    # `exclude_smiles` is the input molecule: it must never be returned as one of its own
    # "generated" outputs. Compared without stereochemistry on purpose -- the decoder mostly
    # emits stereo-free SMILES, so an exact isomeric match would miss most echoes of an input
    # that carries stereocenters.
    exclude_flat = None
    if exclude_smiles:
        exclude_mol = Chem.MolFromSmiles(exclude_smiles)
        if exclude_mol is not None:
            exclude_flat = Chem.MolToSmiles(exclude_mol, isomericSmiles=False)

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
            if mol is None:
                # MoLeR occasionally serializes a ring-fusion carbon as an explicit
                # [CH] that already has four heavy-atom neighbours, giving valence 5;
                # RDKit rejects it. Drop it rather than emit an unusable SMILES.
                continue
            if exclude_flat is not None and Chem.MolToSmiles(mol, isomericSmiles=False) == exclude_flat:
                continue
            key = Chem.MolToSmiles(mol)
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

def sort_by_similarity(input_smiles, outputs):
    # Most similar to the input first: Tanimoto on Morgan fingerprints (radius 2, 2048 bits)
    # between each output and the whole input. Ties break on the SMILES so the order is
    # reproducible. MoLeR returns the molecules in decoding order, which carries no meaning.
    # Empty cells stay at the end.
    filled = [o for o in outputs if o]
    input_mol = Chem.MolFromSmiles(input_smiles)
    if input_mol is None or not filled:
        return outputs
    gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    input_fp = gen.GetFingerprint(input_mol)

    def similarity(smiles):
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return -1.0
        return DataStructs.TanimotoSimilarity(input_fp, gen.GetFingerprint(mol))

    ordered = sorted(filled, key=lambda smiles: (-similarity(smiles), smiles))
    return ordered + [""] * (len(outputs) - len(ordered))


def main() -> None:
    supress_tensorflow_warnings()
    setup_logging()

    blocks_list = read_blocks()

    input_file = sys.argv[1]
    output_file = sys.argv[2]

    smiles_list = read_smiles(input_file=input_file)

    model_directory = os.path.abspath(os.path.join(ROOT, "..", "..", "checkpoints", "MODEL_DIR"))
    # Fresh entropy on every call, so the same input file gives different molecules at each run.
    # The model wrapper also reseeds the global random, numpy and TensorFlow generators to its
    # `seed` (default 0) when it loads, so draw that seed from our own generator as well.
    rng = random.Random()

    R = []
    with load_model_from_directory(model_directory, seed=rng.randrange(2**31)) as model:
        for input_smi in tqdm(smiles_list):
            # One bad input (invalid SMILES, a scaffold the decoder cannot handle) must not stop
            # the other inputs: it gets an empty row and the run goes on.
            try:
                scaff = get_murcko_scaffold(input_smi)
                decoded = scaffold_based_sampling_unique(scaff, blocks_list, model, rng=rng, exclude_smiles=input_smi)
            except Exception as e:
                print(f"[ERROR] {input_smi!r}: {type(e).__name__}: {e}", file=sys.stderr)
                decoded = [""] * N_SAMPLES
            try:
                decoded = sort_by_similarity(input_smi, decoded)
            except Exception as e:
                # an ordering problem must never cost the row itself
                print(f"[WARN] could not sort the outputs of {input_smi!r}: {type(e).__name__}: {e}", file=sys.stderr)
            R += [decoded]

    with open(output_file, "w") as f:
        writer = csv.writer(f)
        header = ["smi_{:02d}".format(i) for i in range(N_SAMPLES)]
        writer.writerow(header)
        for r in R:
            writer.writerow(r)


if __name__ == "__main__":
    main()
