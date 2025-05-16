# Extending molecular scaffolds with fragments

MoLeR is a graph-based generative model that combines fragment-based and atom-by-atom generation of new molecules with scaffold-constrained optimization. It does not depend on generation history and therefore MoLeR is able to complete arbitrary scaffolds. The model has been trained on the GuacaMol dataset. Here we sample a fragment library from Enamine.

This model was incorporated on 2022-12-06.

## Information
### Identifiers
- **Ersilia Identifier:** `eos9taz`
- **Slug:** `moler-enamine-fragments`

### Domain
- **Task:** `Sampling`
- **Subtask:** `Generation`
- **Biomedical Area:** `Any`
- **Target Organism:** `Not Applicable`
- **Tags:** `Chemical graph model`, `Compound generation`

### Input
- **Input:** `Compound`
- **Input Dimension:** `1`

### Output
- **Output Dimension:** `1000`
- **Output Consistency:** `Variable`
- **Interpretation:** 1000 new molecules are sampled for each input molecule, preserving its scaffold. 

Below are the **Output Columns** of the model:
| Name | Type | Direction | Description |
|------|------|-----------|-------------|
| cpd_000 | string |  | Generated molecule index 0 using the MoLeR molecular generator |
| cpd_001 | string |  | Generated molecule index 1 using the MoLeR molecular generator |
| cpd_002 | string |  | Generated molecule index 2 using the MoLeR molecular generator |
| cpd_003 | string |  | Generated molecule index 3 using the MoLeR molecular generator |
| cpd_004 | string |  | Generated molecule index 4 using the MoLeR molecular generator |
| cpd_005 | string |  | Generated molecule index 5 using the MoLeR molecular generator |
| cpd_006 | string |  | Generated molecule index 6 using the MoLeR molecular generator |
| cpd_007 | string |  | Generated molecule index 7 using the MoLeR molecular generator |
| cpd_008 | string |  | Generated molecule index 8 using the MoLeR molecular generator |
| cpd_009 | string |  | Generated molecule index 9 using the MoLeR molecular generator |

_10 of 1000 columns are shown_
### Source and Deployment
- **Source:** `Local`
- **Source Type:** `External`
- **DockerHub**: [https://hub.docker.com/r/ersiliaos/eos9taz](https://hub.docker.com/r/ersiliaos/eos9taz)
- **Docker Architecture:** `AMD64`, `ARM64`
- **S3 Storage**: [https://ersilia-models-zipped.s3.eu-central-1.amazonaws.com/eos9taz.zip](https://ersilia-models-zipped.s3.eu-central-1.amazonaws.com/eos9taz.zip)

### Resource Consumption
- **Model Size (Mb):** `23`
- **Environment Size (Mb):** `2636`
- **Image Size (Mb):** `2581.98`

**Computational Performance (seconds):**
- 10 inputs: `40.44`
- 100 inputs: `400.61`
- 10000 inputs: `593.91`

### References
- **Source Code**: [https://github.com/microsoft/molecule-generation](https://github.com/microsoft/molecule-generation)
- **Publication**: [https://arxiv.org/abs/2103.03864](https://arxiv.org/abs/2103.03864)
- **Publication Type:** `Preprint`
- **Publication Year:** `2021`
- **Ersilia Contributor:** [anamika-yadav99](https://github.com/anamika-yadav99)

### License
This package is licensed under a [GPL-3.0](https://github.com/ersilia-os/ersilia/blob/master/LICENSE) license. The model contained within this package is licensed under a [MIT](LICENSE) license.

**Notice**: Ersilia grants access to models _as is_, directly from the original authors, please refer to the original code repository and/or publication if you use the model in your research.


## Use
To use this model locally, you need to have the [Ersilia CLI](https://github.com/ersilia-os/ersilia) installed.
The model can be **fetched** using the following command:
```bash
# fetch model from the Ersilia Model Hub
ersilia fetch eos9taz
```
Then, you can **serve**, **run** and **close** the model as follows:
```bash
# serve the model
ersilia serve eos9taz
# generate an example file
ersilia example -n 3 -f my_input.csv
# run the model
ersilia run -i my_input.csv -o my_output.csv
# close the model
ersilia close
```

## About Ersilia
The [Ersilia Open Source Initiative](https://ersilia.io) is a tech non-profit organization fueling sustainable research in the Global South.
Please [cite](https://github.com/ersilia-os/ersilia/blob/master/CITATION.cff) the Ersilia Model Hub if you've found this model to be useful. Always [let us know](https://github.com/ersilia-os/ersilia/issues) if you experience any issues while trying to run it.
If you want to contribute to our mission, consider [donating](https://www.ersilia.io/donate) to Ersilia!
