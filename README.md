# Extending molecular scaffolds with structural motifs

MoLeR is a graph-based generative model that combines fragment-based and atom-by-atom generation of new molecules with scaffold-constrained optimization. It does not depend on generation history and therefore MoLeR is able to complete arbitrary scaffolds. The model has been trained on the GuacaMol dataset.

## Identifiers

* EOS model ID: `eos9taz`
* Slug: `moler-generation`

## Characteristics

* Input: `Compound`
* Input Shape: `Single`
* Task: `Generative`
* Output: `Compound`
* Output Type: `String`
* Output Shape: `List`
* Interpretation: 100 new molecules are sampled for each input molecule. 

## References

* [Publication](https://arxiv.org/abs/2103.03864)
* [Source Code](https://github.com/microsoft/molecule-generation)
* Ersilia contributor: [anamika-yadav99](https://github.com/anamika-yadav99)

## Citation

If you use this model, please cite the [original authors](https://arxiv.org/abs/2103.03864) of the model and the [Ersilia Model Hub](https://github.com/ersilia-os/ersilia/blob/master/CITATION.cff).

## License

This package is licensed under a GPL-3.0 license. The model contained within this package is licensed under a MIT license.

Notice: Ersilia grants access to these models 'as is' provided by the original authors, please refer to the original code repository and/or publication if you use the model in your research.

## About Us

The [Ersilia Open Source Initiative](https://ersilia.io) is a Non Profit Organization ([1192266](https://register-of-charities.charitycommission.gov.uk/charity-search/-/charity-details/5170657/full-print)) with the mission is to equip labs, universities and clinics in LMIC with AI/ML tools for infectious disease research.

[Help us](https://www.ersilia.io/donate) achieve our mission!