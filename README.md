# q2-chemistree

A tool to build a tree of MS1 features for use with phylogenetic methods in metabolomics datasets.

## Installation

Once QIIME2 is [installed](https://docs.qiime2.org/2018.2/install/), activate your QIIME2 environment and install q2-chemistree following the steps below:

```bash
git clone https://github.com/biocore/q2-chemistree.git
cd q2-chemistree
pip install .
qiime dev refresh-cache
```

We also need to download SIRIUS, a software-framework for de-novo identification of metabolites. q2-chemistree uses molecular substrucures predicted by this framework to build a hierarchy of the MS1 features in a dataset. SIRIUS is freely available [here](https://bio.informatik.uni-jena.de/software/sirius/).

## Demonstration

`q2-chemistree` ships with the following three methods:

```
qiime chemistree fingerprint
qiime chemistree make-hierarchy
qiime match-table
```

To generate a tree relating the MS1 features in your experiment, we need to pre-process mass-spectrometry data (.mzXML files) using [MZmine2](http://mzmine.github.io/) and produce the following:

Inputs:

1. An MGF file with both MS1 and MS2 information. This file will be imported into QIIME2 as a `MassSpectrometryFeatures` type.
2. A feature table with peak areas of MS1 ions per sample. This table will need to be imported from a CSV file into the [BIOM](http://biom-format.org/documentation/biom_conversion.html) format, and then into QIIME as a `FeatureTable[Frequency]` artifact.

These input files can be obtained following peak detection in MZmine2. [Here](https://github.com/anupriyatripathi/q2-chemistree/demo/batchQE-MZmine-2.33.xml) is an example MZmine2 batch file used to generate these.

You can create a separate folder to store files for this demonstration with:

```bash
mkdir demo-chemistree
cd demo-chemistree
```

Download a small feature table and MGF file using:

```bash
wget https://github.com/biocore/q2-chemistree/demo/feature-table.biom
wget https://github.com/biocore/q2-chemistree/demo/sirius.mgf
```

Once you [import](https://docs.qiime2.org/2018.8/tutorials/importing/) these files into the appropriate QIIME2 artifact format, you can follow these steps to generate your tree. Alternatively, you can proceed by directly downloading the corresponding artifact files as follows:

```bash
wget https://github.com/biocore/q2-chemistree/demo/feature-table.qza
wget https://github.com/biocore/q2-chemistree/demo/sirius.mgf.qza
```

Additionally, for this demo, we download SIRIUS for macOS as follows (for linux the only thing that changes is the URL from which the binary is downloaded):

```bash
wget https://bio.informatik.uni-jena.de/repository/dist-release-local/de/unijena/bioinf/ms/sirius/4.0.1/sirius-4.0.1-osx64-headless.zip
unzip sirius-4.0.1-osx64-headless.zip
```

First, we generate a contingency table with probabilities of presence of various molecular substructures in an MS1 feature (MS1 fingerprint):

```bash
qiime chemistree fingerprint --p-sirius-path 'sirius-osx64-4.0.1/bin' \
  --i-features sirius.mgf.qza \
  --p-profile orbitrap \
  --p-ppm-max 15 \
  --p-n-jobs 1 \
  --o-collated-fingerprints collatedfp.qza
```

This generates a contingency table with probabilities of various molecular substructures (i.e. molecular fingerprints) within an MS1 feature of type `FeatureTable[Frequency]`.

Now, we generate the molecular tree using:

```bash
qiime chemistree make-hierarchy \
  --i-collated-fingerprints collatedfp.qza \
  --p-prob-threshold 0.5 \
  --o-tree demo-chemisTree.qza
```

This generates a tree relating the MS1 features in these data based on molecular substructures predicted for MS1 features. This is of type `Phylogeny[Rooted]`.

Finally, note that SIRIUS predicts molecular substructures for a subset of features (typically for 70-90% of all MS1 features) in your experiment (based on factors such as sample type, the quality MS2 spectra, and used-defined tolerances when using the `fingerprint` command). Thus, we need to remove the MS1 features without fingerprints from the feature table with:

```bash
qiime chemistree match-table --i-tree demo-chemisTree.qza \
  --i-feature-table feature-table.qza \
  --o-filtered-feature-table filtered-feature-table.qza
```

This filters the MS1 table to include only the MS1 features with molecular fingerprints. The resulting table is also of type `FeatureTable[Frequency]`.

Thus, using these steps, we can generate a tree relating MS1 features in mass-spectrometry dataset along with a matched feature table. These can be used as inputs to perform [UniFrac-based](10.1038/ismej.2010.133) [alpha-diversity](https://docs.qiime2.org/2018.8/plugins/available/diversity/alpha-phylogenetic/) and [beta-diversity](https://docs.qiime2.org/2018.8/plugins/available/diversity/beta-phylogenetic/) analyses.
