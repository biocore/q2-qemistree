# q2-chemistree

[![Build Status](https://travis-ci.org/biocore/q2-chemistree.svg?branch=master)](https://travis-ci.org/biocore/q2-chemistree) [![Coverage Status](https://coveralls.io/repos/github/biocore/q2-chemistree/badge.svg?branch=master)](https://coveralls.io/github/biocore/q2-chemistree?branch=master)

A tool to build a tree of MS1 features to compare chemical composition of samples in metabolomics experiments.

## Installation

Once QIIME 2 is [installed](https://docs.qiime2.org/2018.11/install/), activate your QIIME 2 environment and install q2-chemistree following the steps below:

```bash
git clone https://github.com/biocore/q2-chemistree.git
cd q2-chemistree
pip install .
qiime dev refresh-cache
```

q2-chemistree uses [SIRIUS](https://www.nature.com/articles/s41592-019-0344-8), a software-framework developed for de-novo identification of metabolites. We use molecular substrucures predicted by SIRIUS to build a hierarchy of the MS1 features in a dataset. For this demo, please download and unzip the latest version of SIRIUS from [here](https://bio.informatik.uni-jena.de/sirius/). Below, we download SIRIUS for macOS as follows (for linux the only thing that changes is the URL from which the binary is downloaded):

```bash
wget https://bio.informatik.uni-jena.de/repository/dist-release-local/de/unijena/bioinf/ms/sirius/4.0.1/sirius-4.0.1-osx64-headless.zip
unzip sirius-4.0.1-osx64-headless.zip
```

## Demonstration

`q2-chemistree` ships with the following methods:

```
qiime chemistree compute-fragmentation-trees
qiime chemistree rerank-molecular-formulas
qiime chemistree predict-fingerprints
qiime chemistree make-hierarchy
```

To generate a tree that relates the MS1 features in your experiment, we need to pre-process mass-spectrometry data (.mzXML files) using [MZmine2](http://mzmine.github.io) and produce the following inputs:

1. An MGF file with both MS1 and MS2 information. This file will be imported into QIIME 2 as a `MassSpectrometryFeatures` artifact.
2. A feature table with peak areas of MS1 ions per sample. This table will be imported from a CSV file into the [BIOM](http://biom-format.org/documentation/biom_conversion.html) format, and then into QIIME 2 as a `FeatureTable[Frequency]` artifact.

These input files can be obtained following peak detection in MZmine2. [Here](https://raw.githubusercontent.com/biocore/q2-chemistree/master/q2_chemistree/demo/batchQE-MZmine-2.33.xml) is an example MZmine2 batch file used to generate these.

To begin this demonstration, create a separate folder to store all the inputs and outputs:

```bash
mkdir demo-chemistree
cd demo-chemistree
```

Download a small feature table and MGF file using:

```bash
wget https://raw.githubusercontent.com/biocore/q2-chemistree/master/q2_chemistree/demo/feature-table.biom
wget https://raw.githubusercontent.com/biocore/q2-chemistree/master/q2_chemistree/demo/sirius.mgf
```

We [import](https://docs.qiime2.org/2018.11/tutorials/importing/) these files into the appropriate QIIME 2 artifact formats as follows:

```bash
qiime tools import --input-path feature-table.biom --output-path feature-table.qza --type FeatureTable[Frequency]
qiime tools import --input-path sirius.mgf --output-path sirius.mgf.qza --type MassSpectrometryFeatures
```

First, we generate [fragmentation trees](https://www.sciencedirect.com/science/article/pii/S0165993615000916) for molecular peaks detected using MZmine2:

```bash
qiime chemistree compute-fragmentation-trees --p-sirius-path 'sirius-osx64-4.0.1/bin' \
  --i-features sirius.mgf.qza \
  --p-ppm-max 15 \
  --p-profile orbitrap \
  --p-ionization-mode positive \
  --p-java-flags "-Djava.io.tmpdir=/path-to-some-dir/ -Xms16G -Xmx64G" \
  --o-fragmentation-trees fragmentation_trees.qza
```

This generates a QIIME 2 artifact of type `SiriusFolder`. This contains fragmentation trees with candidate molecular formulas for each MS1 feature detected in your experiment.
**Note**: `/path-to-some-dir/` should be a directory where you have write permissions and sufficient storage space. We use -Xms16G and -Xmx64G as the minimum and maximum heap size for Java virtual machine (JVM). If left blank, q2-chemistree will use default JVM flags.

Next, we select top scoring molecular formula as follows:

```bash
qiime chemistree rerank-molecular-formulas --p-sirius-path 'sirius-osx64-4.0.1/bin' \
  --i-features sirius.mgf.qza \
  --i-fragmentation-trees fragmentation_trees.qza \
  --p-zodiac-threshold 0.95 \
  --p-java-flags "-Djava.io.tmpdir=/path-to-some-dir/ -Xms16G -Xmx64G" \
  --o-molecular-formulas molecular_formulas.qza
```

This produces a QIIME 2 artifact of type `ZodiacFolder` with top-ranked molecular formula for MS1 features. Now, we predict molecular substructures in each feature based on the molecular formulas. We use [CSI:FingerID](https://www.pnas.org/content/112/41/12580) for this purpose as follows:

```bash
qiime chemistree predict-fingerprints --p-sirius-path 'sirius-osx64-4.0.1/bin' \
  --i-molecular-formulas molecular_formulas.qza \
  --p-ppm-max 20 \
  --p-java-flags "-Djava.io.tmpdir=/path-to-some-dir/ -Xms16G -Xmx64G" \
  --o-predicted-fingerprints fingerprints.qza
  ```

This gives us a QIIME 2 artifact of type `CSIFolder` that contains probabilities of molecular substructures (total 2936 molecular properties) within in each feature.
Now, we use these predicted molecular substructures to generate a hierarchy of molecules as follows:

```bash
qiime chemistree make-hierarchy \
  --i-csi-result fingerprints.qza \
  --i-feature-table feature-table.qza \
  --o-tree demo-chemisTree.qza \
  --o-merged-feature-table filtered-feature-table.qza
  --o-merged-feature-data feature-data.qza
```

To support meta-analyses, this method is capable of handling one or more datasets i.e pairs of CSI results and feature tables. **Note:** The input CSI results and feature tables should have a one-to-one correspondance i.e csi results and feature tables from all datasets should be provided in the same order.
This method performs the following tasks:
1. Generates a combined feature table by merging all the input feature tables based on fingerprint identity
2. Filters the MS1 features without fingerprints from the merged feature table such that the feature IDs in the feature table and the tree match. This is done because SIRIUS predicts molecular substructures for a subset of features (typically for 70-90% of all MS1 features) in an experiment (based on factors such as sample type, the quality MS2 spectra, and user-defined tolerances such as `--p-ppm-max`, `--p-zodiac-threshold`). The resulting feature table is also of type `FeatureTable[Frequency]`.
3. Generates a tree relating the MS1 features in these data based on molecular substructures predicted for MS1 features. This is of type `Phylogeny[Rooted]`. By default, we only use PubChem fingerprints (total 489 molecular properties). Adding `--p-no-qc-properties` retains all (2936) the molecular properties in the contingency table.
**Note**: The latest release of [SIRIUS](https://www.nature.com/articles/s41592-019-0344-8) uses PubChem version downloaded on 13 August 2017.
4. Generates a combined feature data file that contains information about unique identifiers of each feature (corresponding to the tip names on the tree) and maps it to the original feature identifiers and feature tables that each feature was detected in. This is of type `FeatureData[Molecules]`.


Thus, using these steps, we can generate a tree (`demo-chemisTree.qza`) relating MS1 features in a mass-spectrometry dataset along with a matched feature table (`filtered-feature-table.qza`). These can be used as inputs to perform chemical phylogeny-based [alpha-diversity](https://docs.qiime2.org/2018.11/plugins/available/diversity/alpha-phylogenetic/) and [beta-diversity](https://docs.qiime2.org/2018.11/plugins/available/diversity/beta-phylogenetic/) analyses.
