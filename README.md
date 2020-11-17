# q2-qemistree
##### Canonically pronounced *chemis-tree*.

[![Build Status](https://travis-ci.org/biocore/q2-qemistree.svg?branch=master)](https://travis-ci.org/biocore/q2-qemistree) [![Coverage Status](https://coveralls.io/repos/github/biocore/q2-qemistree/badge.svg?branch=master)](https://coveralls.io/github/biocore/q2-qemistree?branch=master)

A tool to build a tree of mass-spectrometry (LC-MS/MS) features to perform chemically-informed comparison of untargeted metabolomic profiles. The preprint describing q2-qemistree is available [here](https://www.biorxiv.org/content/10.1101/2020.05.04.077636v1)!

## Installation

Once QIIME 2 is [installed](https://docs.qiime2.org/2019.7/install/), activate your QIIME 2 environment and install q2-qemistree following the steps below:

```bash
git clone https://github.com/biocore/q2-qemistree.git
cd q2-qemistree
pip install .
qiime dev refresh-cache
```

q2-qemistree uses [SIRIUS](https://www.nature.com/articles/s41592-019-0344-8), a software-framework developed for de-novo identification of metabolites. We use molecular substructures predicted by SIRIUS to build a hierarchy of the MS1 features in a dataset. For this demo, please download and unzip the latest version of SIRIUS from [here](https://bio.informatik.uni-jena.de/sirius/).

Below, we download SIRIUS for macOS as follows (for linux the only thing that changes is the URL from which the binary is downloaded):

```bash
wget https://bio.informatik.uni-jena.de/repository/dist-release-local/de/unijena/bioinf/ms/sirius/4.0.1/sirius-4.0.1-osx64-headless.zip
unzip sirius-4.0.1-osx64-headless.zip
```

## Demonstration

`q2-qemistree` ships with the following methods:

```
qiime qemistree compute-fragmentation-trees
qiime qemistree rerank-molecular-formulas
qiime qemistree predict-fingerprints
qiime qemistree make-hierarchy
qiime qemistree get-classyfire-taxonomy
qiime qemistree prune-hierarchy
qiime qemistree plot
```

To generate a tree that relates the MS1 features in your experiment, we need to pre-process mass-spectrometry data (.mzXML, .mzML or .mzDATA files) using [MZmine2](http://mzmine.github.io) and produce the following inputs:

1. An MGF file with both MS1 and MS2 information. This file will be imported into QIIME 2 as a `MassSpectrometryFeatures` artifact.
2. A feature table with peak areas of MS1 ions per sample. This table will be imported from a CSV file into the [BIOM](http://biom-format.org/documentation/biom_conversion.html) format, and then into QIIME 2 as a `FeatureTable[Frequency]` artifact.

These input files can be obtained following peak detection in MZmine2. [Here](https://raw.githubusercontent.com/biocore/q2-qemistree/master/q2_qemistree/demo/batchQE-MZmine-2.33.xml) is an example MZmine2 batch file used to generate these.

To begin this demonstration, create a separate folder to store all the inputs and outputs:

```bash
mkdir demo-qemistree
cd demo-qemistree
```

Download a small feature table and MGF file using:

```bash
wget https://raw.githubusercontent.com/biocore/q2-qemistree/master/q2_qemistree/demo/feature-table.biom
wget https://raw.githubusercontent.com/biocore/q2-qemistree/master/q2_qemistree/demo/sirius.mgf
```

We [import](https://docs.qiime2.org/2018.11/tutorials/importing/) these files into the appropriate QIIME 2 artifact formats as follows:

```bash
qiime tools import --input-path feature-table.biom --output-path feature-table.qza --type FeatureTable[Frequency]
qiime tools import --input-path sirius.mgf --output-path sirius.mgf.qza --type MassSpectrometryFeatures
```

**Note:** If the MGF file has formatting errors (eg. no MS1 are included in the MGF, or if an MS1 entry does not have a corresponding MS2 entry), then an appropriate error message will help users troubleshoot this step before proceeding forward.
First, we generate [fragmentation trees](https://www.sciencedirect.com/science/article/pii/S0165993615000916) for molecular peaks detected using MZmine2:

```bash
qiime qemistree compute-fragmentation-trees --p-sirius-path 'sirius-osx64-4.0.1/bin' \
  --i-features sirius.mgf.qza \
  --p-ppm-max 15 \
  --p-profile orbitrap \
  --p-ionization-mode positive \
  --p-java-flags "-Djava.io.tmpdir=/path-to-some-dir/ -Xms16G -Xmx64G" \
  --o-fragmentation-trees fragmentation_trees.qza
```
**Note**: `/path-to-some-dir/` should be a directory where you have write permissions and sufficient storage space. We use -Xms16G and -Xmx64G as the minimum and maximum heap size for Java virtual machine (JVM). If left blank, q2-qemistree will use default JVM flags.

This generates a QIIME 2 artifact of type `SiriusFolder`. This contains fragmentation trees with candidate molecular formulas for each MS1 feature detected in your experiment.

Next, we select top scoring molecular formula as follows:

```bash
qiime qemistree rerank-molecular-formulas --p-sirius-path 'sirius-osx64-4.0.1/bin' \
  --i-features sirius.mgf.qza \
  --i-fragmentation-trees fragmentation_trees.qza \
  --p-zodiac-threshold 0.95 \
  --p-java-flags "-Djava.io.tmpdir=/path-to-some-dir/ -Xms16G -Xmx64G" \
  --o-molecular-formulas molecular_formulas.qza
```

This produces a QIIME 2 artifact of type `ZodiacFolder` with top-ranked molecular formula for MS1 features. Now, we predict molecular substructures in each feature based on the molecular formulas. We use [CSI:FingerID](https://www.pnas.org/content/112/41/12580) for this purpose as follows:

```bash
qiime qemistree predict-fingerprints --p-sirius-path 'sirius-osx64-4.0.1/bin' \
  --i-molecular-formulas molecular_formulas.qza \
  --p-ppm-max 20 \
  --p-java-flags "-Djava.io.tmpdir=/path-to-some-dir/ -Xms16G -Xmx64G" \
  --o-predicted-fingerprints fingerprints.qza
  ```

This gives us a QIIME 2 artifact of type `CSIFolder` that contains probabilities of molecular substructures (total 2936 molecular properties) within in each feature.
We use these predicted molecular substructures to generate a hierarchy of molecules as follows:

```bash
qiime qemistree make-hierarchy \
  --i-csi-results fingerprints.qza \
  --i-feature-tables feature-table.qza \
  --o-tree qemistree.qza \
  --o-feature-table feature-table-hashed.qza \
  --o-feature-data feature-data.qza
```

To support meta-analyses, this method is capable of handling one or more datasets i.e pairs of CSI results and feature tables. You will need to download a new feature table and csi fingerprint result from another experiment to test this functionality as follows:

```bash
wget https://raw.githubusercontent.com/biocore/q2-qemistree/master/q2_qemistree/demo/feature-table2.biom.qza
wget https://raw.githubusercontent.com/biocore/q2-qemistree/master/q2_qemistree/demo/fingerprints2.qza
```

Below is the q2_qemistree command to co-analyze the datasets together:


```bash
qiime qemistree make-hierarchy \
--i-csi-results fingerprints.qza \
--i-csi-results fingerprints2.qza \
--i-feature-tables feature-table.qza \
--i-feature-tables feature-table2.biom.qza \
--o-tree merged-qemistree.qza \
--o-feature-table merged-feature-table-hashed.qza \
--o-feature-data merged-feature-data.qza
```
Additionally, Qemistree also supports the inclusion of structural annotations made using MS/MS spectral library matches for downstream analysis using the optional input `--i-ms2-matches` as follows:

```bash
qiime qemistree make-hierarchy \
  --i-csi-results fingerprints.qza \
  --i-feature-tables feature-table.qza \
  --i-ms2-matches /path-to-MS2-spectral-matches.qza/ \
  --o-tree qemistree.qza \
  --o-feature-table feature-table-hashed.qza \
  --o-feature-data feature-data.qza
```

**Note:**
1. The input to `--i-ms2-matches` can be obtained using [Feature-based molecular networking or FBMN](https://gnps.ucsd.edu/ProteoSAFe/index.jsp?params=%7B%22workflow%22:%22FEATURE-BASED-MOLECULAR-NETWORKING%22,%22library_on_server%22:%22d.speclibs;%22%7D) workflow supported in the web-based mass-spectrometry data analysis platform, [GNPS](https://gnps.ucsd.edu/). To use MS2 matches in Qemistree, please download the results of FBMN workflow and import the tsv file in the folder `clusterinfo_summary` as a QIIME2 artifact of type `FeatureData[Molecules]` as follows:

```bash
qiime tools import \
  --input-path path-to-MS2-spectral-matches.tsv \
  --output-path path-to-MS2-spectral-matches.qza \
  --type FeatureData[Molecules]
```

2. The input CSI results, feature tables and MS2 match tables should have a one-to-one correspondence i.e CSI results, feature tables and MS2 match tables from all datasets should be provided in the same order.

This method generates the following:
1. A combined feature table by merging all the input feature tables; MS1 features without fingerprints are filtered out of this feature table. This is done because SIRIUS predicts molecular substructures for a subset of features (typically for 70-90% of all MS1 features) in an experiment (based on factors such as sample type, the quality MS2 spectra, and user-defined tolerances such as `--p-ppm-max`, `--p-zodiac-threshold`). This output is of type `FeatureTable[Frequency]`.
2. A tree relating the MS1 features in these data based on molecular substructures predicted for MS1 features. This is of type `Phylogeny[Rooted]`. By default, we retain all fingerprint positions i.e. 2936 molecular properties). Adding `--p-qc-properties` filters these properties to keep only PubChem fingerprint positions (489 molecular properties) in the contingency table.
**Note**: The latest release of [SIRIUS](https://www.nature.com/articles/s41592-019-0344-8) uses PubChem version downloaded on 13 August 2017.
3. A combined feature data file that contains unique identifiers of each feature, their corresponding original feature identifier (row ID from Mzmine2), parent mass (`parent_mass`), retention time (`retention_time`), CSI:FingerID structure predictions (`csi_smiles`), MS2 match structure predictions (`ms2_smiles`), and the table(s) (`table_number`) that each feature was detected in. This is of type `FeatureData[Molecules]`. (The renaming of features helps prevent overlap between non-unique feature identifiers in the original feature tables in case of meta-analyses)

These can be used as inputs to perform chemical phylogeny-based [alpha-diversity](https://docs.qiime2.org/2019.1/plugins/available/diversity/alpha-phylogenetic/) and [beta-diversity](https://docs.qiime2.org/2019.1/plugins/available/diversity/beta-phylogenetic/) analyses.

Furthermore, Qemistree supports the classification of molecules into [Classyfire](https://jcheminf.biomedcentral.com/articles/10.1186/s13321-016-0174-y) chemical taxonomy. We generate a feature data table (also of the type `FeatureData[Molecules]`) which includes classification of molecules into chemical 'kingdom', 'superclass', 'class', 'subclass', and 'direct_parent'. We can run Classyfire using Qemistree as follows:

```bash
qiime qemistree get-classyfire-taxonomy \
  --i-feature-data merged-feature-data.qza \
  --o-classified-feature-data classified-merged-feature-data.qza
```
Qemistree will use `ms2_smiles` to make chemical taxonomy assignments, when MS2 matches are available for a feature. Otherwise, `csi_smiles` will be used. The column `structure_source` in `classified-merged-feature-data.qza` records whether taxonomic assignment was done using CSI:FingerID predictions or MS/MS library matches.

Lastly, Qemistree includes some utility functions that are useful to visualize and explore the molecular hierarchy generated above.
Qemistree trees can be visualized using [q2-empress](https://github.com/biocore/empress) [[preprint](https://www.biorxiv.org/content/10.1101/2020.10.06.327080v1)]. Below are the [installation instructions](https://github.com/biocore/empress#installation) that can be run within your qiime2 environment:

```bash
pip uninstall --yes emperor
pip install git+https://github.com/biocore/empress.git
qiime dev refresh-cache
```

1. Prune molecular hierarchy to keep only the molecules with annotations.

```bash
qiime qemistree prune-hierarchy \
  --i-feature-data classified-merged-feature-data.qza \
  --p-column class \
  --i-tree merged-qemistree.qza \
  --o-pruned-tree merged-qemistree-class.qza
```

Users can choose any of the data columns (`--p-column`) that are in the `classified-merged-feature-data.qza` file to prune the hierarchy. For e.g. '#featureID','kingdom', 'superclass', 'class', 'subclass', 'direct_parent', and 'smiles'. All features with no data in this column will be removed from the phylogeny.

2. Generate an annotated qemistree tree in using q2-empress.

```bash
qiime empress community-plot \
    --i-tree merged-qemistree-class.qza \
    --i-feature-table feature-table-hashed.qza \
    --m-sample-metadata-file path-to-sample-metadata.tsv \
    --m-feature-metadata-file classified-merged-feature-data.qza \
    --o-visualization empress-tree.qzv
```

The output empress QZV can be visualized using [Qiime2 Viewer](https://view.qiime2.org); EMPress can be used to interactively modify the tree visualization.
For instance, if the user has sample metadata columns to compare groups of samples, Empress enables them to visualize feature abundance barcharts at the tips of the tree (abundance or relative abundance) of the feature stratified by the sample metadata column of interest (see example below).
![Empress plot](img/gfop-empress-plot.png)
<div align="center">

Please visit [Empress tutorial](https://github.com/biocore/empress) for all the currently supported tree visualization features that can be leveraged to explore the chemical diversity of your metabolomics dataset.
