(gut-to-soil-16S-tutorial)=
# Gut-to-soil axis 16S rRNA analysis tutorial 💩🌱

:::{warning}
This document is a work in progress as of 8 March 2025.
As of this time, it's **not ready for use**, and only live for testing purposes.

Commands/urls/text may be unreliable while in development.
🚜
:::

:::{note}
This guide assumes you are working in a QIIME 2 deployment containing the amplicon distribution and the standalone plugins [q2-boots](https://doi.org/10.12688/f1000research.156295.1) and [q2-kmerizer](https://doi.org/10.1128/msystems.01550-24).
The conda environment file used to create the deployment where the commands in this tutorial were run [can be found in the project's GitHub repository](https://github.com/caporaso-lab/gut-to-soil-tutorial/blob/main/environment-files/readthedocs.yml).
:::

## Background

In this tutorial you'll learn an end-to-end microbiome data science workflow, building on data presented in [Meilander *et al.* (2024): Upcycling Human Excrement: The Gut Microbiome to Soil Microbiome Axis](https://doi.org/10.48550/arXiv.2411.04148).
The data used here is a subset (a single sequencing run) of that generated for the paper, specifically selected so that this tutorial can be run quickly on a personal computer.
The full data set for the paper can be found in [the paper's Artifact Archive](https://doi.org/10.5281/zenodo.13887456).
In the final step, you'll learn how to adapt the workflow for use in analyzing your own data using [Provenance Replay](https://doi.org/10.1371/journal.pcbi.1011676).

The data used in this tutorial was generated using the [Earth Microbiome Project protocol](https://doi.org/10.1038/ismej.2012.8).
Specifically, the hypervariable region 4 (V4) of the 16S rRNA gene was amplified using the F515-R806 primers - a broad-coverage primer pair for Bacteria that also amplifies some Archaea.
Paired-end sequencing was performed on an Illumina MiSeq.
Full details are presented in [Meilander *et al.* (2024)](https://doi.org/10.48550/arXiv.2411.04148).

(gut-to-soil-tutorial:sample-metadata)=
## Sample metadata

Before starting the analysis, explore the sample metadata to familiarize yourself with the samples used in this study.
The following command will download the sample metadata as tab-separated text and save it in the file `sample-metadata.tsv`.
This `sample-metadata.tsv` file is used throughout the rest of the tutorial.

:::{tip}
To learn more about metadata in QIIME 2, including how it should be formatted, refer to [*Using QIIME 2*'s Metadata file format](https://use.qiime2.org/en/latest/references/metadata.html).
:::

:::{describe-usage}
:scope: gut-to-soil

sample_metadata = use.init_metadata_from_url(
   'sample-metadata',
   'https://www.dropbox.com/scl/fi/irosimbb1aud1aa7frzxf/sample-metadata.tsv?rlkey=f45jpxzajjz9xx9vpvfnf1zjx&st=nahafuvy&dl=1')
:::

QIIME 2's metadata plugin provides a Visualizer called `tabulate` that generates a convenient view of a sample metadata file.
To learn more about `metadata tabulate`, you can find it in the plugin reference [here](../plugin-reference/plugins/metadata/3-tabulate).
Let's run this, and then we'll look at the result.
Here's the first QIIME 2 command that you should run in this tutorial:

:::{describe-usage}
sample_metadata_viz, = use.action(
  use.UsageAction(plugin_id='metadata',
                  action_id='tabulate'),
  use.UsageInputs(input=sample_metadata),
  use.UsageOutputNames(visualization='sample_metadata_viz')
)
:::

This will generate a QIIME 2 Visualization.
Visualizations can be viewed by loading them with [QIIME 2 View](https://view.qiime2.org).
Navigate to QIIME 2 View, and drag and drop the visualization that was created to view it.
(You can learn more about viewing Visualizations, including alternatives to QIIME 2 View if you can't use that for any reason, [here](https://use.qiime2.org/en/latest/how-to-guides/view-visualizations.html).)

## Access already-imported QIIME 2 data

You should ask your sequencing center to provide data already in QIIME 2 "demux" artifacts, or provide a QIIME 2 fastq manifest file for your sequencing data.
It should be easy for them to generate.

:::{describe-usage}
:scope: gut-to-soil

demux = use.init_artifact_from_url(
   'demux',
   'https://www.dropbox.com/scl/fi/73f6rmcq7lelug5qbuhl6/demux-10p.qza?rlkey=s0hoh326fes3z2usvgs62tv3c&st=caz1avkn&dl=1')
:::

:::{tip}
Links are included to view and download precomputed QIIME 2 artifacts and visualizations created by commands in the documentation.
For example, the command above created a single `emp-single-end-sequences.qza` file and that file is precomputed and linked following the command.
These are provided so that you can look at Visualizations before you run them, to get a feeling for whether QIIME 2 is going to do what you need before taking the time to install and run QIIME 2.
The linked Artifacts (i.e., `.qza` files) are useful as you branch out and start using other QIIME 2 functionality as they provide example files that you can use to test other commands.
:::

## Summarizing demultiplexed sequences

After demultiplexing, it's useful to generate a summary of the demultiplexing results.
This allows you to determine how many sequences were obtained per sample, and also to get a summary of the distribution of sequence qualities at each position in your sequence data.

:::{describe-usage}

use.action(
    use.UsageAction(plugin_id='demux',
                    action_id='summarize'),
    use.UsageInputs(data=demux),
    use.UsageOutputNames(visualization='demux'))
:::

:::{note}
All QIIME 2 visualizers (i.e., commands that take a `--o-visualization` parameter) will generate a Visualization (i.e., a`.qzv` file).
As discussed above, Visualizations can be viewed by loading them with [QIIME 2 View](https://view.qiime2.org) or using one of the alternative mechanisms discussed [here](https://use.qiime2.org/en/latest/how-to-guides/view-visualizations.html).
When you see instructions to *view a visualization*, you should use one of the approaches on the `.qzv` file that we're referring to.
Navigate to QIIME 2 View, and drag and drop the visualization that was created to view it.
(You can learn more about viewing Visualizations, including alternatives to QIIME 2 View if you can't use that for any reason, .)
:::

## Sequence quality control and feature table construction

QIIME 2 plugins are available for several quality control methods, including [DADA2](https://www.ncbi.nlm.nih.gov/pubmed/27214047), [Deblur](http://msystems.asm.org/content/2/2/e00191-16), and [basic quality-score-based filtering](http://www.nature.com/nmeth/journal/v10/n1/abs/nmeth.2276.html).
In this tutorial we present this step using [DADA2](https://www.ncbi.nlm.nih.gov/pubmed/27214047).
The result of this method will be a `FeatureTable[Frequency]` QIIME 2 artifact, which contains counts (frequencies) of each unique sequence in each sample in the dataset, and a `FeatureData[Sequence]` QIIME 2 artifact, which maps feature identifiers in the `FeatureTable` to the sequences they represent.

[DADA2](https://www.ncbi.nlm.nih.gov/pubmed/27214047) is a pipeline for detecting and correcting (where possible) Illumina amplicon sequence data.
As implemented in the `q2-dada2` plugin, this quality control process will additionally filter any phiX reads (commonly present in marker gene Illumina sequence data) that are identified in the sequencing data, and will filter chimeric sequences.

The `dada2 denoise-paired` method requires four parameters that are used in quality filtering:
- `--p-trim-left-f a`, which trims off the first `a` bases of each forward read
- `--p-trunc-len-f b` which truncates each forward read at position `b`
- `--p-trim-left-r c`, which trims off the first `c` bases of each forward read
- `--p-trunc-len-r d` which truncates each forward read at position `d`
This allows the user to remove low quality regions of the sequences.
To determine what values to pass for these two parameters, you should review the *Interactive Quality Plot* tab in the `demux.qzv` file that was generated by `qiime demux summarize` above.

:::{tip} Question.
Based on the plots you see in `demux.qzv`, what values would you choose for `--p-trim-left-f` and `--p-trunc-len-f` in this case?
:::

In the `demux.qzv` quality plots, we see that the quality of the initial bases seems to be high, so we won't trim any bases from the beginning of the sequences.
The quality seems good all the way out to the end, so we'll trim at 250 bases.
This next command may take up to 10 minutes to run, and is the slowest step in this tutorial.

:::{describe-usage}

rep_seqs, table, stats = use.action(
    use.UsageAction(plugin_id='dada2',
                    action_id='denoise_paired'),
    use.UsageInputs(demultiplexed_seqs=demux,
                    trim_left_f=0,
                    trunc_len_f=250,
                    trim_left_r=0,
                    trunc_len_r=250),
    use.UsageOutputNames(representative_sequences='rep_seqs',
                         table='table',
                         denoising_stats='stats'))
:::

:::{describe-usage}
stats_as_md = use.view_as_metadata('stats_md', stats)

use.action(
    use.UsageAction(plugin_id='metadata',
                    action_id='tabulate'),
    use.UsageInputs(input=stats_as_md),
    use.UsageOutputNames(visualization='stats'))
:::

## FeatureTable and FeatureData summaries

After the quality filtering step completes, you'll want to explore the resulting data.
You can do this using the following two commands, which will create visual summaries of the data.
The `feature-table summarize` command will give you information on how many sequences are associated with each sample and with each feature, histograms of those distributions, and some related summary statistics.
The `feature-table tabulate-seqs` command will provide a mapping of feature IDs to sequences, and provide links to easily BLAST each sequence against the NCBI nt database.
The latter visualization will be very useful later in the tutorial, when you want to learn more about specific features that are important in the data set.

:::{describe-usage}
_, _, feature_frequencies = use.action(
    use.UsageAction(plugin_id='feature_table',
                    action_id='summarize_plus'),
    use.UsageInputs(table=table,
                    metadata=sample_metadata),
    use.UsageOutputNames(summary='table',
                         sample_frequencies='sample_frequencies',
                         feature_frequencies='feature_frequencies'))

feature_frequencies_as_md = use.view_as_metadata('feature_frequencies_md',
                                                 feature_frequencies)

use.action(
    use.UsageAction(plugin_id='feature_table',
                    action_id='tabulate_seqs'),
    use.UsageInputs(data=rep_seqs,
                    metadata=feature_frequencies_as_md),
    use.UsageOutputNames(visualization='rep_seqs'),
)
:::

If you review the tabulated feature sequences, or the feature detail table of the feature table summary, you'll notice that there are many sequences that are observed in only a single sample.
Let's filter those out to reduce the number of sequences we're working with - this will speed up several slower steps that are coming up.
First we filter our feature table, and then we use the new feature table to filter our sequences to only the ones that are contained in the new table.
We'll include `_ms2` in the new output names to remind us that these are filtered to those features in a **m**inimum of **2** **s**amples.

:::{describe-usage}
table_ms2, = use.action(
    use.UsageAction(plugin_id='feature_table',
                    action_id='filter_features'),
    use.UsageInputs(table=table,
                    min_samples=2),
    use.UsageOutputNames(filtered_table='table_ms2'),
)
:::


:::{describe-usage}
rep_seqs_ms2, = use.action(
    use.UsageAction(plugin_id='feature_table',
                    action_id='filter_seqs'),
    use.UsageInputs(data=rep_seqs,
                    table=table_ms2),
    use.UsageOutputNames(filtered_data='rep_seqs_ms2'),
)
:::

:::{exercise}
:label: filtered-feature-table-summary
Now that you have a second (filtered) feature table, create your own command to summarize it, like we did for the original feature table.
How many features did we lose as a result of this filter?
How many total sequences did we lose?
:::

::::{solution} filtered-feature-table-summary
:class: dropdown
Here's the command you would use:

:::{describe-usage}
use.action(
    use.UsageAction(plugin_id='feature_table',
                    action_id='summarize_plus'),
    use.UsageInputs(table=table_ms2,
                    metadata=sample_metadata),
    use.UsageOutputNames(summary='table_ms2',
                         sample_frequencies='sample_frequencies_ms2',
                         feature_frequencies='feature_frequencies_ms2'))
:::
::::


## Generate a tree for phylogenetic diversity analyses

**Replace with q2-kmerizer and q2-boots?**

:::{describe-usage}
kmer_table, = use.action(
    use.UsageAction(plugin_id='kmerizer',
                    action_id='seqs_to_kmers'),
    use.UsageInputs(table=table_ms2,
                    sequences=rep_seqs_ms2),
    use.UsageOutputNames(kmer_table='kmer_table')
)
:::

:::{describe-usage}
use.action(
    use.UsageAction(plugin_id='feature_table',
                    action_id='summarize'),
    use.UsageInputs(table=kmer_table,
                    sample_metadata=sample_metadata),
    use.UsageOutputNames(visualization='kmer_table')
)
:::

QIIME supports several phylogenetic diversity metrics, including Faith's Phylogenetic Diversity and weighted and unweighted UniFrac.
In addition to counts of features per sample (i.e., the data in the `FeatureTable[Frequency]` QIIME 2 artifact), these metrics require a rooted phylogenetic tree relating the features to one another.
This information will be stored in a `Phylogeny[Rooted]` QIIME 2 artifact.
To generate a phylogenetic tree we will use `align-to-tree-mafft-fasttree` pipeline from the `q2-phylogeny` plugin.

First, the pipeline uses the `mafft` program to perform a multiple sequence alignment of the sequences in our `FeatureData[Sequence]` to create a `FeatureData[AlignedSequence]` QIIME 2 artifact.
Next, the pipeline masks (or filters) the alignment to remove positions that are highly variable.
These positions are generally considered to add noise to a resulting phylogenetic tree.
Following that, the pipeline applies FastTree to generate a phylogenetic tree from the masked alignment.
The FastTree program creates an unrooted tree, so in the final step in this section midpoint rooting is applied to place the root of the tree at the midpoint of the longest tip-to-tip distance in the unrooted tree.

## Alpha and beta diversity analysis

**Add note about q2-boots.**
**Integrate some sample filtering, e.g., to just HE and HEC post-roll?**
**Paired difference testing - pre- versus post-roll?**
**Regress samples?**

QIIME 2's diversity analyses are available through the `q2-diversity` plugin, which supports computing alpha and beta diversity metrics, applying related statistical tests, and generating interactive visualizations.
We'll first apply the `core-metrics-phylogenetic` method, which rarefies a `FeatureTable[Frequency]` to a user-specified depth, computes several alpha and beta diversity metrics, and generates principle coordinates analysis (PCoA) plots using Emperor for each of the beta diversity metrics.
The metrics computed by default are:
-   Alpha diversity
    -   Shannon's diversity index (a quantitative measure of community
        richness)
    -   Observed Features (a qualitative measure of community richness)
    -   Faith's Phylogenetic Diversity (a qualitative measure of
        community richness that incorporates phylogenetic relationships
        between the features)
    -   Evenness (or Pielou's Evenness; a measure of community
        evenness)
-   Beta diversity
    -   Jaccard distance (a qualitative measure of community
        dissimilarity)
    -   Bray-Curtis distance (a quantitative measure of community
        dissimilarity)
    -   unweighted UniFrac distance (a qualitative measure of community
        dissimilarity that incorporates phylogenetic relationships
        between the features)
    -   weighted UniFrac distance (a quantitative measure of community
        dissimilarity that incorporates phylogenetic relationships
        between the features)

:::{tip}
You can find additional information on these and other metrics availability in QIIME 2 in [this excellent forum post by a QIIME 2 community member](https://forum.qiime2.org/t/alpha-and-beta-diversity-explanations-and-commands/2282).
When you're ready, we'd love to have your contributions on the Forum as well!
:::

An important parameter that needs to be provided to this script is `--p-sampling-depth`, which is the even sampling (i.e. rarefaction) depth.
Because most diversity metrics are sensitive to different sampling depths across different samples, this script will randomly subsample the counts from each sample to the value provided for this parameter.
For example, if you provide `--p-sampling-depth 500`, this step will subsample the counts in each sample without replacement so that each sample in the resulting table has a total count of 500.
If the total count for any sample(s) are smaller than this value, those samples will be dropped from the diversity analysis.
Choosing this value is tricky.
We recommend making your choice by reviewing the information presented in the `table_ms2.qzv` file that was created above.
Choose a value that is as high as possible (so you retain more sequences per sample) while excluding as few samples as possible.

:::{describe-usage}
core_metrics = use.action(
    use.UsageAction(plugin_id='boots',
                    action_id='core_metrics'),
    use.UsageInputs(table=kmer_table,
                    metadata=sample_metadata,
                    sampling_depth=62000,
                    n=10,
                    replacement=True,
                    alpha_average_method='median',
                    beta_average_method='medoid'),
    use.UsageOutputNames(
        resampled_tables='bootstrap_tables',
        alpha_diversities='bootstrap_alpha_diversities',
        distance_matrices='bootstrap_distance_matrices',
        pcoas='bootstrap_pcoas',
        emperor_plots='bootstrap_emperor_plots')
)

unweighted_dm = use.get_artifact_collection_member(
    'unweighted_dm', core_metrics.distance_matrices, 'jaccard_distance_matrix')
unweighted_pcoa = use.get_artifact_collection_member(
    'unweighted_pcoa', core_metrics.pcoas, 'jaccard')

weighted_dm = use.get_artifact_collection_member(
    'weighted_dm', core_metrics.distance_matrices, 'braycurtis_distance_matrix')
weighted_pcoa = use.get_artifact_collection_member(
    'weighted_pcoa', core_metrics.pcoas, 'braycurtis')

richness_vector = use.get_artifact_collection_member(
    'richness_vector', core_metrics.alpha_diversities, 'observed_features')
evenness_vector = use.get_artifact_collection_member(
    'evenness_vector', core_metrics.alpha_diversities, 'evenness')
:::

:::{tip} Question.
View the `table_ms2.qzv` QIIME 2 artifact, and in particular the *Interactive Sample Detail* tab in that visualization.
What value would you choose to pass for `--p-sampling-depth`?
How many samples will be excluded from your analysis based on this choice?
How many total sequences will you be analyzing in the `core-metrics-phylogenetic` command?
:::

Here we set the `--p-sampling-depth` parameter to ???.
The three samples that have fewer sequences will be dropped from the `core-metrics-phylogenetic` analyses and anything that uses these results.

::::{note}
In many Illumina runs you'll observe a few samples that have very low sequence counts.
You will typically want to exclude those from the analysis by choosing a larger value for the sampling depth at this stage.
::::

After computing diversity metrics, we can begin to explore the microbial composition of the samples in the context of the sample metadata.
This information is present in the [sample metadata](https://data.qiime2.org/2025.4/tutorials/gut-to-soil/sample_metadata) file that was downloaded earlier.

Finally, ordination is a popular approach for exploring microbial community composition in the context of sample metadata.

**Vizard PCoA 1 and 2 versus time.**
The PCoA results that were used in `core-metrics-phylogeny` are also available, making it easy to generate new PCoA-based visualizations.

:::{describe-usage}

unweighted_pcoa_as_md = use.view_as_metadata('unweighted_pcoa_as_md', unweighted_pcoa)
richness_as_md = use.view_as_metadata('richness_as_md', richness_vector)
unweighted_vizard_md = use.merge_metadata('unweighted_vizard_md', sample_metadata, unweighted_pcoa_as_md, richness_as_md)

use.action(
    use.UsageAction(plugin_id='vizard',
                    action_id='scatterplot_2d'),
    use.UsageInputs(metadata=unweighted_vizard_md),
    use.UsageOutputNames(visualization='unweighted_diversity_scatterplot'))
:::

:::{describe-usage}

weighted_pcoa_as_md = use.view_as_metadata('weighted_pcoa_as_md', weighted_pcoa)
weighted_vizard_md = use.merge_metadata('weighted_vizard_md', sample_metadata, weighted_pcoa_as_md, richness_as_md)

use.action(
    use.UsageAction(plugin_id='vizard',
                    action_id='scatterplot_2d'),
    use.UsageInputs(metadata=weighted_vizard_md),
    use.UsageOutputNames(visualization='weighted_diversity_scatterplot'))
:::



:::{exercise} Interpreting ordination plots.
When plotting PCoA axes 1 and 2 and coloring by SampleType, is the HEC more similar to the food compost or HE sample?

What sample type is the Microbe Mix most similar to?
The inside of the toilet pre-use?
The bulking material?

What other interesting relationships do you see when changing the x- and y-axes and sample coloring?

Which sample has the lowest microbiome richness?
:::

## Alpha rarefaction plotting

In this section we'll explore alpha diversity as a function of sampling depth using the `qiime diversity alpha-rarefaction` visualizer.
This visualizer computes one or more alpha diversity metrics at multiple sampling depths, in steps between 1 (optionally controlled with `--p-min-depth`) and the value provided as `--p-max-depth`.
At each sampling depth step, 10 rarefied tables will be generated, and the diversity metrics will be computed for all samples in the tables.
The number of iterations (rarefied tables computed at each sampling depth) can be controlled with `--p-iterations`.
Average diversity values will be plotted for each sample at each even sampling depth, and samples can be grouped based on metadata in the resulting visualization if sample metadata is provided with the `--m-metadata-file` parameter.

:::{describe-usage}
use.action(
    use.UsageAction(plugin_id='diversity',
                    action_id='alpha_rarefaction'),
    use.UsageInputs(table=kmer_table,
                    max_depth=62000,
                    metadata=sample_metadata),
    use.UsageOutputNames(visualization='alpha_rarefaction'))
:::

The visualization will have two plots.
The top plot is an alpha rarefaction plot, and is primarily used to determine if the richness of the samples has been fully observed or sequenced.
If the lines in the plot appear to "level out" (i.e., approach a slope of zero) at some sampling depth along the x-axis, that suggests that collecting additional sequences beyond that sampling depth would not be likely to result in the observation of additional features.
If the lines in a plot don't level out, this may be because the richness of the samples hasn't been fully observed yet (because too few sequences were collected), or it could be an indicator that a lot of sequencing error remains in the data (which is being mistaken for novel diversity).

The bottom plot in this visualization is important when grouping samples by metadata.
It illustrates the number of samples that remain in each group when the feature table is rarefied to each sampling depth.
If a given sampling depth `d` is larger than the total frequency of a sample `s` (i.e., the number of sequences that were obtained for sample `s`), it is not possible to compute the diversity metric for sample `s` at sampling depth `d`.
If many of the samples in a group have lower total frequencies than `d`, the average diversity presented for that group at `d` in the top plot will be unreliable because it will have been computed on relatively few samples.
When grouping samples by metadata, it is therefore essential to look at the bottom plot to ensure that the data presented in the top plot is reliable.

::::{note}

The value that you provide for `--p-max-depth` should be determined by reviewing the "Frequency per sample" information presented in the `table_ms2.qzv` file that was created above.
In general, choosing a value that is somewhere around the median frequency seems to work well, but you may want to increase that value if the lines in the resulting rarefaction plot don't appear to be leveling out, or decrease that value if you seem to be losing many of your samples due to low total frequencies closer to the minimum sampling depth than the maximum sampling depth.
::::

:::{tip} Question.
When grouping samples by "body-site" and viewing the alpha rarefaction plot for the "observed_features" metric, which body sites (if any) appear to exhibit sufficient diversity coverage (i.e., their rarefaction curves level off)?
How many sequence variants appear to be present in those body sites?
:::

:::{tip} Question.
When grouping samples by "body-site" and viewing the alpha rarefaction plot for the "observed_features" metric, the line for the "right palm" samples appears to level out at about 40, but then jumps to about 140.
What do you think is happening here?
(Hint: be sure to look at both the top and bottom plots.)
:::

## Taxonomic analysis

**[Weighted classifier](https://data.qiime2.org/classifiers/sklearn-1.4.2/gtdb/gtdb_diverse_weighted_classifier_r220.qza) is using too much memory for RTD builds.**
**Bump the annotation step up for an "upstream" versus "downstream" split in this document.**

In the next sections we'll begin to explore the taxonomic composition of the samples, and again relate that to sample metadata.
The first step in this process is to assign taxonomy to the sequences in our `FeatureData[Sequence]` QIIME 2 artifact.
We'll do that using a pre-trained Naive Bayes classifier and the `q2-feature-classifier` plugin.
This classifier was trained on the Greengenes 13_8 99% OTUs, where the sequences have been trimmed to only include 250 bases from the region of the 16S that was sequenced in this analysis (the V4 region, bound by the 515F/806R primer pair).
We'll apply this classifier to our sequences, and we can generate a visualization of the resulting mapping from sequence to taxonomy.

::::{note}
The taxonomic classifier used here is very specific to the sample preparation and sequencing protocols used for this study, and [Greengenes 13_8, which it is trained on, is an outdated reference database](https://forum.qiime2.org/t/introducing-greengenes2-2022-10/25291).
The reason we use it here is because the reference data is relatively small, so classification can be run on most modern computers with this classifier.

Generally speaking, we think it's not necessary to use a classifier that was trained on sequences that were trimmed based on the primers you sequenced with.
In practice we notice very minor differences, if any, relative to those trained on full-length sequences.
Far more impactful is the use of environment-weighted classifiers, as described in [Kaehler et al., (2019)](https://doi.org/10.1038/s41467-019-12669-6).

When you're ready to work on your own data, one of the choices you'll need to make is what classifier to use for your data.
You can find pre-trained classifiers on our [*Resources* page](https://resources.qiime2.org/), and [lots of discussion of this topic on the Forum](https://forum.qiime2.org/tag/taxonomy).
::::

:::{describe-usage}

def classifier_factory():
    from urllib import request
    from qiime2 import Artifact
    fp, _ = request.urlretrieve(
        'https://data.qiime2.org/classifiers/sklearn-1.4.2/greengenes/gg-13-8-99-515-806-nb-classifier.qza')

    return Artifact.load(fp)

classifier = use.init_artifact('16S-rRNA-classifier', classifier_factory)
:::

:::{describe-usage}

taxonomy, = use.action(
    use.UsageAction(plugin_id='feature_classifier',
                    action_id='classify_sklearn'),
    use.UsageInputs(classifier=classifier,
                    reads=rep_seqs_ms2),
    use.UsageOutputNames(classification='taxonomy'))

taxonomy_as_md = use.view_as_metadata('taxonomy_as_md', taxonomy)

use.action(
    use.UsageAction(plugin_id='metadata',
                    action_id='tabulate'),
    use.UsageInputs(input=taxonomy_as_md),
    use.UsageOutputNames(visualization='taxonomy'))
:::

:::{tip} Question.
Recall that our `rep-seqs.qzv` visualization allows you to easily BLAST the sequence associated with each feature against the NCBI nt database.
Using that visualization and the `taxonomy.qzv` visualization created here, compare the taxonomic assignments with the taxonomy of the best BLAST hit for a few features.
How similar are the assignments?
If they're dissimilar, at what *taxonomic level* do they begin to differ (e.g., species, genus, family, ...)?
:::

Next, we can view the taxonomic composition of our samples with interactive bar plots.
Generate those plots with the following command and then open the visualization.

:::{describe-usage}

use.action(
    use.UsageAction(plugin_id='taxa',
                    action_id='barplot'),
    use.UsageInputs(table=table_ms2,
                    taxonomy=taxonomy,
                    metadata=sample_metadata),
    use.UsageOutputNames(visualization='taxa_bar_plots'))
:::

:::{tip} Question.
Visualize the samples at *Level 2* (which corresponds to the phylum level in this analysis), and then sort the samples by `body-site`, then by `subject`, and then by `days-since-experiment-start`.
What are the dominant phyla in each in `body-site`? Do you observe any consistent change across the two subjects between `days-since-experiment-start` `0` and the later timepoints?
:::

<!--
## Differential abundance testing with ANCOM-BC

ANCOM-BC can be applied to identify features that are differentially abundant (i.e. present in different abundances) across sample groups.
As with any bioinformatics method, you should be aware of the assumptions and limitations of ANCOM-BC before using it.
We recommend reviewing the [ANCOM-BC paper](https://pubmed.ncbi.nlm.nih.gov/32665548/) before using this method.

::::{note}
Accurately identifying features that are differentially abundant across sample types in microbiome data is a challenging problem and an open area of research.
There is one QIIME 2 plugin that can be used for this: `q2-composition` (used in this section).
In addition to the methods contained in this plugin, new approaches for differential abundance testing are regularly introduced, and it's worth assessing the current state of the field when performing differential abundance testing to see if there are new methods that might be useful for your data.
::::

ANCOM-BC is a compositionally-aware linear regression model that allows for testing differentially abundant features across groups while also implementing bias correction, and is currently implemented in the `q2-composition` plugin.

Because we expect a lot of features to change in abundance across body sites, in this tutorial we'll filter our full feature table to only contain gut samples.
We'll then apply ANCOM-BC to determine which, if any, sequence variants and genera are differentially abundant across the gut samples of our two subjects.

We'll start by creating a feature table that contains only the gut samples.
(To learn more about filtering, see the `Filtering Data <filtering>`{.interpreted-text role="doc"} tutorial.)

:::{describe-usage}

gut_table, = use.action(
    use.UsageAction(plugin_id='feature_table',
                    action_id='filter_samples'),
    use.UsageInputs(table=table,
                    metadata=sample_metadata,
                    where='[body-site]="gut"'),
    use.UsageOutputNames(filtered_table='gut_table'))
:::

ANCOM-BC operates on a `FeatureTable[Frequency]` QIIME 2 artifact.
We can run ANCOM-BC on the subject column to determine what features differ in abundance across gut samples of the two subjects.

:::{describe-usage}

ancombc_subject, = use.action(
    use.UsageAction(plugin_id='composition',
                    action_id='ancombc'),
    use.UsageInputs(table=gut_table,
                    metadata=sample_metadata,
                    formula='subject'),
    use.UsageOutputNames(differentials='ancombc_subject'))

use.action(
    use.UsageAction(plugin_id='composition',
                    action_id='da_barplot'),
    use.UsageInputs(data=ancombc_subject,
                    significance_threshold=0.001),
    use.UsageOutputNames(visualization='da_barplot_subject'))
:::

:::{tip} Question.
1. Which ASV is most enriched, relative to the reference?
   Which is most depleted?
2. What would you expect to change if the `reference-level` was changed from `subject-1` (the default) to `subject-2`?
:::

We're also often interested in performing a differential abundance test at a specific taxonomic level.
To do this, we can collapse the features in our `FeatureTable[Frequency]` at the taxonomic level of interest, and then re-run the above steps.
In this tutorial, we collapse our feature table at the genus level (i.e. level 6 of the Greengenes taxonomy).

:::{describe-usage}
l6_gut_table, = use.action(
    use.UsageAction(plugin_id='taxa',
                    action_id='collapse'),
    use.UsageInputs(table=gut_table,
                    taxonomy=taxonomy,
                    level=6),
    use.UsageOutputNames(collapsed_table='gut_table_l6'))

l6_ancombc_subject, = use.action(
    use.UsageAction(plugin_id='composition',
                    action_id='ancombc'),
    use.UsageInputs(table=l6_gut_table,
                    metadata=sample_metadata,
                    formula='subject'),
    use.UsageOutputNames(differentials='l6_ancombc_subject'))

use.action(
    use.UsageAction(plugin_id='composition',
                    action_id='da_barplot'),
    use.UsageInputs(data=l6_ancombc_subject,
                    significance_threshold=0.001),
    use.UsageOutputNames(visualization='l6_da_barplot_subject'))
:::

:::{tip} Question.
1. Which genus is most enriched?
   Which is most depleted?
2. Do we see more differentially abundant features in the `da-barplot-subject.qzv` visualization, or in the `l6-da-barplot-subject.qzv` visualization?
   Why might you expect this?
:::

-->

## Next steps

You might next want to try to adapt the commands presented in this tutorial to your own data, adjusting parameter settings and metadata column headers as is relevant.
If you need help, head over to the [QIIME 2 Forum](https://forum.qiime2.org).


