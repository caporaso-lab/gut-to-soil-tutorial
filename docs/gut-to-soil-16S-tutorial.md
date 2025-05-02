(gut-to-soil-16S-tutorial)=
# Gut-to-soil axis 16S rRNA analysis tutorial 💩🌱

:::{warning}
This document is a work in progress as of 17 April 2025, and some aspects of the workflow are still in development.
Commands/urls/text may be unreliable while in development.
🚜
:::

:::{note}
This document was built with its own conda environment that includes the amplicon distribution and the standalone plugins [q2-boots](https://doi.org/10.12688/f1000research.156295.1) and [q2-kmerizer](https://doi.org/10.1128/msystems.01550-24).
You can download the environment file that was used from the download link on the top-right of this article.
:::

## Background

In this tutorial you'll learn an end-to-end microbiome marker-gene data science workflow, building on data presented in [Meilander *et al.* (2024): Upcycling Human Excrement: The Gut Microbiome to Soil Microbiome Axis](https://doi.org/10.48550/arXiv.2411.04148).
The data used here is a subset (a single sequencing run) of that generated for the paper, specifically selected so that this tutorial can be run quickly on a personal computer.
The full data set for the paper can be found in [the study's Artifact Repository](https://doi.org/10.5281/zenodo.13887456).
In the final step (**not yet written, as of 17 April 2025**), you'll learn how to adapt the workflow for use in analyzing your own data using [Provenance Replay](https://doi.org/10.1371/journal.pcbi.1011676).

The data used in this tutorial was generated using the [Earth Microbiome Project protocol](https://doi.org/10.1038/ismej.2012.8).
Specifically, the hypervariable region 4 (V4) of the 16S rRNA gene was amplified using the F515-R806 primers - a broad-coverage primer pair for Bacteria that also amplifies some Archaea.
Paired-end sequencing was performed on an Illumina MiSeq.
Full details are presented in [Meilander *et al.* (2024)](https://doi.org/10.48550/arXiv.2411.04148).

::::{margin}
:::{tip} Why study human excrement composting?
Get our take on that [here](https://gut-to-soil-tutorial.readthedocs.io/en/latest/why.html).
:::
::::

(gut-to-soil-tutorial:sample-metadata)=
## Sample metadata

Before starting the analysis, explore the sample metadata to familiarize yourself with the samples used in this study.
The following command will download the sample metadata as tab-separated text and save it in the file `sample-metadata.tsv`.
This `sample-metadata.tsv` file is used throughout the rest of the tutorial.

::::{margin}
:::{tip}
To learn more about metadata in QIIME 2, including how it should be formatted, refer to [*Using QIIME 2*'s Metadata file format](https://use.qiime2.org/en/latest/references/metadata.html).

To learn more about metadata standards, you can refer to [Chloe Herman's video on this topic](https://www.youtube.com/watch?v=erklD1bofzE), which was developed in collaboration with the [National Microbiome Data Collaborative (NMDC)](https://microbiomedata.org/).
:::
::::

:::{describe-usage}
:scope: gut-to-soil

sample_metadata = use.init_metadata_from_url(
   'sample-metadata',
   'https://www.dropbox.com/scl/fi/irosimbb1aud1aa7frzxf/sample-metadata.tsv?rlkey=f45jpxzajjz9xx9vpvfnf1zjx&st=nahafuvy&dl=1')
:::

QIIME 2's [metadata plugin](xref:q2doc-amplicon-target#q2-plugin-metadata) provides a [Visualizer called `tabulate`](xref:q2doc-amplicon-target#q2-action-metadata-tabulate) that generates a convenient view of a sample metadata file.
Let's run this, and then we'll look at the result.
Here's the first QIIME 2 command that you should run in this tutorial:

(sample-metadata-tabulate-viz)=
:::{describe-usage}
use.action(
  use.UsageAction(plugin_id='metadata',
                  action_id='tabulate'),
  use.UsageInputs(input=sample_metadata),
  use.UsageOutputNames(visualization='sample_metadata')
)
:::

::::{margin}
:::{tip}
You can learn more about viewing Visualizations, including alternatives to QIIME 2 View if you can't use that for any reason, [in *Using QIIME 2*](https://use.qiime2.org/en/latest/how-to-guides/view-visualizations.html).
:::
::::

This will generate a QIIME 2 [Visualization](https://amplicon-docs.readthedocs.io/en/latest/explanations/getting-started.html#getting-started-artifacts-and-visualizations).
Visualizations can be viewed by loading them with [QIIME 2 View](https://view.qiime2.org).
Navigate to QIIME 2 View, and drag and drop the visualization that was created to view it.

## Access already-imported QIIME 2 data

This tutorial begins with paired-end read sequencing data that has already been demultiplexed and imported into a QIIME 2 Artifact.
Because sequence data can be delivered to you in many different forms, it's not possible to cover the varieties here.
Instead we refer you to [*How to import data for use with QIIME 2*](https://amplicon-docs.readthedocs.io/en/latest/how-to-guides/how-to-import.html) to learn how to import your data.
If you want to learn why importing is necessary, refer to [Why importing is necessary](https://amplicon-docs.readthedocs.io/en/latest/explanations/why-importing.html).

::::{margin}
:::{tip}
You should ask your sequencing center to provide data already demultiplexed.
In some cases, they may be able to provide a "QIIME 2 demux artifacts".
If not, ask them to provide a [fastq manifest file](https://amplicon-docs.readthedocs.io/en/latest/how-to-guides/how-to-import.html#import-fastq-manifest) for your sequencing data.
It should be easy for them to generate.
:::
::::

:::{describe-usage}
:scope: gut-to-soil

demux = use.init_artifact_from_url(
   'demux',
   'https://www.dropbox.com/scl/fi/73f6rmcq7lelug5qbuhl6/demux-10p.qza?rlkey=s0hoh326fes3z2usvgs62tv3c&st=caz1avkn&dl=1')
:::

## Summarize demultiplexed sequences

When you have demultiplexed sequence data, the next step is typically to generate a visual summary of it.
This allows you to determine how many sequences were obtained per sample, and also to get a summary of the distribution of sequence qualities at each position in your sequence data.

(demux-summary-viz)=
:::{describe-usage}

use.action(
    use.UsageAction(plugin_id='demux',
                    action_id='summarize'),
    use.UsageInputs(data=demux),
    use.UsageOutputNames(visualization='demux'))
:::

:::{exercise} Exploring the "demux summary".
How many samples are represented in this sequencing data?
What is the median number of sequence reads obtained per sample?
What is the median quality score at position 200 of the forward reads?
:::

## Upstream data analysis

Generally, the term "upstream" is used to refer to data analysis pre-feature-asv_table, and "downstream" is used to refer to data analysis post-feature-asv_table.
Let's jump into our upstream analysis.

### Sequence quality control and feature table construction

QIIME 2 plugins are available for several quality control methods, including [DADA2](https://doi.org/10.1038/nmeth.3869), [Deblur](https://doi.org/10.1128/msystems.00191-16), and [basic quality-score-based filtering](https://doi.org/10.1038/nmeth.2276).
In this tutorial we present this step using [DADA2](https://www.ncbi.nlm.nih.gov/pubmed/27214047).
The result of this method will be a `FeatureTable[Frequency]` QIIME 2 artifact, which contains counts (frequencies) of each unique sequence in each sample in the dataset, and a `FeatureData[Sequence]` QIIME 2 artifact, which maps feature identifiers in the `FeatureTable` to the sequences they represent.

DADA2 is a pipeline for detecting and correcting (where possible) Illumina amplicon sequence data.
As implemented in the [dada2 plugin](xref:q2doc-amplicon-target#q2-plugin-dada2), this quality control process will additionally filter any phiX reads (commonly present in marker gene Illumina sequence data) that are identified in the sequencing data, filter chimeric sequences, and merge paired end reads.

The [`denoise-paired` action](xref:q2doc-amplicon-target#q2-action-dada2-denoise-paired), which we'll use here, requires four parameters that are used in quality filtering:
- `trim-left-f a`, which trims off the first `a` bases of each forward read
- `trunc-len-f b` which truncates each forward read at position `b`
- `trim-left-r c`, which trims off the first `c` bases of each forward read
- `trunc-len-r d` which truncates each forward read at position `d`
This allows the user to remove low quality regions of the sequences.
To determine what values to pass for these parameters, you should review the *Interactive Quality Plot* tab in the [`demux.qzv`](#demux-summary-viz) file that was generated above.

:::{exercise} Choosing "trim" and "trunc" parameter values.
:label: dada2-trim-trunc
Based on the plots you see in the [`demux.qzv`](#demux-summary-viz) file that was generated above, what values would you choose for `trim-left-f` and `trunc-len-f` in this case?
What about `trim-left-r` and `trunc-len-r`?
:::

:::{solution} dada2-trim-trunc
:class: dropdown
The quality of the initial bases seems to be high, so I choose not to trim any bases from the beginning of the sequences.
The quality also seems good all the way out to the end, though maybe dropping off after 250 bases.
I'll therefore truncate at 250.
I'll keep these values the same for both the forward and reverse reads, though that is not a requirement.
:::

Now run your DADA2 command.
This step may take up to 10 minutes to complete - it's the longest running step in this tutorial.

:::{describe-usage}

asv_seqs, asv_table, stats = use.action(
    use.UsageAction(plugin_id='dada2',
                    action_id='denoise_paired'),
    use.UsageInputs(demultiplexed_seqs=demux,
                    trim_left_f=0,
                    trunc_len_f=250,
                    trim_left_r=0,
                    trunc_len_r=250),
    use.UsageOutputNames(representative_sequences='asv_seqs',
                         table='asv_table',
                         denoising_stats='stats'))
:::

One of the outputs created by DADA2 is a summary of the denoising run.
That is generated as an [Artifact](https://amplicon-docs.readthedocs.io/en/latest/explanations/getting-started.html#getting-started-artifacts-and-visualizations), so can't be viewed directly.
However this is one of many QIIME 2 types that can be [viewed as Metadata](https://use.qiime2.org/en/latest/how-to-guides/artifacts-as-metadata.html) - a very powerful concept that we'll use again later in this tutorial.
Learning to view artifacts as Metadata creates nearly infinite possibilities for how you can explore your microbiome data with QIIME 2.

Here, we'll again use the [metadata plugins `tabulate` visualizer](xref:q2doc-amplicon-target#q2-action-metadata-tabulate), but this time we'll apply it to the DADA2 statistics.

:::{describe-usage}
stats_as_md = use.view_as_metadata('stats_as_md', stats)

use.action(
    use.UsageAction(plugin_id='metadata',
                    action_id='tabulate'),
    use.UsageInputs(input=stats_as_md),
    use.UsageOutputNames(visualization='stats'))
:::

:::{exercise} Exploring the DADA2 denoising statistics.
Which three samples had the smallest percentage of input reads passing the quality filter?
Refer back to your [tabulated metadata](#sample-metadata-tabulate-viz): what do you know about those samples (e.g., what values do they have in the `SampleType` column)?
(Hint: the `sample-id` is the key that connects data across these two visualizations.)
:::

:::{exercise} Merging metadata.
:label: merge-metadata
If it's annoying to switch back and forth between visualizations to answer the question in the last exercise, you can create a combined tabulation of the metadata.
Try to do that by adapting the instructions in [*How to merge metadata*](https://use.qiime2.org/en/latest/how-to-guides/merge-metadata.html).

This is also useful if you want to create a large tabular summary describing your samples following analysis, as you can include as many different metadata objects as you'd like in these summaries.
:::

::::{solution} merge-metadata
:class: dropdown

:::{describe-usage}
sample_metadata_and_dada2_stats_md = use.merge_metadata('sample_metadata_and_dada2_stats_md', sample_metadata, stats_as_md)

use.action(
    use.UsageAction(plugin_id='metadata',
                    action_id='tabulate'),
    use.UsageInputs(input=sample_metadata_and_dada2_stats_md),
    use.UsageOutputNames(visualization='sample_metadata_w_dada2_stats'))
:::
::::

### Feature table and feature data summaries

After DADA2 completes, you'll want to explore the resulting data.
You can do this using the following two commands, which will create visual summaries of the data.
The [`feature-table summarize` action](xref:q2doc-amplicon-target#q2-action-feature-table-summarize) command will give you information on how many sequences are associated with each sample and with each feature, histograms of those distributions, and some related summary statistics.

:::{describe-usage}
_, _, asv_frequencies = use.action(
    use.UsageAction(plugin_id='feature_table',
                    action_id='summarize_plus'),
    use.UsageInputs(table=asv_table,
                    metadata=sample_metadata),
    use.UsageOutputNames(summary='asv_table',
                         sample_frequencies='sample_frequencies',
                         feature_frequencies='asv_frequencies'))
:::

:::{exercise} Exploring the feature table summary -- part 1.
What is the total number of sequences represented in the feature table?
What is the identifier of the feature that is observed the most number of times (i.e., has the highest frequency)?
:::

The [`feature-table tabulate-seqs` action](xref:q2doc-amplicon-target#q2-action-feature-table-tabulate-seqs) command will provide a mapping of feature IDs to sequences, and provide links to easily BLAST each sequence against the NCBI nt database.
We can also include the feature frequency information in this visualization by passing it as metadata, similar to how we [merged metadata](#merge-metadata) in the exercise above.
In this case, however, we're looking a *feature metadata*, as opposed to *sample metadata*.
As far as QIIME 2 is concerned, there is no difference between these two - in our case, it'll only be the identifiers that differ.

This visualization will be very useful later in the tutorial, when you want to learn more about specific features that are important in the data set.

:::{describe-usage}

asv_frequencies_as_md = use.view_as_metadata('asv_frequencies_md',
                                                 asv_frequencies)

use.action(
    use.UsageAction(plugin_id='feature_table',
                    action_id='tabulate_seqs'),
    use.UsageInputs(data=asv_seqs,
                    metadata=asv_frequencies_as_md),
    use.UsageOutputNames(visualization='asv_seqs'),
)
:::

:::{exercise} Exploring feature data -- part 1.
What is the taxonomy associated with the most frequently observed feature, based on a BLAST search?
:::

### Filtering features from a feature table

If you review the tabulated feature sequences, or the feature detail table of the feature table summary, you'll notice that there are many sequences that are observed in only a single sample.
Let's filter those out to reduce the number of sequences we're working with - this will speed up several slower steps that are coming up.

This is a two-step process.
First we filter our feature table, and then we use the new feature table to filter our sequences to only the ones that are contained in the new table.

::::{margin}
:::{tip}
I include `_ms2` in the new output names here to remind us that these are filtered to those features in a **m**inimum of **2** **s**amples.
Generally speaking, file names are a convenient place to store information like this, but they're unreliable.
File names can easily be changed, and therefore could be modified to contain inaccurate information about the data they contain.

Luckily, [QIIME 2's provenance tracking system](https://amplicon-docs.readthedocs.io/en/latest/explanations/getting-started.html#getting-started-provenance) records all of the information that we need about how results were generated.
We're therefore free to include information like this in file names if it's helpful for us, but we shouldn't ever rely on the file names.
If in doubt -- and always be in doubt 🕵️‍♀️ -- refer to the data provenance.

QIIME 2's data provenance is your source for the truth about how a result was created.
:::
::::

:::{describe-usage}
asv_table_ms2, = use.action(
    use.UsageAction(plugin_id='feature_table',
                    action_id='filter_features'),
    use.UsageInputs(table=asv_table,
                    min_samples=2),
    use.UsageOutputNames(filtered_table='asv_table_ms2'),
)
:::

:::{describe-usage}
asv_seqs_ms2, = use.action(
    use.UsageAction(plugin_id='feature_table',
                    action_id='filter_seqs'),
    use.UsageInputs(data=asv_seqs,
                    table=asv_table_ms2),
    use.UsageOutputNames(filtered_data='asv_seqs_ms2'),
)
:::

:::{exercise} Exploring the feature table summary -- part 2.
:label: summarize-asv-table-ms2
Now that you have a second (filtered) feature table, create your own command to summarize it, like we did for the original feature table.
How many features did we lose as a result of this filter?
How many total sequences did we lose?
:::

::::{solution} summarize-asv-table-ms2
:class: dropdown
Here's the command you would use:

:::{describe-usage}
_, _, asv_frequencies_ms2 = use.action(
    use.UsageAction(plugin_id='feature_table',
                    action_id='summarize_plus'),
    use.UsageInputs(table=asv_table_ms2,
                    metadata=sample_metadata),
    use.UsageOutputNames(summary='asv_table_ms2',
                         sample_frequencies='sample_frequencies_ms2',
                         feature_frequencies='asv_frequencies_ms2'))
:::

Be sure to run this as we're going to use one of the results below.
::::

### Taxonomic annotation

Before we complete our upstream analysis steps, we'll generate taxonomic annotations for our sequences using the [feature-classifier plugin](xref:q2doc-amplicon-target#q2-plugin-feature-classifier).

::::{note}
The taxonomic classifier used here is very specific to the sample preparation and sequencing protocols used for this study, and Greengenes 13_8, which it is trained on, [is an outdated reference database](https://forum.qiime2.org/t/introducing-greengenes2-2022-10/25291).
The reason we use it here is because the reference data is relatively small, so classification can be run on most modern computers with this classifier[^build-requirements-exceed-resources].

When you're ready to work on your own data, one of the choices you'll need to make is what classifier to use for your data.
You can find pre-trained classifiers the QIIME 2 Library [*Resources* page](https://library.qiime2.org/data-resources).
If you don't find a classifier that will work for you there, you may be able to [find one on the Forum](https://forum.qiime2.org/tag/taxonomy), or you can learn how to train your own [by referencing the RESCRIPt documentation](https://github.com/bokulich-lab/RESCRIPt).

We strongly recommend the use of environment-weighted classifiers, as described in [Kaehler et al., (2019)](https://doi.org/10.1038/s41467-019-12669-6), and you can find builds of these on the [*Resources* page](https://library.qiime2.org/data-resources).
::::

First, we'll download a pre-trained classifier artifact.

:::{describe-usage}

def classifier_factory():
    from urllib import request
    from qiime2 import Artifact
    fp, _ = request.urlretrieve(
        'https://data.qiime2.org/classifiers/sklearn-1.4.2/greengenes/gg-13-8-99-515-806-nb-classifier.qza')

    return Artifact.load(fp)

classifier = use.init_artifact('suboptimal-16S-rRNA-classifier', classifier_factory)
:::

Then, we'll apply it to our sequences using [`classify-sklearn`](xref:q2doc-amplicon-target#q2-action-feature-classifier-classify-sklearn).

:::{describe-usage}
taxonomy, = use.action(
    use.UsageAction(plugin_id='feature_classifier',
                    action_id='classify_sklearn'),
    use.UsageInputs(classifier=classifier,
                    reads=asv_seqs_ms2),
    use.UsageOutputNames(classification='taxonomy'))
:::

Then, to get an initial look at our taxonomic classifications, let's integrate taxonomy in the sequence summary, like the one we generated above.

::::{margin}
:::{tip}
If you want to compare taxonomic annotations achieved with different classifiers, you can do that with the [`feature-table tabulate-seqs` action](xref:q2doc-amplicon-target#q2-action-feature-table-tabulate-seqs) by passing in multiple `FeatureData[Taxonomy]` artifacts.
See an example of what that result might look like [here](https://view.qiime2.org/visualization/?src=https://zenodo.org/api/records/13887457/files/asv-seqs-ms10.qzv/content).

While you have that visualization loaded, take a look at the data provenance.
The complexity of that data provenance should give you an idea of why it's helpful to have the computer record all of this information, rather than trying to embed it all in file names or keep track of it in your written notes.

What was used as the DADA2 trim and trunc parameters for the data leading to this visualization?
(Hint: use the provenance search feature).
:::
::::

:::{describe-usage}

asv_frequencies_ms2_as_md = use.view_as_metadata('asv_frequencies',
                                                 asv_frequencies_ms2)

taxonomy_collection = use.construct_artifact_collection(
    'taxonomy_collection', {'Greengenes_13_8': taxonomy}
)

use.action(
    use.UsageAction(plugin_id='feature_table',
                    action_id='tabulate_seqs'),
    use.UsageInputs(data=asv_seqs_ms2,
                    taxonomy=taxonomy_collection,
                    metadata=asv_frequencies_ms2_as_md),
    use.UsageOutputNames(visualization='asv_seqs_ms2'),
)
:::

:::{exercise} Exploring feature data -- part 1.
What is the taxonomy associated with the most frequently observed feature, based on a BLAST search?
How does that compare to the taxonomy assigned by our feature classifier?

In general, we tend to trust the results of our feature classifier over those generated by BLAST, though the BLAST results are good for a quick look or a sanity check.
The reason for this is that the reference databases used with our feature classifiers tend to be specifically designed for this purpose and in some cases all of the sequences included are vetted.
The BLAST databases can contain mis-annotations that may negatively impact the classifications.
:::

:::{tip} Question.
Recall that our `asv-seqs.qzv` visualization allows you to easily BLAST the sequence associated with each feature against the NCBI nt database.
Using that visualization and the `taxonomy.qzv` visualization created here, compare the taxonomic assignments with the taxonomy of the best BLAST hit for a few features.
How similar are the assignments?
If they're dissimilar, at what *taxonomic level* do they begin to differ (e.g., species, genus, family, ...)?
:::

(phylogenetic-tree-building)=
### Building a tree for phylogenetic diversity calculations

QIIME supports several phylogenetic diversity metrics, including Faith's Phylogenetic Diversity and weighted and unweighted UniFrac.
In addition to counts of features per sample (i.e., the data in the `FeatureTable[Frequency]` QIIME 2 artifact), these metrics require a rooted phylogenetic tree relating the features to one another.
The amplicon distribution offers a few ways to build these trees, including a reference-based approach in the [fragment-insertion plugin](xref:q2doc-amplicon-target#q2-plugin-fragment-insertion) and *de novo* (i.e., reference-free) approaches in the [phylogeny plugin](xref:q2doc-amplicon-target#q2-plugin-phylogeny).

The reference based approach, by default, is specific to 16S rRNA marker gene analysis.
We could use that here, but the runtime is too long for our documentation.[^build-requirements-exceed-resources]
If you'd like to see this demonstrated, you can refer to the [*Parkinson's Mouse* tutorial](https://docs.qiime2.org/2024.10/tutorials/pd-mice/).

The *de novo* approach is known to generate low quality trees when very short sequences are used as input, but can be used with any phylogenetically informative marker gene (not just 16S).
If you'd like to see this demonstrated, you can refer to the [*Moving Pictures* tutorial](https://amplicon-docs.readthedocs.io/en/latest/tutorials/moving-pictures.html#generate-a-tree-for-phylogenetic-diversity-analyses).

For those reasons, we're going to skip building phylogenetic trees and instead use an analog of phylogenetic diversity metrics here.

## Downstream data analysis

As mentioned above, we tend to think of "downstream" analysis as beginning with a feature table, taxonomic annotation of our features, and optionally a phylogenetic tree.
Now that we have those (with the exception of the tree, [which we won't use here](#phylogenetic-tree-building)), let's jump in.
This is where it starts to get fun! ⛷️

(gut-to-soil:kmer-based-diversity-analysis)=
### Kmer-based diversity analysis

As mentioned above, we're going to skip building phylogenetic trees and instead use an analog of phylogenetic diversity metrics here.
This will use two [stand-alone QIIME 2 plugins](xref:q2doc-amplicon-target#term-stand-alone-plugin), [q2-boots](https://library.qiime2.org/plugins/caporaso-lab/q2-boots) and [q2-kmerizer](https://library.qiime2.org/plugins/bokulich-lab/q2-kmerizer), which are integrated through the [`kmer-diversity`](xref:q2-boots-target#q2-action-boots-kmer-diversity) action in q2-boots.
To learn about kmerization of features and how this relates to phylogenetic diversity metrics, read the [q2-kmerizer](https://doi.org/10.1128/msystems.01550-24) paper.
[q2-boots](https://library.qiime2.org/plugins/caporaso-lab/q2-boots) provides actions that mirror the interface of the diversity metric calculation actions in the diversity plugin, but generates more robust results because it integrates rarefaction and/or bootstrapping.
You can learn more about this in the [q2-boots paper](https://doi.org/10.12688/f1000research.156295.1).

::::{margin}
:::{note}
When using [q2-kmerizer](https://library.qiime2.org/plugins/bokulich-lab/q2-kmerizer), normalization should be done at the ASV level before kmerization, not at the kmer level.
This is automated by the [`kmer-diversity`](xref:q2-boots-target#q2-action-boots-kmer-diversity) action that we're using here.
This is because you are trying to normalize by sequencing depth, not kmer complexity.
In the end, the difference should not be huge but this distinction could be important for some metrics.
:::
::::

[`kmer-diversity`](xref:q2-boots-target#q2-action-boots-kmer-diversity) is a [`qiime2.Pipeline`](xref:q2doc-amplicon-target#term-pipeline): a type of action that links multiple other QIIME 2 actions together for convenience.
As such, it does a lot of work.
Here are the steps that it takes:
1. Resample the input feature table (i.e., the ASV table) to contain exactly `sampling-depth` sequences per sample, either by bootstrapping or rarefaction, `n` times.
   Samples with fewer than `sampling-depth` sequences will be removed from the feature table and not included in the subsequent steps.
   This will result in `n` feature tables.
2. For each feature table resulting from step 1, using the input sequences (i.e., the ASV sequences), kmerize all sequences into kmers of length `kmer-size`. [^iab-database-searching]
   Use this information to create one kmer table per resampled feature table.
   This will result in `n` feature tables, where the features are kmers (instead of ASVs, as in the input feature table).
3. Compute the user-requested alpha- and beta-diversity metrics on each of the kmer tables resulting from step 2. [^forum-diversity-metrics]
   The metrics computed by default are:
    * Alpha diversity
      * Shannon's diversity index (a quantitative measure of community richness)
      * Observed Features (a qualitative measure of community richness)
      * Evenness (i.e., Pielou's Evenness; a measure of community evenness)
    * Beta diversity
      * Jaccard distance (a qualitative measure of community dissimilarity)
      * Bray-Curtis distance (a quantitative measure of community dissimilarity)
4. For each diversity metric, average the results computed across the `n` kmer tables.
   These results can be used in subsequent analysis steps (e.g., ordination, statistical modeling, machine learning).
5. Perform PCoA ordination on the averaged beta diversity distance matrices resulting from Step 4. [^iab-machine-learning]
6. Generate an interactive [q2-vizard scatter plot](xref:q2doc-amplicon-target#q2-action-vizard-scatterplot-2d) that contains all user-provided sample metadata, all averaged alpha diversity metrics, and the first three ordination axes for all PCoA matrices computed in step 5.

::::{margin}
:::{note}
[q2-kmerizer](https://library.qiime2.org/plugins/bokulich-lab/q2-kmerizer) decomposes each sequence into its constituent kmers.
This should be carefully considered when interpreting alpha-diversity in particular, as the number of observed features (for example) would correspond to the number of unique kmers observed in a sample (representing the genetic diversity), not the number of unique sequences or taxa.
For more information, read the [q2-kmerizer](https://doi.org/10.1128/msystems.01550-24) paper.
:::
::::

A key parameter that needs to be provided to [`kmer-diversity`](xref:q2-boots-target#q2-action-boots-kmer-diversity) is `sampling-depth`, which is the even sampling (i.e., bootstrapping or rarefaction) depth.
Because most diversity metrics are sensitive to different sampling depths (i.e., sequence counts) across different samples, the tables are randomly subsampled such that the total frequency for each sample is the user-specified sampling depth.
For example, if you set `sampling-depth=500`, this step will subsample the sequence counts in each sample so that each sample in the resulting table has a total frequency of 500.
If the total frequency (i.e., the number of sequences observed) for any sample(s) is smaller than this value, those samples will be dropped from the diversity analysis.
Choosing this value is tricky.
We recommend making your choice by reviewing the information presented in the `asv-table-ms2.qzv` file that [you created above](#summarize-asv-table-ms2).

::::{margin}
:::{tip}
[This video](https://www.youtube.com/watch?v=W1W9jcrMFHs) provides guidance on choosing an even sampling depth.
:::
::::

:::{exercise} Choosing a sampling depth.
View the `asv-table-ms2.qzv` QIIME 2 artifact, and in particular the *Interactive Sample Detail* tab in that visualization.
What value would you choose to pass for `sampling-depth`?
How many samples will be excluded from your analysis based on this choice?
How many total sequences will you be analyzing in the `kmer-diversity` command?
:::

::::{margin}
:::{tip}
In many Illumina runs you'll observe a few samples that have very low sequence counts.
You will typically want to exclude those from the analysis by choosing a larger value for the sampling depth at this stage.
:::
::::

I'm going to choose values that is around the first quartile of the sample total frequencies.

:::{describe-usage}
use.action(
    use.UsageAction(plugin_id='boots',
                    action_id='kmer_diversity'),
    use.UsageInputs(table=asv_table_ms2,
                    sequences=asv_seqs_ms2,
                    metadata=sample_metadata,
                    sampling_depth=96,
                    n=10,
                    replacement=True,
                    alpha_average_method='median',
                    beta_average_method='medoid'),
    use.UsageOutputNames(
        resampled_tables='bootstrap_tables',
        kmer_tables='kmer_tables',
        alpha_diversities='bootstrap_alpha_diversities',
        distance_matrices='bootstrap_distance_matrices',
        pcoas='bootstrap_pcoas',
        scatter_plot='kmer_diversity_scatter_plot')
)
:::

::::{margin}
:::{warning}
Remember that we're working with a small subset of a full sequencing run in this tutorial, to keep the runtime short for the tutorial.
As a result the value used for `sampling-depth` here is very low.
Often this might be closer to 10,000 for an Illumina MiSeq run (as of 2025), but this is highly dependent on the sequencing run, the number of samples included, and other factors.
:::
::::

After computing diversity metrics, we can begin to explore the microbial composition of the samples in the context of the sample metadata.
You can review the sample metadata using one of the tabulated views of this file that [we created above](#sample-metadata-tabulate-viz).

:::{exercise} Learn to use the interactive scatter plot.
Open the scatter plot that was generated by the previous command, and plot the first two ordination axes computed from the Bray-Curtis distances.

Experiment with plotting different information on the x and y axes.
Change the *colorBy* and *colorPalette* selections.
Which of the metadata categories results in samples grouping most by color?
:::

:::{exercise} Interpreting ordination plots.
When plotting Bray-Curtis PCoA axes 1 and 2 and coloring by SampleType, is the HEC more similar to the food compost or HE samples?

What sample type is the Microbe Mix most similar to?
The inside of the toilet pre-use?
The bulking material?

What other interesting relationships do you see when changing the x- and y-axes and sample coloring?
Which diversity axes (PCoA or alpha diversity results) appear to be correlated with one another?

Which sample has the lowest microbiome richness?
:::

### Alpha rarefaction plotting

In this section we'll explore alpha diversity as a function of sampling depth using the[`alpha-rarefaction` action](xref:q2doc-amplicon-target#q2-action-diversity-alpha-rarefaction).
This visualizer computes one or more alpha diversity metrics at multiple sampling depths, in steps between 1 (optionally controlled with `min-depth`) and the value provided as `max-depth`.
At each sampling depth step, 10 rarefied tables will be generated, and the diversity metrics will be computed for all samples in the tables.
The number of iterations (rarefied tables computed at each sampling depth) can be controlled with the `iterations` parameter.
Average diversity values will be plotted for each sample at each even sampling depth, and samples can be grouped based on metadata in the resulting visualization if sample metadata is provided.

The value that you provide for `max-depth` should be determined by reviewing the "Frequency per sample" information presented in the `asv-table-ms2.qzv` file.
In general, choosing a value that is somewhere around the median frequency seems to work well, but you may want to increase that value if the lines in the resulting rarefaction plot don't appear to be leveling out, or decrease that value if you seem to be losing many of your samples due to low total frequencies closer to the minimum sampling depth than the maximum sampling depth.

:::{describe-usage}
use.action(
    use.UsageAction(plugin_id='diversity',
                    action_id='alpha_rarefaction'),
    use.UsageInputs(table=asv_table_ms2,
                    max_depth=260,
                    metadata=sample_metadata),
    use.UsageOutputNames(visualization='alpha_rarefaction'))
:::

The visualization will have two plots.
The top plot is an alpha rarefaction plot, and is primarily used to determine if the richness of the samples has been fully observed or sequenced.
If the lines in the plot appear to "level out" (i.e., approach a slope of zero) at some sampling depth along the x-axis, that suggests that collecting additional sequences beyond that sampling depth would not be likely to result in the observation of additional features.
If the lines in the plot don't level out, this may be because the richness of the samples hasn't been fully observed yet (because too few sequences were collected), or it could be an indicator that a lot of sequencing error remains in the data (which is being mistaken for novel diversity).

The bottom plot in this visualization is important when grouping samples by metadata.
It illustrates the number of samples that remain in each group when the feature table is rarefied to each sampling depth.
If a given sampling depth `d` is larger than the total frequency of a sample `s` (i.e., the number of sequences that were obtained for sample `s`), it is not possible to compute the diversity metric for sample `s` at sampling depth `d`.
If many of the samples in a group have lower total frequencies than `d`, the average diversity presented for that group at `d` in the top plot will be unreliable because it will have been computed on relatively few samples.
When grouping samples by metadata, it is therefore essential to look at the bottom plot to ensure that the data presented in the top plot is reliable.

:::{exercise} Relative richness.
When grouping samples by "SampleType" and viewing the alpha rarefaction plot for the "observed_features" metric, which sample types (if any) appear to exhibit sufficient diversity coverage (i.e., their rarefaction curves level out)?
:::

### Taxonomic analysis

Next, we can view the taxonomic composition of our samples with interactive bar plots.
Generate those plots with the following command and then open the visualization.

:::{describe-usage}

use.action(
    use.UsageAction(plugin_id='taxa',
                    action_id='barplot'),
    use.UsageInputs(table=asv_table_ms2,
                    taxonomy=taxonomy,
                    metadata=sample_metadata),
    use.UsageOutputNames(visualization='taxa_bar_plots'))
:::

:::{exercise} Which input table?
:label: why-asv-table
Why are we using the ASV table when generating our taxonomic barplot, rather than the kmer table that we've been using in the last few steps?
:::

:::{solution} why-asv-table
:class: dropdown
We use the ASV table here because the feature ids in that table are the same as the ones used in the `FeatureData[Taxonomy]` artifact.
Additionally, kmerization of our data is a tool used for computing diversity metrics - not something we generally intend to use throughout our analyses.
:::

:::{exercise} Taxa bar plots.
Visualize the samples at *Level 2* (which corresponds to the phylum level in this analysis), and then sort the samples by `SampleType`.
What are the dominant phyla in each in `SampleType`?
:::

### Differential abundance testing with ANCOM-BC2

[ANCOM-BC2](https://doi.org/10.1038/s41592-023-02092-7) is a compositionally-aware linear regression model that allows testing for differentially abundant features across sample groups while also implementing bias correction.
This can be accessed using the [`ancombc2` action](xref:q2doc-amplicon-target#q2-action-composition-ancombc2) in the [composition plugin](xref:q2doc-amplicon-target#q2-plugin-composition).

::::{margin}
:::{warning} Differential abundance testing is easy to get wrong! ☠️
Accurately identifying individual features that are differentially abundant across sample types in microbiome data is a challenging problem and an open area of research, particularly if you don't have an *a priori* hypothesis about which feature(s) are differentially abundant.
A q-value that suggests that you've identified a feature that is differentially abundant across sample groups should be considered a hypothesis, not a conclusion, and you need new samples to test that new hypothesis.

In addition to the methods contained in the [composition plugin](xref:q2doc-amplicon-target#q2-plugin-composition), new approaches for differential abundance testing are regularly introduced.
It's worth assessing the current state of the field when performing differential abundance testing to see if there are new methods that might be useful for your data.
If in doubt, consult a statistician.
:::
::::

We'll perform this analysis in a few steps.

#### Filter samples from the feature table

First, we'll filter our samples from our feature table such that we only have the three groups that we have the most samples for.

:::{describe-usage}

asv_table_ms2_dominant_sample_types, = use.action(
    use.UsageAction(plugin_id='feature_table',
                    action_id='filter_samples'),
    use.UsageInputs(table=asv_table_ms2,
                    metadata=sample_metadata,
                    where='[SampleType] IN ("Human Excrement Compost", "Human Excrement", "Food Compost")'),
    use.UsageOutputNames(filtered_table='asv_table_ms2_dominant_sample_types'))
:::

#### Apply differential abundance testing

Then, we'll apply ANCOM-BC2 to see which ASV are differentially abundant across those sample types.
I specify a reference level here as this defines what each group is compared against.
Since the focus of this study is HEC, I choose that as my reference level.
That will let us see what ASVs are over- or under-represented in the other two sample groups (*Human Excrement* and *Food Compost*) relative to HEC, as HEC defines the "global intercept" that will be measured against.

:::{describe-usage}
ancombc2_results, = use.action(
    use.UsageAction(plugin_id='composition',
                    action_id='ancombc2'),
    use.UsageInputs(table=asv_table_ms2_dominant_sample_types,
                    metadata=sample_metadata,
                    fixed_effects_formula='SampleType',
                    reference_levels=['SampleType::Human Excrement Compost']),
    use.UsageOutputNames(ancombc2_output='ancombc2_results'))
:::

Finally, we'll visualize the results.
Taxonomic annotations can optionally be provided here: this helps make ASVs ids more interpretable.

:::{describe-usage}
use.action(
    use.UsageAction(plugin_id='composition',
                    action_id='ancombc2_visualizer'),
    use.UsageInputs(data=ancombc2_results,
                    taxonomy=taxonomy),
    use.UsageOutputNames(visualization='ancombc2-barplot'))
:::

:::{exercise}
Which genus is most enriched in HEC relative to Food Compost?
Which genus is most enriched in HEC relative to Human Excrement?

Which genus is most depleted in HEC relative to Food Compost?
Which genus is most depleted in HEC relative to Human Excrement?
:::

:::{exercise} Perform ANCOM-BC2 on genera, instead of ASVs.
:label: ancombc2-genera
You might be interested in performing ANCOM-BC2 (or other analyses) with ASVs grouped based on the genus they are derived from, rather than on the ASVs themselves.
This may or may not increase your statistical power -- for example, if all organisms in a genus have a similar impact on their host or environment, you may want to know if that genus is differentially abundance across some category of interest.
Take a look at the documentation for the [`collapse`](xref:q2doc-amplicon-target#q2-action-taxa-collapse) action.
How would you use this to collapse our ASV table at the genus level, and then run ANCOM-BC2 on that table?
:::

::::{solution} ancombc2-genera
:class: dropdown

To collapse our ASVs into genera (i.e. level 6 of the Greengenes taxonomy), we can use the following command.

:::{describe-usage}
genus_table_ms2_dominant_sample_types, = use.action(
    use.UsageAction(plugin_id='taxa',
                    action_id='collapse'),
    use.UsageInputs(table=asv_table_ms2_dominant_sample_types,
                    taxonomy=taxonomy,
                    level=6),
    use.UsageOutputNames(collapsed_table='genus_table_ms2_dominant_sample_types'))
:::

We can then provide the resulting table as the input to ANCOM-BC2.

:::{describe-usage}
genus_ancombc2_results, = use.action(
    use.UsageAction(plugin_id='composition',
                    action_id='ancombc2'),
    use.UsageInputs(table=genus_table_ms2_dominant_sample_types,
                    metadata=sample_metadata,
                    fixed_effects_formula='SampleType',
                    reference_levels=['SampleType::Human Excrement Compost']),
    use.UsageOutputNames(ancombc2_output='genus_ancombc2_results'))
:::

And finally, we can visualize the results.
Notice that in this case we're not providing the taxonomy, because we've already intergrated that information by collapsing at the genus level.

:::{describe-usage}
use.action(
    use.UsageAction(plugin_id='composition',
                    action_id='ancombc2_visualizer'),
    use.UsageInputs(data=genus_ancombc2_results),
    use.UsageOutputNames(visualization='genus-ancombc2-barplot'))
:::

::::

## That's it for now, but more is coming soon!

In the near future (as of 17 April 2025) we plan to integrate analyses using the [sample-classifier plugin](xref:q2doc-amplicon-target#q2-plugin-sample-classifier) and [longitudinal plugin](xref:q2doc-amplicon-target#q2-plugin-longitudinal).
In the meantime, here are some suggestions to continue your learning:
1. Build a machine learning classifier that classifies samples accordining to the three dominant sample types in the feature table that we used with ANCOM-BC.
 (Hint: see [`classify-samples`](xref:q2doc-amplicon-target#q2-action-sample-classifier-classify-samples).)
1. Perform a longitudinal analysis that tracks samples from different buckets over time. Which taxa change most over time? (Hint: see [`feature-volatility`](xref:q2doc-amplicon-target#q2-action-longitudinal-feature-volatility).)
1. Remember that the full data set (five sequencing runs) are available in the [gut-to-soil Artifact Repository](https://doi.org/10.5281/zenodo.13887456).
 Grab one of the larger sequencing runs (we worked with a small sequencing run that was generated as a preliminary test), and adapt the commands in this tutorial to work on a bigger data set.

We're also in the process of refactoring our statistical methods for assessing alpha and beta diversity across groups, using the new [stats plugin](xref:q2doc-amplicon-target#q2-plugin-stats).
We're therefore holding off on integrating statistical analysis until we have that ready.
In the meantime, you can refer to you can refer to the [*Moving Pictures*](https://amplicon-docs.readthedocs.io/en/latest/tutorials/moving-pictures.html) tutorial, as well as the [sample-classifier](https://docs.qiime2.org/2024.10/tutorials/sample-classifier/) and [longitudinal](https://docs.qiime2.org/2024.10/tutorials/longitudinal/) tutorials.

## Replay provenance (work in progress!)

:::{warning}
This section is a work in progress as of 17 April 2025.
Right now this section only demonstrates replaying provenance from a single visualization, but this will be expanded to reflect the full analysis.
Also, note that all commands for training the feature classifier are included in the resulting replay script (in case there are some commands in there that you don't remember running).
More on this coming soon!
:::

You might next want to try to adapt the commands presented in this tutorial to your own data, adjusting parameter settings and metadata column headers as is relevant.
QIIME 2's provenance replay functionality can help with this.
Assuming that you ran all of the steps above in a directory called `gut-to-soil/`, run the following command to generate a template script that you can adapt for your workflow:

```shell
qiime tools replay-provenance \
  --in-fp gut-to-soil/taxa-bar-plots.qzv \
  --out-fp g2s-replayed.bash
```

If you need help, head over to the [QIIME 2 Forum](https://forum.qiime2.org).


[^build-requirements-exceed-resources]: The resource requirements exceed those provided by the [*Read the Docs* (RTD) build system](https://docs.readthedocs.com/platform/stable/builds.html#build-resources), which is used to build the documentation that you're reading.
 RTD provides systems with 7GB of RAM for 30 minutes maximum to build documentation.
 That's a very reasonable (and generous) allocation for building documentation, so we choose to work within those contraints rather than creating our own documentation build system like we've had in the past (e.g., for `https://docs.qiime2.org`).
[^iab-database-searching]: kmerization of biological sequences is described in the [*Database Searching* chapter of *An Introduction to Applied Bioinformatics*](https://readiab.org/database-searching.html#kmer-content).
[^iab-machine-learning]: This process is discussed in the [*Machine Learning in Bioinformatics* chapter of *An Introduction to Applied Bioinformatics*](https://readiab.org/machine-learning.html#unsupervised-learning).
[^forum-diversity-metrics]: Learn more about the available metrics in [this QIIME 2 Forum post](https://forum.qiime2.org/t/alpha-and-beta-diversity-explanations-and-commands/2282).
