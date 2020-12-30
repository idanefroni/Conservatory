**MafFilter** is a program dedicated to the analysis of genome alignments. It parses and manipulates [MAF files](https://genome.ucsc.edu/FAQ/FAQformat.html#format5) as well as more simple fasta files. Despite various filtering options and format conversion tools, **MafFilter** can compute a wide range of statistics (phylogenetic trees, nucleotide diversity, inferrence of selection, etc.). Current version is 1.3.0.


## What can MafFilter do?

**MafFilter** applies a series of "filters" to a MAF file, in order to clean it, extract data and computer statistics while keeping track of the associated meta-data such as genome coordinates and quality scores.

* It can process the alignment to remove low-quality / ambiguous / masked regions.
* It can export data into a single or multiple alignment file in format such as Fasta or Clustal.
* It can read annotation data in GFF or GTF format, and extract the corresponding alignment.
* It can perform sliding windows calculations.
* It can reconstruct phylogeny/genealogy along the genome alignment.
* It can compute population genetics statistics, such as site frequency spectrum, number of fixed/polymorphic sites, etc.

## How can I get it?

**MafFilter** is built using the [Bio++ libraries](http://biopp.univ-montp2.fr), as well as the boost iostream library for handling of compressed files. [Debian](https://packages.debian.org/search?keywords=maffilter&searchon=names&suite=all&section=all) and [RPM](https://download.opensuse.org/repositories/home:/jdutheil:/Bio++2.4.0/) packages are available.

For compiling the programs yourself, from the downloaded sources or from the git repository, please follow the instructions from the [Bio++ website](http://biopp.univ-montp2.fr/wiki/index.php/Installation).

## How do I use it? 

Several example data sets are distributed along with the source code of the package.
PDF and HTML documentation can be generated from the source distribution using `make pdf` and `make html`, respectively.
The HTML reference manual is also browsable as [here](maffilter.html) (single page) and [here](Manual/index.html) (multi-pages).
Questions can be posted on the dedicated forum: [here](https://groups.google.com/forum/?hl=en#!forum/maffilter).

## References

The **MafFilter** program was published in

```
Dutheil JY, Gaillard S, Stukenbrock EH.
BMC Genomics. 2014 Jan 22;15:53.
MafFilter: a highly flexible and extensible multiple genome alignment files processor.
```
Please consider citing this work if you are using the program. 

MafFilter was originally developped in the context of the study of the Gorilla genome sequence

```
Scally A, Dutheil JY, Hillier LW, Jordan GE, Goodhead I, Herrero J, Hobolth A, Lappalainen T, Mailund T, Marques-Bonet T, McCarthy S, Montgomery SH, Schwalie PC, Tang YA, Ward MC, Xue Y, Yngvadottir B, Alkan C, Andersen LN, Ayub Q, Ball EV, Beal K, Bradley BJ, Chen Y, Clee CM, Fitzgerald S, Graves TA, Gu Y, Heath P, Heger A, Karakoc E, Kolb-Kokocinski A, Laird GK, Lunter G, Meader S, Mort M, Mullikin JC, Munch K, O'Connor TD, Phillips AD, Prado-Martinez J, Rogers AS, Sajjadian S, Schmidt D, Shaw K, Simpson JT, Stenson PD, Turner DJ, Vigilant L, Vilella AJ, Whitener W, Zhu B, Cooper DN, de Jong P, Dermitzakis ET, Eichler EE, Flicek P, Goldman N, Mundy NI, Ning Z, Odom DT, Ponting CP, Quail MA, Ryder OA, Searle SM, Warren WC, Wilson RK, Schierup MH, Rogers J, Tyler-Smith C, Durbin R.
Nature. 2012 Mar 7;483(7388):169-75.
Insights into hominid evolution from the gorilla genome sequence.
```
and was further developped in the following studies:

```
Stukenbrock EH, Bataillon T, Dutheil JY, Hansen TT, Li R, Zala M, McDonald BA, Wang J, Schierup MH.
Genome Res. 2011 Dec;21(12):2157-66.
The making of a new pathogen: insights from comparative population genomics of the domesticated wheat pathogen Mycosphaerella graminicola and its wild sister species.

Stukenbrock EH, Christiansen FB, Hansen TT, Dutheil JY, Schierup MH.
Proc Natl Acad Sci U S A. 2012 Jul 3;109(27):10954-9.
Fusion of two divergent fungal individuals led to the recent emergence of a unique widespread pathogen species.

```
