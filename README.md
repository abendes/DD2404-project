# DD2404-project
Project in applied bioinformatics (DD2404)
"Reducing noise in protein multialignments"

Annika Bendes (910427-4288) abendes@kth.se
</br>
Hanna Hassan (920414-0165) hanhas@kth.se

#Structure of the project
/data: Contains all test data
</br>
/doc: Contains all documents including report and notebooks
</br>
/results: Our noise reduced data
</br>
/src: Source code


#Project description:
Multialignments are noisy. Homologous proteins can contain regions that are not inherited and should therefore not be aligned, and other regions may have evolved so fast that the correct multialignment is impossible to infer. In order to get rid of such problematic regions in subsequential analysis, in particular for phylogeny inference, "bad" columns are often removed using tools such as GBlocks or TrimAl.

In this project, you shall implement a simple noise-reduction method and evaluate its impact on phylogeny inference.
Removing noisy columns

In this project, a column is considered noisy if

    there are more than 50% indels,
    at least 50% of amino acids are unique,
    no amino acid appears more than twice.

You will write a Python program that takes a multialignment as input and removes columns that fulfill at least one of the criteria above. If all columns are removed, your program should exit with an error.
Evaluation

Besides validating your program (making sure it computes what it was supposed to do), you must also evaluate its performance. The basic question is whether it is worthwhile to use multialignment noise reduction?

Test data

For testing, you will use data used by the authors of TrimAl in their evaluation. I have downloaded and reduced the dataset to make it easier to work with. The data is available in/info/appbio12/data/noise_project on the school computers, or downloadable from this site. If you download to a Linux computer, you unpack the data using this command:

prompt> tar -zxf testdata.tar.gz

and you probably want to delete the remanining file testdata.tar.gz.

The testdata contains six subdirectories:

$ ls
asymmetric_0.5  asymmetric_1.0  asymmetric_2.0  
symmetric_0.5  symmetric_1.0  symmetric_2.0

and each such directory contains one tree file, containing a reference tree, and 300 alignments that were created by evolving (using a computer program) sequences along the reference tree and then running muscle to align the resulting sequences. The trees are either symmetric ("easy") or assymetric ("less easy") and the number in the filename denotes the average amount of mutations per site in the sequences. Hence, the *_2.0 directories has alignments with on average two mutations per sequence position, from the root to the leaves and this is of course the hardest case.
What to test

You shall, for each gene alignment, infer a tree using the programs fastprot and fnj (see also the FastPhylo package, or use the -h option to get hints) on both the given alignment and the noice-reduced alignment. How often is the reference tree recovered with the two methods?

BioPython contains a module called Bio.Phylo which unfortunately is somewhat limited. However, it might help you to read phylogenetic trees from file!

However, Bio.Phylo has been problematic to at least one student, and the next recommendation, a Python package called DendroPy, is so much better that you are strongly advised to use it.

DendroPy is installed on the CSC computers and which you can download to your computer. If you use it, you can also try a metric called symmetric distance (often referred to as "Robinson-Foulds", but this is claimed to be wrong according to the DendroPy documentation) to quantify the difference between an inferred tree and the reference tree.
Your report

In your report, you are expected to present clear statistics regarding the value of noise reduction (compared to no reduction) for the six reference trees.
