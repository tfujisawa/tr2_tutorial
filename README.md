#An introduction to the tr2 delimitation

The trinomial distribution of triplet (tr2) model is a method for multilocus species delimitation. It measures concordance/discordance of gene trees and finds the best delimitation based on a distribution model of rooted triplets.

Detailed description of the method is is the following paper:

Fujisawa, T., Aswad, A. and Barraclough, T.G. (2016) A rapid and scalable method for multilocus species delimitation using Bayesian model  comparison and rooted triplets. Syst.Biol. 65(5): 759-771. [link](https://academic.oup.com/sysbio/article/65/5/759/2223552)

##Software dependecies
* Python packages: numpy and scipy

* Program for building a guide tree: triplec <http://www.cibiv.at/software/triplec/>

##Installation
The tr2 does not require a special procedure for installation. You can download files of tr2 from its repository. <https://bitbucket.org/tfujisawa/tr2-delimitation-python3>
Then, put the tr2-delimitation directory wherever you want. If you want to run Triplec, download the Triplec.jar from its website <http://www.cibiv.at/software/triplec/>, then create a directory named "bin" in the tr2-delimitation directory and put the Triplec.jar in the created "bin" directory.

There are two version of tr2, python2 and python3. I recommend the python3 version unless you have special reasons to use python2.

If you have an environment with Mercurial installed, you can clone the repository.

```
$hg clone https://tfujisawa@bitbucket.org/tfujisawa/tr2-delimitation-python3
```

##Basic command line

```
$ run_tr2.py -t genetree.tre -g guide.tre -o outputname
```


#####input:
The tr2 requires tree files in Newick format.
* The input indivudual gene trees are in Newick format. **_Trees must be rooted. (Using unrooted trees significantly biases the results.)_** The file must contain one tree per one line.
* The input guide tree are also in Newick. Only first tree is used if the file contains multiple trees. The guide tree must also be rooted. The guide tree does not have to contain all individual samples in gene trees. When individuals are missing, they are simply ignored in the delimitation.

#####output:
The tr2 output two files.
* A tab-delimited text of species delimitation
* A tree annotated with delimitation results in Newick format.

When the "-o" option is given, files, *.table.txt and *.tre, are created. If the "-o" option is omitted, results are output in the console.

##An example with Sistrurus rattle snake data

The Sistrurus data set is a multilocus data set which contains sequence alingments of 19 loci sampled from two species (seven subspecies) of *Sistrurus* rattle snakes. [Kubatko et al.](https://academic.oup.com/sysbio/article/60/4/393/1605022) inferred species-level phylogeny and species status using them.

The alignments and reconstructed trees are available from its TreeBase site <https://www.treebase.org/treebase-web/search/study/summary.html?id=11174>.
1
###Tree reconstruction

First, reconstruct gene trees from alignments with RAxML. We use the "-f d" option to run a rapid search to reduce execusion time.
Outgroups are specified with "-o" option.

```
raxml -m GTRGAMMA -T 2 -f d -p 12345 -s sistrurus.ATP.fasta -n ATP -o Agkistrodon_contortrix_1,Agkistrodon_contortrix_2,Agkistrodon_piscivorus_1,Agkistrodon_piscivorus_2
```

This command works only when outgroups are available. If you do not have infromative outgroups, use a rooting command of RAxML.

```
$ raxml -m GTRGAMMA -T 2 -f a -p 12345 -s sistrurus.ATP.fasta -# 100 -x 100 -n ATP
$ raxml -f I -m GTRGAMMA -t RAxML_bipartitions.ATP -n ATP
```

Rooting trees without outgroups sometimes introduces errors. However, it often has a reasonable accuracy.

Automate the reconstruction with a Shell script to run all 19 reconstruction processes.
```
#! /bin/bash
for sq in "$@"
do
	raxml -m GTRGAMMA -T 2 -f d -p 12345 -s $sq -n ${sq//.fasta} -o  Agkistrodon_contortrix_1,Agkistrodon_contortrix_2,Agkistrodon_piscivorus_1,Agkistrodon_piscivorus_2
done
```

Save these lines in run.raxml.root.sh and run it after giving them a permission for execusion.
```
$ chmod +x run.raxml.sh
$ ./run.raxml.sh sistrurus.*.fasta
```

Just check if all trees are correctly reconstructed and rooted.
Then, concatenate all trees into a single file by the "cat" unix command.
```
$ cat RAxML_bipartitions.* > raxml.all.tre
```

Check the raxml outputs. We only have 17 gene trees because some outgroup samples are missing from two files. For now, we ignore these two loci, but alternative rooting methods can be used for these two loci.

Now, we are ready to run the tr2.

###Delimitation
####Inference with gene trees
```
~/tr2-delimitation-python3/run_tr2.py -t raxml.all.tre
```
This command excute a delimitation from gene trees. First, it creates a consensus tree using triplec program, then find the best delimitation on the consensus tree.

Major outputs are a table of delimitation and an annotated tree.
The table shows assingnments of species membership to individual samples.
```
species	sample
1	Sistrurus_catenatus_edwardsii_CO_1
1	Sistrurus_catenatus_tergeminus_KS3_1
1	Sistrurus_catenatus_tergeminus_KS2_2
1	Sistrurus_catenatus_tergeminus_KS2_1
1	Sistrurus_catenatus_edwardsii_NM1_1
1	Sistrurus_catenatus_edwardsii_NM1_2
... ...
4	Sistrurus_miliarius_barbouri_FL2_1
4	Sistrurus_miliarius_barbouri_FL3_2
4	Sistrurus_miliarius_miliarius_NC_1
4	Sistrurus_miliarius_miliarius_NC_2
5	Agkistrodon_piscivorus_2
5	Agkistrodon_piscivorus_1
6	Agkistrodon_contortrix_2
6	Agkistrodon_contortrix_1
```

The tree shows the best delimitation on the guide tree. It also contains support values of delimitation.
Visualize them by R script.
```
library(ape)

tr <- read.tree("./sistrurus.delimit.tre")

tr$node.label <- substr(tr$node.label, 1,6)
plot(tr, show.node.label=T)
```
The delimitation result is mappned on the guide tree. Nodes indicated by "*"'s define species. Descendants of those nodes are clustered into species. Positive values on nodes indicate that these branching are between-species branching. Negative values are within-species.

Just compare this result with [Fujisawa et al (2016)](https://academic.oup.com/view-large/figure/90603751/syw028f6.png). The figure ... says that with 17 loci the estimated number of species are between 3 and 4. Random resolution of nodes sometimes introduces uncertainty on the result. However, overall patterns are consistent.
