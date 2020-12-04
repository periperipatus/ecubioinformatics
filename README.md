# ECU-UNC Bioinformatics Course

Welcome! This course is developed by [Dr Peri Bolton](#https://peribolton.weebly.com/) and [A/Prof. Chris Balakrishnan](http://www.balalab.com/) at ECU. It is a mixture of tutorials from other ECU faculty, as well as ones from judicious internet scrounging (cited appropriately).

Please refer to the [Course Outline](#course-outline) and Course Schedule for more information on delivery and assessment.

Please **DO NOT jump ahead in the tutorials**, ONLY do the work assigned on the day in the schedule as I will make changes throughout the course to best fit our needs. 

# Table of Contents

* [Software Requirements](#software-requirements)
* [Course Outline](#course-outline)
* [Assignments](#assignments)

# Software requirements

Please ensure the following programs are installed on your desktop **prior to commencing the class**. 

* R https://www.r-project.org/
* R-Studio https://rstudio.com/
* Cisco AnyConnect VPN https://ecu.teamdynamix.com/TDClient/1409/Portal/Requests/ServiceDet?ID=11945 & https://itcs.ecu.edu/2020/03/16/working-remotely/
* Decent text editor software, such as [TextWrangler/BBEdit](https://www.barebones.com/products/textwrangler/) (for Mac) or [Notepad++](https://notepad-plus-plus.org/downloads/) (for Windows).
* Slack https://slack.com/ You can use this on the web browser, OR you can download it. The majority of course communications will be through Slack.

Some method for accessing the ECU computers via unix commandline is required (see below).

## Windows Users

Windows 10 offers a linux subsystem that enables you to use linux scripts. https://docs.microsoft.com/en-us/windows/wsl/install-win10. 
RECOMMENDED: I use the Ubuntu subsystem (Step 6), and many bio servers are based off Ubuntu so this is the best option.

If you are not running Windows 10, you can use PuTTY to access the ECU biology server (https://www.putty.org/). Please contact me if you are struggling to set up this option.

## Mac Users

You already have access to unix scripting through your OSX terminal. You can search for it using "Terminal" or find it in the ```Utilities``` folder in ```Applications```


## Required R packages 

The following packages need to be installed on R on your desktop. Everything else can be installed along the way.

```r
install.packages("BiocManager")
```

<br/>
<div align="right">
    <b><a href="#table-of-contents">^ back to TOC</a></b>
</div>
<br/>

# Course Outline

## Course Goal
 
This course is designed to develop practical bioinformatics skills for beginners, with a particular focus on analysis of **RNAseq datasets**. Although there are many graphical-user-interface and command line tools out there for bioinformatics, we will focus on two commonly used command line scripting languages: R and unix. 

We aim to provide the basic tools necessary for participants to analyse an RNAseq dataset from beginning to end on their own.


## Strategy & delivery
Previously this course has been delivered over a high intensity week and been broadly focussed on analysis of various genomic datasets. However, due to COVID circumistances, this course will be delivered remotely over two weeks, by instructor Dr Peri Bolton, with a mixture of scheduled zoom contact and in-your-own-time practical sessions. To minimize zoom fatigue, the tutorials will be conducted independently with opportunity for help using Slack. The tutorial content may be completed any time during the day prior to the daily zoom-based discussion and check in. 

To maximize student retention of material, we will also have student led discussion of relevant papers and an independent group assignment for presentation at completion of the course. 

## Assessment
There will be student led discussions and a final assignment to complete. These will not be graded, but feedback will be provided. Students will work in pairs on a dataset, and work to produce a final 12 minute presentation + 3 minutes of questions. The final mini “conference” will be held on Friday 18th December. 
Student pairs are encouraged to collaborate and troubleshoot tutorial problems together (on Slack, or whatever medium) and make progress on their independent project as they learn new skills throughout the length of the workshop.


## Schedule
See Excel spreadsheet. Note that tutorials will not take the full length of time allotted in the calendar. 

<br/>
<div align="right">
    <b><a href="#table-of-contents">^ back to TOC</a></b>
</div>
<br/>

# Assignments

You and your partner will work together to analyse one of the three RNAseq datasets and present to the group a 15 minute informal conference style presentation on your findings. These data are real data, but I have selected fewer individuals because I wanted to make analysis time shorter.

All data are held in ```~/Bioinformatics_Workshop/<species>_reads/```
Each presentation should include the following:

* A very brief explanation on the study organism and question (no more than 1 slide)
* Hypothesis and Predictions
* Methods: e.g. what techniques did you use to quantify your reads? What genome did you use? What variables did you compare? What model form did use use to quantify differential gene expression for those variables? 
* Results: 
	* Number of reads in experiment
	* Read Quality?
	* Number of reads counted/mapped
	* Any samples removed? Why?
	* Number of genes differentially expressed for each comparison
	* PCA plot colour coded relevant to your variables
	* Any GO enrichment
* Pick 1 gene that was differentially expressed and plot it. What does it do? 
* Also pick 1 or 2 genes that you believe, *a priori* would be involved in the processes under investigation. Plot those too. Are these significant or not? Why?

* Discuss your results - conclusions, limitations, lessons learned.

## Dataset 1: Wire-tailed manakin brain data

These data were published in Horton, BM, Ryder, TB, Moore, IT, Balakrishnan, CN. (2019). Gene expression in the social behavior network of the wire-tailed manakin (*Pipra filicauda*) brain. *Genes, Brain and Behavior*: e12560

There are data from 4 individuals across 3 brain regions that have been selected for you to analyse: Nucleus Taenia of the medial Amygdala (TnA), Bed Nucleus of the Stria Terminalis (BSTm), and Medial Preoptic Area (POM). The fastqs in ```pfil_reads``` are single end.

Use the NCBI RefSeq genome for *Pipra filicauda* as your reference. You may choose to quantify using the transcriptome ```*_rna.fna```, or to align against the annotated genome ```*_genomic.fna``` and ```*_genomic.gtf```. 

You should think about exploring differences in gene expression between tissues, and differences in gene expression according to social status (Territorial or Floater). Territorial males hold a lek territory, and have subordinate and younger floater males who dance with them.

## Dataset 2: White-throated sparrow morph nestling data
 
These data were published in Newhouse, DJ, Barcelo-Serra, M, Tuttle, EM, Gonser, RA, Balakrishnan, CN. (2019). Parent and offspring genotypes influence gene expression in early life. *Molecular Ecology*, 28: 4166-4180.

These birds have alternate plumages (Tan or White) and behaviors, that is determined by a large inversion on chromosome 2. White birds are aggressive and provide less parental care, and these morphs occur in both sexes. 

You have data 12 individuals from Tan and White morph offspring, from the same parental pairings. You might want to look at differences in gene expression with respect to Morph or Sex. The fastqs in ```wtsp_reads``` are single end. 

## Dataset 3: Stickleback male parental care

These data were published in Bukhari, SA, Saul, MC, James, N, Bensky, MK, Stein, LR, Trapp, R, Bell, AM. (2019). Neurogenomic insights into paternal care and its relation to territorial aggression. *Nature Communications*, 10: 1-11.

20 individuals across non-parental and the egg stages and two brain regions were selected. The fastqs in ```stickleback_reads``` are single end. 
Use the repeat masked 75th Release Ensembl Genome of *Gasterosteus aculeatus* for your analyses. 

<br/>
<div align="right">
    <b><a href="#table-of-contents">^ back to TOC</a></b>
</div>
<br/>

