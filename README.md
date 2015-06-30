# egg
a developing tool by Daliang Ning in IEG, University of Oklahoma, Norman, USA.

2015.5.31 Setup: Now, it is just a version control system of egg package for myself and several friends.

##What is new and what is coming?
New functions are ready:
egg.cor.r:  correlation test between env factors and alpha, beta, each category, respectively, using pearson, kendall, spearman and mantel

Almost done function: main.cca.r: CCA test.

Next, I will write a tool to compare two groups of numbers, considering normality, variance homogenity, independence, effect size, etc., in order to do comparison in a right way. 2015.6.16

## [A] What is Egg?
Egg is an Environmental Microbial Ecology analysis tool based on R. Till now, it only has some basic functions.
###(1) Friendly for beginner
A bit more friendly for people who are not familiar with R. You may just need to collect data files and edit input parameters in an excel file, then you can get all output files you need after run a R file.
###(2) One-run-All-out
The main function of this tool is to finish almost all basic biodiversity analysis with just one action.
Life will be better if you just need to run one file and then have a drink when waiting for all results.
###(3) Growing
It is still developing. I am adding functions every week. You may come back to find something new and download the newest version if it helps. Please let me know if you have new ideas for me to try.
##[B] How to use it now?
Please install R before next step. (http://www.r-project.org/)
###(1) Download
Use "Download ZIP" on the github page of egg, then decompress it (so-called Egg folder).
###(2) Have a try
i> Edit input parameters: Open "input/1.input.csv". You will find the name, description and example of each parameter in the 1st, 3rd, 4th columns. Change the values in 2nd column. Begin with "wd" (work directory) and edit each row of 2nd column as you wish. Please change the "prefix" before next step.

ii> Open R, then click "File-open scripts" to open "Rcode/main.r" in Egg folder.

iii> Edit the first line "setwd("...")", input the path of the folder. (you need to change "\" to "/")

iv> If you do not have required R packages in the following line(s), install them by deleting the "#" and run the line(s) [ctrl+r].

v> Run all [in R, click "Edit-Run All"].

vi> Check the output folder.

###(3) Deal with your own data

i> Creat a folder (in any directory in your computer) for your data. This folder is your work directory ("wd").

ii> Creat an empty folder "output" under the work directory.

iii> Copy "input" folder of Egg into the work directory.

iv> Copy your data files to "input", change them to EXACTLY the same format as the example files in "input", use other names.

v> Edit parameters in "input/1.input.csv" file (only the 2nd column) in your work directory. Remember to change the cells about file names to what your data files named. The "code.wd" could be the path of "Rcode" folder in Egg folder. Check the parameters one by one before next step.

vi> Open R. Open "Rcode/main.r" in the Egg folder.

v> Edit the first line "setwd("...")", input the path of your work directory. (you need to change "\" to "/")

vi> Run all [in R, click "Edit-Run All"], and wait..., and then check the "output" folder :)

## Recommend data analysis pipeline

The two website-based pipelines are quit useful. Some output files of Egg (e.g. "...rawOTU.txt", "...ieg.comm.txt") are compatible with their upload requirement.

IEG sequencing data analysis pipeline: http://zhoulab5.rccc.ou.edu:8080/

IEG microarray pipeline: http://ieg.ou.edu/microarray/

* The microarray pipeline includes various statistic methods and can also be used to analyze sequencing results.

## Citation

Feel free to use it.

If you used it in any publication, you need to cite R and vegan package, and you may acknowledge me (Daliang Ning in IEG, Univeresity of Oklahoma, Norman, USA) or cite Egg as what shows in the last few lines of Rcode/main.r.
