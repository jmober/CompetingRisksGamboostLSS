##############################################################################
#####         Modeling postoperative mortality in older patients         ##### 
#####          by boosting discrete-time competing risks models          #####
##############################################################################
#####                    Electronic Supplement                           #####
#####                    Author: Moritz Berger                           #####
##############################################################################
#####			 Content: Instructions				 ##### 
##############################################################################

#####
This repository contains all the R-Code to conduct the analyses and to 
reproduce the figures and tables presented in Section 5 and Section 6 of 
the manuscript, as well as the figures and tables of the supplementary material.

To execute the R-Code it is useful to load the RStudio project 
'CompetingRisksGamboostLSS.Rproj'. All paths are set relative to this directory. 

The repository contains the following folders: 
- Application: R-Code of the application presented in Section 5 
- Packages_Functions: R-package DCRvalidation as tar.gz (not on CRAN) and 
		      functions for model fitting 
- Simulation: R-Code of the simulation study presented in Section 6 
#####

#####
Application: The folder contains R-Code to reproduce Figure 3, 4, 5 and S7. The 
main analysis is contained in 'main_analysis.R'. Figure S7 of the second part 
('Comparison to alternative methods') can be reproduced by 'evaluationResampling.R'. 
Note: For confidentiality reasons, the data set provided is not identical with 
the data used in the paper. Instead, the R-Code employs an anonymized data set 
that was generated from the original data in order to produce results that are 
similar (but not identical!) to those in the paper. In particular, data lines 
in the anonymized data set do not refer to real individuals.
#####

#####
Simulation: Each R-program 'simulation..R' contains a function 'one_sim()', to 
conduct one single run of the simulation. Corresponding examples are included at the 
end of each program. The whole simulation (all parameter combinations) was executed 
on a Linux Cluster. The original call is inserted as a comment. The proposed GB approach 
was fitted using the R-packages gamboostLSS and the self-implemented family function 
CSHM(). The folder 'raw_results' contains all the results stored in .rda-files. Using 
the results files, Figure1, 2, S1, S2, S3, S4, S5, and S6 can be reproduced running the 
R-programs 'evaluation..R'. 
#####

All the analyses were run a HPC cluster consisting of 35 nodes with a total of 944 cores. 
The cores are a mixture of Nehalem, Opteron and AMD Epyc processores. The total memory of 
the system is 4.6 TB, connections between nodes are using 40 or 56Gbit/s Infiniband links. 

