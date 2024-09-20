# neigborhood-effects-on-dispersal-fitness

R code and data to reproduce all the analyses and plots in: 

Danielle K Barnes, Scott C. Burgess. (in press) Fitness consequences of marine larval dispersal: the role of neighborhood density, arrangement, and genetic relatedness on survival, growth, reproduction, and paternity in a sessile invertebrate. *Journal of Evolutionary Biology* 

1. Title of Dataset: Fitness consequences of marine larval dispersal: the role of neighborhood density, arrangement, and genetic relatedness on survival, growth, reproduction, and paternity in a sessile invertebrate

2. Author Information  
	A. Principal Investigator Contact Information  
		Name: Scott Burgess  
		Institution: Florida State University  
		Address: 319 Stadium Drive, Tallahassee, FL, USA 32306  
		Email: sburgess@bio.fsu.edu  
	B. First Author Contact Information  
		Name: Danielle Barnes. 
		Institution: Florida State University  
		Address: 319 Stadium Drive, Tallahassee, FL, USA 32306  
		Email: barnes.danielle001@gmail.com  

3. Date of data collection: 2021-2022   

4. Location of data collection: Turkey Point, Florida, USA  

5. Funding: National Science Foundation (NSF; OCE-1948788)  



DATA & FILE OVERVIEW

`Experiment_1.csv`  
`Experiment_2.csv`  
`Experiment_3.csv`  
`D0_9_final_Paternity.txt`  
`D0_9_final_BestCluster.csv`  
`D0_9_final_BestFSFamily.csv`  



2. Relationship between files: 

`Experiment 1.Rmd` uses `Experiment_1.csv` is used to make Figure 2

`Experiment 2.Rmd` uses `Experiment_2.csv` is used to make Figure 3

`Experiment 3.Rmd` uses `Experiment_3.csv`, `D0_9_final_Paternity.txt`, `D0_9_final_BestCluster.csv`, `D0_9_final_BestFSFamily.csv` to make Figure 4, 5, and 6.


\

3. Metadata

# Experiment_1.csv

*Unique.ID*: Unique identifier for each F1 colony out planted to the field
*Mother.colony*: Mother identifier for each F1 colony out planted to the field.

*Density*: Number of colonies per replicate (within a 1.5 x 3 inch grid on a petri dish lid; ranges from 2 to 20 colonies in increments of 2)

*Relatedness*: Non=Non-related, F1 colonies from different mothers; Sib=Siblings, F1 colonies from the same mother

*Grid.position*: Numbers correspond to specific positions of colonies on the grid of the petri dish lid. 1 though 5 on top row right to left, 6 through 10 on second row right to left, 11 through 15 on third row right to left, and 16 through 20 on bottom row right to left.

*X*: X coordinate on the grid, 1-5

*Y*: Y coordinate on the grid, 1-4 

*Inside.outside*: Relative position on the grid, positions 1-5, 6, 10, 11, 15, and 16-20 are considered outside, others are inside.

*Date*: Date of data collection

*Age.days*: Age of F1 focal colonies; days since settlement

*Bifurcations*: Number of bifurcations on the focal colony counting longest chain

*Zooids*: Number of zooids on the focal colony counted under the microscope

*Survival*: 1=survived, 0=died

\

# Experiment_2.csv

*Unique.ID*: Unique identifier for each F1 colony out planted to the field

*Mother.colony*: Mother identifier for each F1 colony out planted to the field.

*Relatedness*: Non=Non-related, F1 colonies from different mothers; Sib=Siblings, F1 colonies from the same mother

*Grid.position*: Numbers correspond to specific positions of colonies on the grid of the petri dish lid. 1 though 5 on top row right to left, 6 through 10 on second row right to left, 11 through 15 on third row right to left, and 16 through 20 on bottom row 
right to left.

*X*: X coordinate on the grid, 1-5

*Y*: Y coordinate on the grid, 1-4 

*Position*: Relative position on the grid, positions 1-5, 6, 10, 11, 15, and 16-20 are considered outside, others are inside.

*Date*: Date of data collection

*Age.days*: Age of F1 focal colonies; days since settlement

*Bifurcations*: Number of bifurcations on the focal colony counting longest chain

*Zooids*: Number of zooids on the focal colony counted under the microscope

*Survival*: 1=survived, 0=died

*Fertilized.ovicells*: Number of fertilized ovicells (black) on a focal colony, counted under the microscope

*Ovicells.total*: Number of total ovicells (fertilized and unfertilized) on a focal colony, counted under the microscope

\

# Experiment_3.csv

*Colony*: Unique identifier for each colony out planted to the field (F1 generation). Number corresponds to the parent in the F0 generation.

*Block*: Spatial section of the field site. Block 1 is western most, Block 4 eastern most. Each block separated by 5m. Each block contains one replicate of a treatment.

*Treatment*: The treatment each Colony was placed into in the field. Near=two colonies 0.15m apart, far=two colonies 1m apart, both=three colonies (two 0.15m apart and one 1m away), alone=single colony at least 5m from all other outplanted colonies

*Position*: A=alone on the center transect line; B=0.15m from C on the center transect line; C=0.15m from B; D=1m from F, 0.15m from E, on the center transect line; E=0.15m from D, 0.75m from F; F=1m from D, 0.75m from E; G=1m from H, on the center transect 
line; H=1m from G.

*Direction*: N=North (on the north side of the center transect line, on the side closest to land); S=South (on the south side of the center transect line, on the Gulf side); C=On the center transect line.

*X*: Distance along the X axis transect (parallel to shore, moving west to east) in meters

*Y*: Distance along the Y axis (North-South) in meters, corresponds to direction with positive being toward the land (North) and negative being toward the Gulf (South)

*Time_days*: Number of days since the start of the experiment (start being the first day in the field, not the settlement of F1 colonies)

*Bifurcations*: Number of bifurcations counting longest chain

*Zooids*: Number of zooids counted under the microscope

*rgr40*: Relative growth rate from day 0 to day 40. Zooids per zooid per day. (log(Zooids40) - log(Zooids0)) / (40-0)

*Offspring*: Total number of offspring (F2) produced from the timepoint of previous collection to the current timepoint - larvae and settlers.

*Samples*: The number of samples preserved for genotyping (not the number actually genotyped)

*Survival*: 1=survived, 0=died

\

# D0_9_final_Paternity.txt

Output from the program COLONY when the probability of a father in the candidates provided is 0.9

*OffspringID*: Unique identifier for each offspring

*InferredDad1*: The identification of the inferred father from the list of candidate fathers estimated in COLONY

*ProbDad1*: The probability that the InferredDad1 is the father of each offspring in OffspringID


\

# D0_9_final_BestCluster.csv

Output from the program COLONY when the probability of a father in the candidates provided is 0.9

*ClusterIndex*: Unique ID for each pedigree cluster estimated by the program COLONY


*Probability*: The probability of each Cluster estimated by the program COLONY (not used)

*OffspringID*: Unique identifier for each offspring genotyped

*FatherID*: Unique identifier for each inferred father estimated by the program COLONY

*MotherID*: Unique identifier for each known mother

\

# D0_9_final_BestFSFamily.csv

*FullSibshipIndex*: Unique ID for each full sib family

*Prob(Inc.)*: The inclusion probability. The probability that all individuals of a given full sib family are full sibs 

*Prob(Exc.)*: The exclusion probability. The probability that no other individuals are full-sibs with a given full-sib family

*Member1* through *Member18*: The OffspringID for members of each full sib family

