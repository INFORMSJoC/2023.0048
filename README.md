[![INFORMS Journal on Computing Logo](https://INFORMSJoC.github.io/logos/INFORMS_Journal_on_Computing_Header.jpg)](https://pubsonline.informs.org/journal/ijoc)

# Guideline for Generating Results

In what below we provide a detailed description on how to use our codes to generate each of the figures which are in the results folder. We note that
we have also added all data generated from the codes and needed in other codes (e.g., plotting codes in the data folder). So users might use those
data files for validation purposes.

## Data for Prevalence Rates Used in all code files
Calculations are being done in the excel file named "DataForPrevalenceRates" in the data folder in the sheets "GRP", "Budget(1)" and "Budget(2)".
Columns C, D, and E in Sheet "Budget(2)" are the mean, upper, and lower bounds of the prevalence rates for each disease we have in the code files



## Figure named "Computational Complexity":
1) run the matlab code "numericalAnalysis_resourceAllocation.m". Note that all other matlab files need to be in the same directory.
In the code, you set N to the desired level (25, 50, or 100). Once the code terminates, the results are stored in variable Yexport. 
These codes should result in data similar to "caseStudyNumericalExperimentsminmaxregret" or close, since the code is generating problem instances at random.)
This is the only part where matlab code is required. However for the rest the code is python.
Again this is requried to generate the csv file titled "caseStudyNumericalExperimentsminmaxregret".
2) run the code "Refined_Algo_ComputationalComplexity-Randomized" to generate the numpy data titled "data when N=25", "data when N=50", and "data when N=100". (these are not uploaded,
due to their size, however they can be generated and stored in the same directory)
3) run the code "ChangeFormatComputationalComplextity_random and minmax" to change and fix the format of the above outputs and generate
files "CaseStudyNumericalExperimentsrandomized when N=25", "CaseStudyNumericalExperimentsrandomized when N=50" and 
"CaseStudyNumericalExperimentsrandomized when N=100" used in plotting, as in what next.
4) run the code "PlottingComputationalComplexity_random and minmax" to generate the figures


## Figure named "GroupSizeWithTime":
1) run the code "Regret_GroupTesting" which generates the following data files 
 a) "F_constrained_matrix_of_ks_RegretAPP1"
 b) "solutionovertime_RegretAPP1"
 c) "group_size_path_of_ksRegretAPP1"
 d) "path_matrix_of_ksRegretAPP1"
Those files are used to get the numbers inside the code file "PlottingGroupSizewithTime"
2) run the code "PlottingGroupSizewithTime" to generate the figure


## Figure named "PerTimeCost_GRP":
1) run the code "Regret_GroupTesting" which generates the following data files for all K values
 a) "F_constrained_matrix_of_ks_RegretAPP1"
 b) "solutionovertime_RegretAPP1"
 c) "group_size_path_of_ksRegretAPP1"
 d) "path_matrix_of_ksRegretAPP1"
2) run the code "PureRobust_GroupTesting" which generates the following data files 
 a) "F_constrained_matrix_of_ks_RobustAPP1"
 b) "solutionovertime_RobustAPP1"
 c) "group_size_path_of_ksRobustAPP1"
 d) "path_matrix_of_ksRobustAPP1"
Those files are used to extract the numbers used in "PlottingFunction_GroupTesting_ActualData-with Robust", note that one could also read those files,
but in this code we copied and pasted the values for those data files inside the code, from the console where the numbers for steps 1 and 2 were printed.
Note that for the data files 1a), 1b), 1c), and 1d), and 2a), 2b), 2c), and 2d), 
we copied the data for K=0, 1, 2, 3, only i.e., the first 4 elements of each data structure: this is according to the following:
F_of_n_star_continuous in "PlottingFunction_GroupTesting_ActualData-with Robust" is F_of_n_star_continuous in "Regret_GroupTesting"
F_of_n_equal_16 in "PlottingFunction_GroupTesting_ActualData-with Robust" is F_of_n_equal_16 in "Regret_GroupTesting"
F_constrained_matrix_of_ks in "PlottingFunction_GroupTesting_ActualData-with Robust" is F_constrained_matrix_of_ks in "Regret_GroupTesting"
F_constrained_matrix_of_ks_Robust in "PlottingFunction_GroupTesting_ActualData-with Robust" is F_constrained_matrix_of_ks in "PureRobust_GroupTesting"
path_matrix_of_ks in "PlottingFunction_GroupTesting_ActualData-with Robust" is path_matrix_of_ks in "Regret_GroupTesting"
group_size_path_of_ks in "PlottingFunction_GroupTesting_ActualData-with Robust" is group_size_path_of_ks in "Regret_GroupTesting"

Similarly, numbers for MinMaxDistance, Total_Average_distance, MinMaxDistance_from16, and Total_Average_distance_from16_Robust which are in 
the code file "PlottingFunction_GroupTesting_ActualData-with Robust" are also copied and pasted from the console 
for the following variables: MinMaxDistance, Total_Average_distance, MinMaxDistance_from16, Total_Average_distance_from16
from 1) and Total_Average_distance_from16 from 2)(resp.), (again first 4 elements)
Now that we have all values needed in "PlottingFunction_GroupTesting_ActualData-with Robust" we,
2) run the code "PlottingFunction_GroupTesting_ActualData-with Robust" to generate the figure



## Figure named "PerTimeCost_Budget":
1) run the code "Regret_Budget" which generates the following data files for all K values
 a) "F_constrained_matrix_of_ks_RegretAPP2"
 b) "solutionovertime_RegretAPP2"
 c) "group_size_path_of_ksRegretAPP2"
 d) "path_matrix_of_ksRegretAPP2"
2) run the code "PureRobust_Budget" which generates the following data files 
 a) "F_constrained_matrix_of_ks_RobustAPP2"
 b) "solutionovertime_RobustAPP2"
 c) "group_size_path_of_ksRobustAPP2"
 d) "path_matrix_of_ksRobustAPP2"
Those files are used to extract the numbers used in "PlottingFunction_Budget_ActualData-with Robust", note that one could also read those files,
but in this code we copied and pasted the values for those data files inside the code, from the console where the numbers for steps 1 and 2 were printed.
Note that for the data files 1a), 1b), 1c), and 1d), and 2a), 2b), 2c), and 2d), 
we copied only the data for K=0, 1, 2, 3, only i.e., the first 4 elements of each data structure: this is according to the following:
F_continuous_array in "PlottingFunction_Budget_ActualData-with Robust" is function_continuous_case in "Regret_Budget"
F_constrained_matrix_of_ks in "PlottingFunction_Budget_ActualData-with Robust" is F_constrained_matrix_of_ks in "Regret_Budget"
F_constrained_matrix_of_ks_Robust in "PlottingFunction_Budget_ActualData-with Robust" is F_constrained_matrix_of_ks in "PureRobust_Budget"
path_matrix_of_ks in "PlottingFunction_Budget_ActualData-with Robust" is path_matrix_of_ks in "Regret_Budget"
budget_path_3d_matrix_of_ks in "PlottingFunction_Budget_ActualData-with Robust" is budget_path_3d_matrix_of_ks in "Regret_Budget"

Similarly, numbers for MinMaxDistance, Total_Average_distance, and Total_Average_distance_Rob which are in 
the code file "PlottingFunction_Budget_ActualData-with Robust" are also copied and pasted from the console 
for the following variables: MinMaxDistance, Total_Average_distance from 1) and Total_Average_distance from 2)(resp.), (again first 4 elements)
Now that we have all values needed in "PlottingFunction_Budget_ActualData-with Robust" we,
2) run the code "PlottingFunction_Budget_ActualData-with Robust" to generate the figure


## Figure named "Prev_Babesiosis_WNV":
1) run the code "AdditionalLaTeXFigures", lines 1-5 (to import packages)and lines 249-319 (to generate the figure)


## Figure named "Prev_WNV_Reported Numbers":
1) run the code "AdditionalLaTeXFigures", lines 1-5 (to import packages)and lines 137-183 (to generate the figure)


## Figure named "UncertaintySet_CI":
1) run the code "AdditionalLaTeXFigures", lines 1-26 (to import packages and define function)and lines 326-359 (to generate the figure)


## Figure named "YearlyBudget":
1) run the code "YearlyBudget" and copy from the console the following data and paste them in the file "PlottingYearlyBudget_K=3"
2) run code "PerTimeCost_Budget" to generate 1)a, 1)b, 1)c, and 1)d and copy from the console the following data and paste them in the file "PlottingYearlyBudget_K=3"
Step 1) and 2) here are as follows:
Path_avg_B in "PlottingYearlyBudget_K=3" is avg_bd_solution_path in "YearlyBudget"
avg_bd_solution_path in "PlottingYearlyBudget_K=3" is avg_bd_solution in "YearlyBudget"
X_path_avg_Bnotrounded in "PlottingYearlyBudget_K=3" is X_allocation_path in "YearlyBudget"
f_avg_B in "PlottingYearlyBudget_K=3" is function_constrained_case in "YearlyBudget"

Now for the code "PerTimeCost_Budget", consider only the 4rth element, i.e., for K=3:
Path_constant_B in "PlottingYearlyBudget_K=3" is path_matrix_of_ks[3] in "PerTimeCost_Budget"
constant_bd_solution_path in "PlottingYearlyBudget_K=3" is 45
X_path_constant_Bnotrounded in "PlottingYearlyBudget_K=3" is budget_path_3d_matrix_of_ks[3] in "PerTimeCost_Budget"
f_constant_B in "PlottingYearlyBudget_K=3" is F_constrained_matrix_of_ks[3] in "PerTimeCost_Budget"

3) run "PlottingYearlyBudget_K=3" to generate the figure


## Figures named "SimulationBasedAnalysis_GRP" and "TotalCostwrtK_GRP":
1) run the code "Regret_GroupTesting", which generates the following data files for all K values
 a) "F_constrained_matrix_of_ks_RegretAPP1"
 b) "solutionovertime_RegretAPP1"
 c) "group_size_path_of_ksRegretAPP1"
 d) "path_matrix_of_ksRegretAPP1"
2) run the code "PureRobust_GroupTesting" which generates the following data files 
 a) "F_constrained_matrix_of_ks_RobustAPP1"
 b) "solutionovertime_RobustAPP1"
 c) "group_size_path_of_ksRobustAPP1"
 d) "path_matrix_of_ksRobustAPP1"
3) run code "GroupTesting_SamplingApproach" which reads those solutions and then performs the simulation based analysis. 
Lines 1-231 generate figure "SimulationBasedAnalysis_GRP" and lines 247-294 generate figure "TotalCostwrtK_GRP"



## Figures named "SimulationBasedAnalysis_Budget" and "TotalCostwrtK_Budget":
1) run the code "Regret_Budget", this will generate the following data files for all K values
 a) "F_constrained_matrix_of_ks_RegretAPP2"
 b) "solutionovertime_RegretAPP2"
 c) "group_size_path_of_ksRegretAPP2"
 d) "path_matrix_of_ksRegretAPP2"
2) run the code "PureRobust_Budget", this will generate the following data files 
 a) "F_constrained_matrix_of_ks_RobustAPP2"
 b) "solutionovertime_RobustAPP2"
 c) "group_size_path_of_ksRobustAPP2"
 d) "path_matrix_of_ksRobustAPP2"
3) run code "Budget_SamplingApproach" which reads those solutions and then perform the simulation based analysis.
Lines 1-275 generate figure "SimulationBasedAnalysis_Budget" and lines 290-336 generate figure "TotalCostwrtK_Budget"



# Cite
To cite the contents of this repository, please cite both the paper and this repo, using their respective DOIs.

https://doi.org/10.1287/ijoc.2023.0048

https://doi.org/10.1287/ijoc.2023.0048.cd

Below is the BibTex for citing this snapshot of the repository.

@misc{Shams2025Optimal,
  author =        {Shams Eddin, Marwan and El-Amine, Hadi and Aprahamian, Hrayer},
  publisher =     {INFORMS Journal on Computing},
  title =         {{Optimal System Adjustment Under Operational Constraints with Applications to Infectious Disease Screening},
  year =          {2025},
  doi =           {10.1287/ijoc.2023.0048.cd},
  url =           {https://github.com/INFORMSJoC/2023.0048},
  note =          {Available for download at https://github.com/INFORMSJoC/2023.0048},
}
