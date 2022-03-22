# Comparison of COVID-19 vaccine prioritization strategies

This repository contains code and data for the analyses in 'Risk factor targeting for vaccine prioritization during the COVID-19 pandemic' [1]. It includes code for running simulations of a static model of a number of different COVID-19 vaccine prioritization strategies. It also includes code used to simulate the entire population of California (N=39,148,760) with data on demographic factors (age, sex, race/ethnicity, and county of residence), 'special population' status, and comorbidity status for each individual derived from the 2018 American Community Survey [2], California Health Interview Survey (CHIS) [3], and US Bureau of Labor Statistics [4]. The model is calibrated to data on 28,175 COVID-19 deaths in California up to December 30, 2020, which was provided by the California Department of Public Health. We use the model to predict the impact of the different vaccine prioritization strategies on key health outcomes (COVID-19 clinical cases, deaths and disability-adjusted life years) in California over the first 6 months of 2021.

## Prerequisites

* [R version 4.0.3](https://www.r-project.org/)

* The following R packages are required to run the code:
  * triangle
  * dplyr
  * tidycensus
  * acs
  * censusapi
  * devtools
  * hashmap
  * foreach
  * doParallel
  * readxl
  * truncnorm
  * lubridate
  * reshape2
  * data.table
  * matrixStats
  * abind
  * broom
  * sf
  * ggplot2
  * RColorBrewer
  * viridis
  * cowplot
  * patchwork

## Data

The California Department of Public Health data required to calibrate the model contains personally identifiable information and therefore cannot be made publicly available. Individuals interested in accessing the data should contact the California Department of Public Health. The simulated population data required to run the vaccine prioritization simulations is available [here](https://doi.org/10.5281/zenodo.4516526), and the rest of the data required for running the simulations is in the [Data](Data) subfolder.

The CHIS data used to simulate comorbidity statuses for each individual cannot be made freely available. Individuals interested in accessing the data should apply [here](https://healthpolicy.ucla.edu/chis/data/Pages/GetCHISData.aspx). Dummy data for individuals' comorbidity statuses is simulated in [simulate_all_CA_counties.R](Code/simulate_all_CA_counties.R). To use the comorbidity statuses simulated from the CHIS data, the simulated population data available [here](https://doi.org/10.5281/zenodo.4516526) should be downloaded into the [Data](Data) subfolder as described below.

## Installing

Clone/download this project onto your machine using the green button at the top right of this page.

### Running the code

Download the simulated population data from [here](https://doi.org/10.5281/zenodo.4516526) into the [Data](Data) subfolder from the cloned/downloaded repository.

Install the required R packages if they are not already installed by running the following line of code in R

```R
install.packages(c("triangle","dplyr","tidycensus","acs","censusapi","devtools","foreach","doParallel","readxl","truncnorm","lubridate","reshape2","data.table","matrixStats","abind","broom","sf","ggplot2","RColorBrewer","viridis","cowplot","patchwork"))
```

The 'hashmap' package used in [simulate_all_CA_counties.R](Code/simulate_all_CA_counties.R) is no longer available on CRAN, so is automatically installed from GitHub in the code if it is not already installed.

The vaccine prioritization simulations can then be run in R by entering

```R
setwd("<path>/<to>/Code/")
source("run_vaccine_impact_prediction_death_IFR_model.R")
```

at the command prompt, where `<path>/<to>` is the path to the [Code](Code) subfolder on your machine, or by navigating to the [Code](Code) subfolder in the downloaded repository in a terminal window on Mac/Linux and entering

```
Rscript run_vaccine_impact_prediction_death_IFR.R
```
 
or in Windows command line by entering

```
"C:\<path>\<to>\Rscript.exe" C:<path>\<to>\run_vaccine_impact_prediction_death_IFR_model.R
```

where `<path>\<to>` represents the paths to Rscript.exe and 'run_vaccine_prediction_death_IFR_model.R' on your machine.

## Built With

* [R version 4.0.3 (2020-10-10)](https://www.r-project.org/)

## Authors

* Lloyd Chapman: <lloyd.chapman1@lshtm.ac.uk>
* Poojan Shukla: <pshukla@berkeley.edu>

## License

This project is licensed under the GNU General Public License v3.0 - see the [LICENSE.txt](LICENSE.txt) file for details

## References
1. L.A.C. Chapman, P. Shukla, I. Rodriguez-Barraquer, P.B. Shete, T.M. Leon, K. Bibbins-Domingo, G.W. Rutherford, R. Schechter, N.C. Lo. Risk factor targeting for vaccine prioritization during the COVID-19 pandemic. Scientific Reports. 2022. [https://doi.org/10.1038/s41598-022-06971-5](https://doi.org/10.1038/s41598-022-06971-5).

2. US Census Bureau, American Community Survey. 2018 American Community Survey 5-year Estimates. (2018) Available at: [https://www.census.gov/programs-surveys/acs/data.html](https://www.census.gov/programs-surveys/acs/data.html) [Accessed January 19, 2021].

3. UCLA Center for Health Policy Research, California Health Interview Survey. Available
at: [https://healthpolicy.ucla.edu/chis/about/Pages/about.aspx](https://healthpolicy.ucla.edu/chis/about/Pages/about.aspx) [Accessed January 19,
2021].

4. U.S. Bureau of Labor Statistics, Employed persons by detailed occupation and age (2019). Available at: [https://www.bls.gov/cps/cpsaat11b.htm](https://www.bls.gov/cps/cpsaat11b.htm) [Accessed January 19, 2021].
