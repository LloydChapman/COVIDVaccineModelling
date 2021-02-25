# Comparison of COVID-19 vaccine prioritization strategies in the United States

This repository contains code and data for the analyses in 'Comparison of COVID-19 vaccination prioritization strategies in the United States'. It includes code for running simulations of a static model of a number of different COVID-19 vaccine prioritization strategies. It also includes code used to simulate the entire population of California (N=39,148,760) with data on demographic factors (age, sex, race/ethnicity, and county of residence), 'special population' status, and comorbidity status for each individual derived from the 2018 American Community Survey [1], California Health Interview Survey (CHIS) [2], and US Bureau of Labor Statistics [3]. The model is calibrated to data on 28,175 COVID-19 deaths in California up to December 30, 2020, which was provided by the California Department of Public Health. We use the model to predict the impact of the different vaccine prioritization strategies on key health outcomes (COVID-19 clinical cases, deaths and disability-adjusted life years) in California over the first 6 months of 2021.

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

## Data

The California Department of Public Health data required to calibrate the model contains personally identifiable information and therefore cannot be made publicly available. Individuals interested in accessing the data should contact the California Department of Public Health. The simulated population data required to run the vaccine prioritization simulations is available [here](https://doi.org/10.5281/zenodo.4516526), and the rest of the data required for running the simulations is in the [Data](Data) subfolder.

The CHIS data used to simulate comorbidity statuses for each individual cannot be made freely available. Individuals interested in accessing the data should apply [here](https://healthpolicy.ucla.edu/chis/data/Pages/GetCHISData.aspx). Dummy data for individuals' comorbidity statuses is simulated in 'simulate_all_CA_counties.R'. To use the comorbidity statuses simulated from the CHIS data, the simulated population data available [here](https://doi.org/10.5281/zenodo.4516526) should be downloaded into the 'Data' subfolder as described below.

## Installing

Clone/download this project onto your machine using the green button at the top right of this page.

### Running the code

Download the simulated population data from [here](https://doi.org/10.5281/zenodo.4516526) into the 'Data' subfolder from the cloned/downloaded repository.

Install the required R packages if they are not already installed by running the following line of code in R

```R
> install.packages(c("triangle","dplyr","tidycensus","acs","censusapi","devtools","foreach","doParallel","readxl","truncnorm","lubridate","reshape2","data.table","matrixStats","abind","broom","sf","ggplot2","RColorBrewer","viridis"))
```

The 'hashmap' package used in 'simulate_all_CA_counties.R' is no longer available on CRAN, so is automatically installed from GitHub in the code if it is not already installed.

The vaccine prioritization simulations can then be run in R by entering

```R
> source("run_vaccine_impact_prediction_death_IFR_model.R")
```

at the command prompt, or by navigating to the downloaded code folder in a terminal window on Mac/Linux and entering

```
% Rscript run_vaccine_impact_prediction_death_IFR.R
```
 
or in Windows command line by entering

```
C:\>"C:\<path>\<to>\Rscript.exe" C:<path>\<to>\run_vaccine_impact_prediction_death_IFR_model.R
```

where `<path>\<to>` represents the paths to Rscript.exe and 'run_vaccine_prediction_death_IFR_model.R' on your machine.

## Built With

* [R version 4.0.3 (2020-10-10)](https://www.r-project.org/)

## Authors

* Lloyd Chapman: <lloyd.chapman@ucsf.edu>
* Poojan Shukla: <pshukla@berkeley.edu>

## License

This project is licensed under the GNU General Public License v3.0 - see the [LICENSE.txt](LICENSE.txt) file for details

## References
1. US Census Bureau, American Community Survey. 2018 American Community Survey 5-year Estimates. (2018) Available at: [https://www.census.gov/programs-surveys/acs/data.html](https://www.census.gov/programs-surveys/acs/data.html) [Accessed January 19, 2021].

2. UCLA Center for Health Policy Research, California Health Interview Survey. Available
at: [https://healthpolicy.ucla.edu/chis/about/Pages/about.aspx](https://healthpolicy.ucla.edu/chis/about/Pages/about.aspx) [Accessed January 19,
2021].

3. U.S. Bureau of Labor Statistics, Employed persons by detailed occupation and age (2019). Available at: [https://www.bls.gov/cps/cpsaat11b.htm](https://www.bls.gov/cps/cpsaat11b.htm) [Accessed January 19, 2021].
