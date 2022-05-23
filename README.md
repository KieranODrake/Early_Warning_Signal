# Early_Warning_Signal
 Early warning signals in relation to Covid-19 waves in the UK

The aim is to define Covid-19 wave start dates in the UK and investigate a range of early warning signals,
comparing their lead time.
___________

**Data inputs**

UK Covid-19 case numbers (https://coronavirus.data.gov.uk/details/cases)
data_2022-May-09 - UK Covid-19 cases by sample date.csv
UK Covid-19 hospitalisations (https://coronavirus.data.gov.uk/details/healthcare)
data_2022-May-09 - UK Covid-19 hospital admissions.csv


**Scripts**

data_load.R
- Used to load in data files 

gam_fitting.R
- Used to fit generalised additive model (GAM) to Covid-19 cases or hospitalisations

wave_define_growth_model.R
- Used to define wave dates based on growth of the log of the smoothed cases or hospitalisations as 
	produced from the GAM

wave_define_derivative_model.R
- Used to define wave dates based on the first derivative of the smoothed cases or hospitalisations as 
	produced from the GAM as described in O'Brien & Clements (2021) DOI: 10.1098/rsbl.2021.0487 and 
	in supplementary methods

wave_define_overlay.R
- Main script used to control parameters of interest and to call functions described above
- Also used to plot outputs (although data is not returned from functions and so must be generated from within the relevant functions)

wave_define_hosp_growth_overlay.R
- As above but focused on hospitalisation data and varying the growth threshold