This folder contains two python scripts and one .mat data file:

1. VDA_functions.py:
	This file includes all functions needed to complete the VDA13 model
	all functions include a comprehensive description and the associated
	author's name is listed. 
	@staticmethod is used so that these functions may be called into other
	files.
2. VDA13_SD97_Model_Validation.py:
	This file employs the VDA_functions and applies them to a case study
	based on data collected during the '97 Sandy Duck experiment. The code
	is well commented and explains what portion of the project each function
	is addressing.
3. A_intrawave_velocity_timeseries_output_file.mat:
	This file includes the wave data obtained from the gauge deployed nearest
	to the shore during the '97 Sandy Duck experiment, and its data are used
	as inputs to the VDA functions.