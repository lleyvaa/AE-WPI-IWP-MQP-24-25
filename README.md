# AE-WPI-IWP-MQP-24-25
Reopsitory to store matlab scripts and data for the Worcester Polytecnic Institute Ionic Wind Propulsion Major Qualifying Project (AY 2024-2025)

INSTRUCTIONS FOR USE OF HVPC DESIGN TOOL:
  *This design tool is adapted from code used in [1] functional for the design space pertaining to the isolation stage of the high voltage power converter. Current cores used in analysis are ETD 3f35 cores.

  1. Download files in "MATLAB and .xls Files" folder.
  2. Open "Inverter_Transformer_Analysis.m" file. This code takes data from the "CoreLossData.xls" file as input to the design space (the "Ecore_actual_EEER_xfmer_LCC_V2.m" matlab script) and generates an excel file containing a selection of transfmromer configurations. "3_3_25_xxx.xls" is an example of the output data
  3. Adjust electrical parameters in "Inverter_Transformer_Analysis.m" and "Ecore_actual_EEER_xfmer_LCC_V2.m" code to suit needs. Refer to MQP report section 3.2.3 and relevant references used for deatils.
  4. Run the code in to generate excel sheet with output parameters.

[1] Y. He, Towards Lightweight High-voltage Power Conversion, Massachusetts Institute of Technology, Department of Electrical Engineering and Computer Science, 2020.
