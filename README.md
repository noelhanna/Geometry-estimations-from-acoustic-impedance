# Geometry-estimations-from-acoustic-impedance
Matlab fitting algorithm to estimate vocal tract area functions from acoustic impedance spectra

- Method described in: "Estimation of vocal tract and trachea area functions from impedance spectra measured through the lips", A Rodriguez, N Hanna, A Almeida, J Smith, J Wolfe, Australasian International Conference on Speech Science and Technology, 77-80, 2018

Uses simplified version of the repository Vocal-tract-acoustic-impedance-simulations (currently private - undergoing revisions)
- Vocal_Tract_test_GENERAL.m

## Features
- Uses a GUI to select a 3 microphone acoustic impedance measurement, then makes a fit based on impedance magnitude
- If no data is loaded (user presses 'cancel') a default file is loaded with closed glottis, inhalation and exhalation data for a male subject.
- Currently, user can select 1-10 VT cylinders for a closed glottis measurement
- Currently, for open glottis measurements, the VT must have 7 cylinders. The sugglottal tract can have 1-4 cylinders.

## Instructions
- run GUI by typing GUI_fit_v0_2
- Alternatively, Vocal_Tract_test_GENERAL.m can be run separately if the Airway structure is specified, e.g. Airway.VT.length = 0.17; Airway.VT.radius = 0.01;, Airway.alpha_multiplier = 5; Airway.SG.radius = ...

## Significant contributions
- Anne Rodriguez

### Not currently working
- Choose to fit based on impedance magnitude, phase or a combination of both. Currently only impedance magnitude is fitted. This can be changed by changing the output of the fit_general_cylinder_model_v0_3.m function.
- Non-rigid parameters are used by default. For a rigid model need to change output of fit_general_cylinder_model_v0_3.m to Z rather than Z_flex.
