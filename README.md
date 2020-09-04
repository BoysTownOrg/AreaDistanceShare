# AreaDistanceShare

[Version 09/04/2020]

This repository AreaDistanceShare contains MATLAB-readable software to generate the figures contained in the article:

"Sound-field estimation near the tympanic membrane using area-distance measurements in the ear canal," Douglas H. Keefe. *J. Acoust. Soc. Am.* 148, 1193-1214 (2020).   

The software was tested using MATLAB Release R2020a. Please refer to the above published article for physical explanations.

The correspondence between MATLAB scripts and figure numbers is:  

| Figure # | MATLAB script file      |
| -------- | ----------------------- |
| 1        | `PlotSpaceTime_Share.m` |
| 2	| `DoOneWayModel_Share.m` |
| 3 | `DoCalib_Share.m` |
| 4 | `DriverCone_Share.m` |
| 5, 7-9 | `DoNearTM_Share.m` |
| 6 | `DoAreas_Share.m` |

<div style="page-break-after: always"></div>

### Notes on MATLAB code

1. The variable `Base.altitude` holds the altitude of the test site in meters, 300 m for Omaha, NE in the published figures. This value may be adjusted to zero to use thermodynamic constants of air evaluated at sea level. These air constants are calculated along with other general code in the classdef file `BaseRFShare.m`. There is a comment in the code of each MATLAB script at the relevant line at which `Base.altitude` can be adjusted.

2. The code for Figs. 5-9 uses as input files the data files (MAT) that were written out by other MATLAB software used to measure ear responses in D. H.  Keefe, "Causality-constrained measurements of aural acoustic reflectance and reflection functions,"  *J. Acoust. Soc. Am.* 147, 300-324 (2020). The code for Figs. 2-3 also includes software to calculate the decimation filtered model reflection function of a closed cylindrical tube. This is adapted from software used in the 2020 paper cited in this note.

3. Running the script `DoNearTM_Share` over-writes and saves intermediate data in the file `LosslessResults.MAT`, which is subsequently loaded back in by the script. This part of the code compares loss-less and lossy outputs near the tympanic membrane.




### License information	

The repository AreaDistanceShare holds files in the main folder that are copyrighted and licensed under the MIT License (see License.txt in the main folder).

`PlotSpaceTime_Share.m` calls two other external files to make Fig. 1 that are not covered by the above license. These files are contained in subfolders of the main folder as follows:  

- `arrow.m` with its `license.txt` file in subfolder `arrow` downloaded from https://www.mathworks.com/matlabcentral/fileexchange/278-arrow.

- `drawaxis.m` with its `license.txt` in subfolder `drawaxis`was obtained via a download link at this MathWorks website address: https://www.mathworks.com/matlabcentral/answers/100506-how-do-i-move-the-x-axis-so-that-it-always-intercepts-the-y-axis-at-y-0-on-a-two-dimensional-plot. 

  



