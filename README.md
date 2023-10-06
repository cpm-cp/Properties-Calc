# Properties-Calc
This repository contains a project code in Python language to calculate the properties (enthalpy &amp; entropy) for a unique current. This makes it using a data set by the handbook: Introduction to chemical engineering thermodynamics by J.M. Smith, Hendrick Van Ness, Michael Abbott, Mark Swihart

## Why this project?
The finished to this project is realize a complete program to calculate the different thermodynamic propertis by the a specific process. The first task is realize a code for solvinng a turbine system with a inner current and a out current, this system work with six substances and have a heat losses at surrounding.

## How to use?
This code is condensed in a file that import by a _"[tools](https://github.com/cpm-cp/Properties-Calc/tree/main/__init__/tools)"_ module. The first step is neccesary that the user configure the susbtances that involve in the process and also configure the molar fraction for each one substances. For realize this is necceary to change some parameters, as the inner temperature, outer temperature, inner pressure, outer pressure, heat losses rate and the reference parameters to pressure and temoperature.

## What's the project goal?
- Config to solve this speciic case.
- Add a GUI.
- Add config to have different unit process as, turbine, compressor, exchanger heat, pump and other common unit process.
- Connfig to the different unit process can be intereact with other unit process to solving or simulate the process.
- Custom the GUI for a friendly by the users.