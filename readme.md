## Estimating horizontal englacial velocities from (simulated) repeat-pass ice-penetrating radar data

This code accompanies a talk at AGU 2022. It is an early snapshot version of the work. Please cite as:

> Thomas Teisberg, Dustin Schroeder, Paul Summers (2022), Methods for Constraining Englacial Velocity Fields using Airborne Ice-penetrating Radar Data, Abstract C45B-05 presented at 2022 AGU Fall Meeting, 12-16 Dec.

## Running the code

There are two components to the code. One is a set of scripts that use [ISSM](https://issm.jpl.nasa.gov/) (Larour et al., 2012) to simulate 3D englacial velocity fields to create synthetic input measurements. These scripts are located in `issm/`. The main file to run is `issm/generate_simulation_results.m` with your choice of parameters in the top section.

Two example output files have also been included in this repository if you'd prefer to just run the inversion file.

The second component is a Julia script (`invert_velocity.jl`) that processes the results of the ISSM simulation, solves the ODE described in the talk, and produces the figures used. The input data file can be selected at the top of the script.
