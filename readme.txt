Readme file:

This is an implementation of the ion-channel noise method of Fox and Lu, as described 
by Goldwyn and Shea-Brown, applied to the Hodgkin-Huxley model and to the 
Rothman-Manis model of VCN neurons, as examples.

Reference: Brett A. Schmerl and Mark D. McDonnell, "Channel noise induced stochastic 
facilitation in an auditory brainstem neuron model", Physical Review E 88:052722, 2013
(http://link.aps.org/doi/10.1103/PhysRevE.88.052722). Preprint available at http:/
arxiv.org/abs/1311.2643

The following files have been uploaded to ModelDB:

SchmerlMcDonnell_Driver.m -- this script can be modified by the user to control the 
inputs and simulation conditions, and select from amongst the four example neuron models 
implemented, i,e. Hodgkin-Huxley, Rothman-Manis Type I-II, Rothman-Manis Type I-C,
Rothman-Manis Type II.

EulerMaruyama.m -- function that controls the updating of all variables for the next time step

UpdateEquations.m -- function that calculates the RHS of differential equations

UpdateDiffusionMatrixSquareRoot.m -- function that calculates the matrix square root
needed to obtain the diffusion matrix

UpdateDriftAndOccupancies.m -- function that updates the transition matrix and
occupancies

Params_HodgkinHuxley.m -- parameters for the Hodgkin-Huxley model

Params_RothmanManisTypeII.m  -- parameters for the Rothman-Manis Type II model

Params_RothmanManisTypeI_C.m  -- parameters for the Rothman-Manis Type I-C model

Params_RothmanManisTypeI_II.m  -- parameters for the Rothman-Manis Type I-II model


Usage:

Running SchmerlMcDonnell_Driver.m from the command line will, by default, solve the
deterministic version of the Hodgkin-Huxley model with an input current of 7.2 microAmps,
for 10000 channels, over a simulation duration of 400 ms, with initial conditions given in 
Params_HodgkinHuxley.m. The solution for the membrane potential will be plotted against
time. The user can edit SchmerlMcDonnell_Driver.m to change the model, to use the SSE
stochastic model, and to change the input current, simulation duration, and number of
channels

