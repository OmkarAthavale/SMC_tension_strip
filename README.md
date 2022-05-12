# minimum_project
This project uses a small 2D mesh to check the computation of slow wave electrophysiology in simplified Imtiaz ICC models (Du et al., 2013, PMID: 24111242)

## Milestones (branch)
For 2D and 3D meshes.

1) Simplified Imtiaz model (master)
2) With checkpointing (chkpt)
3) With neural data (inp_neural)
4) With checkpointing and neural data (neural_chkpt)

## Progress
|***Passing TestMinimal for milestone number***|Series|Parallel|
|--|--|--|
|Laptop|1<sub>2D,3D</sub>, 2<sub>2D,3D</sub>||
|HPC|1<sub>2D</sub>, 2<sub>2D</sub>|1<sub>2D,3D</sub>, 2<sub>2D,3D</sub>|

## Planned neural input implementation
This approach can avoid issues with additional functions in derived cell classes and serialisation issues with the neural data. Should be able to serialise the singleton ParamConfig class separate from the Chaste checkpoint BidomainProblem function.

- *Python neural preprocessing script* takes histogram-containing txt file from NEURON output
  - Converts firing rates to parameter values using calibration equations
  - Identifies regions that change at each histogrammed time bin
  - Outputs txt file containing a tidy table where each row is (time, ctrlRegionNum, parameterName, parameterValue)
- *ParamConfig singleton class* is set up to load and match neural data to mesh before solving Chaste bidomain problem
  - Read neural data tidy table
  - Match neural control regions to node global ids
  - Generate list of (globalNodeIndex, parameterName, parameterValue) triplets for a given simulation time. Should only include those that have changed between since the last sim time step to be efficient. The simulation timestep is shorter than the time dimension bins from NEURON histogram output (0.1 vs 2 ms). 
  - Be serialisable to store this data alongside checkpoint
- *BidomainProblemNeural* is a child class of BidomainProblem
  - Redefines AtBeginningOfTimestep function to update cell paramters via SetParameter based on triplets from ParamConfig for a given time
  - Also calls parent's BidomainProblem::AtBeginningOfTimestep to update electrodes.

## Notes
- Laptop has Chaste 2021.1 in Docker on [WSL](https://docs.microsoft.com/en-us/windows/wsl/install) using standard [Chaste Docker](https://github.com/Chaste/chaste-docker) instructions. Passes all tests, including parallel tests, including after rebuild using build_chaste.sh script supplied in Docker image. HPC has Chaste 2019.1 and does not pass all tests in the Parallel test pack. However, it does successfully run some custom project tests in parallel, with checkpointing (ICC3D_Longit).
- chaste_codegen is used to generate Du2013_neural derived cell classes. This seems to (so far!) be backward compatible with Chaste 2019.1 on HPC, and the 2021.1 release notes do not that that chaste_codegen and PyCml outputs are not backward+forward compatible.
- Du2013_neural is generated from Du2013_neural.cellml (.txt in repo to avoid compilation conflicts). This file has certain parameters annotated as modifiable-parameters, derived-output or special Chaste-recognised parameters (eg V_m, C_m, i_Ca). Parameters match Imtiaz2002d_noTStart.cellml

## Future work
- Have *Python neural preprocessing script* work as part of NEURON output procedures OR call the preprocessing script at runtime in Chaste
  - This would need information about the control region grid to be passed between scripts/programmes automatically not by manual entry
