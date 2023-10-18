# Parameter estimation under quantum jump unravelling
This code accompanies the manuscript [TO APPEAR]. If you are not interested in the metrological aspect, but only in the simulation of the quantum jump process with Gillespie, you might want to check out our Gillespie library [here](https://github.com/marcoradaelli/GillespieQuantumJumps), and the associated preprint [here](https://arxiv.org/abs/2303.15405).

## Installation
At the moment, the package is not (yet) a registered Julia package. To install it, in the Package REPL give the command

    add "https://github.com/marcoradaelli/WaitingTimeMetrology/"
The package manager should automatically find and install the required dependencies.

To include the package in your Julia file, write

    using FisherGillespie
All the functions detailed below should then be called by using the syntax

    FisherGillespie.[function name here]([arguments here])

## Pure states evolution
If all the jumps are rank-1 projectors, the initial state is pure and the monitoring is complete, then the conditional quantum state of the system will always be pure. 

### Evolving the monitoring operator
The function `gillespie_fisher` allows to evolve the monitoring operator, along with the conditional quantum state, on quantum jump trajectories. The function's prototype is the following:

    function gillespie_fisher(
        H::Matrix{ComplexF64},
        Hp::Matrix{ComplexF64},
        Hm::Matrix{ComplexF64},
        M_l::Vector{Matrix{ComplexF64}},
        Mp_l::Vector{Matrix{ComplexF64}},
        Mm_l::Vector{Matrix{ComplexF64}},
        dθ::Float64,
        ψ0::Vector{ComplexF64},
        t_final::Float64,
        dt::Float64,
        number_trajectories::Int64,
        verbose::Bool=true)

#### Returns 
The function returns a tuple of four elements. 
- `trajectories_results`, is a list of lists of dictionaries. Each of the elements of the outer list corresponds to a trajectory, and each trajectory contains a list of dictionaries; each dictionary corresponds to a quantum jump. For simplicity of implementation, a first jump in channel `nothing` is always added at time zero.
-  `V` is the list of the no-jump evolution operators, pre-computed at `dt` intervals from zero to `t_final`. 
- `V_dot` is the list of derivatives of such `V` operators with respect to the parameter $\theta$. 
- `t_range` is the Julia range `0:dt:t_final`.

#### Arguments
- `H` is the Hamiltonian of the system, calculated for the true value of the parameter $\theta$. 
-  `Hp` is the Hamiltonian of the system, calculated for the shifted value $\theta + d\theta$.
- `Hm` is the Hamiltonian of the system, calculated for the shifted value $\theta - d\theta$.
- `M_l` is the list of the jump operators, calculated for the true value of $\theta$.
- `Mp_l` is the list of the jump operators, calculated for the shifted value $\theta + d\theta$.
- `Mm_l` is the list of the jump operators, calculated for the shifted value $\theta - d\theta$.
- `dθ` is the increment in $\theta$, used to compute the numerical derivatives.
- `ψ0` is the initial state, given in the form of a state vector (ket).
- `t_final` is the final time of the evolution.
- `dt` is the evolution's time step.
- `number_trajectories` is the total number of trajectories that have to be pre-computed.
- `verbose` is an optional parameter, defaulting at `false`. When `verbose` is set to `true`, some more output may be shown during the evolution. **Warning!** When a large number of trajectories or a long evolution time is employed, setting `verbose` to true may potentially crash the program.

### Obtaining the Fisher information
The Fisher information is obtained by averaging over the 