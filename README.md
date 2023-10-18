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
The Fisher information is obtained by averaging over many realisations the quantity $\text{Tr}[\xi_t]^2$. In order to do this, it is first necessary to fill in the gaps left by the Gillespie evolution. The function `fisher_at_time_on_trajectory` is responsible for that.

    function fisher_at_time_on_trajectory(
        t_range::StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}, Int64},
        relevant_times::StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}, Int64},
        V::Vector{Matrix{ComplexF64}},
        Vdot::Vector{Matrix{ComplexF64}},
        trajectory_data::Vector{Dict{String, Any}},
        ψ0::Vector{ComplexF64})

#### Returns
The function returns a vector of floats, containing the (stochastic) Fisher information at all times specified in the `relevant_times` argument, on a specific trajectory.

#### Arguments
- `t_range` is the range of the times for which the `V` and `Vdot` arguments are given. It can coincide with `relevant_times`, and it does in most cases, but it is not compulsory.
- `relevant_times` is the range of times at which the calculation of the Fisher information is requested. 
- `V` is the list of no-jump evolution operators at all times specified in `t_range`.
- `Vdot` is the list of the derivatives of the no-jump evolution operators with respect to $\theta$ at all times specified in `t_range`.
- `trajectory_data` is a list of dictionaries of the same form as the output of the `gillespie_fisher` function.
- `ψ0` is the initial state, given as a state vector (ket).

### Calculating the monitoring operator for different parameter values
When dealing with MLE estimation, it is important to be able to evolve the monitoring operator along the trajectory for different values of $\theta$, in order to then be able to minimise over $\text{Tr}[\xi_t]^2$. In this case, the trajectory is fixed and given as an input.

    function evolve_ξ_on_trajectory(θ::Float64, 
        trajectory::Vector{Dict{String, Any}}, 
        func_H_θ::Function, 
        func_M_l_θ::Function, 
        ψ0::Vector{ComplexF64},
        t_final::Float64, 
        dt::Float64, 
        dθ::Float64)

#### Returns
The monitoring operator `ξ` at the final time `t_final`, computed for the parameter `θ`.

#### Arguments
- `θ` is the value of the parameter for which the evolution has to be calculated.
- `trajectory` is the list of dictionaries that specifies a trajectory, of the same format of one of the components of the output of the `gillespie_fisher` function.
- `func_H_θ` is a function that gives the Hamiltonian, when given the parameter $\theta$ as input.
- `func_M_l_θ` is a function that gives the list of jump operators, when given the parameter $\theta$ a input.
- `ψ0` is the initial state, given as a state vector (ket).
- `t_final` is the final time of the evolution, at which the monitoring operator will be calculated.
- `dt` is the timestep for the evolution.
- `dθ` is the step in $\theta$ to compute the derivatives.

## Mixed states, partial monitoring, merging channels
If the initial state is not pure, the jump operators are not rank-one projectors or partial monitoring is considered, then the complete formalism of superoperators has to be deployed.

### Evolving the monitoring operator
