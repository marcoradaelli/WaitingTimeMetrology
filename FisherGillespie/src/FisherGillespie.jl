module FisherGillespie

using LinearAlgebra
using Random
using StatsBase
using Plots
using Polynomials
using ProgressMeter

function verify_working()
    println("The package has been imported correctly")
    return true
end

function find_nearest(a,x)
    length(a) > 0 || return 0:-1
    r = searchsorted(a,x)
    length(r) > 0 && return r
    last(r) < 1 && return searchsorted(a,a[first(r)])
    first(r) > length(a) && return searchsorted(a,a[last(r)])
    x-a[last(r)] < a[first(r)]-x && return searchsorted(a,a[last(r)])
    x-a[last(r)] > a[first(r)]-x && return searchsorted(a,a[first(r)])
    return first(searchsorted(a,a[last(r)])):last(searchsorted(a,a[first(r)]))
end

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
    verbose::Bool=True)

    t_range = 0.:dt:t_final

    # Constructs the overall jump operator.
    J = zero(M_l[1])
    for M in M_l
        J += M' * M
    end
    # Effective (non-Hermitian) Hamiltonian.
    He = H - 1im/2. * J

    # Constructs the no-jump evolution operators for all the relevant times.
    V = Matrix{ComplexF64}[] # List of the no-jump evolution operators.
    Qs = Matrix{ComplexF64}[] # List of the non-state-dependent part of the waiting time distribution.
    for t in t_range
        ev_op = exp(-1im * He * t)
        push!(V, ev_op)
        nsd_wtd = ev_op' * J * ev_op
        push!(Qs, nsd_wtd)
    end

    # Prints the matrix norm for the latest Qs.
    error = norm(last(Qs))
    println("-> Truncation error given by norm of latest Qs matrix: " * string(error))

    # Displaced version of the effective Hamiltonian.
    Jp = zero(Mp_l[1])
    for Mp in Mp_l
        Jp += Mp' * Mp
    end
    Hep = Hp - 1im/2. * Jp
    Jm = zero(Mm_l[1])
    for Mm in Mm_l
        Jm += Mm' * Mm
    end
    Hem = Hm - 1im/2. * Jm

    # Vector of the derivatives of the evolution operator wrt the parameter.
    Vdot = Matrix{ComplexF64}[]
    for t in t_range
        vd = (exp(-1im * Hep * t) - exp(-1im * Hem * t)) / (2 * dθ)
        push!(Vdot, vd)
    end
    # Derivatives of all the jump operators.
    Mdot = Matrix{ComplexF64}[]
    for n_M in eachindex(M_l)
        push!(Mdot, (Mp_l[n_M] - Mm_l[n_M]) / (2 * dθ))
    end
    # Precomputation of all the derivatives of MkV(t).
    Δ = Array{Matrix{ComplexF64}}[]
    for n_M in eachindex(M_l)
        fixed_M = Matrix{ComplexF64}[]
        for n_t in eachindex(t_range)
            obj = Mdot[n_M] * V[n_t] + M_l[n_M] * Vdot[n_t]
            push!(fixed_M, obj)
        end
        push!(Δ, fixed_M)
    end
        
    # List for the results.
    trajectories_results = Array{Dict{String, Any}}[]

    # Cycle over the trajectories.
    @showprogress 1 "Fisher-Gillespie evolution..." for trajectory in 1:number_trajectories
        
        # Initial state.
        ψ = ψ0
        # Absolute time.
        τ = 0
        # Initial ξ matrix (same size and type as the Hamiltonian).
        ξ = zero(H)
        
        results = Dict{String, Any}[]
        dict_initial = Dict("AbsTime" => 0,
            "TimeSinceLast" => 0,
            "JumpChannel" => nothing,
            "ξAfter" => ξ,
            "ψAfter" => ψ0)
        push!(results, dict_initial)
        
        while τ < t_final
            dict_jump = Dict()
            
            # Pass to density matrix formalism.
            ρ = ψ * ψ'
            
            # Compute the waiting time distribution, exploiting the pre-computed part.
            Ps = Float64[]
            for Q in Qs
                wtd = real(ψ' * Q * ψ)
                push!(Ps, wtd)
            end
            
            # Sample from the waiting time distribution.
            n_T = sample(1:length(t_range), Weights(Ps))
            
            # Increase the absolute time.
            τ += t_range[n_T]
            merge!(dict_jump, Dict("AbsTime" => τ, "TimeSinceLast" => t_range[n_T]))
            
            # Update the state.
            ψ = V[n_T] * ψ
            # Chooses where to jump.
            weights = Float64[]
            for M in M_l
                weight = real(ψ' * M' * M * ψ)
                push!(weights, weight)
            end
            n_jump = sample(1:length(M_l), Weights(weights))
            merge!(dict_jump, Dict("JumpChannel" => n_jump))
            # Update the state after the jump.
            ψ = M_l[n_jump] * ψ
            norm_state = norm(ψ)
            # Renormalize the state.
            ψ = ψ / norm_state
            
            if verbose
                println(string(dict_jump))
            end
            
            # Compute the ξ matrix.
            ξ = 1/(norm_state^2) * (M_l[n_jump] * V[n_T] * ξ * (V[n_T])' * (M_l[n_jump])')
            ξ += 1/(norm_state^2) * (Δ[n_jump][n_T] * ρ * (V[n_T])' * (M_l[n_jump])')
            ξ += 1/(norm_state^2) * (M_l[n_jump] * V[n_T] * ρ * (Δ[n_jump][n_T])')
            trξ2 = real(tr(ξ))^2
            merge!(dict_jump, Dict("ξAfter" => ξ, "ψAfter" => ψ, "trξ2" => trξ2))
            push!(results, dict_jump)
        end
        
        push!(trajectories_results, results)
    end

    return trajectories_results, V, Vdot, t_range
end

function fisher_at_time_on_trajectory(
    t_range::StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}, Int64},
    relevant_times::StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}, Int64},
    V::Vector{Matrix{ComplexF64}},
    Vdot::Vector{Matrix{ComplexF64}},
    trajectory_data::Vector{Dict{String, Any}},
    ψ0::Vector{ComplexF64})

    v_fisher = Float64[]

    # Creates an array of jump times.
    jump_times = [trajectory_data[i]["AbsTime"] for i in eachindex(trajectory_data)]
    # Creates an array of states after the jumps.
    ψ_after_jumps = [trajectory_data[i]["ψAfter"] for i in eachindex(trajectory_data)]
    # Creates an array of ξ after the jumps.
    ξ_after_jumps = [trajectory_data[i]["ξAfter"] for i in eachindex(trajectory_data)]

    # Cycles over the jumps times.
    for n_jump in 1:length(jump_times)-1
        next_jump_time = jump_times[n_jump + 1]
        # Determines the set of relevant times between this jump and the following one.
        relevant_times_in_interval = [t for t in t_range if jump_times[n_jump] <= t < next_jump_time]
        # Cycles over the relevant times.
        for t_abs in relevant_times_in_interval
            ψ = ψ_after_jumps[n_jump]
            ρ = ψ * ψ'
            ξ = ξ_after_jumps[n_jump]
            n_t = find_nearest(t_range, t_abs - jump_times[n_jump])[1]
            trace = real(tr(V[n_t] * ρ * (V[n_t])'))
            ξ = Vdot[n_t] * ρ * (V[n_t])' + V[n_t] * ρ * (Vdot[n_t])' + V[n_t] * ξ * (V[n_t])'
            ξ = ξ / trace
            fisher = real(tr(ξ))^2
            push!(v_fisher, fisher)
        end
    end

    # Now computes the ξ for all times after the latest jump.
    last_jump_absolute_time = last(jump_times)
    relevant_times_after_last_jump = [t for t in t_range if t >= last_jump_absolute_time]
    for t_abs in relevant_times_after_last_jump
        ψ = last(ψ_after_jumps)
        ρ = ψ * ψ'
        ξ = last(ξ_after_jumps)
        n_t = find_nearest(t_range, t_abs - last_jump_absolute_time)[1]
        trace = real(tr(V[n_t] * ρ * (V[n_t])'))
        ξ = Vdot[n_t] * ρ * (V[n_t])' + V[n_t] * ρ * (Vdot[n_t])' + V[n_t] * ξ * (V[n_t])'
        ξ = ξ / trace
        fisher = real(tr(ξ))^2
        push!(v_fisher, fisher)
    end

    return v_fisher
end

function compute_fisher_information(
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
    verbose::Bool=True)

    trajectories_results, V, Vdot, t_range = gillespie_fisher(H, Hp, Hm, M_l, Mp_l, Mm_l, dθ, ψ0, t_final, dt, number_trajectories, verbose)
    average_fisher = zeros(length(t_range))
    println()
    @showprogress 1 "Filling in the gaps..." for n_trajectory in eachindex(trajectories_results)
        v_fisher = fisher_at_time_on_trajectory(t_range, t_range, V, Vdot, trajectories_results[n_trajectory], ψ0)
        average_fisher += v_fisher / number_trajectories
    end

    return t_range, average_fisher
end

end # module
