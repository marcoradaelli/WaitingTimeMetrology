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

function vectorise(x)
    return vec(x)
end

function unvectorise(x)
    d = trunc(Int, sqrt(length(x)))
    return reshape(x, (d,d))
end

function gillespie_fisher_mixed(
        H::Matrix{ComplexF64},
        Hp::Matrix{ComplexF64},
        Hm::Matrix{ComplexF64},
        M_l_l::Vector{Vector{Matrix{ComplexF64}}},
        Mp_l_l::Vector{Vector{Matrix{ComplexF64}}},
        Mm_l_l::Vector{Vector{Matrix{ComplexF64}}},
        S_l::Vector{Matrix{ComplexF64}},
        Sp_l::Vector{Matrix{ComplexF64}},
        Sm_l::Vector{Matrix{ComplexF64}},
        ρ0::Matrix{ComplexF64},
        t_final::Float64,
        dt::Float64,
        dθ::Float64,
        number_trajectories::Int64,
        verbose::Bool=false,
        ξ0=nothing)
    
    display("Starting Gillespie evolution")
    
    # Range of times.
    t_range = 0.:dt:t_final
    
    # M_l_l is a list of lists. Each sub-list represents jump operators that have to be merged.
    
    # Creates the appropriate density matrix.
    d = size(H)[1]
    ide = Matrix{Float64}(I, d, d)
    
    # Creates the J operator. The overall J operator is independent on how the jumps are grouped together.
    J = zero(ide)
    for M_l in M_l_l
        for M in M_l
            J += M' * M
        end
    end
    
    # Vectorized version of J.
    vect_J = vectorise(J)
    
    # Vectorized version of L_0 (no-jump evolution).
    # Hamiltonian part.
    vect_L0 = - 1im * kron(ide, H) + 1im * kron(transpose(H), ide)
    # Non-Hamiltonian part: unmonitored operators.
    for S in S_l
        vect_L0 += kron(conj.(S), S) - 0.5 * kron(ide, S' * S) - 0.5 * kron(transpose(S' * S), ide)
    end
    # Non-Hamiltonian part: monitored operators.
    for M_l in M_l_l
        for M in M_l
            vect_L0 += - 0.5 * kron(ide, M' * M) - 0.5 * kron(transpose(M' * M), ide)
        end
    end
    
    # Vectorized version of L_0_dagger.
    # Hamiltonian part.
    vect_L0_dagger = 1im * kron(ide, H) - 1im * kron(transpose(H), ide)
    # Non-Hamiltonian part: unmonitored operators.
    for S in S_l
        vect_L0_dagger += kron(transpose(S), S') - 0.5 * kron(transpose(S' * S), ide) - 0.5 * kron(ide, S' * S)
    end
    # Non-Hamiltonian part: monitored operators.
    for M_l in M_l_l
        for M in M_l
            vect_L0_dagger += - 0.5 * kron(transpose(M' * M), ide) - 0.5 * kron(ide, M' * M)
        end
    end
    
    # Pre-computation stage.
    # Construct the no-jump evolution for all times in the range.
    V = Matrix{ComplexF64}[]
    Qs = Matrix{ComplexF64}[]
    @showprogress "Precomputing no-jumps..." for t in t_range
        ev_op = exp(vect_L0 * t)
        push!(V, ev_op)
        nsd_wtd = unvectorise(exp(vect_L0_dagger * t) * vect_J)
        push!(Qs, nsd_wtd)
    end
    
    # Construct the derivative of V w.r.t. the parameter.
    vect_L0p = - 1im * kron(ide, Hp) + 1im * kron(transpose(Hp), ide)
    # Non-Hamiltonian part: unmonitored operators.
    for Sp in Sp_l
        vect_L0p += kron(conj.(Sp), Sp) - 0.5 * kron(transpose(Sp' * Sp), ide) - 0.5 * kron(ide, Sp' * Sp)
    end
    # Non-Hamiltonian part: monitored operators.
    for Mp_l in Mp_l_l
        for Mp in Mp_l
            vect_L0p += - 0.5 * kron(ide, Mp' * Mp) - 0.5 * kron(transpose(Mp' * Mp), ide)
        end
    end
    
    vect_L0m = - 1im * kron(ide, Hm) + 1im * kron(transpose(Hm), ide)
    # Non-Hamiltonian part: unmonitored operators.
    for Sm in Sm_l
        vect_L0m += kron(conj.(Sm), Sm) - 0.5 * kron(transpose(Sm' * Sm), ide) - 0.5 * kron(ide, Sm' * Sm)
    end
    # Non-Hamiltonian part: monitored operators.
    for Mm_l in Mm_l_l
        for Mm in Mm_l
            vect_L0m += - 0.5 * kron(ide, Mm' * Mm) - 0.5 * kron(transpose(Mm' * Mm), ide)
        end
    end
    Vdot = Matrix{ComplexF64}[]
    @showprogress "Precomputing derivatives..." for t in t_range
        ev_op_p = exp(vect_L0p * t)
        ev_op_m = exp(vect_L0m * t)
        derivative = (ev_op_p - ev_op_m) / (2 * dθ)
        push!(Vdot, derivative)
    end
        
    # Computes all the J maps corresponding to the distinct jumps.
    J_l = Matrix{ComplexF64}[]
    for M_l in M_l_l
        J_dist = zero(kron(ide, ide))
        for M in M_l
            J_dist += kron(conj(M), M)
        end
        push!(J_l, J_dist)
    end
    
    display("Number of created jump maps: " * string(length(J_l)))
        
    # Computes also all the derivatives w.r.t. θ of the J maps.
    Jdot_l = Matrix{ComplexF64}[]
    for i in eachindex(M_l_l)
        Jdot_dist = zero(kron(ide, ide))
        for j in eachindex(M_l_l[i])
            Jp_dist = kron(conj(Mp_l_l[i][j]), Mp_l_l[i][j])
            Jm_dist = kron(conj(Mm_l_l[i][j]), Mm_l_l[i][j])
            Jdot_dist += (Jp_dist - Jm_dist) / (2 * dθ)
        end
        push!(Jdot_l, Jdot_dist)
    end
    
    # Array for final results.
    trajectories_results = Array{Dict{String, Any}}[]
        
    # Proper Gillespie evolution.
    @showprogress 1 "Gillespie evolution..." for trajectory in 1:number_trajectories 
        # Initial state.
        ρ = ρ0
        # Initial monitoring operator.
        if ξ0 == nothing
            ξ = zero(ρ)
        else
            ξ = ξ0
        end
        # Absolute time.
        τ = 0
        
        # Creates the array of results for the single trajectory, and pushes the initial state as a fictitious first jump.
        results = Dict{String, Any}[]
        dict_initial = Dict("AbsTime" => 0,
            "TimeSinceLast" => 0,
            "JumpChannel" => nothing,
            "ρAfter" => ρ0,
            "ξAfter" => ξ)
        push!(results, dict_initial)
        
        while τ < t_final
            dict_jump = Dict()
            
            # Compute the waiting time distribution, exploiting the pre-computed part.
            Ps = Float64[]
            for Q in Qs
                wtd = real(tr(Q * ρ))
                push!(Ps, wtd)
            end
            
            # Samples from the waiting time distribution the next jump time.
            n_T = sample(1:length(t_range), Weights(Ps))
            
            # Increases the absolute time, and save in dictionary.
            τ += t_range[n_T]
            merge!(dict_jump, Dict("AbsTime" => τ, "TimeSinceLast" => t_range[n_T]))
            
            # Stores the old value of ρ, before starting to edit it.
            pre_ρ = ρ
            vect_pre_ρ = vectorise(ρ)
            
            # Updates the state.
            vect_ρ = V[n_T] * vectorise(ρ)
            ρ = unvectorise(vect_ρ)
                
            # Chooses where to jump.
            weights = Float64[]
            for J_dist in J_l
                weight = real(tr(unvectorise(J_dist * vect_ρ)))
                push!(weights, weight)
            end
            n_jump = sample(1:length(J_l), Weights(weights))
            merge!(dict_jump, Dict("JumpChannel" => n_jump))
            
            # Updates the state after the jump.
            vect_ρ = J_l[n_jump] * vect_ρ
            ρ = unvectorise(vect_ρ)
            norm_state = real(tr(ρ))
            # Renormalises the state.
            ρ = ρ / norm_state
            vect_ρ = vectorise(ρ)
            merge!(dict_jump, Dict("ρAfter" => ρ))
            
            # Evolves the monitoring operator.
            vect_ξ = vectorise(ξ)
            vect_ξ =  (Jdot_l[n_jump] * V[n_T] * vect_pre_ρ + J_l[n_jump] * Vdot[n_T] * vect_pre_ρ + J_l[n_jump] * V[n_T] * vect_ξ)            
            
            ξ = unvectorise(vect_ξ) / norm_state
            
            merge!(dict_jump, Dict("ξAfter" => ξ))
            
            if verbose
                println(string(dict_jump))
            end
            
            push!(results, dict_jump)
        end
        
        push!(trajectories_results, results)
    end
    
    display("All good");
    
    return trajectories_results, V, Vdot, t_range
end

function fisher_at_time_on_trajectory_mixed(
    t_range::StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}, Int64},
    relevant_times::StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}, Int64},
    V::Vector{Matrix{ComplexF64}},
    Vdot::Vector{Matrix{ComplexF64}},
    trajectory_data::Vector{Dict{String, Any}},
    ψ0::Matrix{ComplexF64})

    v_fisher = Float64[]

    # Creates an array of jump times.
    jump_times = [trajectory_data[i]["AbsTime"] for i in eachindex(trajectory_data)]
    # Creates an array of states after the jumps.
    ρ_after_jumps = [trajectory_data[i]["ρAfter"] for i in eachindex(trajectory_data)]
    # Creates an array of ξ after the jumps.
    ξ_after_jumps = [trajectory_data[i]["ξAfter"] for i in eachindex(trajectory_data)]

    # Cycles over the jumps times.
    for n_jump in 1:length(jump_times)-1
        next_jump_time = jump_times[n_jump + 1]
        # Determines the set of relevant times between this jump and the following one.
        relevant_times_in_interval = [t for t in t_range if jump_times[n_jump] <= t < next_jump_time]
        # Cycles over the relevant times.
        for t_abs in relevant_times_in_interval
            ρ = ρ_after_jumps[n_jump]
            ξ = ξ_after_jumps[n_jump]
            vect_ρ = vectorise(ρ)
            vect_ξ = vectorise(ξ)
            ξ = ξ_after_jumps[n_jump]
            n_t = find_nearest(t_range, t_abs - jump_times[n_jump])[1]
            trace = real(tr(unvectorise(V[n_t] * vect_ρ)))
            vect_ξ = (Vdot[n_t] * vect_ρ + V[n_t] * vect_ξ) / trace
            ξ = unvectorise(vect_ξ)
            fisher = real(tr(ξ))^2
            push!(v_fisher, fisher)
        end
    end

    # Now computes the ξ for all times after the latest jump.
    last_jump_absolute_time = last(jump_times)
    relevant_times_after_last_jump = [t for t in t_range if t >= last_jump_absolute_time]
    for t_abs in relevant_times_after_last_jump
        ρ = last(ρ_after_jumps)
        ξ = last(ξ_after_jumps)
        vect_ρ = vectorise(ρ)
        vect_ξ = vectorise(ξ)
        n_t = find_nearest(t_range, t_abs - last_jump_absolute_time)[1]
        trace = real(tr(unvectorise(V[n_t] * vect_ρ)))
        vect_ξ = (Vdot[n_t] * vect_ρ + V[n_t] * vect_ξ) / trace
        ξ = unvectorise(vect_ξ)
        fisher = real(tr(ξ))^2
        push!(v_fisher, fisher)
    end

    return v_fisher
end
