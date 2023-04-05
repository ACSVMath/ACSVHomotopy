export find_min_crits_comb

"""
  find_min_crits_comb(h; r=Nothing, show_progress=false)
Find the minimal critical points of the variety defined by the vanishing of the polynomial h in the variables vars, or fail.
Optionally, r may be given to specify a diagonal other than the main diagonal.
When show_progress is true, display progress of HomotopyContinuation solvers.

Example:
julia> using ACSVHomotopy

julia> @polyvar x y
(x, y)

julia> h = 1-x-y
-x - y + 1

julia> minimal = find_min_crits_comb(h)
1-element Array{Array{Float64,1},1}:
 [0.5, 0.5]

"""
function find_min_crits_comb(h::DynamicPolynomials.Polynomial; r=Nothing, show_progress=false)

    z = variables(h)
    n = length(z)

    if r == Nothing
      r = ones(n)
    else
      if !(typeof(r) <: Vector{<:Number}) || length(r) != n
        println("r should be a vector with length equal to the number of variables in h")
        return
      end
    end

    if h(z => zeros(n)) == 0
      println("No series expansion exists when h(0) = 0")
      return
    end

    crits = get_crits(h, r)
    crit_intervals = certified_solution_interval_after_krawczyk.(crits)

    # HC.jl gives certified real solutions
    certified_real = certified_solution_interval_after_krawczyk.(filter(is_real, crits))

    if length(certified_real) == 0
      println("No Real Positive Critical Points.")
      return
    end

    # Build the CombCase system & solve
    λ, s, t = @polyvar λ s t
    tidx = n+3 # get the real solutions to sys with t ∉ (0, 1)
    vars =[z; λ; s; t]
    eqs = [z .* differentiate(h, z) .- r*λ; h; h(z => t.*z); s*(1-t) - 1]
    sys = System(eqs, variables=vars);
    certs = get_certified_sols(sys; show_progress=show_progress)

    # each interval I₀ ∈ `sol_intervals` contains a distinct real solution to `sys`
    sol_intervals = certified_solution_interval_after_krawczyk.(certs)

    # Check that solution intervals are able to separate all roots
    for I₀ in sol_intervals
      # if a certified solution is close to multiple critical points, then HC.jl was unable to seperate the solutions
      if length(filter(crit -> compare_acb(crit[1:n], I₀[1:n]), crit_intervals)) > 1
        println("Fail! Unable to separate solutions!")
        return
      end
    end

    # refine sol_intervals and check for sols where t ∈ (0, 1)
    sol_intervals = map(I₀ -> refine_t(eqs, vars, I₀, tidx), sol_intervals)
    t_in_zero_to_one = filter(I₀ -> zero_to_one(real(I₀[tidx])), sol_intervals)

    # Filter those critical points which are real & positive (the n+1'th coord is λ)
    real_positive_crit_intervals = filter(z₀ -> all(0 .< real.(z₀[1:n])), certified_real)

    # ζ is the unique minimal point, there are no solutions σ of sys such that t ∈ (0, 1) and σ agrees with ζ on first n coords
    ζ = filter(crit -> !any([compare_acb(crit[1:n], sol[1:n]) for sol in t_in_zero_to_one]), real_positive_crit_intervals)

    if length(ζ) == 0
      println("Fail! No Minimal Points.")
      return
    end

    if length(ζ) > 1
      println("Fail! There are multiple real minimal points.")
      return
    end

    ζ = ζ[1] # there is exactly one real positive minimal point
    crit_intervals = filter(I -> !Arblib.overlaps(I, ζ), crit_intervals)

    # TODO: figure out what number to put here
    maxcoeff = maximum(coefficients(h))
    deg = degree(leadingmonomial(h))
    eta = log2(maxcoeff) + deg * log2(1+n)
    Hn = deg^n + eta*deg^(n-1)
    Cn = deg^n
    B = Hn + (log2(max(deg^n*(deg^n-1), 1)) + 4*log2(deg+3))*Cn
    println("B", B)

    minimal = [ζ]
    for I in crit_intervals
      if Arblib.overlaps(ArbMatrix(abs.(I)), ArbMatrix(abs.(ζ)))
        refine([h; z .* differentiate(h, z) .- r*λ], [z; λ], AcbMatrix(I), 2^(-B))
        refine([h; z .* differentiate(h, z) .- r*λ], [z; λ], AcbMatrix(ζ), 2^(-B))
      end

      if Arblib.overlaps(ArbMatrix(abs.(I)), ArbMatrix(abs.(ζ)))
        push!(minimal, I)
      end
    end
    minimal = filter(crit -> compare_acb(abs.(certified_solution_interval_after_krawczyk(crit)[1:n]), ζ[1:n]), crits) # the other minimal points have same coordwise modulus

    return [is_real(w) ? real(solution_approximation(w))[1:n] : solution_approximation(w)[1:n] for w in minimal]
end
