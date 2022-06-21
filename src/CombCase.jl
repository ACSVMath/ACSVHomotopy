export find_min_crits_comb

"""
    find_min_crits_comb(h)
Find the minimal critical points of the variety defined by the vanishing of
the polynomial h in the variables vars, or fail.
"""
function find_min_crits_comb(h; show_progress=false)

    z = variables(h)
    n = length(z)

    crits = get_crits(h, n, z)
    if length(crits) == 0
      println("No Critical Points.")
      return
    end

    @polyvar λ μ t
    tidx = n+3 # get the real solutions to sys with t ∉ (0, 1)
    sys = System([z .* differentiate(h, z) .- λ; h; h(z => t.*z); μ*(1-t) - 1], variables=[z; λ; μ; t]);
    certs = get_certified_sols(sys, show_progress)
    certs = filter(cert -> all(real(certified_solution_interval(cert)) .> 0), certs) # every coord real and +tive
    certs = filter(cert -> zero_to_one(real(certified_solution_interval(cert)[tidx])), certs) # t ∈ (0, 1)

    # Check that solution intervals are able to separate all roots
    real_positive = filter(crit -> is_real(crit) && all(real(certified_solution_interval(crit)[1:n]) .> 0), crits)
    for cert in certs
      if length(filter(crit -> compare_acb(certified_solution_interval(crit)[1:n], certified_solution_interval(cert)[1:n]), real_positive)) > 1
        println("Fail! Unable to separate solutions!")
        return
      end
    end

    # ζ is the minimal point, there are no solutions σ of sys such that t ∈ (0, 1) and σ agrees with ζ on first n coords
    ζ = filter(crit -> !any([compare_acb(certified_solution_interval(crit)[1:n], certified_solution_interval(cert)[1:n]) for cert in certs]), real_positive)

    if length(ζ) == 0
      println("No Minimal Points.")
      return
    end

    if length(ζ) > 1
      println("Fail! There are multiple real minimal points.")
      return
    end

    ζ = certified_solution_interval(ζ[1])[1:n] # there is exactly one real positive minimal point
    minimal = solution_candidate.(filter(w -> compare_acb(abs.(certified_solution_interval(w)[1:n]), ζ), crits))
    return [w[1:n] for w in minimal]
end

