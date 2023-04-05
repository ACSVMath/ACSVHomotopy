export find_min_crits

include("ApproxCrit.jl")

"""
  find_min_crits(h; show_progress=false, approx_crit=false, monodromy=false)
Find the minimal critical points of the variety defined by the vanishing of the polynomial h in the variables vars, or fail.
Optionally, r may be given to specify a diagonal other than the main diagonal.
When show_progress is true, display progress of HomotopyContinuation solvers.
When approx_crit is true, use an approximation to the critical points in place of the algebraic representations to solve smaller systems.
When monodromy is true, use the monodromy method.

Example:
julia> using ACSVHomotopy

julia> @polyvar x y
(x, y)

julia> h = 1-x-y
-x - y + 1

julia> minimal = find_min_crits(h)
1-element Array{Array{Float64,1},1}:
 [0.5, 0.5]

"""
function find_min_crits(h::DynamicPolynomials.Polynomial; r=Nothing, show_progress=false, approx_crit=false, monodromy=false)

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
  end

  # solves a simpler system to get the critical points
  crits = get_crits(h, r)
  crit_intervals = certified_solution_interval_after_krawczyk.(crits)

  if length(crits) == 0
    println("No Critical Points")
    return []
  end

  # use alternative method for minimal points
  if approx_crit || monodromy
    return find_min_crits_approx(h, r; show_progress=show_progress, monodromy=monodromy)
  end

  # split h = u + v*im in the variables x, y
  u, v, x, y = split(h)
  @assert (h(z => x + im.*y) == u + im*v)

  a, b, λᴿ, λᴵ, s, t, ν, μ, ν₁, ν₂ = @polyvar a[1:n] b[1:n] λᴿ λᴵ s t ν μ ν₁ ν₂
  criteqs = [ u(x => a, y => b); v(x => a, y => b);
              a.*differentiate(u(x => a, y => b), a) + b.*differentiate(u(x => a, y => b), b) .- r*λᴿ;
              a.*differentiate(v(x => a, y => b), a) + b.*differentiate(v(x => a, y => b), b) .- r*λᴵ ]
  circeqs = x.^2 + y.^2 - t.*(a.^2 + b.^2)
  J2eqs = (ν₁.*y - ν₂.*x) .* differentiate(u, x) - (ν₁.*x + ν₂.*y) .* differentiate(u, y)

  # equations, systems and solutions for solving
  tidx = 4*n+3
  eqs = [
    [u; v; circeqs; criteqs; subs(J2eqs, [ν₁, ν₂] => [1, ν]); s*(1-t) - 1],
    [u; v; circeqs; criteqs; subs(J2eqs[2:n], [ν₁, ν₂] => [0, 1]); s*(1-t) - 1],
  ]

  vars = [
    [a; b; x; y; λᴿ; λᴵ; t; s; ν],
    [a; b; x; y; λᴿ; λᴵ; t; s]
  ]
  system = [
    System(eqs[1], variables=vars[1]),
    System(eqs[2], variables=vars[2]),
  ]

  certs = [
    filter(cert -> is_real(cert), get_certified_sols(system[1], show_progress=show_progress)),
    filter(cert -> is_real(cert), get_certified_sols(system[2], show_progress=show_progress)),
  ]

  intervals = [
    certified_solution_interval_after_krawczyk.(certs[1]),
    certified_solution_interval_after_krawczyk.(certs[2])
  ]

  # assume every point is minimal then:
  # for each critical point, check for corresponding solutions of sys1, sys2, sys3 which have t ∈ (0, 1)
  minimal = [true for _ in 1:length(crits)]
  for idx in 1:length(crits)

    # check each of the 2 systems
    for j in 1:2
      # solutions matching crit[idx] on first n coords
      matched = filter(I₀ -> compare_acb(
        real.(I₀)[1:n] + im.*real.(I₀)[n+1:2*n], # solution a, b vals
        crit_intervals[idx][1:n] # critical point
      ), intervals[j])

      # tighten the intervals so none have a 1 as a possible t value, then check t ∈ (0, 1)
      matched = map(I₀ -> refine_t(eqs[j], vars[j], I₀, tidx), matched)
      t_in_zero_to_one = filter(I₀ -> zero_to_one(real(I₀[tidx])), matched)

      if length(t_in_zero_to_one) > 0
        minimal[idx] = false
      end
    end
  end

  return [is_real(crits[i]) ? real(solution_approximation(crits[i]))[1:n] : solution_approximation(crits[i])[1:n] for i in 1:length(crits) if minimal[i]]
end
