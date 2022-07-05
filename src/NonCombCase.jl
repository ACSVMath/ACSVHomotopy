export find_min_crits

include("ApproxCrit.jl")

"""
  find_min_crits(h; certificates=false, show_progress=false, approx_crit=false, monodromy=false)
Find the minimal critical points of the variety defined by the vanishing of the polynomial h in the variables vars, or fail.
When certificates is true, return results of HomotopyContinuation::certify.
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
function find_min_crits(h; certificates=false, show_progress=false, approx_crit=false, monodromy=false)

  vars = variables(h)
  n = length(vars)

  if h(vars => zeros(n)) == 0
    println("No series expansion exists when h(0) = 0")
  end

  if approx_crit
    return find_min_crits_approx(h; certificates=certificates, show_progress=show_progress, monodromy=monodromy)
  end

  if monodromy
    return find_min_crits_approx(h; certificates=certificates, show_progress=show_progress, monodromy=monodromy)
  end

  # solves a simpler system to get the critical points
  crits = get_crits(h, n, vars)
  if length(crits) == 0
    println("No Critical Points")
    return []
  end

  # split h = u + v*im in the variables x, y
  u, v, x, y = split(h, n, vars)
  @assert (h(vars => x + im.*y) == u + im*v)

  @polyvar a[1:n] b[1:n] λᴿ λᴵ t ν μ ν₁ ν₂
  criteqs = [ u(x => a, y => b); v(x => a, y => b);
              a.*differentiate(u(x => a, y => b), a) + b.*differentiate(u(x => a, y => b), b) .- λᴿ;
              a.*differentiate(v(x => a, y => b), a) + b.*differentiate(v(x => a, y => b), b) .- λᴵ ]
  circeqs = x.^2 + y.^2 - t.*(a.^2 + b.^2)
  J2eqs = (ν₁.*y - ν₂.*x) .* differentiate(u, x) - (ν₁.*x + ν₂.*y) .* differentiate(u, y)

  sys1 = System([u; v; circeqs; criteqs; subs(J2eqs[2:n], [ν₁, ν₂] => [0, 1])], variables=[a; b; x; y; λᴿ; λᴵ; t])
  sys2 = System([u; v; circeqs; criteqs; subs(J2eqs[2:n], [ν₁, ν₂] => [1, 0])], variables=[a; b; x; y; λᴿ; λᴵ; t])
  sys3 = System([u; v; circeqs; criteqs; subs(J2eqs, [ν₁, ν₂] => [1, ν]); 1 - ν*μ], variables=[a; b; x; y; λᴿ; λᴵ; t; ν; μ])
  certs = [get_certified_sols(sys1, show_progress); get_certified_sols(sys2, show_progress); get_certified_sols(sys3, show_progress)]

  # assume every point is minimal
  minimal = [true for _ in 1:length(crits)]

  for (idx, crit) in enumerate(crits)
    sols = filter(cert -> compare_acb(
      real.(certified_solution_interval_after_krawczyk(cert))[1:n] + im.*real.(certified_solution_interval_after_krawczyk(cert))[n+1:2*n], # solution a, b vals
      certified_solution_interval_after_krawczyk(crit)[1:n] # critical point
    ), certs)
    tidx = 4*n+3
    if any([zero_to_one(real(certified_solution_interval_after_krawczyk(sol)[tidx])) for sol in sols])
      minimal[idx] = false
    end
  end

  if certificates
    return [crits[i] for i in 1:length(crits) if minimal[i]]
  else
    return [is_real(crits[i]) ? real(solution_approximation(crits[i]))[1:n] : solution_approximation(crits[i])[1:n] for i in 1:length(crits) if minimal[i]]
  end
end

