"""
  find_min_crits_approx(h; show_progress=false, monodromy=false)
Find the minimal critical points of the variety defined by the vanishing of
the polynomial h in the variables vars, or fail.
"""
function find_min_crits_approx(h; show_progress=false, monodromy=false)
  vars = variables(h)
  n = length(vars)

  # solves a simpler system to get the critical points
  crits = get_crits(h, n, vars)

  # split h = u + v*im in the variables x, y
  u, v, x, y = split(h, n, vars)
  @assert (h(vars => x + im.*y) == u + im*v)

  # build the system & solve
  @polyvar a[1:n] b[1:n] λᴿ λᴵ t ν μ ν₁ ν₂
  circeqs = x.^2 + y.^2 - t.*(a.^2 + b.^2)
  J2eqs = (ν₁.*y - ν₂.*x) .* differentiate(u, x) - (ν₁.*x + ν₂.*y) .* differentiate(u, y)

  # assume every point is minimal
  minimal = [true for _ in 1:length(crits)]

  # run monodromy on each crit to determine minimality
  for (idx, crit) in enumerate(crits)
    abvals = [real(solution_candidate(crit)[1:n]); imag(solution_candidate(crit)[1:n])]

    # sys1: ν₁ = 0 in J2eqs
    sys1 = System([u; v; circeqs; subs(J2eqs[2:n], [ν₁, ν₂] => [0, 1])], variables=[x; y; t], parameters=[a; b])
    minimal[idx] &= minimality_check(sys1, abvals, 2*n+1; show_progress=show_progress, monodromy=monodromy)

    # sys2: ν₂ = 0 in J2eqs
    sys2 = System([u; v; circeqs; subs(J2eqs[2:n], [ν₁, ν₂] => [1, 0])], variables=[x; y; t], parameters=[a; b])
    minimal[idx] &= minimality_check(sys2, abvals, 2*n+1; show_progress=show_progress, monodromy=monodromy)

    # sys3: ν₁ and ν₂ are non-zero in J2eqs
    sys3 = System([u; v; circeqs; subs(J2eqs, [ν₁, ν₂] => [1, ν]); 1 - ν*μ], variables=[x; y; t; ν; μ], parameters=[a; b])
    minimal[idx] &= minimality_check(sys3, abvals, 2*n+1; show_progress=show_progress, monodromy=monodromy)
  end

  return [solution_candidate(crits[i])[1:n] for i in 1:length(crits) if minimal[i]]
end

