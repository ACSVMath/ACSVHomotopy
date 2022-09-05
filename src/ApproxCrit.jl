# TRUE if any of the solutions in `sols` has t ∈ (0, 1)
function check_solutions(
  sys::System,
  eqs::Vector{<:DynamicPolynomials.Polynomial{true, <:Number}},
  vars::Vector{DynamicPolynomials.PolyVar{true}},
  params::Vector{DynamicPolynomials.PolyVar{true}},
  target_params::Vector{<:Number},
  sols::Vector{<:Vector{<:Number}},
  tidx::Integer
)

eqs = subs(eqs, params => target_params)
certs = distinct_certificates(certify(sys, sols, target_params))
certs = filter(cert -> is_certified(cert), certs)
intervals = certified_solution_interval_after_krawczyk.(certs)
intervals = map(I₀ -> refine_t(eqs, vars, I₀, tidx), intervals)
t_in_zero_to_one = filter(I₀ -> zero_to_one(real(I₀[tidx])), intervals)

return length(t_in_zero_to_one) > 0 # at least one sol has t ∈ (0, 1)
end

# Helper function for performing minimality check in ApproxCrit method
# FALSE if the system has no real solutions with t ∈ (0, 1)
function minimality_check(
  eqs::Vector{<:DynamicPolynomials.Polynomial{true, <:Number}},
  vars::Vector{DynamicPolynomials.PolyVar{true}},
  params::Vector{DynamicPolynomials.PolyVar{true}},
  target_params::Vector{<:Number},
  tidx::Integer;
  show_progress=false,
  monodromy=false
)

# solve the system using the appropriate method, then check for solutions with t ∈ (0, 1)
sys = System(eqs, variables=vars, parameters=params)
if monodromy
  init = get_start_solution(sys, target_params)
  if !isnothing(init)
    res = results(monodromy_solve(sys, init, target_params, show_progress=show_progress))
    sols = real_solutions(res)
    if length(sols) > 0
      # if there exists a solution with t ∈ (0, 1), then the point is not minimal
      return !check_solutions(sys, eqs, vars, params, target_params, sols, tidx)
    end
  end
else
  res = results(solve(sys; target_parameters=target_params, show_progress=show_progress))
  sols = real_solutions(res)
  if length(sols) > 0
    # if there exists a solution with t ∈ (0, 1), then the point is not minimal
    return !check_solutions(sys, eqs, vars, params, target_params, sols, tidx)
  end
end

# the point still may or may not be minimal (false positives are possible, the minimality check is inconclusive in this case)
# because we can not verify that HC.jl has found all the solutions
return true
end

function find_min_crits_approx(h::DynamicPolynomials.Polynomial, r::Vector{<:Number}; show_progress=false, monodromy=false)

vars = variables(h)
n = length(vars)

# solves a simpler system to get the critical points
crits = get_crits(h, r)

# split h = u + v*im in the variables x, y
u, v, x, y = split(h)
@assert (h(vars => x + im.*y) == u + im*v)

# build the system & solve
a, b, λᴿ, λᴵ, s, t, ν, μ, ν₁, ν₂ = @polyvar a[1:n] b[1:n] λᴿ λᴵ s t ν μ ν₁ ν₂
circeqs = x.^2 + y.^2 - t.*(a.^2 + b.^2)
J2eqs = (ν₁.*y - ν₂.*x) .* differentiate(u, x) - (ν₁.*x + ν₂.*y) .* differentiate(u, y)

tidx = 2n+1
eqs = [
  [u; v; circeqs; subs(J2eqs[2:n], [ν₁, ν₂] => [0, 1]); s*(1-t) - 1],
  [u; v; circeqs; subs(J2eqs[2:n], [ν₁, ν₂] => [1, 0]); s*(1-t) - 1],
  [u; v; circeqs; subs(J2eqs, [ν₁, ν₂] => [1, ν]); 1 - ν*μ; s*(1-t) - 1]
]

vars = [
  [x; y; t; s],
  [x; y; t; s],
  [x; y; t; s; ν; μ]
]

# assume every point is minimal
minimal = [true for _ in 1:length(crits)]

# run monodromy on each crit to determine minimality
for (idx, crit) in enumerate(crits)
  abvals = [real(solution_approximation(crit)[1:n]); imag(solution_approximation(crit)[1:n])]

  # run the minimality check on each system
  minimal[idx] &= minimality_check(eqs[1], vars[1], [a; b], abvals, tidx; show_progress=show_progress, monodromy=monodromy)
  minimal[idx] &= minimality_check(eqs[2], vars[2], [a; b], abvals, tidx; show_progress=show_progress, monodromy=monodromy)
  minimal[idx] &= minimality_check(eqs[3], vars[3], [a; b], abvals, tidx; show_progress=show_progress, monodromy=monodromy)
end

return [is_real(crits[i]) ? real(solution_approximation(crits[i]))[1:n] : solution_approximation(crits[i])[1:n] for i in 1:length(crits) if minimal[i]]
end
