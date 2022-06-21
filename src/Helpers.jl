export leading_asymptotics

# split h into real and imagary parts, h = u + im*v
function split(h, n, vars)
  @polyvar x[1:n] y[1:n] I
  h = h(vars => x + I.*y)
  u = sum(filter(t -> iseven(degree(t, I)), terms(h)))
  v = sum(filter(t -> isodd(degree(t, I)), terms(h)))
  return subs(u, I => im), -im*subs(v, I => im), x, y
end

# determine the critical points
function get_crits(h, n, vars)
  @polyvar λ
  sys = System([h; vars .* differentiate(h, vars) - λ], variables=[vars; λ])
  sols = solutions(solve(sys; show_progress=false))
  if length(sols) > 0
    certs = filter(cert -> is_certified(cert), certificates(certify(sys, sols)))
    return filter(cert -> abs.(certified_solution_interval(cert)[n+1]) < 0 || 0 < abs.(certified_solution_interval(cert)[n+1]), certs)
  else
    return []
  end
end

# check if 0 ∈ abs(a - b) for arbitrary prec complex ball matrix, true if a == b
function compare_acb(a, b)
  return !any((abs.(a - b) .< 0) .| (0 .< abs.(a - b)))
end

function get_certified_sols(sys, show_progress)
  res = solve(sys; show_progress=show_progress)
  if length(real_solutions(res)) > 0
    certs = certificates(certify(sys, real_solutions(res)))
    return filter(cert -> is_certified(cert), certs)
  else
    return []
  end
end

# determine if x ∈ (0, 1)
function zero_to_one(x)
  return (0 < x) && (x < 1)
end

# get a start solution
function get_start_solution(sys, params; n_iters=1000)
  for _ in 1:n_iters
    z₀ = randn(ComplexF64, nvariables(sys))
    res = newton(sys, z₀, params)
    if is_success(res)
      return res.x
    end
  end
end

# TRUE if any solution has t ∈ (0, 1)
function check_solutions(certs, tidx)
  return length(certs) > 0 && any([zero_to_one(abs(certified_solution_interval(cert)[tidx])) for cert in certs])
end

# FALSE if the system has no real solutions with t ∈ (0, 1)
function minimality_check(sys, params, tidx; show_progress=false, monodromy=false)
  if monodromy
    init = get_start_solution(sys, params)
    if !isnothing(init)
      res = monodromy_solve(sys, init, params, show_progress=show_progress)
      if length(real_solutions(results(res))) != 0
        certs = distinct_certificates(certify(sys, real_solutions(results(res)), params))
        certs = filter(cert -> is_certified(cert), certs)
        # if there exists a solution with t ∈ (0, 1), then the point is not minimal
        return !check_solutions(certs, tidx)
      end
    end
  else
    res = solve(sys; target_parameters=params, show_progress=show_progress)
    if length(real_solutions(results(res))) != 0
      certs = distinct_certificates(certify(sys, real_solutions(results(res)), params))
      certs = filter(cert -> is_certified(cert), certs)
      # if there exists a solution with t ∈ (0, 1), then the point is not minimal
      return !check_solutions(certs, tidx)
    end
  end
  # the point still may or may not be minimal (false positives are possible)
  return true
end

# Check assumption A1
function checkA1(h, vars, crits; ϵ=1e-8)
  ∇h = differentiate(h, vars)
  return !any([∇h(vars => crit) < ϵ for crit in crits])
end

function leading_asymptotics(g, h, crits)
  vars = variables(h)
  n = length(vars)

  res = ""
  for ζ in crits
    λ = ζ[1] * subs(differentiate(h, vars[1]), vars => ζ)
    U = [[convert(ComplexF64, ζ[k]*ζ[l]*subs(differentiate(differentiate(h, vars[k]), vars[l]), vars => ζ)) for k in 1:n] for l in 1:n]
    ℋ = [[i == j ? prod(ζ)/convert(ComplexF64, λ*ζ[i]ζ[j])*(U[i][n] + U[j][n] - U[i][j] - U[n][n] - 2*λ) : prod(ζ)/convert(ComplexF64, λ*ζ[i]ζ[j])*(U[i][n] + U[j][n] - U[i][j] - U[n][n] - λ) for i in 1:n] for j in 1:n]
    # map(r -> convert(Vector{ComplexF64}, r), ℋ )
    ℋ = reduce(hcat, ℋ )
    # ℋ = convert(Matrix{ComplexF64}, ℋ )

    coeff = (2π)^((1-n)/2) / sqrt(convert(ComplexF64, prod(ζ.^(3-n)) / ζ[n]^2 * det(ℋ ))) * (-subs(g, vars => ζ) / ζ[n]*subs(differentiate(h, vars[n]), vars => ζ))
    convert(ComplexF64, coeff)
    if res != ""
      res += " + "
    end
    p = '(' * string(prod(ζ)) * ')' * "^(-n)"
    res *= p * ' ' * "n^" * '(' * string((1-n)/2) * ')' * ' ' * '(' * string(coeff) * ')'
  end

  return res
end

