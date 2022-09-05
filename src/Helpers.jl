export leading_asymptotics

# split h into real and imagary parts, h = u + im*v
function split(h::DynamicPolynomials.Polynomial{true, <:Number})
  z = variables(h)
  n = length(z)

  x, y, I = @polyvar x[1:n] y[1:n] I
  h = h(z => x + I.*y)
  u = sum(filter(t -> iseven(degree(t, I)), terms(h)))
  v = sum(filter(t -> isodd(degree(t, I)), terms(h)))
  return subs(u, I => im), -im*subs(v, I => im), x, y
end

# determine the critical points of h
function get_crits(h::DynamicPolynomials.Polynomial{true, <:Number}, r::Vector{<:Number})
  z = variables(h)
  n = length(z)

  # solve critical point system
  λ, = @polyvar λ
  sys = System([h; z .* differentiate(h, z) .- r*λ], variables=[z; λ])
  sols = solutions(solve(sys; show_progress=false))

  # certify found solutions
  if length(sols) > 0
    certs = filter(cert -> is_certified(cert), certificates(certify(sys, sols)))
    # remove solutions with λ = 0
    return filter(
      cert -> abs.(certified_solution_interval_after_krawczyk(cert)[n+1]) < 0 || 0 < abs.(certified_solution_interval_after_krawczyk(cert)[n+1]),
      certs
    )
  else
    return []
  end
end

# check if 0 ∈ abs(a - b) for arbitrary prec complex ball vectors. (true if a == b is possible)
function compare_acb(a::AbstractVector{<:Number}, b::AbstractVector{<:Number})
  return !any((abs.(a - b) .< 0) .| (0 .< abs.(a - b)))
end

# get certified solutions to sys with no parameters
function get_certified_sols(sys::System; show_progress=false)
  res = solve(sys; show_progress=show_progress)
  if length(real_solutions(res)) > 0
    certs = certificates(certify(sys, real_solutions(res)))
    return filter(cert -> is_certified(cert), certs)
  else
    return []
  end
end

# determine if x ∈ (0, 1)
function zero_to_one(x::Real)
  return (0 < x) && (x < 1)
end

# get a start solution to sys
function get_start_solution(sys::System, params::Vector{<:Number}; n_iters=1000)
  for _ in 1:n_iters
    z₀ = randn(ComplexF64, nvariables(sys))
    res = newton(sys, z₀, params)
    if is_success(res)
      return res.x
    end
  end
end

function krawczyk(
    F::Vector{<:DynamicPolynomials.Polynomial{true, <:Number}},
    z::Vector{DynamicPolynomials.PolyVar{true}},
    I₀::Vector{Arblib.Acb},
    JF::Arblib.AcbMatrix,
    Y::Arblib.AcbMatrix
  )

  # get mid point of I₀
  x₀ = AcbMatrix(I₀)
  Arblib.get_mid!(x₀, AcbMatrix(I₀))

  # for substitution into F
  x₀ = convert(Array{Acb}, x₀)
  x₀ = vcat(x₀...)
  Fx = subs(F, z => x₀)
  Fx = Arblib.AcbMatrix(convert(Array{Acb, 1}, Fx))

  return vcat(x₀ - Y*Fx + (I - Y*JF)*(I₀ - x₀)...)
end

# Given a system of polynomials `F` in the variables `z` and a certified solution interval I₀ as returned by `certified_solution_interval_after_krawczyk`,
# this function returns a smaller certified solution interval J₀ so that 1 is not contained in the t coordinate of J₀.
function refine_t(
    F::Vector{<:DynamicPolynomials.Polynomial{true, <:Number}},
    z::Vector{DynamicPolynomials.PolyVar{true}},
    I₀::Arblib.AcbMatrix,
    tidx::Integer;
    max_iters=100
  )

  # convert matrix to vector for substitutions to make sense
  I₀ = vcat(I₀...)
  n = length(I₀)

  # Jacobian of system
  J = differentiate(F, z)
  J = subs(J, z => I₀)

  # Convert to AcbMatrix - Arblib internal type & take mid points
  J = AcbMatrix(convert(Array{Acb, 2}, J))
  Arblib.get_mid!(J, J)

  # Y matrix for Krawczyk
  Y = inv(J)
  Y = AcbMatrix(Y)

  for iter in 1:max_iters
    if !Arblib.contains_zero(I₀[tidx] - 1)
      return I₀ # interval does not contain 1
    end

    nxt = krawczyk(F, z, I₀, J, Y)
    if !all(Arblib.contains(I₀[j], nxt[j]) for j in 1:n)
      # The interval is not contracting
      println("WARN: Intervals not contracting, false positive results possible.")
      return I₀
    end
    I₀ = nxt
  end

  # max iterations reached
  println("WARN: Max iterations reached, false positive results possible.")
  return I₀
end

function refine(
    F::Vector{<:DynamicPolynomials.Polynomial{true, <:Number}},
    z::Vector{DynamicPolynomials.PolyVar{true}},
    I₀::Arblib.AcbMatrix,
    width::Real,
    max_iters=100
  )

  # convert matrix to vector for substitutions to make sense
  I₀ = vcat(I₀...)
  n = length(I₀)

  # Jacobian of system
  J = differentiate(F, z)
  J = subs(J, z => I₀)

  # Convert to AcbMatrix - Arblib internal type & take mid points
  J = AcbMatrix(convert(Array{Acb, 2}, J))
  Arblib.get_mid!(J, J)

  # Y matrix for Krawczyk
  Y = inv(J)
  Y = AcbMatrix(Y)

  for iter in 1:max_iters
    nxt = krawczyk(F, z, I₀, J, Y)
    if !all(Arblib.contains(I₀[j], nxt[j]) for j in 1:n)
      # The interval is not contracting
      println("WARN: Intervals not contracting, false positive results possible.")
      return I₀
    end

    I₀ = nxt
    rerad = maximum(radius.(real(I₀)))
    imrad = maximum(radius.(imag(I₀)))
    if max(rerad, imrad) < width
      # interval is sufficiently tight
      return I₀
    end
  end

  # max iterations reached
  println("WARN: Max iterations reached, false positive results possible.")
  return I₀
end

"""
  leading_asymptotics(g, h, crits)
Given the minimal critical points, leading_asymptotics prints the leading term of the asymptotic formula approximating the coefficients of the power series expansion of g/h around 0.
Optionally, r may be given to specify a diagonal other than the main diagonal.

Example:
julia> using ACSVHomotopy

julia> @polyvar x y
(x, y)

julia> h = 1-x-y
-x - y + 1

julia> minimal = find_min_crits_comb(h)
1-element Array{Array{Float64,1},1}:
 [0.5, 0.5]

julia> leading_asymptotics(1, h, minimal)
"(0.25)^(-n) n^(-0.5) (0.5641895835477563 + 0.0im)"

"""
function leading_asymptotics(g, h::DynamicPolynomials.Polynomial, crits::Vector{<:Vector{<:Number}}; r=Nothing)

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

  res = ""
  for ζ in crits
    λ = convert(ComplexF64, ζ[1] * subs(differentiate(h, z[1]), z => ζ))
    U = [[convert(ComplexF64, ζ[k]*ζ[l]*subs(differentiate(differentiate(h, z[k]), z[l]), z => ζ))/λ for k in 1:n] for l in 1:n]
    ℋ = [[i == j ? convert(ComplexF64, 2 + U[i][j] - U[i][n] - U[j][n] + U[n][n]) : convert(ComplexF64, 1 + U[i][j] - U[i][n] - U[j][n] + U[n][n]) for i in 1:n-1] for j in 1:n-1]
    ℋ = reduce(hcat, ℋ )

    coeff = (2π)^((1-n)/2) / sqrt(convert(ComplexF64, det(r[n]*ℋ ))) * (-sign(r[n])*subs(g, z => ζ) / convert(ComplexF64, ζ[n]*subs(differentiate(h, z[n]), z => ζ)))
    convert(ComplexF64, coeff)

    if res != ""
      res += " + "
    end
    p = '(' * string(prod(ζ)) * ')' * "^(-nr)"
    res *= p * ' ' * "n^" * '(' * string((1-n)/2) * ')' * ' ' * '(' * string(coeff) * ')'
  end

  return res
end
