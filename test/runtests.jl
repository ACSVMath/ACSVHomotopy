using Test, ACSVHomotopy

function test_comb_case()

  @polyvar x y z w

  # test monodromy method
  # 2 variables
  h = 1 - x - y
  println("Running CombCase h = ", h)
  crits = find_min_crits_comb(h)
  @assert length(crits) == 1
  @assert all(abs.(crits[1] - [1/2, 1/2]) .< 1e-5)

  h = (1-x-y)*(20-x-40*y)-1
  println("Running CombCase h = ", h)
  crits = find_min_crits_comb(h)
  @assert length(crits) == 1
  @assert all(abs.(crits[1] - [0.549, 0.309]) .< 0.001)

  # Apéry's sequence
  h = 1 - w*(1 + x)*(1 + y)*(1 + z)*(x*y*z + y*z + y + z + 1)
  println("Running CombCase h = ", h)
  crits = find_min_crits_comb(h)
  @assert length(crits) == 1
  @assert all(abs.(crits[1] - [1+sqrt(2), sqrt(2)/2, sqrt(2)/2, 58*sqrt(2)-82]) .< 1e-5)

  println("comb case tests completed successfully")
  return true
end

function test_non_comb_case()

  @polyvar x y z w

  # test monodromy
  h = 1 - x - y
  println("Running Approx Crit + Monodromy: h = ", h)
  crits = find_min_crits(h; approx_crit=true, monodromy=true)

  h = (1-x-y)*(20-x-40*y)-1
  println("Running Approx Crit + Monodromy: h = ", h)
  crits = find_min_crits(h; approx_crit=true, monodromy=true)

  # no critical points
  h = 2 + y - x*(1+y)^2
  println("Running Approx Crit + Monodromy: h = ", h)
  crits = find_min_crits(h; approx_crit=true, monodromy=true)

  # test approx crit + solve
  h = (1-x-y)*(20-x-40*y)-1
  println("Running Approx Crit + Solve: h = ", h)
  crits = find_min_crits(h; approx_crit=true)

  # no critical points
  h = 2 + y - x*(1+y)^2
  println("Running Approx Crit + Solve: h = ", h)
  crits = find_min_crits(h; approx_crit=true)

  # 3 variable examples
  h = 1 - (x+y+z) + (81/8)*x*y*z
  println("Running Approx Crit + Solve: h = ", h)
  crits = find_min_crits(h; approx_crit=true, monodromy=false)

  h = 1-(x^2+x*y+y^2+z^2)
  println("Running Full Solve: h = ", h)
  crits = find_min_crits(h; approx_crit=true)

  # Apéry's sequence (It may take a while)
  # h = 1 - w*(1 + x)*(1 + y)*(1 + z)*(x*y*z + y*z + y + z + 1)
  # println("Running Approx Crit + Solve: h = ", h)
  # crits = find_min_crits(h; approx_crit=true, show_progress=true)

  # ----- test non-comb method ----- #
  h = 1 - x - y
  println("Running Full Solve: h = ", h)
  crits = find_min_crits(h)

  h = (1-x-y)*(20-x-40*y)-1
  println("Running Full Solve: h = ", h)
  crits = find_min_crits(h)

  # 3 variables
  h = 1 - (x+y+z) + (81/8)*x*y*z
  println("Running Full Solve: h = ", h)
  crits = find_min_crits(h)

  println("non-comb case tests completed successfully")
  return true
end

@test test_comb_case()
@test test_non_comb_case()

