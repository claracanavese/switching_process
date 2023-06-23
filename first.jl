using Catalyst
using DifferentialEquations
using Plots

# Basic model
rn = @reaction_network begin
    b, 0 --> X
    d, X --> 0
end

u0 = [:X => 1.0]
tspan = (0.0, 10.0)
params = [:b => 1.0, :d => 0.2]

oprob = ODEProblem(rn, u0, tspan, params)

sol = solve(oprob)
plot(sol)

# SIR
sir_model = @reaction_network begin
    b, S + I --> 2I
    k, I --> R
end 

u0 = [:S => 50, :I => 1, :R => 0.0]
tspan = (0.0, 10.0)
params = [:b => 0.2, :k => 1.0]

dprob = DiscreteProblem(sir_model, u0, tspan, params)
jprob = JumpProblem(sir_model, dprob, Direct())

sol = solve(jprob, SSAStepper())
plot(sol)

# Reaction DSL
# syntax
rn = @reaction_network begin
    2.0, X + Y --> XY
    1.0, XY --> Z1 + Z2
end

osys = convert(ODESystem, rn)

# define symbolic variables
@variables t
@species X(t) Y(t) Z(t) XY(t) Z1(t) Z2(t)

# create the mapping
u0 = [X => 1.0, Y => 1.0, XY => 1.0, Z1 => 1.0, Z2 => 1.0]
# alternative
u0 = symmap_to_varmap(rn, [:X => 1.0, :Y => 1.0, :XY => 1.0, :Z1 => 1.0, :Z2 => 1.0])

tspan = (0.0, 1.0)
oprob = ODEProblem(osys, u0, tspan, [])
oprob = ODEProblem(rn, u0, tspan, [])

# defining parameters and species
rn = @reaction_network begin
    k1, X --> Y
    k2, Y --> X
end

# production, destruction, stoichiometry
rn = @reaction_network begin
    2.0, 0 --> X
    1.0, X --> 0
end

rn = @reaction_network begin
    1.0, 2X --> Y
  end

convert(ODESystem, rn)

# arrow variants
rn = @reaction_network begin
    1.0, X + Y --> XY
    1.0, X + Y → XY
    1.0, XY ← X + Y
    1.0, XY <-- X + Y
  end

rn2 = @reaction_network begin
(2.0,2.0), X + Y <--> XY
end

# variable reaction rates
rn = @reaction_network begin
    1.0, X --> ∅
    k*X, Y --> ∅
end

convert(ODESystem, rn)

# my sistem
rn = @reaction_network begin
    @species begin
        Z1(t)=1.0
        Z2(t)=0.0
    end
    @parameters begin
        a1 = 15.0
        b1 = 0.0
        a2 = 10.0
        b2 = 0.0
        s1 = 0.1
        s2 = 0.01
    end
    a1, Z1 --> 2Z1
    b1, Z1 --> 0
    a2, Z2 --> 2Z2
    b2, Z2 --> 0
    s1, Z1 --> Z1 + Z2
    s2, Z2 --> Z1 + Z2  
end
u0 = [:Z1 => 1.0, :Z2 => 0.0]
tspan = (0.0, 1.0)
oprob = ODEProblem(rn, u0, tspan, [])
sol = solve(oprob)
plot(sol)