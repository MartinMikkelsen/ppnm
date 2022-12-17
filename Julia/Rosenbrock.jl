# Import the package
import Pkg; Pkg.add("OptimizationNLopt")
import Pkg; Pkg.add("Optim")

using Optimization, OptimizationNLopt, Optim

# Define the problem to optimize
L(u,p) =  (p[1] - u[1])^2 + p[2] * (u[2] - u[1]^2)^2
u0 = zeros(2)
p  = [1.0,100.0]
prob = OptimizationProblem(L, u0, p, lb = [-1.0,-1.0], ub = [1.0,1.0])

# Solve the optimization problem
sol = solve(prob,NLopt.LN_NELDERMEAD())

# Analyze the solution
@show sol.u, L(sol.u,p)
