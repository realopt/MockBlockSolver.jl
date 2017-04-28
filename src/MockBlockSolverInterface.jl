type MockSolver <: AbstractMathProgSolver
  vars_decomposition
  cstrs_decomposition
  sp_mult
  sp_prio
  obj_magnitude
  obj_lb
  obj_ub
  vars_branching_priorities
  oracles
  genvarcallbacks
  gencstrcallbacks
end
export MockSolver

function MockSolver()
  MockSolver(nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing)
end

type MockMathProgModel <: AbstractMathProgModel
  solver::MockSolver
end

LinearQuadraticModel(s::MockSolver) = MockMathProgModel(s)
ConicModel(s::MockSolver) = LPQPtoConicBridge(LinearQuadraticModel(s))
supportedcones(s::MockSolver) = [:Free,:Zero,:NonNeg,:NonPos]

################################################################################
######################   BlockDecomposition.BlockSolverInterface ###############
################################################################################
set_constrs_decomposition!(s::MockSolver, cstrs_decomposition_list) =
  (s.cstrs_decomposition = cstrs_decomposition_list)

set_vars_decomposition!(s::MockSolver, vars_decomposition_list) =
  (s.vars_decomposition = vars_decomposition_list)

set_sp_mult!(s::MockSolver, sp_multiplities) =
  (s.sp_mult = sp_multiplities)

set_var_branching_prio!(s::MockSolver, var_priorities) =
  (s.vars_branching_priorities = var_priorities)

set_sp_prio!(s::MockSolver, sp_priorities) =
  (s.sp_prio = sp_priorities)

set_var_generic!(s::MockSolver, genvarsfcts) =
  (s.genvarcallbacks = genvarsfcts)

set_cstr_generic!(s::MockSolver, gencstrcbs) =
  (s.gencstrcallbacks = gencstrcbs)

function set_objective_bounds_and_magnitude!(s::MockSolver, magnitude, lb, ub)
  s.obj_magnitude = (lb > magnitude || ub < magnitude)? NaN : magnitude
  s.obj_lb = (lb > -Inf)? lb : NaN
  s.obj_ub = (ub < Inf)? ub : NaN
end

getcurrentcost(m::MockMathProgModel, varidx) = error("GetCurrentCost not implemented")

set_oracles!(s::MockSolver, oracles) = (s.oracles = oracles)

################################################################################
##################### MathProgBase.SolverInterface #############################
## Loads problem data from the given file
function loadproblem!(m::MockMathProgModel, filename::AbstractString)
  println("loadproblem!(error) not implemented")
end

## Loads the provided problem data to set up the linear programming problem
##  loadproblem!(m::MockMathProgModel, A, l, u, c, lb, ub, sense)
##  which is the following LP
##        minₓ cᵀx
##        s.c. lb ≤ Ax ≤ ub
##              l ≤  x ≤ u
## sense specifies the direction of optimization problem (:Min or :Max)
function loadproblem!(m::MockMathProgModel, A, l, u, c, lb, ub, sense)
  println("Fake optimization.")
end


## Returns a vector containing the lower bounds l on the variables.
getvarLB(m::MockMathProgModel) = println("getvarLB! not implemented")

## Sets the lower bounds on the variables.
setvarLB!(m::MockMathProgModel, l) = println("setvarLB! not implemented")

##  Returns a vector containing the upper bounds u on the variables.
getvarUB(m::MockMathProgModel) = println("getvarUB not implemented")

## Sets the upper bounds on the variables.
setvarUB!(m::MockMathProgModel, u) = println("setvarUB! not implemented")

## Returns a vector containing the lower bounds lb on the linear constraints.
getconstrLB(m::MockMathProgModel) = println("getconstrLB not implemented")

## Sets the lower bounds on the linear constraints.
setconstrLB!(m::MockMathProgModel, lb) = println("setconstrLB! not implemented")

## Returns a vector containing the upper bounds ub on the linear constraints.
getconstrUB(m::MockMathProgModel) = println("getconstrUB not implemented")

## Sets the upper bounds on the linear constraints.
setconstrUB!(m::MockMathProgModel, ub) = println("setconstrUB! not implemented")

## Returns a vector containing the linear objective coefficients c.
getobj(m::MockMathProgModel) = println("getobj not implemented")

## Sets the linear objective coefficients.
setobj!(m::MockMathProgModel, c) = println("setobj! not implemented")

## Returns the full linear constraint matrix A, typically as a SparseMatrixCSC.
getconstrmatrix(m::MockMathProgModel) = error("getconstrmatrix! not implemented")

## Adds a new variable to the model, with lower bound l (-Inf if none), upper
## bound u (Inf if none), and objective coefficient objcoef. Constraint
## coefficients for this new variable are specified in a sparse format: the
## constrcoef vector contains the nonzero coefficients, and the constridx vector
## contains the indices of the corresponding linear constraints.
function addvar!(m::MockMathProgModel, constridx, constrcoef, l, u, objcoef)
  println("addvar! not implemented")
end

## Adds a new variable to the model, with lower bound l (-Inf if none),
## upper bound u (Inf if none), and objective coefficient objcoef. This is
## equivalent to calling the above method with empty arrays for the constraint
## coefficients.
function addvar!(m::MockMathProgModel, l, u, objcoef)
  println("addvar! not implemented")
end

## Adds a new linear constraint to the model, with lower bound lb (-Inf if none)
## and upper bound ub (Inf if none). Coefficients for this new constraint are
## specified in a sparse format: the coef vector contains the nonzero
## coefficients, and the varidx vector contains the indices of the corresponding
## variables.
# function addconstr!(m::MockMathProgModel, varidx, coef, lb, ub)
#   println("addconstr! not implemented")
# end
function addconstr!(m::MockMathProgModel, varidx, coef, lb, ub)
  println("addconstr! not implemented")
end

## Returns the number of linear constraints in the model.
numlinconstr(m::MockMathProgModel) = println("numlinconstr not implemented")

## Returns a vector containing the values of the linear constraints at the
## solution. This is the vector Ax.
getconstrsolution(m::MockMathProgModel) =
  println("getconstrsolution not implemented")

## Returns the dual solution vector corresponding to the variable bounds,
## known as the reduced costs. Not available when integer variables are present.
getreducedcosts(m::MockMathProgModel) =
  println("getreducedcosts not implemented")


######################## SOLVE INTERFACE ######################################

## Returns the solution vector found by the solver.
getsolution(m::MockMathProgModel) = println("getsolution not implemented")

## Returns the objective value of the solution found by the solver. In
## particular, this may be the objective value for the best feasible solution if
## optimality is not attained or proven.
getobjval(m::MockMathProgModel) = println("getobjval not implemented")

## Solves the optimization problem.
function optimize!(m::MockMathProgModel)
  println("Fake optimize!")
end

## Returns the termination status after solving. Possible values include
## :Optimal, :Infeasible, :Unbounded, :UserLimit (iteration limit or timeout),
## and :error. Solvers may return other statuses, for example, when presolve
## indicates that the model is either infeasible or unbounded, but did not
## determine which.
function status(m::MockMathProgModel)
  println("Fake status")
  return :Optimal
end

## Returns the best known bound on the optimal objective value. This is used,
## for example, when a branch-and-bound method is stopped before finishing.
getobjbound(m::AbstractMathProgModel) = println("getobjbound not implemented")

## Returns the final relative optimality gap as optimization terminated.
getobjgap(m::AbstractMathProgModel) = println("getobjgap not implemented")

## Returns an object that may be used to access a solver-specific API for this
## model.
getrawsolver(m::AbstractMathProgModel) = println("getrawsolver not implemented")

##  Returns the total elapsed solution time as reported by the solver.
getsolvetime(m::AbstractMathProgModel) = println("getsolvetime not implemented")

## Sets the optimization sense of the model. Accepted values are :Min and :Max.
setsense!(m::MockMathProgModel, sense) = println("setsense! not implemented")

## Returns the optimization sense of the model.
getsense(m::MockMathProgModel) = println("getsense not implemented")

## Returns the number of variables in the model.
numvar(m::MockMathProgModel) = println("numvar not implemented")

## Returns the total number of constraints in the model.
numconstr(m::MockMathProgModel) = println("numconstr not implemented")

## Release any resources and memory used by the model. Note that the Julia
## garbage collector takes care of this automatically, but automatic collection
## cannot always be forced. This method is useful for more precise control of
## resources, especially in the case of commercial solvers with licensing
## restrictions on the number of concurrent runs. Users must discard the model
## object after this method is invoked.
freemodel!(m::AbstractMathProgModel) = println("freemodel! not implemented")

## Sets the types of the variables to those indicated by the vector v. Valid
## types are :Int for integer, :Cont for continuous, :Bin for binary, :SemiCont
## for semicontinuous, and :SemiInt for semi-integer.
setvartype!(m::MockMathProgModel, typ::Vector{Symbol}) =
    println("setvartype! not implemented")

## Returns a vector indicating the types of each variable, with values described
## above.
getvartype(m::MockMathProgModel) = println("getvartype not implemented")

## It is the philosophy of MathProgBase to not abstract over most solver
## parameters, because very few are universal, and we do not want to make it
## difficult for users who already familiar with the parameters of their solver
## of choice. However, in certain situations an abstraction over parameters is
## needed, for example, in meta-solvers which must adjust parameters of a
## subsolver inside a loop. Parameters set using these methods should override
## any parameters provided in a solver-dependent way. Solvers/models may chose
## to implement the following method:

## Sets solver-independent parameters via keyword arguments. Curent valid
## parameters are TimeLimit and Silent. TimeLimit (Float64): If the solve is not
## completed to optimality tolerances within the given number of seconds, the
## solver should return immediately with status :UserLimit. Silent (Bool): If
## set to true, the solver should disable any output to the console. If set to
## false, the parameter has no effect.
function setparameters!(m::Union{AbstractMathProgSolver, AbstractMathProgModel}; kwargs...)
    println("setparameters! not yet implemented")
end

## If these parameter-setting methods are called on an AbstractMathProgSolver,
## then they should apply to all new models created from the solver
## (but not existing models). If they are called on an AbstractMathProgModel,
## they should apply to that model only. Unrecognized parameters
## (those not listed above) should be ignored or trigger a warning message.
