module MockBlockSolver
  import Base.convert, Base.show, Base.copy, Base.pointer

  using JuMP, BlockDecomposition

  importall MathProgBase.SolverInterface
  importall BlockDecomposition.BlockSolverInterface

  include("MockBlockSolverInterface.jl")

end # module
