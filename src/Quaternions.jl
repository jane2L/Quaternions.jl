__precompile__()

module Quaternions

  import Base: +, -, *, /, ^
  import Base: abs, abs2, angle, conj, cos, exp, inv, isfinite, log, real, sin, sqrt
  import Base: convert, promote_rule
  import LinearAlgebra: norm, normalize,dot,cross

  include("Quaternion.jl")
  include("Octonion.jl")
  include("DualQuaternion.jl")
  include("ComplexQuaternion.jl")

  export Quaternion
  export quat
  export ComplexQuaternion
  export Cquat
  export snorm
  export snorm2
  export Octonion
  export octo
  export DualQuaternion
  export dualquat
  export angleaxis
  export angle
  export axis
  export linpol
  export normalize
  export normalizea
  export dconj
  export quatrand
  export nquatrand
  export octorand
  export dualquatrand
  export ndualquatrand
  export qrotation
  export rotationmatrix
  export slerp
end
