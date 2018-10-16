__precompile__()

module Quaternions

  import Base: +, -, *, /, ^
  import Base: abs, abs2, angle, conj, cos, exp, inv, isfinite, log, real,imag, sin, sqrt, atan, tan
  import Base: convert, promote_rule
  import LinearAlgebra: norm, normalize,dot,cross

  include("Quaternion.jl")
  include("Octonion.jl")
  include("DualQuaternion.jl")
  include("ComplexQuaternion.jl")

  export Quaternion
  export quat

  # ** ComplexQuaternion **
  export ComplexQuaternion
  export Cquat
  export snorm
  export snorm2
  export conj
  export complexconj
  export reversioninvolution
  export inv
  export scalar
  export vector
  export complexForm
  export complexFormAxis
  export complexPolarForm
  export hamiltonPolarForm
  export hamiltonPolarAmplitude
  export snormalize
  export arg
  export I



# ** Octonion **
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
