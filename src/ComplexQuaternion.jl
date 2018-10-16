#using ComplexNumbers

struct ComplexQuaternion{T<:Real} <: Number
  qr::Quaternion{T}
  qi::Quaternion{T}
  norm::Bool
end



ComplexQuaternion(qr::Quaternion, qi::Quaternion, n::Bool = false) =
  ComplexQuaternion(promote(qr, qi)..., n)

ComplexQuaternion(d1::Complex, d2::Complex, d3::Complex, d4::Complex, n::Bool = false) =
  ComplexQuaternion(Quaternion(d1.re, d2.re, d3.re, d4.re, n),
                  Quaternion(d1.im, d2.im, d3.im, d4.im), n)

ComplexQuaternion(x::Real) = ComplexQuaternion(Quaternion(x), Quaternion(zero(x)), abs(x) == one(x))

ComplexQuaternion(c::Complex) = ComplexQuaternion(c, zero(c), zero(c), zero(c), abs(c) == one(c.re))

ComplexQuaternion(q::Quaternion) = ComplexQuaternion(q, zero(q), q.norm)

ComplexQuaternion(a::Vector) = ComplexQuaternion(zero(Quaternion{typeof(a[1])}), Quaternion(a))

# Creating the pseudo-scalar I
const I = ComplexQuaternion(Quaternion(0,0,0,0),Quaternion(1,0,0,0))

convert(::Type{ComplexQuaternion{T}}, x::Real) where {T} =
    ComplexQuaternion(convert(Quaternion{T}, x), convert(Quaternion{T}, 0))

convert(::Type{ComplexQuaternion{T}}, c::Complex) where {T} =
    ComplexQuaternion(convert(Complex{T}, c), convert(Complex{T}, 0), convert(Complex{T}, 0), convert(Complex{T}, 0))

convert(::Type{ComplexQuaternion{T}}, q::Quaternion) where {T} =
    ComplexQuaternion(convert(Quaternion{T}, q), convert(Quaternion{T}, 0), q.norm)

convert(::Type{ComplexQuaternion{T}}, q::ComplexQuaternion{T}) where {T <: Real} = q

convert(::Type{ComplexQuaternion{T}}, cq::ComplexQuaternion) where {T} =
    ComplexQuaternion(convert(Quaternion{T}, cq.qr), convert(Quaternion{T}, cq.qi), cq.norm)

promote_rule(::Type{ComplexQuaternion{T}}, ::Type{T}) where {T <: Real} = ComplexQuaternion{T}
promote_rule(::Type{ComplexQuaternion}, ::Type{T}) where {T <: Real} = ComplexQuaternion
promote_rule(::Type{ComplexQuaternion{T}}, ::Type{S}) where {T <: Real, S <: Real} = ComplexQuaternion{promote_type(T, S)}
promote_rule(::Type{Quaternion{T}}, ::Type{ComplexQuaternion{S}}) where {T <: Real, S <: Real} = ComplexQuaternion{promote_type(T, S)}
promote_rule(::Type{ComplexQuaternion{T}}, ::Type{ComplexQuaternion{S}}) where {T <: Real, S <: Real} = ComplexQuaternion{promote_type(T, S)}

Cquat(q1, q2, n=false) = ComplexQuaternion(q1, q2, n)
Cquat(d1, d2, d3, d4, n=false) = ComplexQuaternion(d1, d2, d3, d4, n)
Cquat(x) = ComplexQuaternion(x)

function Base.show(io::IO, cq::ComplexQuaternion)
  print(io,cq.qr.s ," + ", cq.qr.v1, "i + ", cq.qr.v2 ,"j + " ,cq.qr.v3,"k + I(",cq.qi.s, " + ", cq.qi.v1,"i + " ,cq.qi.v2,"j + ", cq.qi.v3,"k ) \n" )
end

qr(cq::ComplexQuaternion) = cq.qr
qi(cq::ComplexQuaternion) = cq.qi
c1(cq::ComplexQuaternion) = Complex(cq.qr.s,cq.qi.s)
c2(cq::ComplexQuaternion) = Complex(cq.qr.v1,cq.qi.v1)
c3(cq::ComplexQuaternion) = Complex(cq.qr.v2,cq.qi.v2)
c4(cq::ComplexQuaternion) = Complex(cq.qr.v3,cq.qi.v3)

(/)(cq::ComplexQuaternion, x::Real) = ComplexQuaternion(cq.qr/x, cq.qi/x)

(/)(cq::ComplexQuaternion, d::Complex) =
  ComplexQuaternion(c1(cq) / d,
                  c2(cq) / d,
                  c3(cq) / d,
                  c3(cq) / d)

(/)(cq::ComplexQuaternion,q::Quaternion) = ComplexQuaternion(cq.qr/q,cq.qi/q)

(*)(cq::ComplexQuaternion, cr::ComplexQuaternion) = ComplexQuaternion(cq.qr*cr.qr-cq.qi*cr.qi,cq.qr*cr.qi+ cq.qi*cr.qr)
(*)(cq::ComplexQuaternion,c::Complex) = ComplexQuaternion(c1(cq)*c,c2(cq)*c,c3(cq)*c,c4(cq)*c)
(+)(cq::ComplexQuaternion,c::Complex) = ComplexQuaternion(c1(cq)+c,c2(cq)+c,c3(cq)+c,c4(cq)+c)
(-)(cq::ComplexQuaternion,c::Complex) = ComplexQuaternion(c1(cq)-c,c2(cq)-c,c3(cq)-c,c4(cq)-c)


abs2(cq::ComplexQuaternion) = cq.norm ? one(cq.qr.s) :
  abs2(cq.qr)+abs2(cq.qi)

abs(cq::ComplexQuaternion) = cq.norm ? one(cq.qr.s) : sqrt(abs2(cq))

conj(cq::ComplexQuaternion) = ComplexQuaternion(conj(cq.qr), conj(cq.qi), cq.norm)
complexconj(cq::ComplexQuaternion) = ComplexQuaternion(cq.qr, -cq.qi, cq.norm)
reversioninvolution(cq::ComplexQuaternion) = ComplexQuaternion(cq.c1,-cq.c2,-cq.c3,cq.c4) # TO CHECK # TO RENAME

function snorm2(cq::ComplexQuaternion) # Returns a complex
    if cq.norm
        return Complex(one(cq.qr.s))
    else
        sn = cq*conj(cq)
        return (Complex(sn.qr.s,sn.qi.s))
    end
end

snorm(cq::ComplexQuaternion) = cq.norm ? Complex(one(cq.qr.s)) : Complex(sqrt(snorm2(cq)))

inv(cq::ComplexQuaternion) = (snorm2(cq) != 0) ? conj(cq) / snorm2(cq) : error("Not inversible")

scalar(cq::ComplexQuaternion) = (ComplexQuaternion(c1(cq),Complex(0),Complex(0),Complex(0)))
vector(cq::ComplexQuaternion) = (ComplexQuaternion(Complex(0),c1(cq),c2(cq),c3(cq)))


function complexForm(cq::ComplexQuaternion)
    a = scalar(cq)
    V = vector(cq)
    b = snorm(V)
    xi = V/b
    return (a,xi,b)
end

function complexFormAxis(cq::ComplexQuaternion)
    V = vector(cq)
    b = snorm(V)
    xi = V/b
    return (xi)
end

function hamiltonPolarForm(cq::ComplexQuaternion)
    a,xi,b = complexForm(cq)
    a = Complex(c1(a))
    if abs(a) > 0
        theta = atan(b/a)
    else
        theta = pi/2  # TODO chek the notation pi (works in the terminal)
    end
    R = snorm(cq)
    return (R,xi,theta)
end

function hamiltonPolarAmplitude(cq::ComplexQuaternion)
    R,xi,theta = hamiltonPolarForm(cq)
    cq = R
end



function complexPolarForm(cq::ComplexQuaternion)
    if abs2(cq.qr) > 0
        angle = 1/cq.qr*cq.qi
        psi = atan(Quaternion(angle))
        Q = (cq.qr +cq.qi)/(cos(psi)+sin(psi))
        return (Q,psi)
    else
        cq = cq#TODO
        return (0,0)
    end

end




function normalize(cq::ComplexQuaternion)
  if (cq.norm)
    return cq
  end
  a = abs(cq)
  if abs(a) > 0
    qa = cq / a
    return ComplexQuaternion(qa.qr, qa.qi, true)
  else
    return cq
  end
end

function snormalize(cq::ComplexQuaternion)
    a = snorm(cq)
    if abs2(a) > 0
      qa = cq / a
      return (ComplexQuaternion(qa.qr, qa.qi))
    else
      cq
    end
end

function normalizea(cq::ComplexQuaternion)
  if (cq.norm)
    return (cq, one(Complex))
  end
  a = abs(cq)
  if abs(a) > 0
    qa = cq / a
    return (ComplexQuaternion(qa.qr, qa.qi, true), a)
  else
    return (cq, zero(Complex))
  end
end

function snormalizea(cq::ComplexQuaternion)
    a = snorm(cq)
    if abs2(a) > 0
        qa = cq/a
        return (ComplexQuaternion(qa.qr,qa.qi),a)
    else
        return (cq,zero(Complex))
    end
end



(-)(cq::ComplexQuaternion) = ComplexQuaternion(-cq.qr, -cq.qi, cq.norm)

(+)(cq::ComplexQuaternion, cw::ComplexQuaternion) = ComplexQuaternion(cq.qr + cw.qr, cq.qi + cw.qi)
(-)(cq::ComplexQuaternion, cw::ComplexQuaternion) = ComplexQuaternion(cq.qr - cw.qr, cq.qi - cw.qi)
(*)(cq::ComplexQuaternion, cw::ComplexQuaternion) = ComplexQuaternion(cq.qr * cw.qr - cq.qi*cw.qi,
                                                              cq.qi * cw.qr + cq.qr * cw.qi,
                                                              cq.norm && cw.norm)
(/)(cq::ComplexQuaternion, cw::ComplexQuaternion) = cq * inv(cw)


function arg(cq::ComplexQuaternion) # argument according to Chappell definition
    Z = scalar(cq)
    F = vector(cq)
    phi = atan(Complex(snorm(F))/Complex(c1(Z)))
    return phi
end

# function angleaxis(cq::ComplexQuaternion)
#   tq = cq.qi * conj(cq.qr)
#   t = [2.0 * tq.v1, 2.0 * tq.v2, 2.0 * tq.v3]
#   qrs = cq.qr.s
#   th0, s0 = angleaxis(cq.qr)
#   sqr = quat(0.0, s0)
#   if abs(abs(qrs) - one(qrs)) == 0
#     th = Complex(th0, 0.5 * abs(quat(0, t)))
#     th, ComplexQuaternion(sqr)
#   else
#     th = Complex(th0, 0.5 * dot(t, s0))
#     s0c1 = cross(s0, t)
#     tanth = tan(th0)
#     s0c2 = (s0c1 / tanth + t) * 0.5
#     sqiv = cross(s0c2, s0)
#     th, ComplexQuaternion(sqr, quat(0.0, sqiv))
#   end
# end
#
# function angle(cq::ComplexQuaternion)
#   th, ax = angleaxis(cq)
#   th
# end
#
# function axis(cq::ComplexQuaternion)
#   th, ax = angleaxis(cq)
#   ax
# end

# Implementation of exp/log and trigonometric functions
# Based on the work in: Chappell 2014, Functions of multivector variables

function exp(cq::ComplexQuaternion)
  se = scalar(cq)
  se = exp(se)
  cq = vector(cq)
  sq = snorm(cq)
  nq = snormalize(cq)
  if abs2(sq) > 0
    ComplexQuaternion(se) * (ComplexQuaternion(cos(sq)) + nq * ComplexQuaternion(sin(sq)))
  else
    ComplexQuaternion(se)*(1+cq)
  end
end


function log(cq::ComplexQuaternion)
  Z = scalar(cq)
  F = vector(cq)
  N = snorm(cq)
  lgN = log(N)
  phi = arg(cq)
  Quaternion(N,0,0,0) + phi*snormalize(F)
end

(^)(cq::ComplexQuaternion, cw::ComplexQuaternion) = exp(cw * log(cq))

function sqrt(cq::ComplexQuaternion)
  exp(0.5 * log(cq))
end

# Cosine , sine functions ill-defined




Complexquatrand() = ComplexQuaternion(quatrand(), quatrand())
nComplexquatrand() = normalize(Complexquatrand())
