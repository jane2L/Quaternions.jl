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

function show(io::IO, cq::ComplexQuaternion)
  show(io, cq.qr)
  print(io, " + I ")
  show(io, cq.qi)
  print(io)
end

qr(cq::ComplexQuaternion) = cq.qr
qi(cq::ComplexQuaternion) = cq.qi
c1(cq::ComplexQuaternion) = Complex(cq.qr.s,cq.qi.s)
c2(cq::ComplexQuaternion) = Complex(cq.qr.v1,cq.qi.v1)
c3(cq::ComplexQuaternion) = Complex(cq.qr.v2,cq.qi.v2)
c4(cq::ComplexQuaternion) = Complex(cq.qr.v3,cq.qi.v3)

(/)(cq::ComplexQuaternion, x::Real) = ComplexQuaternion(cq.qr / x, cq.qi / x)

(/)(cq::ComplexQuaternion, d::Complex) =
  ComplexQuaternion(Complex(cq.qr.s, cq.qi.s) / d,
                  Complex(cq.qr.v1, cq.qi.v1) / d,
                  Complex(cq.qr.v2, cq.qi.v2) / d,
                  Complex(cq.qr.v3, cq.qi.v3) / d)

abs2(cq::ComplexQuaternion) = cq.norm ? one(cq.qr.s) :
  abs2(cq.qr)+abs2(cq.qi)

abs(cq::ComplexQuaternion) = cq.norm ? one(cq.qr.s) : sqrt(abs2(cq))

conj(cq::ComplexQuaternion) = ComplexQuaternion(conj(cq.qr), conj(cq.qi), cq.norm)
complexconj(cq::ComplexQuaternion) = ComplexQuaternion(cq.qr, -cq.qi, cq.norm)
reversioninvolution(cq::ComplexQuaternion) = ComplexQuaternion(cq.c1,-cq.c2,-cq.c3,cq.c4) # TO CHECK # TO RENAME

snorm2(cq::ComplexQuaternion) = cq.norm ? Complex(one(cq.qr.s)) : cq*conj(cq)
snorm(cq::ComplexQuaternion) = cq.norm ? Complex(one(cq.qr.s)) : sqrt(snorm2(cq))


inv(cq::ComplexQuaternion) = (snorm2(cq) != 0) ? conj(cq) / snorm2(cq) : error("Not inversible")

function normalize(cq::ComplexQuaternion)
  if (cq.norm)
    return cq
  end
  a = abs(cq)
  if abs(a) > 0
    qa = cq / a
    Cquat(qa.qr, qa.qi, true)
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
    Cquat(qa.qr, qa.qi, true), a
  else
    cq, zero(Complex)
  end
end

(-)(cq::ComplexQuaternion) = ComplexQuaternion(-cq.qr, -cq.qi, cq.norm)

(+)(cq::ComplexQuaternion, cw::ComplexQuaternion) = ComplexQuaternion(cq.qr + cw.qr, cq.qi + cw.qi)
(-)(cq::ComplexQuaternion, cw::ComplexQuaternion) = ComplexQuaternion(cq.qr - cw.qr, cq.qi - cw.qi)
(*)(cq::ComplexQuaternion, cw::ComplexQuaternion) = ComplexQuaternion(cq.qr * cw.qr - cq.qi*cw.qi,
                                                              cq.qi * cw.qr + cq.qr * cw.qi,
                                                              cq.norm && cw.norm)   # TO CHECK: condition for norm of the product to be one
(/)(cq::ComplexQuaternion, cw::ComplexQuaternion) = cq * inv(cw)

function exp(cq::ComplexQuaternion)
  se = Complex(cq.qr.s, cq.qi.s)
  se = exp(se)
  cq = Cquat(quat(0.0, imag(cq.qr)), quat(0.0, imag(cq.qi)))
  cq, th = normalizea(cq)
  if cq.norm
    Cquat(se) * (Cquat(cos(th)) + cq * Cquat(sin(th)))
  else
    Cquat(se)
  end
end

function angleaxis(cq::ComplexQuaternion)
  tq = cq.qi * conj(cq.qr)
  t = [2.0 * tq.v1, 2.0 * tq.v2, 2.0 * tq.v3]
  qrs = cq.qr.s
  th0, s0 = angleaxis(cq.qr)
  sqr = quat(0.0, s0)
  if abs(abs(qrs) - one(qrs)) == 0
    th = Complex(th0, 0.5 * abs(quat(0, t)))
    th, Cquat(sqr)
  else
    th = Complex(th0, 0.5 * dot(t, s0))
    s0c1 = cross(s0, t)
    tanth = tan(th0)
    s0c2 = (s0c1 / tanth + t) * 0.5
    sqiv = cross(s0c2, s0)
    th, Cquat(sqr, quat(0.0, sqiv))
  end
end

function angle(cq::ComplexQuaternion)
  th, ax = angleaxis(cq)
  th
end

function axis(cq::ComplexQuaternion)
  th, ax = angleaxis(cq)
  ax
end

function exp(cq::ComplexQuaternion)
  se = Complex(cq.qr.s, cq.qi.s)
  se = exp(se)
  cq = Cquat(quat(0.0, imag(cq.qr)), quat(0.0, imag(cq.qi)))
  cq, th = normalizea(cq)
  if cq.norm
    Cquat(se) * (Cquat(cos(th)) + cq * Cquat(sin(th)))
  else
    Cquat(se)
  end
end

function log(cq::ComplexQuaternion)
  cq, a = normalizea(cq)
  sl = log(a)
  th, s = angleaxis(cq)
  s * Cquat(th) + Cquat(sl)
end

(^)(cq::ComplexQuaternion, cw::ComplexQuaternion) = exp(cw * log(cq))

function sqrt(cq::ComplexQuaternion)
  exp(0.5 * log(cq))
end

Cquatrand() = Cquat(quatrand(), quatrand())
nCquatrand() = normalize(Cquatrand())
