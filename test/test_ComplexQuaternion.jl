using ..Quaternions

a = ComplexQuaternion(2+2im,3+1im,4+2im,5+0im)
println("show a: ")
show(a)
println("\n")

println("qr et qi: ")
show(a.qr)
println("\n")

show(a.qi)
println("\n")

println("norm et snorm: ")
show(abs2(a))
println("\n")

show(snorm2(a))
show(abs2(snorm(a)))
println("\n")

println("scalar et vector: ")
show(scalar(a))
println("\n")

show(vector(a))
println("\n")

println("complex form: ")
show(complexForm(a))
println("\n")

show(complexFormAxis(a))
println("\n")

println("hamilton polar form: ")
show(hamiltonPolarForm(a))
println("\n")

println("arg: ")
show(arg(a))
println("\n")


#@test typeof(hamiltonPolarAmplitude(a)) = Complex

a/snorm(a)
