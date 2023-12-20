using SpecialFunctions
using LegendrePolynomials
using Plots
using Plots.PlotMeasures
using LaTeXStrings

global fact = zeros(typeof(big(1)*big(1)), 10000)
fact[1] = 1 

###
# @func
# factorial: Funcao para programacao dinamica de fatoriais dado o vetor global fact,
# proporcionando um tradeoff entre memoria e desempenho
# @params
# x: O fatorial a ser calculado
###
function factorial(x::T) where T<:BigInt
    if x == 0
        return 1
    end
    if fact[x] == zero(BigInt)
        i::typeof(x*x) = x
        while fact[i] == zero(typeof(x*x))
            i-=1
        end
        i+=1
        while i != x + 1
            fact[i] = fact[i-1] * i
            i+=1
        end
    end
    return fact[x]
end

###
# @func
# besbeam: Funcao para calculo do feixe de bessel ideal e de maneira direta
# @param
# psiamp: Amplitude do feixe
# axiconang: Angulo de axicon do feixe
# order: Ordem do feixe
# rho, phi, z: Coordenadas espaciais cilíndricas do ponto a ser calculado
###
function besbeam(psiamp, axiconang, order, rho, phi, z) 
    k = 1000 * 2 * big(pi) / big(1.54) #fixed wavelength 1.54mm
    kz = k * big(cos(axiconang))
    krho = k * big(sin(axiconang))
    return big(psiamp * besselj(order, krho * rho) * cis(order * phi) * cis(kz * z))
end

###
# @func
# bsccalc: Funcao de calculo dos fatores de forma da expansao em ondas parciais
# @param
# n,m: numeros inteiros de acoplamento do feixe
# axiconang: angulo de axicon do feixe expandido
###
function bsccalc(n, m, axiconang)
    fract = (im ^ (n - m)) * (2 * n + 1) * factorial(big(n - m)) / factorial(big(n+m))
    special = Plm(cos(axiconang), n, m)
    return fract * special 
end

###
# @func
# erro: Erro percentual logaritmico da aproximacao
# @param
# x,y,z: Coordenadas espaciais cartesianas do ponto analisado
###
function erro(x, y, z)
    amp = 1
    axang = deg2rad(1)
    ord = 1
    rho =  sqrt(big(x)^2 + big(y)^2)
    phi = big(atan(y, x))
    x0 = 0
    y0 = 0
    z0 = 0
    rho0 = sqrt(x0^2 + y0^2)
    phi0 = big(atan(y0, x0))
    appr = abs(partialwavexp(amp, big(axang), ord, 0, 0, 0, x - x0, y - y0, z - z0)) ^ 2
    exac = abs(besbeam(amp, big(axang), ord, sqrt(rho^2 + rho0^2 - 2 * rho * rho0 * cos(phi - phi0)), phi, z)) ^ 2
    return log10(abs(exac - appr) / exac) 
end

###
# @func
# partialwavexp: Expansao em ondas parciais do feixe analisado
# @param
# psiamp: Amplitude do feixe
# axiconang: Angulo de axicon do feixe
# order: Ordem do feixe
# r, theta, phi: Coordenadas esfericas do ponto a ser analisado
# x0, y0, z0: Coordenadas de origem do feixe
###
function partialwavexp(psiamp, axiconang, order, r, theta, phi, x0, y0, z0)
    k = 1000 * 2 * big(pi) / big(1.54)
    kr = k * r
    krho = k * sin(axiconang)
    kz = k * cos(axiconang)
    phi0 = Base.atan(y0,x0)
    rho0 = (x0^2+y0^2)^(1/2)
    nmax = Int64(ceil(kr + (big(405) / 100) * (kr^(1/3)) + 2))
    psi = 0
    besselm = []
    cisvalm = []
    for m in -nmax:nmax
        push!(cisvalm,  cis(-(m - order) * phi0) * cis(-kz * z0))
        push!(besselm, besselj(m - order, Float64(krho * rho0)))
    end
    for n in 0:nmax
        spher = sphericalbesselj(Int64(n), Float64(kr))
        for m in -n:n
            BSC = bsccalc(n, m, axiconang)
            psi += cisvalm[m + nmax + 1] * besselm[m + nmax + 1] * BSC * spher * Plm(cos(theta), n, m) * cis(m * phi)
        end
    end
    return psi * psiamp
end







