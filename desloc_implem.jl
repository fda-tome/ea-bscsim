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


# @func
# erro: Erro percentual logaritmico da aproximacao
# @param
# x,y,z: Coordenadas espaciais cartesianas do ponto analisado
function error_function(x, y)
    amp = 1
    axang = deg2rad(1)
    ord = 1
    rho =  sqrt(big(x)^2 + big(y)^2)
    phi = big(atan(y, x))
    x0 = 0
    y0 = 0
    z0 = 0
    z = 0
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

using StatsBase
using Random


# Initialize population
function init_population(npop, bounds)
    pop = []
    for _ in 1:npop
        x = rand() * (bounds[1, 2] - bounds[1, 1]) + bounds[1, 1]
        y = rand() * (bounds[2, 2] - bounds[2, 1]) + bounds[2, 1]
        push!(pop, (x,y))
    end
    pop
end

# Crossover operator
function crossover(parent1, parent2)
    return ((parent1[1] + parent2[1]) / 2, (parent1[2] + parent2[2]) / 2)
end

# Mutation operator
function mutate(individual, indpb, bounds)
    newtuple = [individual...]
    for i in 1:length(individual)
        if rand() < indpb
            newtuple[i] += randn() * (bounds[i, 2] - bounds[i, 1]) / 10000
            newtuple[i] = max(min(newtuple[i], bounds[i, 2]), bounds[i, 1])
        end
    end
    return (newtuple...,)
end

# Evaluate fitness
function evaluate(population)
    return [error_function(ind...) for ind in population]
end

function select(population, fitness, k, elitism_percentage)
    population_size = length(population)
    elite_size = round(Int, elitism_percentage * population_size)
    
    # Select the top individuals based on fitness for elitism
    elite_indices = sortperm(fitness)[1:elite_size]
    elite = [population[i] for i in elite_indices]
    
    # Select the remaining individuals using tournament selection
    selected_indices = []
    for _ in 1:(k - elite_size)
        tournament_indices = sample(1:population_size, 2, replace=false)
        winner_index = argmin(fitness[tournament_indices])
        push!(selected_indices, tournament_indices[winner_index])
    end
    
    return vcat(elite, [population[i] for i in selected_indices])
end


# Genetic algorithm function
function genetic_algorithm(error_function, ngen, npop, cxpb, indpb, bounds, elitism_percentage)
    # Set up initial population
    population = init_population(npop, bounds)

    bestFit = []

    avgFit = []

    offspring = []
    
    fitness = evaluate(population)

    best_index = 0

    for gen in 1:ngen
        offspring = []
        # Evaluate fitness
        fitness = evaluate(population)
        
        # Select parents
        parents = select(population, fitness, 2, elitism_percentage)
        
        # Apply crossover and mutation to create offspring
        for i in 1:npop
            if rand() < cxpb 
                push!(offspring, mutate(crossover(parents[1], parents[2]), indpb, bounds))
            else
                push!(offspring, parents[1])
            end
        end
        
        # Replace old population with offspring
        population = offspring

        best_index = argmin(fitness)

        push!(bestFit, fitness[best_index])
        
        push!(avgFit, mean(fitness))
    end
    
    # Return the best individual and its fitness
    return population[best_index], fitness[best_index], bestFit, avgFit
end

# Set genetic algorithm parameters
ngen = 100
npop = 1000
cxpb = 0.6
indpb = 0.1
elitism_percentage = 0.05
bounds = [[-0.2, -0.2] [0.2, 0.2]]

# Run the genetic algorithm
best_individual, best_fitness, bestFit, avgFit = genetic_algorithm(error_function, ngen, npop, cxpb, indpb, bounds, elitism_percentage)

p = plot(1:ngen, [bestFit, avgFit], label = ["Best Fitness" "Mean Fitness"], xlabel = "Generation", ylabel = "Fitness")

savefig(p, "plot.png")

# Display the results
println("Best Individual: ", best_individual)
println("Best Fitness: ", best_fitness)

