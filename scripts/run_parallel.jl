using Distributed
addprocs()

# To avoid problems in future files remove this: from julia settings
# @everywhere include("dome_final.jl")
@everywhere include("C:\\Users\\iicpe\\Desktop\\papers\\Dome_Truss_Optimization\\scripts\\RandomDomeTruss.jl")

run_algorithms(
    vcat(metaheuristics(maxevals, n_iterations),
         model_based(maxevals)))



# vcat(metaheuristics(10, n_iterations),
#      model_based(maxevals))[1:2])
