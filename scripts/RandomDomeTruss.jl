using KhepriFrame3DD
using KhepriAutoCAD

truss_nodes(ps) =
  map(truss_node, ps)

truss_bars(ps, qs) =
  map(truss_bar, ps, qs)

function fixedpoint(f, x0; residualnorm = (x -> norm(x,Inf)), tol = 1E-10, maxiter=100)
  residual = Inf
  iter = 1
  xold = x0
  while residual > tol && iter < maxiter
    xnew = f(xold)
    residual = residualnorm(xold - xnew)
    xold = xnew
    iter += 1
  end
  return xold
end

catenary_points(p, h, r, ϕ, n) =
  let a = -fixedpoint(a->r/acosh(h/a+1), 1.0)
    map_division(ρ -> add_cyl(p, ρ, ϕ, h + a*(cosh(ρ/a) - 1)),
                 0, r, n)
  end

# Rib

truss_rib(as, bs) =
  let fixed_truss_node_family =
        truss_node_family_element(default_truss_node_family(),
                                  support=truss_node_support(ux=true, uy=true, uz=true)),
      c = intermediate_loc(as[end], bs[end]),
      as = as[1:end-1],
      bs = bs[1:end-1]
    truss_nodes(as)
    truss_nodes(bs)
    truss_node(c, fixed_truss_node_family)
    truss_bars([c,c], [as[end], bs[end]])
    truss_bars(as, bs)
    truss_bars(as[1:end-1], as[2:end])
    truss_bars(bs[1:end-1], bs[2:end])
    truss_bars(as[1:2:end], bs[2:2:end])
    truss_bars(bs[2:2:end], as[3:2:end])
  end

dome_truss_ribs(p, h, Δ, qs, n) =
  let
    for q in qs
      v = q - p
      ρ = pol_rho(v)
      ϕ = pol_phi(v)
      truss_rib(catenary_points(p, h, ρ, ϕ, n),
                catenary_points(p, h - Δ, ρ - Δ, ϕ, n))
    end
  end

#
random_frontier(p, rmin, rmax, n) =
  map_division(ϕ->add_pol(p, random_range(rmin, rmax), ϕ), 0, 2π, n, false)

random_seed(12345)
frontier = random_frontier(u0(), 10, 14, 20)
# delete_all_shapes()
# slab(closed_polygonal_path(frontier))
# dome_truss_ribs(u0(), 9, 1, frontier, 20)
# @backend frame3dd view_truss_deformation(truss_analysis(vz(-1e4)), autocad)

#=
@time let d = 0.1,
    e = d/10,
    load = vz(-5e4),
    center = xy(1,2),
    visualizer = autocad #rhino #autocad
  with_truss_node_family(radius=d*1.2) do
    with_truss_bar_family(radius=d/2, inner_radius=d/2-e) do
      with(current_backends, (visualizer, #=robot,=# frame3dd)) do
        delete_all_shapes()
        dome_truss_ribs(center, 9, 1, frontier, 10)
        with(current_backend, visualizer) do
          with(current_layer, create_layer("Frame3DD", true, rgb(1,0,0))) do
            @backend frame3dd view_truss_deformation(truss_analysis(load), visualizer)
          end
        end
      end
    end
  end
end
=#

# -------------------------------------------------------------------------
#                      Multi-Objective Optimization
# -------------------------------------------------------------------------
using ADOPT
using ADOPT.Sampling
using ADOPT.Platypus
using ADOPT.ScikitLearnModels

using DelimitedFiles


# res_dir = "C:\\Users\\iicpe\\Desktop\\BS_Optimization_truss\\results"
res_dir = "C:\\Users\\iicpe\\Desktop\\papers\\Dome_Truss_Optimization\\results"

# train_dir = "C:\\Users\\iicpe\\Desktop\\BS_Optimization_truss\\scripts\\training_solutions.csv"

## Initial Population -----------------------------------------------------
unscale(value, nmin, nmax, omin=0, omax=1) =
      (value .- omin) ./ (omax - omin) .* (nmax - nmin) .+ nmin

get_initial_population(variables, pop_ini) = begin
   ADOPT.Sampling.set_seed!(12345)
   X = ADOPT.latinHypercube(length(variables), pop_ini)
   unscaled_X = zeros(size(X, 2), size(X, 1))

   for (i, var) in (enumerate(variables))
       unscaled_X[:, i] = unscale(X[i, :], ADOPT.lower_bound(var), ADOPT.upper_bound(var))
   end
   map(1:size(unscaled_X, 1)) do i
       round.(unscaled_X[i, :])
   end
end

## Materials Young's Modulus and Cost -------------------------------------
#=
fonts:
   https://www.engineeringtoolbox.com/young-modulus-d_417.html
   http://www.metalary.com/titanium-price/
   https://www.alibaba.com/showroom/astm-a36-steel-price.html
   http://www.matweb.com/search/QuickText.aspx?SearchText=Carbon%20Steel
   http://aesteironsteelpipes.com/american-petroleum-institute-5l-astm-a106-astm-a53-grade-b-seamless-steel-pipe-price-list-412.html

   https://mechanicalc.com/reference/material-properties-tables
   https://www.justintools.com/unit-conversion/pressure.php?k1=psi&k2=gigapascals
=#

materials_e=[164090000000.0, # Cast Iron ASTM A536
             180000000000.0, # Steel, stainless AISI 302
             200000000000.0, # Carbon Steel, Structural ASTM A516
             206840000000.0, # Alloy Steel, ASTM A242
             204770000000.0, # Alloy Steel, AISI 4140
             193000000000.0, # Stainless Steel AISI 201
            ]

materials_cost=[460.00, # Cast Iron ASTM A536
                940.00, # Steel, stainless AISI 302
                650.00, # Carbon Steel, Structural ASTM A516
                750.00, # Alloy Steel, ASTM A242
                1750.00, # Alloy Steel, AISI 4140
                1225.00, # Stainless Steel AISI 201
               ]

## Objectives -------------------------------------------------------------
#=
The objectives for the optimization are minimizing the maximum displacement,
while minimizing the material cost of the truss structure.
=#

map_to_int(x, min, step=0.01) = min + step*x

truss_volume = 0

function displacement(material, diameter, x_center, y_center)
  let d = map_to_int(diameter, 0.05),
      e = d/10,
      x_c = map_to_int(x_center, -3)
      y_c = map_to_int(y_center, -3)
      load = vxz(1e3, -1e4)
      set_backend_family(default_truss_bar_family(),
           frame3dd,
           frame3dd_truss_bar_family(
             E=materials_e[Int(material)], # (Young's modulus)
             G=81000000000.0, # (Kirchoff's or Shear modulus)
             p=0.0, # Roll angle.
             d=77010.0)) # Density
    with_truss_node_family(radius=d*1.2) do
      with_truss_bar_family(radius=d/2, inner_radius=d/2-e) do
        KhepriFrame3DD.with(current_backend, #=robot,=# frame3dd) do
          delete_all_shapes()
          dome_truss_ribs(xy(x_c, y_c), 9, 1, frontier, 10)
          global truss_volume = truss_bars_volume()
          max_displacement(truss_analysis(load))
        end
      end
    end
  end
end

function displacement_obj((material, diameter, x_center, y_center))
  let d_max_aceitavel = 0.2,
      peso_d_max = 0.1,
      d = displacement(material, diameter, x_center, y_center)

      d > d_max_aceitavel ?
        d - d_max_aceitavel :
        peso_d_max*(d - d_max_aceitavel)
    end
  end

# cost((x_center, y_center)) = 0.1 * materials_cost[5]
cost((material, diameter, x_center, y_center))=truss_volume * materials_cost[Int(material)]

objs = [Objective(displacement_obj, :MIN),
        Objective(cost, :MIN)]

## Variables --------------------------------------------------------------
material = IntVariable(1, 6)
diameter = IntVariable(0, 20) # min= 0.05, max=0.2, step=0.01
x_center = IntVariable(0, 600) # min= -3, max=3, step=0.01
y_center = IntVariable(0, 600) # min= -3, max=3, step=0.01

vars = [material, diameter, x_center, y_center]

problem = Model(vars, objs)

# shared parameters
nruns = 1
maxevals = 600
# maxevals = 600

# metaheuristics parameters
n_iterations = 25

## Problem Definition -----------------------------------------------------
### for METAHEURISTICS ------------------------------
function metaheuristics(maxevals, n_iterations)
  nparticles = div(maxevals, n_iterations)
  problem() = Model(vars, objs)
  initial_population() = get_initial_population(vars, nparticles)
  generator() = Dict(:name => ADOPT.Platypus.InjectedPopulation,
                   :solutions => initial_population(),
                   :problem => problem())

  EAs_params() = Dict(:population_size => nparticles,
                    :generator => generator())

  PSOs_params() = Dict(:leader_size => nparticles,
                     :swarm_size => nparticles,
                     :max_iterations => nparticles * 2,
                     :generator => generator())

  OMOPSO_params() = Dict(:leader_size => nparticles,
                       :swarm_size => nparticles,
                       :max_iterations => nparticles * 2,
                       :epsilons => [1, 1], # os epsilons sao referentes aos objectivos
                       :generator => generator())
  [(NSGAII, EAs_params()),
   (MOEAD, EAs_params()),
   (PESA2, EAs_params()),
   # (SPEA2, EAs_params()),
   (SMPSO, PSOs_params()),
   (OMOPSO, OMOPSO_params())]
end

### for MODEL-BASED ---------------------------------
### Test 1 - GPR
# Linear models
ADOPT.ScikitLearnModels.@sk_import linear_model: (LinearRegression, LogisticRegression, BayesianRidge, ARDRegression)
ADOPT.ScikitLearnModels.@sk_import tree: (DecisionTreeRegressor, ExtraTreeRegressor)
ADOPT.ScikitLearnModels.@sk_import ensemble: (RandomForestRegressor, ExtraTreesRegressor)
ADOPT.ScikitLearnModels.@sk_import neural_network: (MLPRegressor, BernoulliRBM)
ADOPT.ScikitLearnModels.@sk_import svm: (SVR, NuSVR, LinearSVR)
ADOPT.ScikitLearnModels.@sk_import gaussian_process: (GaussianProcessRegressor,)
ADOPT.ScikitLearnModels.@sk_import gaussian_process.kernels: (ConstantKernel, DotProduct, Matern, RBF, RationalQuadratic, WhiteKernel)

function model_based(maxevals)
  ### Metaheuristic Solvers
  nparticles = 100
  iter = 20
  nevals_mtsolver = nparticles * iter * 2
  params = Dict(:population_size => nparticles)
  pso_params = Dict(:leader_size => nparticles,
                    :swarm_size => nparticles,
                    :mutation_probability => 0.5,
                    :mutation_perturbation => 0.5)

  solvers =
    [PlatypusSolver(NSGAII, max_eval=nevals_mtsolver, algorithm_params=params, nondominated_only=true),
     PlatypusSolver(SPEA2, max_eval=nevals_mtsolver, algorithm_params=params, nondominated_only=true),
     PlatypusSolver(SMPSO, max_eval=nevals_mtsolver, algorithm_params=pso_params, nondominated_only=true)]
  ### Model-based Solver
  #=
  Below are two options for the meta solver:
     1. The first where the training solutions are generated and evaluated by the algorithm.
     2. The second where the training solutions are given to the algorithm.
  =#

  #### Option 1
  n_samples = Int(maxevals * 0.7)
  meta_solver(metamodel, solver) = let
   params = Dict(:sampling_function => randomMC, :nsamples => n_samples)
   surrogate = Surrogate(metamodel, objectives=objs, creation_f=sk_fit!,
                         update_f=sk_fit!, evaluation_f=sk_predict)
   MetaSolver(solver; surrogates=[surrogate], max_eval=maxevals, sampling_params=params, nondominated_only=true)
  end

  #### Option 2
  #=
  meta_solver(metamodel, solver) = let
   data = open(train_dir, "r") do io
       readline(io); # To skip the header in the file. If you don't have an header, remove this line
       DelimitedFiles.readdlm(io, ',', Float64, '\n')
   end
   X, y = data[1:420, 4:7]', data[1:420, 8:9]'
   params = Dict(:X => X, :y => y)
   surrogate = Surrogate(metamodel, objectives=objs, creation_f=sk_fit!,
                         update_f=sk_fit!, evaluation_f=sk_predict)
   MetaSolver(solver; surrogates=[surrogate], max_eval=maxevals, sampling_params=params, nondominated_only=true)
  end
  =#
  regressors = [GaussianProcessRegressor(), RandomForestRegressor(), ExtraTreesRegressor()]

  [meta_solver(regressor, solver)
   for regressor in regressors
   for solver in solvers]
end


# -------------------------------------------------------------------------
#                           PARALLELIZATION
# -------------------------------------------------------------------------
using Distributed

run_algorithm(algorithm) =
  ADOPT.with(ADOPT.results_dir, res_dir) do
     benchmark(nruns=nruns,
               algorithms=[algorithm],
               problem=Model(vars, objs),
               max_evals=maxevals)
      end

run_algorithms(algorithms) =
  pmap(run_algorithm, algorithms)


algorithms = vcat(metaheuristics(maxevals, n_iterations),
                  model_based(maxevals))

# run_algorithm(metaheuristics(maxevals, n_iterations)[1])
