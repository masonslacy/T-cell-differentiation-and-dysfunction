## download required packages (only need if the packages are not already downloaded)
using Pkg
Pkg.add("Plots")
Pkg.add("CairoMakie")
Pkg.add("Colors")
Pkg.add("LinearAlgebra")
Pkg.add("ProgressMeter")
Pkg.add("JLD")
Pkg.add("Statistics")
Pkg.add("Random")
Pkg.add("Distributions")
Pkg.add("LsqFit")
Pkg.add("Combinatorics")
Pkg.add("Optim")
Pkg.add("HypothesisTests")
Pkg.add("LatinHypercubeSampling")


## import packages and functions, define variables
using Plots: palette
using CairoMakie
using Colors
using LinearAlgebra
using ProgressMeter
using JLD
using Statistics
using Random
using Distributions
using LsqFit
using Combinatorics
using Optim
using HypothesisTests
using LatinHypercubeSampling

include("functions.jl") # functions used for simulating and fitting the models

# Terminology: Note that when I refer to exhausted/dysfunctional cells, I am referring to cells that upregulate inhibitory receptors associated with exhaustion (and progenitor vs terminally exhausted cells are theoretically differentiatied by their level of inhibitory receptor expression)
#              Similarly, the terms "memory" and "effector" refer to cells from the model compartments containing memory-like (stem memory and central memory) and effector-like (effector memory and terminal effector) cells
#              Otherwise inaccurate terminology is correct in the paper, and are only inaccurately used here for brevity


# define/initialise variables

# model choices that we are not changing (not considered in this study, so the code to handle these cases may not be up to date if they are toggled):
include_diff = true; #  true to include differentiation pathways between memory and effector cell types (otherwise only activation)
act_diff = false; #  true to include activation before differentiation (naive -> activated -> differentiated cell types)
prog_to_eff = false; #  true to allow progenitor exhausted cells to differentiate to effector cells
include_exh = true; #  true to include exhaustion pathways from activated/differentiated cells to exhausted cell types
prog_exh = false; #  true to include progenitor exhausted cells (step before terminal exhaustion)
no_prolif = false; #  true to remove proliferation from memory, effector and progenitor exhausted cells (but not activated cells)
cell_sec = false; #  true to assume that activated cells secrete IL-2
const_IL2 = true; #  true to use a model which assumes that IL-2 is constant ("large enough" for all IL-2 dependent rates to be at maximum)

# reaction rates that do not change as a result (no need to consider these)
f_n_a(ρ) = 0; # transition from naive to activated
f_n_p(ρ, I) = 0; # transition from naive to progenitor exhausted 
f_a_n_prolif(I) = 0; # asymmetric activated cell proliferation creating naive
f_p_n_prolif(I) = 0; # asymmetric progenitor exhausted cell proliferation creating naive
f_a_m(ρ) = 0; # transition from activated to memory
f_a_e(ρ) = 0; # transition from activated to effector
f_a_p(ρ, I) = 0; # transition from activated to progenitor exhausted
f_a_d(ρ, I) = 0; # transition from activated to terminally exhausted
f_a_a_prolif(I) = 0; # symmetric activated cell proliferation
f_a_death = 0; # activated cell death
f_m_p(ρ, I) = 0; # transition from memory to progenitor exhausted 
f_e_p(ρ, I) = 0; # transition from effector to progenitor exhausted 
f_p_a(ρ) = 0; # transition from progenitor exhausted to activated
f_p_e(ρ) = 0; # transition from progenitor exhausted to effector
f_p_p_prolif(I) = 0; # symmetric progenitor exhausted cell proliferation 
f_I_ρ(ρ, t) = 0; # IL-2 secretion from micro-rods
f_I_sec = 0; # IL-2 secretion by activated cells 
f_I_n_uptake(I) = 0; # IL-2 uptake by naive cells 
f_I_a_uptake(I) = 0; # IL-2 uptake by activated cells 
f_I_m_uptake(I) = 0; # IL-2 uptake by memory cells 
f_I_e_uptake(I) = 0; # IL-2 uptake by effector cells 
f_I_p_uptake(I) = 0; # IL-2 uptake by progenitor exhausted cells 


# get experimental data, input "true" to plot data
init_num, mean_init_num, std_init_num, final_num, mean_final_num, std_final_num = get_data(false); # init_num is at day 0 [naive/activated/memory/effector/TPEX/TEX, healthy/patient],  final_num is at day 8 [naive/activated/memory/effector/TPEX/TEX, 7 stimulus amounts, healthy/patient]


# generate 2D scaffold so that we can use an average scaffold density (in reality, the scaffold generation/selection will not affect fitting results, only parameter values) - this code was used/explored more in our previous work (DOI: 10.1098/rsos.251979)
area_scale = 0.02; # scales down the domain area from the total surface area of the well / dish
mass_scale = 0.005; # scales down the mass / number of micro-rods we consider for simulation
mean_MSR_len = 83.0; sd_MSR_len = 32.3/2; # mean and standard deviation for MSR lengths (μm)                (Zhang, et al., 2023)
mean_MSR_diam = 14.6; sd_MSR_diam = 3.3/2; # mean and standard deviation for MSR diameter (or width) (μm)           ^
MSR_SA = 570*1000000; # mean surface area of MSR per gram (μm^2/μg)                                               ^
MSR_m = (π*mean_MSR_diam*mean_MSR_len+1/4*pi*mean_MSR_diam^2)/MSR_SA; # mean mass of individual MSR (μg)
h = 5; # distance between nodes in grid (μm)
use_MvPDF = false; # true to use a multivariate PDF to define micro-rod x and y positions
A_well = 0.32; # total surface area of dish or well used for cell culture (cm^2)    -    well of a 96 well plate (Zhang et al., 2023), surface area from (ThermoFisher, “Useful Numbers for Cell Culture")
APC_ms_conc = 93.75; # concentration of inputted APC-ms (μg/cm^2)   (Zhang et al., 2023)         or 333 μg/mL (Cheung et al., 2018)
APC_ms_vol = 100e-3; # input volume of APC-ms (mL)   (Cheung et al., 2018)
m_input = APC_ms_conc*A_well; # input mass of MSRs (μg) assuming majority of APC-ms mass is MSRs         or  APC_ms_conc*APC_ms_vol  if using volume
max_d = Inf; # maximum number of micro-rods overlapped at any lattice site in the domain
N_rod_2D = mass_scale/(mean_MSR_diam/(APC_ms_vol/A_well*10000)); # number of rods stacked vertically (number of rod diameters) to approximate a 3D domain in 2D
# probability density functions (for x and y) that govern the x and y positions of the centre of micro-rods (functions of the domain width and height)
xpos_PDF(max) = Uniform(1,max); # options:  Uniform(1,max), Normal(max/2,max/5), SkewNormal(max/2,max/5,-max/10), Arcsine(1,max), ...
ypos_PDF(max) = Uniform(1,max); # options cont        Exponential(max/5), Laplace(max/2), SymTriangularDist(max/2,max/5), TriangularDist(1,max,1)
# probability density function (of θ, in rad ∈ [0, π)) that governs the rotation of micro-rods    -    check plot with:   plot(x->pdf(rot_PDF,x),xlab="θ (rad)",ylab="Probability density",legend=false)
rot_PDF = Uniform(0,π); # options:  Uniform(0,π), Normal(π/2,π/20) 

act_stim = 1; # set activating stimulus mass ratio 1 so that it can be scaled by other stimuli scales
seed_ρ = 1; # seed for ρ generation 
ρ_avg_base = gen_scaffold(); # generate scaffold density field with other scaffold-dependent quantities (detailed in gen_scaffold function)

stim_vals = [0.8, 2, 4, 6, 8, 10, 12]; # values for the activating stimulus ratio in experiments
T = 8*24*60; # simulate for 8 days (matches data)

r_d = 0.1/(24*60); # death rate (1/min), estimated   (same as λ in the paper)

# parameters for plotting 
col_nai = RGB(44/255, 143/255, 163/255); # blue for naive cells
col_mem = RGB(115/255, 37/255, 97/255); # purple for memory cells
col_eff = RGB(237/255, 62/255, 62/255); # red for effector cells
col_exh = RGB(89/255, 68/255, 19/255); # brown for exhausted cells
alpha_val = 0.25; # alpha value (opacity) for displaying data standard deviations
data_names = ["Healthy", "Patient"]; # names of datasets for plotting

set_theme!(fonts = (; regular = "Times New Roman", bold = "Times New Roman Bold")); # set fonts





## code to fit models to data 

# choices for model fitting
stim_vals_fitting = [0.8, 2, 4, 6, 8, 10, 12]; # values for the activating stimuli ratio to fit models to data - must be a subset of {0.8, 2, 4, 6, 8, 10, 12} but not necessarily all (may fit the models to less data)
donors = [2]; # indices for donor data to fit models to data - 1: healthy, 2: patient  (only select 1 to fit to the healthy model data for example, but can select both to fit models to both data sets at the same time)
solve_alg = BFGS(); # optimisation problem solver algorithm   -   Fminbox(BFGS()), BFGS(), LBFGS(), OptimizationOptimisers.Adam(0.01), LFBGS(linesearch = BackTracking()).       I find BFGS to work best
scaler = 10^10; # scale simulation results and data by this amount (estimated to be sufficient through testing) to avoid large gradients
fit_mich_consts = true; # true to fit Michaelis constants alongside maximum rates for differentiation and exhaustion rate functions
fit_prolif = true; # true to fit the maximum proliferation rate alongside maximum rates for differentiation and exhaustion rate functions 
plot_fits = false; # true to plot the fits against the data during fitting (usually just for testing/debugging - this will be thousands of plots if fitting all models)
save_data = true; # true to save all data (will not save if errors are encountered or the fitting is stopped early)
save_name = "test_name"; # file name for saved data

# parameter choices: note that parameter names here are slightly different to in the paper - r_dm, r_de and r_te are maximum rates of differentiation and dysfunction (or upregulation of inhibitory markers),  sρ_dm, sρ_de and sρ_te are the Michaelis constants for these rates, and r_p is the proliferation rate (μ)
limit_search = true; # true to limit the search in parameter space to the lower and upper bounds defined below
num_params = 7; # total number of parameters that will be fitted to data
lower_all = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]; # lower bounds for all parameters [r_dm, r_de, r_te, sρ_dm, sρ_de, sρ_te, r_p]
upper_all = [20/24/60, 20/24/60, 20/24/60, 30*ρ_avg_base, 30*ρ_avg_base, 30*ρ_avg_base, 10/24/60]; # upper bounds for all parameters [r_dm, r_de, r_te, sρ_dm, sρ_de, sρ_te, r_p]
init_set = [1/24/60, 2/24/60, 2/24/60, 10*ρ_avg_base, 10*ρ_avg_base, 13*ρ_avg_base, 2/24/60]; # initial values for all parameters, used if not using multi-start algorithm [r_dm, r_de, r_te, ρ_dm, ρ_de, ρ_te, r_p]
search_time = 10; # maximum time (s) to search for an individual model fit

multi_start = true; # true to use a multi-start approach for model fitting 
num_starts = 10; # number of initial guesses (evenly spaced in parameter space using Latin hypercube sampling)
LHS_iters = 20; # number of iterations to optimise the space-filling criterion for Latin hypercube sampling

if !multi_start
    num_starts = 1; # changes the number of simulations per model to 1 if not using a multi-start approach
end

mass_action = false; # true to use mass action terms to define differentiation and exhaustion rates (only one parameter per rate)


# further specify model feature selections (leave as [true,false] to check both)
eff_to_mem_options = [false, true]; # (1a) switch linear differentiation from effector to memory
div_diff_options = [false, true]; # (1b) divergent/branching differentiation from naive to memory and effector
bi_diff_options = [false, true]; # (1c) naive -> effector then bidirectional from effector <-> memory
no_diff_options = [false, true]; # (1d) no differentiation between memory and effector (makes sense for a divergent differentiation model)
sig_stren_options = [false, true]; # (1e) alter differentiation rates such that memory cells are produced more by less antigen stimulation
no_stim_options = [false, true]; # (1f) turn off antigen-dependent differentiation after the first differentiation event
back_diff_options = [false, true]; # (1g) back-differentiation from exhausted to effector
naive_to_exh_options = [false, true]; # (2a) allow naive cells to become exhausted 
mem_to_exh_options = [false, true]; # (2b) allow memory cells to become exhausted 
eff_to_exh_options = [false, true]; # (2c) effector cells become exhausted
asymm_div_options = [false, true]; # (3a) proliferating cells produce a new naive-like cell
eff_prolif_options = [false, true]; # (3b) effector cells proliferate
naive_prolif_options = [false, true]; # (3c) naive cells proliferate
eff_death_options = [false, true]; # (3d) effector cells die
exh_death_options = [false, true]; # (3e) exhausted cells die


# define acceptable model choice combinations (not redundant/nonsensical)
accept_model(eff_to_mem,div_diff,bi_diff,no_diff,sig_stren,no_stim,back_diff,naive_to_exh,mem_to_exh,eff_to_exh,asymm_div,eff_prolif,naive_prolif,eff_death,exh_death) = !((no_diff && (!div_diff || eff_to_mem || bi_diff || no_stim)) || (sig_stren && no_stim) || (!naive_to_exh && !mem_to_exh && !eff_to_exh) || (!asymm_div && eff_prolif && eff_death) || (eff_to_mem && bi_diff))# || (back_diff && !naive_to_exh && !mem_to_exh && eff_to_exh)); # remove specific combinations of model choices, defines the whole set of reasonable models
# can define your own set of accepted models (as a true/false expression on feature selections) if you wish to fit a smaller subset of models:
#  accept_model(eff_to_mem,div_diff,bi_diff,no_diff,sig_stren,no_stim,back_diff,naive_to_exh,mem_to_exh,eff_to_exh,asymm_div,eff_prolif,naive_prolif,eff_death,exh_death) = ...

# determine number of models and record model feature selections
models = []; # array to contain model feature selections for each model, each model listed as a Boolean array of true and false for each model choice
num_models = 0; # number of unique models
for eff_to_mem = eff_to_mem_options # switch linear differentiation from effector to memory
    for div_diff = div_diff_options # divergent differentiation from naive to memory and effector
        for bi_diff = bi_diff_options # naive -> effector then bidirectional from effector <-> memory
            for no_diff = no_diff_options # no differentiation between memory and effector (makes sense for a divergent differentiation model)
                for sig_stren = sig_stren_options # alter differentiation rates such that memory cells are produced more by less antigen stimulation
                    for no_stim = no_stim_options # turn off antigen-dependent differentiation after the first differentiation event
                        for back_diff = back_diff_options # exhausted cells back-differentiate to effector
                            for naive_to_exh = naive_to_exh_options # allow naive cells to become terminally exhausted 
                                for mem_to_exh = mem_to_exh_options # allow memory cells to become terminally exhausted 
                                    for eff_to_exh = eff_to_exh_options # effector cells become terminally exhausted
                                        for asymm_div = asymm_div_options # proliferating cells produce a new naive cell
                                            for eff_prolif = eff_prolif_options # effector cells proliferate
                                                for naive_prolif = naive_prolif_options # naive cells proliferate
                                                    for eff_death = eff_death_options # effector cells die
                                                        for exh_death = exh_death_options # terminally exhausted cells die
                                                            if accept_model(eff_to_mem,div_diff,bi_diff,no_diff,sig_stren,no_stim,back_diff,naive_to_exh,mem_to_exh,eff_to_exh,asymm_div,eff_prolif,naive_prolif,eff_death,exh_death) # ignore redundant model choices
                                                                num_models += 1; # add to the number of unique models
                                                                append!(models, [[eff_to_mem,div_diff,bi_diff,no_diff,sig_stren,no_stim,back_diff,naive_to_exh,mem_to_exh,eff_to_exh,asymm_div,eff_prolif,naive_prolif,eff_death,exh_death]]); # append model choices to array
                                                            end
                                                        end
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end
println("Fitting $(num_models) models...")

# set up for model fitting
num_stim = length(stim_vals_fitting); # number of stimulus ratios to simulate for fitting

healthy_IC = mean_init_num[[1,3,4,6],1]; # initial numbers of cells for healthy donor data 
patient_IC = mean_init_num[[1,3,4,6],2]; # initial numbers of cells for patient data
ICs = [healthy_IC, patient_IC]; # initial conditions for both healthy and patient-derived samples

if mass_action 
    fit_mich_consts = false; # don't fit Michaelis constants if there are no extra parameters for the differentiation and exhaustion rates
    upper_all[1:3] = upper_all[1:3]/ρ_avg_base; # rescale bounds and initial guesses on mass action parameters
    lower_all[1:3] = lower_all[1:3]/ρ_avg_base;
    init_set[1:3] = init_set[1:3]/ρ_avg_base; 
end
if fit_mich_consts # if fitting Michaelis constants alongside maximum rates
    if fit_prolif # if also fitting proliferation rate
        lower = lower_all; # lower bounds for parameters [r_dm, r_de, r_te, sρ_dm, sρ_de, sρ_te, r_p]
        upper = upper_all; # upper bounds for parameters [r_dm, r_de, r_te, sρ_dm, sρ_de, sρ_te, r_p]
        init_params = init_set; # initial values for parameters (initialising fitting) [r_dm, r_de, r_te, ρ_dm, ρ_de, ρ_te, r_p]
    else # only maximum rates and Michaelis constants
        num_params = 6; # reduced total number of parameters
        lower = lower_all[1:6]; # lower bounds for parameters [r_dm, r_de, r_te, sρ_dm, sρ_de, sρ_te]
        upper = upper_all[1:6]; # upper bounds for parameters [r_dm, r_de, r_te, sρ_dm, sρ_de, sρ_te]
        init_params = init_set[1:6]; # initial values for parameters (initialising fitting) [r_dm, r_de, r_te, ρ_dm, ρ_de, ρ_te]
    end
else 
    if fit_prolif # fitting maximum rates for diff and exh, and proliferation rate
        num_params = 4;
        lower = lower_all[[1:3;7]]; # lower bounds for parameters [r_dm, r_de, r_te, r_p]
        upper = upper_all[[1:3;7]]; # upper bounds for parameters [r_dm, r_de, r_te, r_p]
        init_params = init_set[[1:3;7]]; # initial values for parameters (initialising fitting) [r_dm, r_de, r_te, r_p]
    else # otherwise only fitting maximum rates for diff and exh
        num_params = 3; 
        lower = lower_all[1:3]; # lower bounds for parameters [r_dm, r_de, r_te]
        upper = upper_all[1:3]; # upper bounds for parameters [r_dm, r_de, r_te]
        init_params = init_set[1:3]; # initial values for parameters (initialising fitting) [r_dm, r_de, r_te]
    end
end
if !limit_search # if the search is not limited by any upper bound
    upper .= Inf; # replace all limits in upper by infinity
end


# initialise output arrays
model_err = zeros(num_models)*NaN; # error / loss for each fitted model to the dataset/s that it was fitted to (initialised at NaN)

nai_num = [zeros(length(stim_vals),2)*NaN for i in 1:num_models]; # array to contain the number of naive cells for each donor and stimulus ratio for each combination of model choices, initialised at NaN
mem_num = [zeros(length(stim_vals),2)*NaN for i in 1:num_models]; # memory cells 
eff_num = [zeros(length(stim_vals),2)*NaN for i in 1:num_models]; # effector cells 
exh_num = [zeros(length(stim_vals),2)*NaN for i in 1:num_models]; # terminally exhausted cells

param_fits = [zeros(length(init_params))*NaN for i in 1:num_models]; # parameter fits [r_dm, r_de, r_te, sρ_dm, sρ_de, sρ_te], initialised at NaN

num_converged = zeros(num_models); # number of times Optim converged to a fit for each model (up to num_starts)


# simulate and fit models to data
time_elapsed = 0; # total time elapsed 
for model_n = eachindex(models) # for each model to be fitted
    start_time = time(); # start time for current loop 
    # extract all feature selections of the current model
    eff_to_mem, div_diff, bi_diff, no_diff, sig_stren, no_stim, back_diff, naive_to_exh, mem_to_exh, eff_to_exh, asymm_div, eff_prolif, naive_prolif, eff_death, exh_death = models[model_n]; 

    # define rates of differentiation, exhaustion and proliferation based on current feature selections and as functions of each fitted parameter values
    r_dm_fit(ρ, r_dm_max, sρ_dm) = mass_action ? ((eff_to_mem || bi_diff || div_diff) && no_stim ? r_dm_max*ρ_avg_base : sig_stren ? r_dm_max*max(12*ρ_avg_base-ρ,0) : r_dm_max*ρ) : ((eff_to_mem || bi_diff || div_diff) && no_stim ? 1/2*r_dm_max : sig_stren ? r_dm_max*(1-ρ/(sρ_dm + ρ)) : r_dm_max*ρ/(sρ_dm + ρ)); # if removing stimulus dependence, set differentiation rate to that at the average stimulus concentration for a stimulus ratio of 1 (mass action) or half its maximum (Michaelis-Menten), otherwise if using signal strength model, change differentiation rate to memory to decrease in ρ
    r_de_fit(ρ, r_de_max, sρ_de) = mass_action ? ((!eff_to_mem && !bi_diff || div_diff) && no_stim ? r_de_max*ρ_avg_base : r_de_max*ρ) : ((!eff_to_mem && !bi_diff || div_diff) && no_stim ? 1/2*r_de_max : r_de_max*ρ/(sρ_de + ρ)); # if removing stimulus dependence, set differentiation rate to that at the average stimulus concentration for a stimulus ratio of 1 (mass action) or half its maximum (Michaelis-Menten)
    r_te_fit(ρ, r_te_max, sρ_te) = mass_action ? r_te_max*ρ : r_te_max*ρ/(sρ_te + ρ); # rate of terminal exhaustion
    r_p_fit(r_p) = r_p; # proliferation rate

    # functions fed into ODE model, dependent on feature selections and parameter values  (ignore call error warnings)
    f_n_m_fit(ρ, r_dm_max, sρ_dm) = !eff_to_mem && !bi_diff || div_diff ? r_dm_fit(ρ, r_dm_max, sρ_dm) : 0; # transition rate from naive to memory
    f_n_e_fit(ρ, r_de_max, sρ_de) = eff_to_mem || bi_diff || div_diff ? r_de_fit(ρ, r_de_max, sρ_de) : 0; # transition rate from naive to effector
    f_m_e_fit(ρ, r_de_max, sρ_de) = (!eff_to_mem || bi_diff) && !no_diff ? r_de_fit(ρ, r_de_max, sρ_de) : 0; # transition rate from memory to effector 
    f_e_m_fit(ρ, r_dm_max, sρ_dm) = (eff_to_mem || bi_diff) && !no_diff ? r_dm_fit(ρ, r_dm_max, sρ_dm) : 0; # transition rate from effector to memory
    f_n_d_fit(ρ, r_te_max, sρ_te) = naive_to_exh ? r_te_fit(ρ, r_te_max, sρ_te) : 0; # transition rate from naive to terminally exhausted
    f_m_d_fit(ρ, r_te_max, sρ_te) = mem_to_exh ? r_te_fit(ρ, r_te_max, sρ_te) : 0; # transition rate from memory to terminally exhausted 
    f_e_d_fit(ρ, r_te_max, sρ_te) = eff_to_exh ? r_te_fit(ρ, r_te_max, sρ_te) : 0; # transition rate from effector to terminally exhausted 
    f_d_e_fit(ρ, r_de_max, sρ_de) = back_diff ? r_de_fit(ρ, r_de_max, sρ_de) : 0; # transition rate from exhausted to effector 
    f_n_n_prolif_fit(r_p) = naive_prolif ? r_p_fit(r_p) : 0; # symmetric naive cell proliferation rate
    f_m_n_prolif_fit(r_p) = asymm_div ? r_p_fit(r_p) : 0; # asymmetric memory cell proliferation rate
    f_m_m_prolif_fit(r_p) = asymm_div ? 0 : r_p_fit(r_p); # symmetric memory cell proliferation rate
    f_e_n_prolif_fit(r_p) = asymm_div && eff_prolif ? r_p_fit(r_p) : 0; # asymmetric effector cell proliferation rate
    f_e_e_prolif_fit(r_p) = !asymm_div && eff_prolif ? r_p_fit(r_p) : 0; # symmetric effector cell proliferation rate
    f_e_death_fit = eff_death ? r_d : 0; # effector cell death rate
    f_d_death_fit = exh_death ? r_d : 0; # terminally exhausted cell death rate

    R_nn_fit(ρ, r_dm_max, r_de_max, r_te_max, sρ_dm, sρ_de, sρ_te, r_p) = -f_n_m_fit(ρ, r_dm_max, sρ_dm) - f_n_e_fit(ρ, r_de_max, sρ_de) - f_n_d_fit(ρ, r_te_max, sρ_te) + f_n_n_prolif_fit(r_p); # stimulus-dependent reaction term driven by naive cells for n (naive cells)
    R_nm_fit(r_p) = f_m_n_prolif_fit(r_p); # reaction term driven by memory cells for n
    R_ne_fit(r_p) = f_e_n_prolif_fit(r_p); # reaction term driven by effector cells for n
    R_mn_fit(ρ, r_dm_max, sρ_dm) = f_n_m_fit(ρ, r_dm_max, sρ_dm); # stimulus-dependent reaction term driven by naive cells for m (memory cells)
    R_mm_fit(ρ, r_de_max, r_te_max, sρ_de, sρ_te, r_p) = -f_m_e_fit(ρ, r_de_max, sρ_de) - f_m_d_fit(ρ, r_te_max, sρ_te) + f_m_m_prolif_fit(r_p); # stimulus-dependent reaction term driven by memory cells for m
    R_me_fit(ρ, r_dm_max, sρ_dm) = f_e_m_fit(ρ, r_dm_max, sρ_dm); # stimulus-dependent reaction term driven by effector cells for m
    R_en_fit(ρ, r_de_max, sρ_de) = f_n_e_fit(ρ, r_de_max, sρ_de); # stimulus-dependent reaction term driven by naive cells for e (effector cells)
    R_em_fit(ρ, r_de_max, sρ_de) = f_m_e_fit(ρ, r_de_max, sρ_de); # stimulus-dependent reaction term driven by memory cells for e
    R_ee_fit(ρ, r_dm_max, r_te_max, sρ_dm, sρ_te, r_p) = -f_e_m_fit(ρ, r_dm_max, sρ_dm) - f_e_d_fit(ρ, r_te_max, sρ_te) + f_e_e_prolif_fit(r_p) - f_e_death_fit; # stimulus-dependent reaction term driven by effector cells for e 
    R_ed_fit(ρ, r_de_max, sρ_de) = f_d_e_fit(ρ, r_de_max, sρ_de); # stimulus-dependent reaction term driven by exhausted cells for e
    R_dn_fit(ρ, r_te_max, sρ_te) = f_n_d_fit(ρ, r_te_max, sρ_te); # stimulus-dependent reaction term driven by naive cells for d (terminally exhausted/dysfunction-associated cells)
    R_dm_fit(ρ, r_te_max, sρ_te) = f_m_d_fit(ρ, r_te_max, sρ_te); # stimulus-dependent reaction term driven by memory cells for d
    R_de_fit(ρ, r_te_max, sρ_te) = f_e_d_fit(ρ, r_te_max, sρ_te); # stimulus-dependent reaction term driven by effector cells for d
    #=UPDATED, WAS CONSTANT=#R_dd_fit(ρ, r_de_max, sρ_de) = -f_d_e_fit(ρ, r_de_max, sρ_de) - f_d_death_fit; # reaction term driven by exhausted cells for d

    ODE_funcs = [R_nn_fit, R_nm_fit, R_ne_fit, R_mn_fit, R_mm_fit, R_me_fit, R_en_fit, R_em_fit, R_ee_fit, R_ed_fit, R_dn_fit, R_dm_fit, R_de_fit, R_dd_fit]; # initialise functions for ODE system


    # define all gradients of each function with respect to each parameter being fitted
    ∂r_dm_∂r_dm_max(ρ, sρ_dm) = mass_action ? ((eff_to_mem || bi_diff || div_diff) && no_stim ? ρ_avg_base : sig_stren ? max(12*ρ_avg_base-ρ,0) : ρ) : ((eff_to_mem || bi_diff || div_diff) && no_stim ? 1/2 : sig_stren ? (1-ρ/(sρ_dm + ρ)) : ρ/(sρ_dm + ρ)); # derivative of r_dm w.r.t. r_dm_max
    ∂r_dm_∂sρ_dm(ρ, r_dm_max, sρ_dm) = mass_action ? 0 : ((eff_to_mem || bi_diff || div_diff) && no_stim ? 0 : sig_stren ? r_dm_max*ρ/(sρ_dm + ρ)^2 : -r_dm_max*ρ/(sρ_dm + ρ)^2); # r_dm w.r.t. sρ_dm
    ∂r_de_∂r_de_max(ρ, sρ_de) = mass_action ? ((!eff_to_mem && !bi_diff || div_diff) && no_stim ? ρ_avg_base : ρ) : ((!eff_to_mem && !bi_diff || div_diff) && no_stim ? 1/2 : ρ/(sρ_de + ρ)); # r_de w.r.t. r_de_max
    ∂r_de_∂sρ_de(ρ, r_de_max, sρ_de) = mass_action ? 0 : ((!eff_to_mem && !bi_diff || div_diff) && no_stim ? 0 : -r_de_max*ρ/(sρ_de + ρ)^2); # r_de w.r.t. sρ_de 
    ∂r_te_∂r_te_max(ρ, sρ_te) = mass_action ? ρ : ρ/(sρ_te + ρ); # etc...
    ∂r_te_∂sρ_te(ρ, r_te_max, sρ_te) = mass_action ? 0 : -r_te_max*ρ/(sρ_te + ρ)^2;
    ∂r_p_∂r_p = 1;

    ∂f_n_m_∂r_dm_max(ρ, sρ_dm) = !eff_to_mem && !bi_diff || div_diff ? ∂r_dm_∂r_dm_max(ρ, sρ_dm) : 0; 
    ∂f_n_m_∂sρ_dm(ρ, r_dm_max, sρ_dm) = !eff_to_mem && !bi_diff || div_diff ? ∂r_dm_∂sρ_dm(ρ, r_dm_max, sρ_dm) : 0;
    ∂f_n_e_∂r_de_max(ρ, sρ_de) = eff_to_mem || bi_diff || div_diff ? ∂r_de_∂r_de_max(ρ, sρ_de) : 0;
    ∂f_n_e_∂sρ_de(ρ, r_de_max, sρ_de) = eff_to_mem || bi_diff || div_diff ? ∂r_de_∂sρ_de(ρ, r_de_max, sρ_de) : 0;
    ∂f_m_e_∂r_de_max(ρ, sρ_de) = (!eff_to_mem || bi_diff) && !no_diff ? ∂r_de_∂r_de_max(ρ, sρ_de) : 0;
    ∂f_m_e_∂sρ_de(ρ, r_de_max, sρ_de) = (!eff_to_mem || bi_diff) && !no_diff ? ∂r_de_∂sρ_de(ρ, r_de_max, sρ_de) : 0;
    ∂f_e_m_∂r_dm_max(ρ, sρ_dm) = (eff_to_mem || bi_diff) && !no_diff ? ∂r_dm_∂r_dm_max(ρ, sρ_dm) : 0;
    ∂f_e_m_∂sρ_dm(ρ, r_dm_max, sρ_dm) = (eff_to_mem || bi_diff) && !no_diff ? ∂r_dm_∂sρ_dm(ρ, r_dm_max, sρ_dm) : 0;
    ∂f_n_d_∂r_te_max(ρ, sρ_te) = naive_to_exh ? ∂r_te_∂r_te_max(ρ, sρ_te) : 0; 
    ∂f_n_d_∂sρ_te(ρ, r_te_max, sρ_te) = naive_to_exh ? ∂r_te_∂sρ_te(ρ, r_te_max, sρ_te) : 0; 
    ∂f_m_d_∂r_te_max(ρ, sρ_te) = mem_to_exh ? ∂r_te_∂r_te_max(ρ, sρ_te) : 0;
    ∂f_m_d_∂sρ_te(ρ, r_te_max, sρ_te) = mem_to_exh ? ∂r_te_∂sρ_te(ρ, r_te_max, sρ_te) : 0;
    ∂f_e_d_∂r_te_max(ρ, sρ_te) = eff_to_exh ? ∂r_te_∂r_te_max(ρ, sρ_te) : 0;
    ∂f_e_d_∂sρ_te(ρ, r_te_max, sρ_te) = eff_to_exh ? ∂r_te_∂sρ_te(ρ, r_te_max, sρ_te) : 0;
    ∂f_d_e_∂r_de_max(ρ, sρ_de) = back_diff ? ∂r_de_∂r_de_max(ρ, sρ_de) : 0;
    ∂f_d_e_∂sρ_de(ρ, r_de_max, sρ_de) = back_diff ? ∂r_de_∂sρ_de(ρ, r_de_max, sρ_de) : 0;
    ∂f_n_n_prolif_∂r_p = naive_prolif ? ∂r_p_∂r_p : 0;
    ∂f_m_n_prolif_∂r_p = asymm_div ? ∂r_p_∂r_p : 0
    ∂f_m_m_prolif_∂r_p = asymm_div ? 0 : ∂r_p_∂r_p; 
    ∂f_e_n_prolif_∂r_p = asymm_div && eff_prolif ? ∂r_p_∂r_p : 0; 
    ∂f_e_e_prolif_∂r_p = !asymm_div && eff_prolif ? ∂r_p_∂r_p : 0;

    ∂R_nn_∂r_dm_max(ρ, sρ_dm) = -∂f_n_m_∂r_dm_max(ρ, sρ_dm);
    ∂R_nn_∂sρ_dm(ρ, r_dm_max, sρ_dm) = -∂f_n_m_∂sρ_dm(ρ, r_dm_max, sρ_dm);
    ∂R_nn_∂r_de_max(ρ, sρ_de) = -∂f_n_e_∂r_de_max(ρ, sρ_de);
    ∂R_nn_∂sρ_de(ρ, r_de_max, sρ_de) = -∂f_n_e_∂sρ_de(ρ, r_de_max, sρ_de);
    ∂R_nn_∂r_te_max(ρ, sρ_te) = -∂f_n_d_∂r_te_max(ρ, sρ_te);
    ∂R_nn_∂sρ_te(ρ, r_te_max, sρ_te) = -∂f_n_d_∂sρ_te(ρ, r_te_max, sρ_te);
    ∂R_nn_∂r_p = ∂f_n_n_prolif_∂r_p; 
    ∂R_nm_∂r_p = ∂f_m_n_prolif_∂r_p; 
    ∂R_ne_∂r_p = ∂f_e_n_prolif_∂r_p; 
    ∂R_mn_∂r_dm_max(ρ, sρ_dm) = ∂f_n_m_∂r_dm_max(ρ, sρ_dm);
    ∂R_mn_∂sρ_dm(ρ, r_dm_max, sρ_dm) = ∂f_n_m_∂sρ_dm(ρ, r_dm_max, sρ_dm);
    ∂R_mm_∂r_de_max(ρ, sρ_de) = -∂f_m_e_∂r_de_max(ρ, sρ_de);
    ∂R_mm_∂sρ_de(ρ, r_de_max, sρ_de) = -∂f_m_e_∂sρ_de(ρ, r_de_max, sρ_de);
    ∂R_mm_∂r_te_max(ρ, sρ_te) = -∂f_m_d_∂r_te_max(ρ, sρ_te);
    ∂R_mm_∂sρ_te(ρ, r_te_max, sρ_te) = -∂f_m_d_∂sρ_te(ρ, r_te_max, sρ_te);
    ∂R_mm_∂r_p = ∂f_m_m_prolif_∂r_p; 
    ∂R_me_∂r_dm_max(ρ, sρ_dm) = ∂f_e_m_∂r_dm_max(ρ, sρ_dm);
    ∂R_me_∂sρ_dm(ρ, r_dm_max, sρ_dm) = ∂f_e_m_∂sρ_dm(ρ, r_dm_max, sρ_dm);
    ∂R_en_∂r_de_max(ρ, sρ_de) = ∂f_n_e_∂r_de_max(ρ, sρ_de);
    ∂R_en_∂sρ_de(ρ, r_de_max, sρ_de) = ∂f_n_e_∂sρ_de(ρ, r_de_max, sρ_de);
    ∂R_em_∂r_de_max(ρ, sρ_de) = ∂f_m_e_∂r_de_max(ρ, sρ_de);
    ∂R_em_∂sρ_de(ρ, r_de_max, sρ_de) = ∂f_m_e_∂sρ_de(ρ, r_de_max, sρ_de);
    ∂R_ee_∂r_dm_max(ρ, sρ_dm) = -∂f_e_m_∂r_dm_max(ρ, sρ_dm);
    ∂R_ee_∂sρ_dm(ρ, r_dm_max, sρ_dm) = -∂f_e_m_∂sρ_dm(ρ, r_dm_max, sρ_dm);
    ∂R_ee_∂r_te_max(ρ, sρ_te) = -∂f_e_d_∂r_te_max(ρ, sρ_te);
    ∂R_ee_∂sρ_te(ρ, r_te_max, sρ_te) = -∂f_e_d_∂sρ_te(ρ, r_te_max, sρ_te);
    ∂R_ee_∂r_p = ∂f_e_e_prolif_∂r_p; 
    ∂R_ed_∂r_de_max(ρ, sρ_de) = ∂f_d_e_∂r_de_max(ρ, sρ_de);
    ∂R_ed_∂sρ_de(ρ, r_de_max, sρ_de) = ∂f_d_e_∂sρ_de(ρ, r_de_max, sρ_de);
    ∂R_dn_∂r_te_max(ρ, sρ_te) = ∂f_n_d_∂r_te_max(ρ, sρ_te);
    ∂R_dn_∂sρ_te(ρ, r_te_max, sρ_te) = ∂f_n_d_∂sρ_te(ρ, r_te_max, sρ_te);
    ∂R_dm_∂r_te_max(ρ, sρ_te) = ∂f_m_d_∂r_te_max(ρ, sρ_te); 
    ∂R_dm_∂sρ_te(ρ, r_te_max, sρ_te) = ∂f_m_d_∂sρ_te(ρ, r_te_max, sρ_te); 
    ∂R_de_∂r_te_max(ρ, sρ_te) = ∂f_e_d_∂r_te_max(ρ, sρ_te);
    ∂R_de_∂sρ_te(ρ, r_te_max, sρ_te) = ∂f_e_d_∂sρ_te(ρ, r_te_max, sρ_te);
    ∂R_dd_∂r_de_max(ρ, sρ_de) = -∂f_d_e_∂r_de_max(ρ, sρ_de);
    ∂R_dd_∂sρ_de(ρ, r_de_max, sρ_de) = -∂f_d_e_∂sρ_de(ρ, r_de_max, sρ_de);

    # store all partial derivatives of all ODE functions 
    if fit_mich_consts # if also fitting Michaelis constants to data
        if fit_prolif # also fitting proliferation rate
            ODE_funcs = [R_nn_fit, R_nm_fit, R_ne_fit, R_mn_fit, R_mm_fit, R_me_fit, R_en_fit, R_em_fit, R_ee_fit, R_ed_fit, R_dn_fit, R_dm_fit, R_de_fit, R_dd_fit]; # store all functions for ODE inside array 

            ODE_func_derivs = [∂R_nn_∂r_dm_max 0 0 ∂R_mn_∂r_dm_max 0 ∂R_me_∂r_dm_max 0 0 ∂R_ee_∂r_dm_max 0 0 0 0 0;
                                ∂R_nn_∂r_de_max 0 0 0 ∂R_mm_∂r_de_max 0 ∂R_en_∂r_de_max ∂R_em_∂r_de_max 0 ∂R_ed_∂r_de_max 0 0 0 ∂R_dd_∂r_de_max;
                                ∂R_nn_∂r_te_max 0 0 0 ∂R_mm_∂r_te_max 0 0 0 ∂R_ee_∂r_te_max 0 ∂R_dn_∂r_te_max ∂R_dm_∂r_te_max ∂R_de_∂r_te_max 0;
                                ∂R_nn_∂sρ_dm 0 0 ∂R_mn_∂sρ_dm 0 ∂R_me_∂sρ_dm 0 0 ∂R_ee_∂sρ_dm 0 0 0 0 0;
                                ∂R_nn_∂sρ_de 0 0 0 ∂R_mm_∂sρ_de 0 ∂R_en_∂sρ_de ∂R_em_∂sρ_de 0 ∂R_ed_∂sρ_de 0 0 0 ∂R_dd_∂sρ_de;
                                ∂R_nn_∂sρ_te 0 0 0 ∂R_mm_∂sρ_te 0 0 0 ∂R_ee_∂sρ_te 0 ∂R_dn_∂sρ_te ∂R_dm_∂sρ_te ∂R_de_∂sρ_te 0
                                ∂R_nn_∂r_p ∂R_nm_∂r_p ∂R_ne_∂r_p 0 ∂R_mm_∂r_p 0 0 0 ∂R_ee_∂r_p 0 0 0 0 0]; 
        else # fitting params for diff and exh rates
            ODE_funcs = [(ρ,r_dm,r_de,r_te,sρ_dm,sρ_de,sρ_te)->R_nn_fit(ρ,r_dm,r_de,r_te,sρ_dm,sρ_de,sρ_te,rp_max), R_nm_fit(rp_max), R_ne_fit(rp_max), R_mn_fit, (ρ,r_de,r_te,sρ_de,sρ_te)->R_mm_fit(ρ,r_de,r_te,sρ_de,sρ_te,rp_max), R_me_fit, R_en_fit, R_em_fit, (ρ,r_dm,r_te,sρ_dm,sρ_te)->R_ee_fit(ρ,r_dm,r_te,sρ_dm,sρ_te,rp_max), R_ed_fit, R_dn_fit, R_dm_fit, R_de_fit, R_dd_fit]; # store all functions for ODE inside array 

            ODE_func_derivs = [∂R_nn_∂r_dm_max 0 0 ∂R_mn_∂r_dm_max 0 ∂R_me_∂r_dm_max 0 0 ∂R_ee_∂r_dm_max 0 0 0 0 0;
                                ∂R_nn_∂r_de_max 0 0 0 ∂R_mm_∂r_de_max 0 ∂R_en_∂r_de_max ∂R_em_∂r_de_max 0 ∂R_ed_∂r_de_max 0 0 0 ∂R_dd_∂r_de_max;
                                ∂R_nn_∂r_te_max 0 0 0 ∂R_mm_∂r_te_max 0 0 0 ∂R_ee_∂r_te_max 0 ∂R_dn_∂r_te_max ∂R_dm_∂r_te_max ∂R_de_∂r_te_max 0;
                                ∂R_nn_∂sρ_dm 0 0 ∂R_mn_∂sρ_dm 0 ∂R_me_∂sρ_dm 0 0 ∂R_ee_∂sρ_dm 0 0 0 0 0;
                                ∂R_nn_∂sρ_de 0 0 0 ∂R_mm_∂sρ_de 0 ∂R_en_∂sρ_de ∂R_em_∂sρ_de 0 ∂R_ed_∂sρ_de 0 0 0 ∂R_dd_∂sρ_de;
                                ∂R_nn_∂sρ_te 0 0 0 ∂R_mm_∂sρ_te 0 0 0 ∂R_ee_∂sρ_te 0 ∂R_dn_∂sρ_te ∂R_dm_∂sρ_te ∂R_de_∂sρ_te 0]; 
        end
    else 
        if fit_prolif # fitting maximum rates for diff and exh, and proliferation rate
            ODE_funcs = [(ρ,r_dm,r_de,r_te,r_p)->R_nn_fit(ρ,r_dm,r_de,r_te,ρ_dm,ρ_de,ρ_te,r_p), R_nm_fit, R_ne_fit, (ρ,r_dm)->R_mn_fit(ρ,r_dm,ρ_dm), (ρ,r_de,r_te,r_p)->R_mm_fit(ρ,r_de,r_te,ρ_de,ρ_te,r_p), (ρ,r_dm)->R_me_fit(ρ,r_dm,ρ_dm), (ρ,r_de)->R_en_fit(ρ,r_de,ρ_de), (ρ,r_de)->R_em_fit(ρ,r_de,ρ_de), (ρ,r_dm,r_te,r_p)->R_ee_fit(ρ,r_dm,r_te,ρ_dm,ρ_te,r_p), (ρ,r_de)->R_ed_fit(ρ,r_de,ρ_de), (ρ,r_te)->R_dn_fit(ρ,r_te,ρ_te), (ρ,r_te)->R_dm_fit(ρ,r_te,ρ_te), (ρ,r_te)->R_de_fit(ρ,r_te,ρ_te), (ρ,r_de)->R_dd_fit(ρ,r_de,ρ_de)]; # re-define all functions for ODE inside array (remove dependence on Michaelis constants)

            ODE_func_derivs = [ρ->∂R_nn_∂r_dm_max(ρ,ρ_dm) 0 0 ρ->∂R_mn_∂r_dm_max(ρ,ρ_dm) 0 ρ->∂R_me_∂r_dm_max(ρ,ρ_dm) 0 0 ρ->∂R_ee_∂r_dm_max(ρ,ρ_dm) 0 0 0 0 0;
                                ρ->∂R_nn_∂r_de_max(ρ,ρ_de) 0 0 0 ρ->∂R_mm_∂r_de_max(ρ,ρ_de) 0 ρ->∂R_en_∂r_de_max(ρ,ρ_de) ρ->∂R_em_∂r_de_max(ρ,ρ_de) 0 ρ->∂R_ed_∂r_de_max(ρ,ρ_de) 0 0 0 ρ->∂R_dd_∂r_de_max(ρ,ρ_de);
                                (ρ)->∂R_nn_∂r_te_max(ρ,ρ_te) 0 0 0 (ρ)->∂R_mm_∂r_te_max(ρ,ρ_te) 0 0 0 (ρ)->∂R_ee_∂r_te_max(ρ,ρ_te) 0 (ρ)->∂R_dn_∂r_te_max(ρ,ρ_te) (ρ)->∂R_dm_∂r_te_max(ρ,ρ_te) (ρ)->∂R_de_∂r_te_max(ρ,ρ_te) 0
                                ∂R_nn_∂r_p ∂R_nm_∂r_p ∂R_ne_∂r_p 0 ∂R_mm_∂r_p 0 0 0 ∂R_ee_∂r_p 0 0 0 0 0]; 
        else # only fitting maximum rates for diff and exh
            ODE_funcs = [(ρ,r_dm,r_de,r_te)->R_nn_fit(ρ,r_dm,r_de,r_te,ρ_dm,ρ_de,ρ_te,rp_max), R_nm_fit(rp_max), R_ne_fit(rp_max), (ρ,r_dm)->R_mn_fit(ρ,r_dm,ρ_dm), (ρ,r_de,r_te)->R_mm_fit(ρ,r_de,r_te,ρ_de,ρ_te,rp_max), (ρ,r_dm)->R_me_fit(ρ,r_dm,ρ_dm), (ρ,r_de)->R_en_fit(ρ,r_de,ρ_de), (ρ,r_de)->R_em_fit(ρ,r_de,ρ_de), (ρ,r_dm,r_te)->R_ee_fit(ρ,r_dm,r_te,ρ_dm,ρ_te,rp_max), (ρ,r_de)->R_ed_fit(ρ,r_de,ρ_de), (ρ,r_te)->R_dn_fit(ρ,r_te,ρ_te), (ρ,r_te)->R_dm_fit(ρ,r_te,ρ_te), (ρ,r_te)->R_de_fit(ρ,r_te,ρ_te), (ρ,r_de)->R_dd_fit(ρ,r_de,ρ_de)]; # re-define all functions for ODE inside array (remove dependence on Michaelis constants)

            ODE_func_derivs = [ρ->∂R_nn_∂r_dm_max(ρ,ρ_dm) 0 0 ρ->∂R_mn_∂r_dm_max(ρ,ρ_dm) 0 ρ->∂R_me_∂r_dm_max(ρ,ρ_dm) 0 0 ρ->∂R_ee_∂r_dm_max(ρ,ρ_dm) 0 0 0 0 0;
                                ρ->∂R_nn_∂r_de_max(ρ,ρ_de) 0 0 0 ρ->∂R_mm_∂r_de_max(ρ,ρ_de) 0 ρ->∂R_en_∂r_de_max(ρ,ρ_de) ρ->∂R_em_∂r_de_max(ρ,ρ_de) 0 ρ->∂R_ed_∂r_de_max(ρ,ρ_de) 0 0 0 ρ->∂R_dd_∂r_de_max(ρ,ρ_de);
                                (ρ)->∂R_nn_∂r_te_max(ρ,ρ_te) 0 0 0 (ρ)->∂R_mm_∂r_te_max(ρ,ρ_te) 0 0 0 (ρ)->∂R_ee_∂r_te_max(ρ,ρ_te) 0 (ρ)->∂R_dn_∂r_te_max(ρ,ρ_te) (ρ)->∂R_dm_∂r_te_max(ρ,ρ_te) (ρ)->∂R_de_∂r_te_max(ρ,ρ_te) 0]; 
        end
    end
    

    try # try to fit the model to data for each stimulus value and for healthy and patient samples, and simulate expansion with fitted model
        
        sol_and_grad_output(s, u0, p) = sol_and_grad(s, u0, T, ρ_avg_base, p, ODE_funcs, ODE_func_derivs, scaler); # get solution at time T, and gradients of the solution with respect to each parameter, as a function of the stimulus ratio s, initial condition u0, and parameter values p
        u_i(s, u0, p) = sol_and_grad_output(s, u0, p)[1]; # extract solution to ODE model
        ∂u_i(s, u0, p) = sol_and_grad_output(s, u0, p)[2]; # extract gradient of the solution 

        if multi_start # if using a multi-start approach
            plan, ~ = LHCoptim(num_starts, num_params, LHS_iters); # sample num_starts sets of parameters in num_params-dimensional space, with LHS_iters iterations to optimise the space-filling criterion
            if num_params == 7 || num_params == 6 || num_params == 3 
                param_set = scaleLHC(plan, [(lower_all[i]+upper_all[i]*1e-5,upper_all[i]-upper_all[i]*1e-5) for i in 1:num_params]); # determine the set of initial parameters
            elseif num_params == 4
                param_set = scaleLHC(plan, [(lower_all[i]+upper_all[i]*1e-5,upper_all[i]-upper_all[i]*1e-5) for i in [1:3;7]]);
            end
        end 

        lowest_loss = Inf; # initialise lowest loss for multi-start approach
        nai_num_best = []; # initialise the best model predictions
        mem_num_best = [];
        eff_num_best = [];
        exh_num_best = [];
        params_best = []; # initialise the best parameter fits
        for start_n = 1:num_starts # for each starting parameter set
            if multi_start # if using multi-start
                init_params = param_set[start_n,:]; # get current set of initial parameters
            end # otherwise just use the set parameter values

            f_optim(p) = loss_and_grad!(zeros(size(init_params)), p, mean_final_num[[1,3,4,6],:,:]/scaler, u_i, ∂u_i, stim_vals, donors); # loss function for optimiser
            g_optim!(grad, p) = loss_and_grad!(grad, p, mean_final_num[[1,3,4,6],:,:]/scaler, u_i, ∂u_i, stim_vals, donors); # gradient of the loss function for optimiser

            res = optimize(f_optim, g_optim!, lower, upper, init_params, Fminbox(solve_alg), Optim.Options(time_limit=search_time)); # results from optimisation
            fitted_params = Optim.minimizer(res); # extract fitted parameters from optimisation results

            num_converged[model_n] += Optim.converged(res); # increment counter for number of fits converged if it did converge


            # simulate ODE model with fitted parameters for each stimulus ratio 
            nai_num_curr = zeros(length(stim_vals),2); # final number of naive cells for each stimulus ratio and donor type
            mem_num_curr = zeros(length(stim_vals),2); # memory cells
            eff_num_curr = zeros(length(stim_vals),2); # effector cells 
            exh_num_curr = zeros(length(stim_vals),2); # exhausted cells
            for donor_n = 1:2 # for each donor type (1: healthy, 2: patient)
                for stim_n = eachindex(stim_vals) # for each stimulus ratio in the data (not just those used to fit parameters)
                    act_stim = stim_vals[stim_n];
                    nai_num_curr[stim_n, donor_n], mem_num_curr[stim_n, donor_n], eff_num_curr[stim_n, donor_n], exh_num_curr[stim_n, donor_n] = u_i(act_stim, ICs[donor_n], fitted_params)*scaler; # store the solutions for current stimulus ratio and initial condition 
                end 
            end 

            curr_loss = mean((nai_num_curr[:,donors]-mean_final_num[1,:,donors]).^2 + (mem_num_curr[:,donors]-mean_final_num[3,:,donors]).^2 + (eff_num_curr[:,donors]-mean_final_num[4,:,donors]).^2 + (exh_num_curr[:,donors]-mean_final_num[6,:,donors]).^2); # get current error / loss for current model against only the dataset/s (healthy/patient) it was fitted to

            if curr_loss < lowest_loss # if this loss is lower than the previous lowest, store the new best fits and parameter values 
                lowest_loss = curr_loss; # overwrite the lowest loss
                nai_num_best = nai_num_curr;
                mem_num_best = mem_num_curr;
                eff_num_best = eff_num_curr;
                exh_num_best = exh_num_curr;
                params_best = fitted_params;
            end
        end # end for each start in the multi-start approach

        param_fits[model_n] = params_best; # store fitted parameters in array
        nai_num[model_n] = nai_num_best; # record final numbers of naive cells in complete array
        mem_num[model_n] = mem_num_best; # memory cells 
        eff_num[model_n] = eff_num_best; # effector cells 
        exh_num[model_n] = exh_num_best; # exhausted cells 
        model_err[model_n] = lowest_loss; # store the best (lowest) loss

        if plot_fits # if plotting the fitted models
            for donor_n = donors # for each donor dataset used for fitting
                pop_plt = CairoMakie.Figure(fontsize=25, size=(700,600)); # figure for plotting

                if fit_mich_consts # depending on the parameters fitted, alter the title
                    if fit_prolif
                        ax = CairoMakie.Axis(pop_plt[1, 1], xlabel="Stimulus-to-micro-rod mass ratio (ng/μg)", ylabel="Cells at day $(T÷24÷60) (millions)", limits = (0, maximum(stim_vals), 0, maximum(mean_final_num[[1,3,4,6],:,donor_n]+std_final_num[[1,3,4,6],:,donor_n])/1e6), title="$(data_names[donor_n]) model $(model_n): $(models[model_n]*1)\nr_dm ≈ $(round(params_best[1]*24*60,sigdigits=3)) day⁻¹, r_de ≈ $(round(params_best[2]*24*60,sigdigits=3)) day⁻¹, r_te ≈ $(round(params_best[3]*24*60,sigdigits=3)) day⁻¹\ns_dm ≈ $(round(params_best[4]/ρ_avg_base,sigdigits=2)) ng/μg, s_de ≈ $(round(params_best[5]/ρ_avg_base,sigdigits=2)) ng/μg, s_te ≈ $(round(params_best[6]/ρ_avg_base,sigdigits=2)) ng/μg\nr_p ≈ $(round(params_best[7]*24*60,sigdigits=3)) day⁻¹");
                    else
                        ax = CairoMakie.Axis(pop_plt[1, 1], xlabel="Stimulus-to-micro-rod mass ratio (ng/μg)", ylabel="Cells at day $(T÷24÷60) (millions)", limits = (0, maximum(stim_vals), 0, maximum(mean_final_num[[1,3,4,6],:,donor_n]+std_final_num[[1,3,4,6],:,donor_n])/1e6), title="$(data_names[donor_n]) model $(model_n): $(models[model_n]*1)\nr_dm ≈ $(round(params_best[1]*24*60,sigdigits=3)) day⁻¹, r_de ≈ $(round(params_best[2]*24*60,sigdigits=3)) day⁻¹, r_te ≈ $(round(params_best[3]*24*60,sigdigits=3)) day⁻¹\ns_dm ≈ $(round(params_best[4]/ρ_avg_base,sigdigits=2)) ng/μg, s_de ≈ $(round(params_best[5]/ρ_avg_base,sigdigits=2)) ng/μg, s_te ≈ $(round(params_best[6]/ρ_avg_base,sigdigits=2)) ng/μg");
                    end
                else
                    if fit_prolif
                        ax = CairoMakie.Axis(pop_plt[1, 1], xlabel="Stimulus-to-micro-rod mass ratio (ng/μg)", ylabel="Cells at day $(T÷24÷60) (millions)", limits = (0, maximum(stim_vals), 0, maximum(mean_final_num[[1,3,4,6],:,donor_n]+std_final_num[[1,3,4,6],:,donor_n])/1e6), title="$(data_names[donor_n]) model $(model_n): $(models[model_n]*1)\nr_dm ≈ $(round(params_best[1]*24*60,sigdigits=3)) day⁻¹, r_de ≈ $(round(params_best[2]*24*60,sigdigits=3)) day⁻¹, r_te ≈ $(round(params_best[3]*24*60,sigdigits=3)) day⁻¹\nr_p ≈ $(round(params_best[4]*24*60,sigdigits=3)) day⁻¹");
                    else
                        ax = CairoMakie.Axis(pop_plt[1, 1], xlabel="Stimulus-to-micro-rod mass ratio (ng/μg)", ylabel="Cells at day $(T÷24÷60) (millions)", limits = (0, maximum(stim_vals), 0, maximum(mean_final_num[[1,3,4,6],:,donor_n]+std_final_num[[1,3,4,6],:,donor_n])/1e6), title="$(data_names[donor_n]) model $(model_n): $(models[model_n]*1)\nr_dm ≈ $(round(params_best[1]*24*60,sigdigits=3)) day⁻¹, r_de ≈ $(round(params_best[2]*24*60,sigdigits=3)) day⁻¹, r_te ≈ $(round(params_best[3]*24*60,sigdigits=3)) day⁻¹");
                    end
                end
                
                CairoMakie.lines!(stim_vals, mean_final_num[1,:,donor_n]/1e6; color=col_nai); # plot naive cells  (simulated)
                CairoMakie.band!(stim_vals, (mean_final_num[1,:,donor_n]-std_final_num[1,:,donor_n])/1e6, (mean_final_num[1,:,donor_n]+std_final_num[1,:,donor_n])/1e6; color=(col_nai,alpha_val)); # naive cells (data)
                CairoMakie.lines!(stim_vals, nai_num_best[:,donor_n]/1e6, color=col_nai, linestyle=:dash, linewidth=3) ;
                CairoMakie.lines!(stim_vals, mean_final_num[3,:,donor_n]/1e6; color=col_mem); # memory cells  (simulated)
                CairoMakie.band!(stim_vals, (mean_final_num[3,:,donor_n]-std_final_num[3,:,donor_n])/1e6, (mean_final_num[3,:,donor_n]+std_final_num[3,:,donor_n])/1e6; color=(col_mem,alpha_val)); # memory cells (data)
                CairoMakie.lines!(stim_vals, mem_num_best[:,donor_n]/1e6, color=col_mem, linestyle=:dash, linewidth=3);
                CairoMakie.lines!(stim_vals, mean_final_num[4,:,donor_n]/1e6; color=col_eff); # effector cells  (simulated)
                CairoMakie.band!(stim_vals, (mean_final_num[4,:,donor_n]-std_final_num[4,:,donor_n])/1e6, (mean_final_num[4,:,donor_n]+std_final_num[4,:,donor_n])/1e6; color=(col_eff,alpha_val)); # effector cells (data)
                CairoMakie.lines!(stim_vals, eff_num_best[:,donor_n]/1e6, color=col_eff, linestyle=:dash, linewidth=3);
                CairoMakie.lines!(stim_vals, mean_final_num[6,:,donor_n]/1e6; color=col_exh); # exhausted cells (simulated)
                CairoMakie.band!(stim_vals, (mean_final_num[6,:,donor_n]-std_final_num[6,:,donor_n])/1e6, (mean_final_num[6,:,donor_n]+std_final_num[6,:,donor_n])/1e6; color=(col_exh,alpha_val)); # exhausted cells (data)
                CairoMakie.lines!(stim_vals, exh_num_best[:,donor_n]/1e6, color=col_exh, linestyle=:dash, linewidth=3);
                
                display(pop_plt); # display current plot
            end
        end

    catch e # if an error was encountered somewhere during fitting or simulation
        if isa(e, InterruptException) # if code was interrupted 
            println("Stopping fitting early.")
            save_data = false; # don't save data to file
            break; # exit loop over all models
        else # otherwise another error occurred during fitting or plotting 
            println(e)
            if plot_fits
                println("Could not fit parameters, or error occurred during plotting.")
            else
                println("Could not fit parameters.")
            end
        end
    end

    time_elapsed += time() - start_time; # add to the elapsed time
    time_left = time_elapsed*(num_models/model_n-1); # estimated remaining runtime (s)
    println("$(round(model_n/num_models*100,digits=2))% complete.  Time left: $(Int(time_left÷86400))d, $(Int((time_left%86400)÷3600))h, $(Int((time_left%3600)÷60))m, $(Int(round(time_left%60)))s") # compute and print estimated remaining computation time
end # end for each model
if all(isnan.(model_err))
    save_data = false; # if no models were successfully fitted, do not save data
end

println("Total time: $(Int(time_elapsed÷86400))d, $(Int((time_elapsed%86400)÷3600))h, $(Int((time_elapsed%3600)÷60))m, $(Int(round(time_elapsed%60)))s") # print total time elapsed

# save key results to files if the code completed successfully (was not interrupted)
if save_data
    if any(isnan.(model_err)) # warning for any NaN values in outputs
        println("WARNING: At least one model was not fitted.")
        # in this case, may need to fit the unfitted model/s separately and replace them in the saved file
    end
    try 
        save("$(save_name).jld","models",models,"model_err",model_err,"param_fits",param_fits,"nai_num",nai_num,"mem_num",mem_num,"eff_num",eff_num,"exh_num",exh_num,"num_converged",num_converged) # save all relevant data to file
    catch 
        println("WARNING: Could not save data.")
    end
end







## code to analyse and plot results in the main text  (except for Figs 3 and 4 which are in the section below this one)
#  code to generate SI Figs and Videos is at the bottom

use_arrows2d = true; # decide which function to use for plotting arrows in pathways (depends on version of CairoMakie - if code throws an error then toggle this variable):  true to use "arrows2d" function, false to use "arrows"

# The code below loads in data that I generated using the code above (the data was saved and is being loaded in without relying on the above code as it takes a significant amount of time to run on typical computers)
# load data on models fitted to healthy cell data set (except for models with exhausted cell back-differentiation)
data = load("healthy_10starts.jld"); # load models fitted to healthy dataset
models = data["models"]; # model choices:  [eff_to_mem, div_diff, bi_diff, no_diff, sig_stren, no_stim, back_diff, naive_to_exh, mem_to_exh, eff_to_exh, asymm_div, eff_prolif, naive_prolif, eff_death, exh_death]
param_fits_h = data["param_fits"]; # parameter fits:  [r_dm_max, r_de_max, r_te_max, s_dm, s_de, s_te, r_p]
nai_num_h = data["nai_num"]; # number of each cell subtype:  (stimulus ratio, healthy/patient)
mem_num_h = data["mem_num"];
eff_num_h = data["eff_num"];
exh_num_h = data["exh_num"];
model_err_h_nonnormal = data["model_err"]; # non-normalised healthy loss
model_err_h = model_err_h_nonnormal/mean(model_err_h_nonnormal); # normalised healthy loss
num_converged_h = data["num_converged"]; # number of starts that converged to a fit for each model

# load data on models fitted to patient cell data set (except for models with exhausted cell back-differentiation)
data = load("patient_10starts.jld"); # load models fitted to patient dataset
param_fits_p = data["param_fits"]; # parameter fits:  [r_dm_max, r_de_max, r_te_max, s_dm, s_de, s_te, r_p]
nai_num_p = data["nai_num"]; # number of each cell subtype:  (stimulus ratio, healthy/patient)
mem_num_p = data["mem_num"];
eff_num_p = data["eff_num"];
exh_num_p = data["exh_num"];
model_err_p_nonnormal = data["model_err"]; # non-normalised patient loss
model_err_p = model_err_p_nonnormal/mean(model_err_p_nonnormal); # normalised patient loss
num_converged_p = data["num_converged"]; # number of starts that converged to a fit for each model


# load in models including exhausted cell back-differentiation and add them to the set of all models (these were originally fitted separately)
data = load("healthy_10starts_backdiffmodels.jld"); # load models fitted to healthy dataset
append!(models, data["models"]); # model choices:  [eff_to_mem, div_diff, bi_diff, no_diff, sig_stren, no_stim, back_diff, naive_to_exh, mem_to_exh, eff_to_exh, asymm_div, eff_prolif, naive_prolif, eff_death, exh_death]
append!(param_fits_h, data["param_fits"]); # parameter fits:  [r_dm_max, r_de_max, r_te_max, s_dm, s_de, s_te, r_p]
append!(nai_num_h, data["nai_num"]); # number of each cell subtype:  (stimulus ratio, healthy/patient)
append!(mem_num_h, data["mem_num"]);
append!(eff_num_h, data["eff_num"]);
append!(exh_num_h, data["exh_num"]);
append!(model_err_h_nonnormal, data["model_err"]); # non-normalised healthy loss
append!(model_err_h, data["model_err"]/mean(data["model_err"])); # normalised healthy loss
append!(num_converged_h, data["num_converged"]); # number of starts that converged to a fit for each model

data = load("patient_10starts_backdiffmodels.jld"); # load models fitted to patient dataset
append!(param_fits_p, data["param_fits"]); # parameter fits:  [r_dm_max, r_de_max, r_te_max, s_dm, s_de, s_te, r_p]
append!(nai_num_p, data["nai_num"]); # number of each cell subtype:  (stimulus ratio, healthy/patient)
append!(mem_num_p, data["mem_num"]);
append!(eff_num_p, data["eff_num"]);
append!(exh_num_p, data["exh_num"]);
append!(model_err_p_nonnormal, data["model_err"]); # non-normalised patient loss
append!(model_err_p, data["model_err"]/mean(data["model_err"])); # normalised patient loss
append!(num_converged_p, data["num_converged"]); # number of starts that converged to a fit for each model



num_models = length(models); # number of models
num_choices = length(models[1]); # number of model choices per model
choice_names = ["1a","1b","1c","1d","1e","1f","1g","2a","2b","2c","3a","3b","3c","3d","3e"]; # names of each feature for plotting

ignore_1a_func(model) = model[3] || model[4]; # boolean function for any model to determine whether it has features 1c or 1d, which are features that overwrite 1a (so do not consider these models to say anything about feature 1a)
ignore_1a = findall(x->x==true, ignore_1a_func.(models)); # indices of models to ignore for feature 1a


err_order_h = sortperm(model_err_h); # model indices in order from lowest to highest loss (healthy)
err_order_p = sortperm(model_err_p); #  "  (patient)

A = hcat(ones(length(model_err_h)), model_err_h); # construct matrix equation to compute coefficients in linear model for healthy vs patient loss
coeffs = A \ model_err_p; # coefficients [c, m] in linear model y=mx+c

test = CorrelationTest(model_err_h, model_err_p); # define test to check the significance of the correlation between the healthy and patient model loss
p = pvalue(test); # find the p-value for the hypothesis test (p < 0.05 indicates that there is strong evidence that the data is correlated)
println("Correlation: $(cor(model_err_h,model_err_p)),  p-value: $p"); # print correlation between healthy and patient model error


model_err_avg = (model_err_h + model_err_p)/2; # average normalised loss from healthy and patient donors
err_order_avg = sortperm(model_err_avg); # model indices in order from lowest to highest average normalised loss


# compute normalised loss of the standard deviation lines to the mean (normalised average variance) for each dataset
std_loss_h = mean(std_final_num[1,:,1].^2 + std_final_num[3,:,1].^2 + std_final_num[4,:,1].^2 + std_final_num[6,:,1].^2)/mean(model_err_h_nonnormal); # healthy 
std_loss_p = mean(std_final_num[1,:,2].^2 + std_final_num[3,:,2].^2 + std_final_num[4,:,2].^2 + std_final_num[6,:,2].^2)/mean(model_err_p_nonnormal); # patient

inbounds_h = trues(num_models); # true if a model is within the mean +- std for all datapoints and for each of memory, effector and exhausted populations
inbounds_p = trues(num_models); 
for model_n = eachindex(models) # for each model
    if !(all(mean_final_num[3,:,1]-std_final_num[3,:,1].<mem_num_h[model_n][:,1].<mean_final_num[3,:,1]+std_final_num[3,:,1]) && all(mean_final_num[4,:,1]-std_final_num[4,:,1].<eff_num_h[model_n][:,1].<mean_final_num[4,:,1]+std_final_num[4,:,1]) && all(mean_final_num[6,:,1]-std_final_num[6,:,1].<exh_num_h[model_n][:,1].<mean_final_num[6,:,1]+std_final_num[6,:,1])) # if the model is not within the bounds for each cell type in the healthy dataset
        inbounds_h[model_n] = false; # that model fails the test for the healthy dataset
    end 
    if !(all(mean_final_num[3,:,2]-std_final_num[3,:,2].<mem_num_p[model_n][:,2].<mean_final_num[3,:,2]+std_final_num[3,:,2]) && all(mean_final_num[4,:,2]-std_final_num[4,:,2].<eff_num_p[model_n][:,2].<mean_final_num[4,:,2]+std_final_num[4,:,2]) && all(mean_final_num[6,:,2]-std_final_num[6,:,2].<exh_num_p[model_n][:,2].<mean_final_num[6,:,2]+std_final_num[6,:,2])) # if the model is not within the bounds for each cell type in the healthy dataset
        inbounds_p[model_n] = false; # that model fails the test for the patient dataset
    end 
end

# define the top models to be those that lie within the bounds for both healthy and patient dataset
inbounds_perc = 1; # define percentage of top models that may lie outside of the bounds
top_model_num = 0; # number of models to take from the "top" models
not_in_bounds = 0; # number of models in the top models that do not lie within all data bounds
for model_n in err_order_avg # for each model, ordered by increasing average loss 
    if !(inbounds_h[model_n] && inbounds_p[model_n]) # if the current model is not within the bounds somewhere in either dataset
        not_in_bounds += 1; # increment the number of models that are not in all bounds
    end
    if not_in_bounds/(top_model_num+1)*100 > inbounds_perc # if the number of models not in the bounds has exceeded the limit
        not_in_bounds -= 1; # remove the last model
        break; # exit the loop and stop adding to top models
    end
    top_model_num += 1; # increment the number of top models 
end

top_model_num_h = 0; # number of top healthy models
not_in_bounds = 0; # number of healthy models in the top models that do not lie within all healthy data bounds
for model_n in err_order_h # for each model, ordered by increasing healthy loss 
    if !inbounds_h[model_n] # if the current model is not within the bounds somewhere in the healthy dataset
        not_in_bounds += 1; # increment the number of models that are not in the bounds
    end
    if not_in_bounds/(top_model_num_h+1)*100 > inbounds_perc # if the number of models not in the bounds has exceeded the limit
        not_in_bounds -= 1; # remove the last model
        break; # exit the loop and stop adding to top models
    end
    top_model_num_h += 1; # increment the number of top models 
end

top_model_num_p = 0; # number of top patient models
not_in_bounds = 0; # number of patient models in the top models that do not lie within all patient data bounds
for model_n in err_order_p # for each model, ordered by increasing patient loss 
    if !inbounds_p[model_n] # if the current model is not within the bounds somewhere in the patient dataset
        not_in_bounds += 1; # increment the number of models that are not in the bounds
    end
    if not_in_bounds/(top_model_num_p+1)*100 > inbounds_perc # if the number of models not in the bounds has exceeded the limit
        not_in_bounds -= 1; # remove the last model
        break; # exit the loop and stop adding to top models
    end
    top_model_num_p += 1; # increment the number of top models 
end




# Fig 2a: model loss of healthy and patient models separately
# healthy model loss (left)
colours = ifelse.(inbounds_h[err_order_h], :blue, :red); # define colours to differentiate between models that lie within data bounds (blue) and those that don't (red)
plt = CairoMakie.Figure(fontsize=28); 
ax = CairoMakie.Axis(plt[1, 1], xlabel="Models ordered by increasing healthy loss", ylabel="Normalised healthy model loss");
CairoMakie.barplot!(eachindex(err_order_h), model_err_h[err_order_h], color=colours)
CairoMakie.barplot!([0],[0], color=:blue, label="Within one std dev") 
CairoMakie.barplot!([0],[0], color=:red, label="Outside one std dev")
CairoMakie.axislegend(position = :lt)
display(plt)
save("figs_and_videos\\Fig_2a_left.png", plt)

# patient model loss (right)
colours = ifelse.(inbounds_p[err_order_p], :blue, :red); # define colours to differentiate between models that lie within data bounds (blue) and those that don't (red)
plt = CairoMakie.Figure(fontsize=28); 
ax = CairoMakie.Axis(plt[1, 1], xlabel="Models ordered by increasing patient loss", ylabel="Normalised patient model loss");
CairoMakie.barplot!(eachindex(err_order_p), model_err_p[err_order_p], color=colours)
CairoMakie.barplot!([0],[0], color=:blue, label="Within one std dev") 
CairoMakie.barplot!([0],[0], color=:red, label="Outside one std dev")
display(plt)
save("figs_and_videos\\Fig_2a_right.png", plt)


case_study_inds = Int.(zeros(4)); # indices for models used in the case studies
case_study_inds[1] = findfirst(x->x==[true, false, false, false, false, false, true,  false, false, true,  false, false, false, false, true], models); # index for base pathway in case study 1
case_study_inds[2] = findfirst(x->x==[false, false, true, false, false, false, true,  false, false, true,  false, false, false, false, true], models); # index for altered pathway in case study 1
case_study_inds[3] = findfirst(x->x==[true, true, false, false, false, false, false,  false, true, false,  false, false, false, true, false], models); # index for base pathway in case study 2
case_study_inds[4] = findfirst(x->x==[false, true, false, false, false, false, false,  false, true, false,  false, false, false, true, false], models); # index for altered pathway in case study 2
not_case_study = .~in.(1:num_models, Ref(case_study_inds)); # returns 1 if the model is not a case study model
case_study_colours = [(:springgreen2, 0.8), (:springgreen2, 0.8), (:darkorange1, 0.8), (:darkorange1, 0.8)]; # marker colours for case study models 
case_study_markers = [:xcross, :circle, :xcross, :circle]; # marker styles for case study models 
case_study_markersize = 25; # markersize for case study models

# Fig 2b: healthy model loss against patient model loss
plt = CairoMakie.Figure(fontsize=26); 
ax = CairoMakie.Axis(plt[1, 1], xlabel="Normalised healthy model loss", ylabel="Normalised patient model loss", aspect = DataAspect());
CairoMakie.scatter!(model_err_h[not_case_study], model_err_p[not_case_study], color=(:black, 0.25)) 
CairoMakie.scatter!(model_err_h[case_study_inds], model_err_p[case_study_inds], color=case_study_colours, marker=case_study_markers, markersize=case_study_markersize, strokecolor=:black, strokewidth=1.5) # add in case study models
display(plt)
save("figs_and_videos\\Fig_2b.png", plt)

plt = CairoMakie.Figure(fontsize=25); # zoomed in plot
ax = CairoMakie.Axis(plt[1, 1], xlabel="Normalised healthy model loss", ylabel="Normalised patient model loss", aspect = DataAspect(), limits=(0,1,0,1));
CairoMakie.scatter!(model_err_h[not_case_study][model_err_h[not_case_study].<=1 .&& model_err_p[not_case_study].<=1], model_err_p[not_case_study][model_err_h[not_case_study].<=1 .&& model_err_p[not_case_study].<=1], color=(:black, 0.25))
for i in eachindex(case_study_inds) # for each case study model 
    if model_err_h[case_study_inds[i]]<=1 && model_err_p[case_study_inds[i]]<=1 # if the model should be plotted 
        CairoMakie.scatter!(model_err_h[case_study_inds[i]], model_err_p[case_study_inds[i]], color=case_study_colours[i], marker=case_study_markers[i], markersize=case_study_markersize, strokecolor=:black, strokewidth=1.5) # add in case study model
    end
end
display(plt)
save("figs_and_videos\\Fig_2b_zoomed.png", plt)



# Fig 5a: averaged model loss
colours = ifelse.(inbounds_h[err_order_avg] .&& inbounds_p[err_order_avg], :blue, :red); # define colours to differentiate between models that lie within data bounds (blue) and those that don't (red)
plt = CairoMakie.Figure(fontsize=28); 
ax = CairoMakie.Axis(plt[1, 1], xlabel="Models ordered by increasing average loss", ylabel="Average normalised model loss", limits=(-400,num_models+400,-0.1,maximum(model_err_avg)+0.2));
CairoMakie.barplot!(eachindex(err_order_avg), model_err_avg[err_order_avg], color=colours)
CairoMakie.vlines!(top_model_num, linestyle=:dash, color=:black, linewidth=3)
CairoMakie.barplot!([0],[0], color=:blue, label="Within one std dev") # dummy plots for making the legend
CairoMakie.barplot!([0],[0], color=:red, label="Outside one std dev")
#CairoMakie.scatter!([findfirst(x->x==i, err_order_avg) for i in case_study_inds], model_err_avg[case_study_inds].+0.15, color=case_study_colours, marker=case_study_markers, markersize=case_study_markersize, strokecolor=:black, strokewidth=1.5)
CairoMakie.axislegend(position = :lt)
display(plt)
save("figs_and_videos\\Fig_5a.png", plt)



pair_num_tot = zeros(num_choices,num_choices); # number of times a pair of model choices appears in all models 
pair_num_h = zeros(num_choices,num_choices); # number of times a pair of model choices appears in the top healthy models 
pair_num_p = zeros(num_choices,num_choices); # number of times a pair of model choices appears in the top patient models 
pair_num_avg = zeros(num_choices,num_choices); # number of times a pair of model choices appears in the top models according to average loss
pair_num_tot_1aoff = zeros(num_choices+1,num_choices+1); # number of times a pair of model choices appears in all models (also considering 1a off as a feature selection)
pair_num_h_1aoff = zeros(num_choices+1,num_choices+1); # number of times a pair of model choices appears in the top healthy models (also considering 1a off as a feature selection)
pair_num_p_1aoff = zeros(num_choices+1,num_choices+1); # number of times a pair of model choices appears in the top patient models (also considering 1a off as a feature selection)
pair_num_avg_1aoff = zeros(num_choices+1,num_choices+1); # number of times a pair of model choices appears in the top models according to average loss (also considering 1a off as a feature selection)
for model_n = eachindex(models) # for each model
    curr_model = models[model_n]; # current model
    for choice_n = findall(x->x==1,curr_model) # for each individual model choice in the current model 
        pair_num_tot[choice_n,choice_n] += 1; # add to the number of times this model choice was used
    end
    for comb = combinations(findall(x->x==1,curr_model),2) # for each combination of model choices in the current model
        pair_num_tot[comb[1],comb[2]] += 1; # add to the number of times this pair of model choices is used 
    end
    if !curr_model[1] && !(model_n in ignore_1a) # if 1a is truly off (not overwritten)
        pair_num_tot_1aoff[1,1] += 1; # increment the number of times feature 1a is off
        for feat_n = 2:num_choices # for each other feature
            if curr_model[feat_n] # if that feature is on (alongside 1a off)
                pair_num_tot_1aoff[1,feat_n+1] += 1; # increment that pair
            end
        end
    end

    if model_n in err_order_h[1:top_model_num]#or top_model_num_h  # if the current model is in the top healthy models
        for choice_n = findall(x->x==1,curr_model) # for each individual model choice in the current model 
            pair_num_h[choice_n,choice_n] += 1; # add to the number of times this model choice was used
        end
        for comb = combinations(findall(x->x==1,curr_model),2) # for each combination of model choices in the current model
            pair_num_h[comb[1],comb[2]] += 1; # add to the number of times this pair of model choices is used 
        end
        if !curr_model[1] && !(model_n in ignore_1a) # if 1a is truly off (not overwritten)
            pair_num_h_1aoff[1,1] += 1; # increment the number of times feature 1a is off
            for feat_n = 2:num_choices # for each other feature
                if curr_model[feat_n] # if that feature is on (alongside 1a off)
                    pair_num_h_1aoff[1,feat_n+1] += 1; # increment that pair
                end
            end
        end
    end

    if model_n in err_order_p[1:top_model_num]#or top_model_num_p  # if the current model is in the top patient models
        for choice_n = findall(x->x==1,curr_model) # for each individual model choice in the current model 
            pair_num_p[choice_n,choice_n] += 1; # add to the number of times this model choice was used
        end
        for comb = combinations(findall(x->x==1,curr_model),2) # for each combination of model choices in the current model
            pair_num_p[comb[1],comb[2]] += 1; # add to the number of times this pair of model choices is used 
        end
        if !curr_model[1] && !(model_n in ignore_1a) # if 1a is truly off (not overwritten)
            pair_num_p_1aoff[1,1] += 1; # increment the number of times feature 1a is off
            for feat_n = 2:num_choices # for each other feature
                if curr_model[feat_n] # if that feature is on (alongside 1a off)
                    pair_num_p_1aoff[1,feat_n+1] += 1; # increment that pair
                end
            end
        end
    end

    if model_n in err_order_avg[1:top_model_num] # if the current model is in the top models according to the average loss
        for choice_n = findall(x->x==1,curr_model) # for each individual model choice in the current model 
            pair_num_avg[choice_n,choice_n] += 1; # add to the number of times this model choice was used
        end
        for comb = combinations(findall(x->x==1,curr_model),2) # for each combination of model choices in the current model
            pair_num_avg[comb[1],comb[2]] += 1; # add to the number of times this pair of model choices is used 
        end
        if !curr_model[1] && !(model_n in ignore_1a) # if 1a is truly off (not overwritten)
            pair_num_avg_1aoff[1,1] += 1; # increment the number of times feature 1a is off
            for feat_n = 2:num_choices # for each other feature
                if curr_model[feat_n] # if that feature is on (alongside 1a off)
                    pair_num_avg_1aoff[1,feat_n+1] += 1; # increment that pair
                end
            end
        end
    end
end
pair_num_tot_1aoff[2:num_choices+1,2:num_choices+1] = pair_num_tot; # add in results for all feature selections other than 1a off
pair_num_h_1aoff[2:num_choices+1,2:num_choices+1] = pair_num_h;
pair_num_p_1aoff[2:num_choices+1,2:num_choices+1] = pair_num_p;
pair_num_avg_1aoff[2:num_choices+1,2:num_choices+1] = pair_num_avg;
for choice_n = 1:num_choices 
    pair_num_avg[choice_n, 1:choice_n-1] .= NaN; # remove the lower triangular portion for plotting
    pair_num_h[choice_n, 1:choice_n-1] .= NaN;
    pair_num_p[choice_n, 1:choice_n-1] .= NaN;
    pair_num_avg_1aoff[choice_n, 1:choice_n-1] .= NaN; # remove the lower triangular portion for plotting
    pair_num_h_1aoff[choice_n, 1:choice_n-1] .= NaN;
    pair_num_p_1aoff[choice_n, 1:choice_n-1] .= NaN;
end
pair_num_avg_1aoff[num_choices+1, 1:num_choices] .= NaN; # remove the lower triangular portion for plotting
pair_num_h_1aoff[num_choices+1, 1:num_choices] .= NaN;
pair_num_p_1aoff[num_choices+1, 1:num_choices] .= NaN;


choice_names_plot = ["1a off","1a on","1b on","1c on","1d on","1e on","1f on","1g on","2a on","2b on","2c on","3a on","3b on","3c on","3d on","3e on"];

# for events X: "feature X is turned on" and A: "model is accepted as a good model (by healthy, patient or average model loss)"
PrX = pair_num_tot_1aoff/num_models; # proportion of all models that have feature X turned on   Pr(X),  independent of model fitting/dataset used 
PrXA_h = pair_num_h_1aoff/top_model_num_h; # proportion of the accepted models for the healthy dataset that have feature X turned on   Pr(X|A)
PrA_h = top_model_num_h/num_models; # proportion of models that were accepted out of all models under criteria involving healthy model loss   Pr(A)
PrXA_p = pair_num_p_1aoff/top_model_num_p; # Pr(X|A) for patient models 
PrA_p = top_model_num_p/num_models; # Pr(A) for patient models
PrXA_avg = pair_num_avg_1aoff/top_model_num; # Pr(X|A) for models ordered by average loss 
PrA_avg = top_model_num/num_models; # Pr(A) for models ordered by average loss

Q_h = (PrXA_h - PrX)./(PrXA_h - 2*PrXA_h.*(PrX .+ PrA_h - PrA_h.*PrXA_h) + PrX);#(OR_h.-1)./(OR_h.+1); # Yule's Q (coefficient of association) for features in healthy models
Q_p = (PrXA_p - PrX)./(PrXA_p - 2*PrXA_p.*(PrX .+ PrA_p - PrA_p.*PrXA_p) + PrX)#(OR_p.-1)./(OR_p.+1); # Yule's Q for features in patient models
Q_avg = (PrXA_avg - PrX)./(PrXA_avg - 2*PrXA_avg.*(PrX .+ PrA_avg - PrA_avg.*PrXA_avg) + PrX)#(OR_avg.-1)./(OR_avg.+1); # Yule's Q for features in models ordered by average loss


# Fig 5b: Yule's Q for each combination of model choices in the best models
plt = CairoMakie.Figure(fontsize=32, size=(800,600)); 
ax = CairoMakie.Axis(plt[1, 1], xlabel="Model feature selection", ylabel="Model feature selection", xticks = (1:length(choice_names_plot), choice_names_plot), yticks = (1:length(choice_names_plot), choice_names_plot), xticklabelrotation=pi/2);
hmap = CairoMakie.heatmap!(1:length(choice_names_plot), 1:length(choice_names_plot), Q_avg; colormap=cgrad(:vik, rev=true), colorrange=(-1,1))
CairoMakie.Colorbar(plt[1, 2], hmap; label="Yule's Q", width=15, ticksize=5, tickalign=1);
CairoMakie.colsize!(plt.layout, 1, Aspect(1, 1.0));
CairoMakie.colgap!(plt.layout, 10)
display(plt)
save("figs_and_videos\\Fig_5b.png", plt)



lit_support = [2 8 7 3 17; 14 2 12 6 5; 1 2 5 26 0; 2 0 7 9 20; 3 5 6 0 0; 4 2 3 2 5; 6 9 1 7 2;
               8 2 1 9 0; 1 6 1 10 1; 10 4 2 6 0; 
               3 4 2 7 13; 4 3 1 14 3; 10 3 2 7 2; 7 13 2 0 0; 7 1 2 2 0]; # number of articles from Table 1 in main text, [feature (in order), likely/implicitly likely/neutral/implicitly unlikely/unlikely]

lit_support_ratio = zeros(num_choices+1); # ratio between -1 (against) and 1 (for) for each feature selection including 1a off,  calculated as (1*likely + 0.5*implicitly likely + 0*neutral - 0.5*implicitly unlikely - 1*unlikely)/(likely + implicitly likely + neutral + implicitly unlikely + unlikely)
lit_support_ratio[2:end] = lit_support*[1;0.5;0;-0.5;-1]./sum(lit_support,dims=2)[:]; # add in features 1a on onwards
lit_support_ratio[1] = -lit_support_ratio[2]; # likeliness for 1a off

# Fig 5c: Yule's Q for individual features compared with the weighted average of values from Table 1
plt = CairoMakie.Figure(fontsize=26.5, size=(600,400)); ax = CairoMakie.Axis(plt[1, 1], xticks = (1:16, choice_names_plot), xticklabelrotation=pi/2, xlabel="Model feature selection", ylabel="Yule's Q / Weighted average", limits=(0,17,-1.1,1.1));
plot_cols = [(:gray,0.9), (col_nai,0.75)]; # colours for barplot
CairoMakie.barplot!([i for i in 1:num_choices+1 for j in 1:2], [lit_support_ratio Q_avg[1:length(choice_names_plot)+1:length(choice_names_plot)^2]]'[:], color=plot_cols[[j for i in 1:num_choices+1 for j in 1:2]])
CairoMakie.barplot!([0], color=plot_cols[1], label="Literature")
CairoMakie.barplot!([0], color=plot_cols[2], label="Top models")
CairoMakie.axislegend(halign=:center, valign=:top);
display(plt)
save("figs_and_videos\\Fig_5c.png", plt)



feat_ignore = ["1b on","2a on","2b on","2c on","3d on","3e on"]; # list of feature selections to ignore (leave empty to plot all)

tot_models = 2000;#top_model_num # total number of models (ordered from best to worst)
sweep_num = 500; # number of models to "sweep" over when checking proportion of models that contain a model choice 

choice_freq_sweep = zeros(num_choices+1, tot_models-sweep_num); # frequency of model choices as we sweep through the top tot_models models
choice_on = [false for i in 1:num_choices+1]; # true if the feature is on for any models in the top tot_models models
avg_avg_loss = zeros(tot_models-sweep_num); # average of the average normalised model loss for all models in each sweep
for start_ind = 1:tot_models-sweep_num # for each starting index (start of the current sweep)
    for model_n = err_order_avg[start_ind:start_ind+sweep_num-1] # for each model index in the current sweep
        choice_freq_sweep[2:end, start_ind] += models[model_n]; # increment the number of times each model feature was used for this sweep
        if !models[model_n][1] && !(model_n in ignore_1a)# if 1a is truly off (not excluding this model)
            choice_freq_sweep[1, start_ind] += 1; # increment number of 1a off
            choice_on[1] = true; # 1a is off in at least one model
        end
        for choice_n = 1:num_choices # for each model feature 
            if models[model_n][choice_n] # if it is on in the current model 
                choice_on[choice_n+1] = true; # feature is on in at least one model
            end
        end
        avg_avg_loss[start_ind] += model_err_avg[model_n]; # add error from the current model
    end
end
avg_avg_loss = avg_avg_loss/sweep_num; # divide by number of entries

PrX = pair_num_tot_1aoff[1:17:16^2]/num_models; # proportion of all models that have feature X turned on   Pr(X)
PrXA_sweep = choice_freq_sweep/sweep_num; # proportion of the accepted models that have feature X turned on   Pr(X|A)
PrA_sweep = sweep_num/num_models; # proportion of models that were accepted out of all models   Pr(A)

# Fig 5d: Yule's Q for models with increasing model loss
loss_x_axis = true; # true to make the x axis the model loss, otherwise the model ID (of the central model in each window)
plt = CairoMakie.Figure(fontsize=32, size=(750,500)); 
if loss_x_axis
    ax = CairoMakie.Axis(plt[1, 1], xlabel="Average normalised model loss", ylabel="Yule's Q");
else
    ax = CairoMakie.Axis(plt[1, 1], xlabel="Models ordered by increasing loss", ylabel="Yule's Q");
end
plot_cols = cgrad(:jet, num_choices+1-length(feat_ignore), categorical = true);
for choice_n_plot = 1:num_choices+1-length(feat_ignore) # for each feature that we are not ignoring
    choice_n = findall(x->!(x in feat_ignore),choice_names_plot)[choice_n_plot];
    if choice_on[choice_n] # if at least one model has the current feature
        Q_sweep = (PrXA_sweep[choice_n,:] .- PrX[choice_n])./(PrXA_sweep[choice_n,:] .- 2*PrXA_sweep[choice_n,:].*(PrX[choice_n] .+ PrA_sweep .- PrA_sweep.*PrXA_sweep[choice_n,:]) .+ PrX[choice_n]); # Yule's Q (coefficient of association) for features in models
        if loss_x_axis
            CairoMakie.lines!(avg_avg_loss, Q_sweep, label="$(choice_names_plot[choice_n])", color=plot_cols[choice_n_plot], linewidth=3);#model_err_avg[err_order_avg[(1:tot_models-sweep_num).+sweep_num÷2]] # plot proportion of models that contain the current model feature for all sweeps
        else
            CairoMakie.lines!((1:tot_models-sweep_num).+sweep_num÷2, Q_sweep, label="$(choice_names_plot[choice_n])", color=plot_cols[choice_n_plot], linewidth=3); # plot proportion of models that contain the current model feature for all sweeps
        end
    end
end
leg = CairoMakie.axislegend(); plt[1, 2] = leg; # make legend and place it next to plot
display(plt)
save("figs_and_videos\\Fig_5d.png", plt)




donor_cols = palette(:viridis, 2); # colours for healthy and patient model values 
param_scales = [24*60, 24*60, 24*60, 1/ρ_avg_base, 1/ρ_avg_base, 1/ρ_avg_base, 24*60]; # scaling for parameter values 
param_names_short = ["r_m", "r_e", "r_d", "s_m", "s_e", "s_d", "r_p"]; # shorter versions of parameter names

# Fig 5e: boxplots for parameter values (in models_inds) in healthy and patient models
# plot for parameters with units 1/day
plt = CairoMakie.Figure(fontsize=28, size=(350,500)); ax = CairoMakie.Axis(plt[1, 1], xticks = (1.5:3:3*4-0.5, param_names_short[[1,2,3,7]]), xticklabelrotation=pi/4, xlabel="Parameter", ylabel="Parameter value (day⁻¹)", limits=(0,nothing,0,nothing));
for param_n_plot = 1:4
    param_n = [1,2,3,7][param_n_plot];
    CairoMakie.boxplot!((3*param_n_plot-2)*ones(top_model_num), [param_fits_h[i][param_n] for i in err_order_avg[1:top_model_num]]*param_scales[param_n]; whiskerwidth = 1, width = 1, color=donor_cols[1]);#, show_outliers=false
    CairoMakie.boxplot!((3*param_n_plot-1)*ones(top_model_num), [param_fits_p[i][param_n] for i in err_order_avg[1:top_model_num]]*param_scales[param_n]; whiskerwidth = 1, width = 1, color=donor_cols[2]);
end
CairoMakie.boxplot!([-1], [0]; whiskerwidth = 1, width = 1, color=donor_cols[1], label="Healthy");
CairoMakie.boxplot!([-1], [0]; whiskerwidth = 1, width = 1, color=donor_cols[2], label="Patient");
CairoMakie.axislegend(position = :rt)
display(plt)
save("figs_and_videos\\Fig_5e_left.png", plt)

# plot for parameters with units ng/μg
plt = CairoMakie.Figure(fontsize=28, size=(280,500)); ax = CairoMakie.Axis(plt[1, 1], xticks = (1.5:3:3*3-0.5, param_names_short[[4,5,6]]), xticklabelrotation=pi/4, xlabel="Parameter", ylabel="Parameter value (ng/μg)", limits=(0,nothing,0,nothing));
for param_n_plot = 1:3
    param_n = [4,5,6][param_n_plot];
    CairoMakie.boxplot!((3*param_n_plot-2)*ones(top_model_num), [param_fits_h[i][param_n] for i in err_order_avg[1:top_model_num]]*param_scales[param_n]; whiskerwidth = 1, width = 1, color=donor_cols[1]);#, show_outliers=false
    CairoMakie.boxplot!((3*param_n_plot-1)*ones(top_model_num), [param_fits_p[i][param_n] for i in err_order_avg[1:top_model_num]]*param_scales[param_n]; whiskerwidth = 1, width = 1, color=donor_cols[2]);
end
CairoMakie.boxplot!([-1], [0]; whiskerwidth = 1, width = 1, color=donor_cols[1], label="Healthy");
CairoMakie.boxplot!([-1], [0]; whiskerwidth = 1, width = 1, color=donor_cols[2], label="Patient");
display(plt)
save("figs_and_videos\\Fig_5e_right.png", plt)


# Fig 5f: boxplots for the differences in parameter values between healthy and patient with specific models selected out of the selected models
# plot for parameters with units 1/day
plt = CairoMakie.Figure(fontsize=28, size=(300,500)); ax = CairoMakie.Axis(plt[1, 1], xticks = (1:4, param_names_short[[1,2,3,7]]), xticklabelrotation=pi/4, xlabel="Parameter", ylabel="Parameter difference, h-p (day⁻¹)", limits=(nothing, nothing, -10, 10));
for param_n_plot = 1:4
    param_n = [1,2,3,7][param_n_plot];
    CairoMakie.boxplot!(param_n_plot*ones(top_model_num), ([param_fits_h[i][param_n] for i in err_order_avg[1:top_model_num]]-[param_fits_p[i][param_n] for i in err_order_avg[1:top_model_num]])*param_scales[param_n]; whiskerwidth = 1, width = 1, color=col_nai);#, show_outliers=false
end
display(plt)
save("figs_and_videos\\Fig_5f_left.png", plt)

# plot for parameters with units ng/μg
plt = CairoMakie.Figure(fontsize=28, size=(250,500)); ax = CairoMakie.Axis(plt[1, 1], xticks = (1:3, param_names_short[[4,5,6]]), xticklabelrotation=pi/4, xlabel="Parameter", ylabel="Parameter difference, h-p (ng/μg)", limits=(nothing, nothing, -20, 20));
for param_n_plot = 1:3
    param_n = [4,5,6][param_n_plot];
    CairoMakie.boxplot!(param_n_plot*ones(top_model_num), ([param_fits_h[i][param_n] for i in err_order_avg[1:top_model_num]]-[param_fits_p[i][param_n] for i in err_order_avg[1:top_model_num]])*param_scales[param_n]; whiskerwidth = 1, width = 1, color=col_nai);#, show_outliers=false
end
display(plt)
save("figs_and_videos\\Fig_5f_right.png", plt)




# plot results under a constraint on feature selections, for example, without bidirectional differentiation
model_assumptions(eff_to_mem,div_diff,bi_diff,no_diff,sig_stren,no_stim,back_diff,naive_to_exh,mem_to_exh,eff_to_exh,asymm_div,eff_prolif,naive_prolif,eff_death,exh_death) = !bi_diff; # condition to exclude models

choice_names_assump = ["1a off","1a on","1b on","1d on","1e on","1f on","1g on","2a on","2b on","2c on","3a on","3b on","3c on","3d on","3e on"];

model_inds = []; # array of indices of the top models that match the assumptions
for model_n = err_order_avg[1:top_model_num] # for each of the top models 
    eff_to_mem, div_diff, bi_diff, no_diff, sig_stren, no_stim, back_diff, naive_to_exh, mem_to_exh, eff_to_exh, asymm_div, eff_prolif, naive_prolif, eff_death, exh_death = models[model_n]; # extract model features
    if model_assumptions(eff_to_mem,div_diff,bi_diff,no_diff,sig_stren,no_stim,back_diff,naive_to_exh,mem_to_exh,eff_to_exh,asymm_div,eff_prolif,naive_prolif,eff_death,exh_death) # if the model satisfies the assumptions 
        append!(model_inds, [model_n]); # add this to the list of models that satisfy the assumption 
    end
end
#println("$(length(model_inds)) models satisfy the assumption")

pair_num_assump = zeros(num_choices, num_choices); # number of times a pair of model features is used in the top models that satisfy the assumptions
pair_num_assump_1aoff = zeros(num_choices+1, num_choices+1); # number of times a pair of model features (including 1a off) is used in the top models that satisfy the assumptions
for model_n in model_inds # for each of the top models
    curr_model = models[model_n]; # current model
    for choice_n = findall(x->x==1,curr_model) # for each individual model choice in the current model 
        pair_num_assump[choice_n,choice_n] += 1; # add to the number of times this model choice was used
    end
    for comb = combinations(findall(x->x==1,curr_model),2) # for each combination of model choices in the current model
        pair_num_assump[comb[1],comb[2]] += 1; # add to the number of times this pair of model choices is used 
    end
    if !curr_model[1] && !(model_n in ignore_1a) # if 1a is truly off (not overwritten)
        pair_num_assump_1aoff[1,1] += 1; # increment the number of times feature 1a is off
        for feat_n = 2:num_choices # for each other feature
            if curr_model[feat_n] # if that feature is on (alongside 1a off)
                pair_num_assump_1aoff[1,feat_n+1] += 1; # increment that pair
            end
        end
    end
end
pair_num_assump_1aoff[2:num_choices+1,2:num_choices+1] = pair_num_assump; # add in results for all feature selections other than 1a off
for choice_n = 1:num_choices 
    pair_num_assump[choice_n, 1:choice_n-1] .= NaN; # remove the lower triangular portion for plotting
    pair_num_assump_1aoff[choice_n, 1:choice_n-1] .= NaN;
end
pair_num_assump_1aoff[num_choices+1, 1:num_choices] .= NaN; # remove the lower triangular portion for plotting

pair_num_tot_assump = zeros(num_choices,num_choices); # number of times a pair of model choices appears in all models under the assumption
pair_num_tot_1aoff_assump = zeros(num_choices+1,num_choices+1); # number of times a pair of model choices appears in all models (also considering 1a off as a feature selection) under the assumption
num_models_assump = 0; # number of all models that satisfy the assumption
for model_n = eachindex(models) # for each model
    curr_model = models[model_n]; # current model
    eff_to_mem, div_diff, bi_diff, no_diff, sig_stren, no_stim, back_diff, naive_to_exh, mem_to_exh, eff_to_exh, asymm_div, eff_prolif, naive_prolif, eff_death, exh_death = curr_model; # extract model features
    if model_assumptions(eff_to_mem,div_diff,bi_diff,no_diff,sig_stren,no_stim,back_diff,naive_to_exh,mem_to_exh,eff_to_exh,asymm_div,eff_prolif,naive_prolif,eff_death,exh_death) # if the model satisfies the assumptions 
        for choice_n = findall(x->x==1,curr_model) # for each individual model choice in the current model 
            pair_num_tot_assump[choice_n,choice_n] += 1; # add to the number of times this model choice was used
        end
        for comb = combinations(findall(x->x==1,curr_model),2) # for each combination of model choices in the current model
            pair_num_tot_assump[comb[1],comb[2]] += 1; # add to the number of times this pair of model choices is used 
        end
        if !curr_model[1] && !(model_n in ignore_1a) # if 1a is truly off (not overwritten)
            pair_num_tot_1aoff_assump[1,1] += 1; # increment the number of times feature 1a is off
            for feat_n = 2:num_choices # for each other feature
                if curr_model[feat_n] # if that feature is on (alongside 1a off)
                    pair_num_tot_1aoff_assump[1,feat_n+1] += 1; # increment that pair
                end
            end
        end
        num_models_assump += 1; # increment number of models that satisfy the assumption
    end
end
pair_num_tot_1aoff_assump[2:num_choices+1,2:num_choices+1] = pair_num_tot_assump; # add in results for all feature selections other than 1a off

PrX = pair_num_tot_1aoff_assump[[1:3;5:end],[1:3;5:end]]/num_models_assump; # proportion of all models (with 1c off) that have feature X turned on   Pr(X),  independent of model fitting/dataset used 
PrXA_assump = pair_num_assump_1aoff[[1:3;5:end],[1:3;5:end]]/length(model_inds); # proportion of the accepted models under the assumption that have feature X turned on   Pr(X|A)
PrA_assump = length(model_inds)/num_models_assump; # proportion of models that were accepted out of all models (with 1c off)   Pr(A)
Q_assump = (PrXA_assump - PrX)./(PrXA_assump - 2*PrXA_assump.*(PrX .+ PrA_assump - PrA_assump.*PrXA_assump) + PrX);#(OR_assump.-1)./(OR_assump.+1); # Yule's Q (coefficient of association) for features in top models that satisfy the assumption


# Fig 6a: Yule's Q for pairs of feature selections in the best models that do not contain bidirectional differentiation
plt = CairoMakie.Figure(fontsize=32, size=(800,600)); 
ax = CairoMakie.Axis(plt[1, 1], xlabel="Model feature selection", ylabel="Model feature", xticks = (1:length(choice_names_assump), choice_names_assump), yticks = (1:length(choice_names_assump), choice_names_assump), xticklabelrotation=pi/2);
hmap = CairoMakie.heatmap!(1:length(choice_names_assump), 1:length(choice_names_assump), Q_assump; colormap=cgrad(:vik, rev=true), colorrange=(-1,1))
CairoMakie.Colorbar(plt[1, 2], hmap; label="Yule's Q", width=15, ticksize=5, tickalign=1);
CairoMakie.colsize!(plt.layout, 1, Aspect(1, 1.0));
CairoMakie.colgap!(plt.layout, 10)
display(plt)
save("figs_and_videos\\Fig_6a.png", plt)

# Fig 6b: Yule's Q for individual feature selections in the best models that do not contain bidirectional differentiation compared to the weighted average of values in Table 1
plt = CairoMakie.Figure(fontsize=25, size=(600,400)); ax = CairoMakie.Axis(plt[1, 1], xticks = (1:length(choice_names_assump), choice_names_assump), xticklabelrotation=pi/2, xlabel="Model feature selection", ylabel="Yule's Q / Weighted average", limits=(0,16,-1.1,1.1));
plot_cols = [(:gray,0.9), (col_nai,0.75)]; # colours for barplot
CairoMakie.barplot!([i for i in 1:length(choice_names_assump) for j in 1:2], [lit_support_ratio[[1:3;5:end]] Q_assump[1:length(choice_names_assump)+1:length(choice_names_assump)^2]]'[:], color=plot_cols[[j for i in 1:length(choice_names_assump) for j in 1:2]])
CairoMakie.barplot!([0], color=plot_cols[1], label="Literature")
CairoMakie.barplot!([0], color=plot_cols[2], label="Top models")
CairoMakie.axislegend(halign=:left, valign=:top);
display(plt)
save("figs_and_videos\\Fig_6b.png", plt)






## code for generating Figs 3 and 4, and also testing other pathways in literature

# select figure to generate: "Fig_3" or "Fig_4"
fig_name = "Fig_3";

if fig_name == "Fig_3"
    # define model feature selections to test from literature, with first being the base pathway and second being an altered pathway.  Feature selections are  [1a, 1b, 1c, 1d, 1e, 1f, 1g,  2a, 2b, 2c,  3a, 3b, 3c, 3d, 3e]
    # for Fig 3:
    model_tests = [[true, false, false, false, false, false, true,  false, false, true,  false, false, false, false, true],   [false, false, true, false, false, false, true,  false, false, true,  false, false, false, false, true]]; # base model from "T cell exhaustion", Wherry (2011),  and altered pathway is by switching bidirectional differentiation on
    solid_name = "1c off"; dash_name = "1c on"; # names for solid (base) and dashed (altered) models in legend
elseif fig_name == "Fig_4"
    # for Fig 4:
    model_tests = [[true, true, false, false, false, false, false,  false, true, false,  false, false, false, true, false],   [false, true, false, false, false, false, false,  false, true, false,  false, false, false, true, false]]; # base model from "T cells in health and disease", Sun et al. (2023),  and altered pathway is by switching the direction of linear differentiation
    solid_name = "1a on"; dash_name = "1a off"; # names for solid (base) and dashed (altered) models in legend
end

# alternatively define your own model: 
#model_tests = [[true, true, false, false, false, false, false,  false, true, false,  false, false, false, true, false],   [false, true, false, false, false, false, false,  false, true, false,  false, false, false, true, false]]; 
#solid_name = "1a on"; dash_name = "1a off"; # names for solid (base) and dashed (altered) models in legend

# Some example pathways (estimated from literature) and their general performance in this framework
#  "Defining ‘T cell exhaustion’", Blank et al. (2019);   [false, true, false, false, false, false, true (?),   true, false, false (?),   false, false, false, true, false]          #7001 (1g on, 2c off), #6983 (1g and 2c off), #1177 (1g and 2c on), #1016 (1g off, 2c on)
#  "T cell exhaustion", Wherry (2011);   [true, false, false, false, false, false, true (?),   false, false, true,   false, false, false, false, true]          #7494 (1g on), #7694 (1g off).     with 1g on: #1 on adding pathway from memory back to effector (1c on), #1268 on adding a pathway from memory to dysfunction-associated (2b on), #668 on adding asymmetric-like division (3a on)
#  "T cells in health and disease", Sun et al. (2023);   [true, true, false, false, false, false, false,   false, true, false,   false, false, false, true, false]          #6922.    #2200 on switching to memory->effector (1a off), #128 on adding cyclic differentiation (1c on), #964 on adding dysfunction-associated cell back-differentiation (1g on), #316 on adding asymmetric-like division (3a on)
#  "Epigenetic control of CD8+ T cell differentiation", Henning et al. (2018);   [false, false, false (?), false, true (?), false, true,   false, true, true,   true (?), false, false, true, true]          linear models: #5786 (asymmetric: 1c off, 1e on, 3a on), #4644 (symmetric: 1c off, 1e on, 3a off),  bidirectional models: #704 (asymmetric: 1c on, 1e off, 3a on), #121 (symmetric: 1c on, 1e off, 3a off)
#  "Effector-like and memory T-cell differentiation: implications for vaccine development", Kaech et al. (2002);   [true (?), false (?), false, false (?), false, true, false,   false, false, true,   false, false, false, true, true]          #7784 (1a and 1e on, 1b and 1d off), #7723 (1a and 1e off, 1b and 1d on)


T = 8*24*60; # final time for expansion (min)
s_range = range(0,15,2001); # values for stimulus ratio s
start_ind = 30; # starting index to avoid plotting numerically unstable values at the start

param_fits_both = [param_fits_h, param_fits_p]; # parameter fits for both healthy and patient datasets
healthy_IC = mean_init_num[[1,3,4,6],1];#[0, 0.1, 0.4, 0.5]*sum(mean_init_num[[1,3,4,6],1]); # initial numbers of cells for healthy donor data 
patient_IC = mean_init_num[[1,3,4,6],2];#[0, 0.1, 0.4, 0.5]*sum(mean_init_num[[1,3,4,6],1]); # initial numbers of cells for patient data
ICs = [healthy_IC, patient_IC]; # initial conditions for both healthy and patient-derived samples

model_outputs = zeros(2, length(s_range), 4, 2); # [model, s, cell, donor]:  cell count predictions for each cell in each model for each donor
for model_n = 1:2 # for both models being tested (i.e. a base pathway, 1, and an altered pathway, 2)
    curr_model = model_tests[model_n]; # get the current model (base or altered)
    eff_to_mem, div_diff, bi_diff, no_diff, sig_stren, no_stim, back_diff, naive_to_exh, mem_to_exh, eff_to_exh, asymm_div, eff_prolif, naive_prolif, eff_death, exh_death = curr_model; # extract model feature selections

    pathway_plt = plot_pathway(model_n, eff_to_mem, div_diff, bi_diff, no_diff, sig_stren, no_stim, back_diff, naive_to_exh, mem_to_exh, eff_to_exh, asymm_div, eff_prolif, naive_prolif, eff_death, exh_death); # plot current model pathway

    if curr_model in models
        model_ind = findfirst(x->x==curr_model,models); # get index for the model in models
        println("Model $model_n is #$(findfirst(x->x==model_ind,err_order_h)) for healthy models (L=$(round(model_err_h[model_ind],digits=3))), #$(findfirst(x->x==model_ind,err_order_p)) for patient models (L=$(round(model_err_p[model_ind],digits=3))), and #$(findfirst(x->x==model_ind,err_order_avg)) on average (L=$(round(model_err_avg[model_ind],digits=3))).")
    else
        error("Model $model_n was not simulated.")
    end

    r_dm_fit(ρ, r_dm_max, sρ_dm) = (eff_to_mem || bi_diff || div_diff) && no_stim ? 1/2*r_dm_max : sig_stren ? r_dm_max*(1-ρ/(sρ_dm + ρ)) : r_dm_max*ρ/(sρ_dm + ρ); # if removing stimulus dependence, set differentiation rate to that at the average stimulus concentration for a stimulus ratio of 1 (mass action) or half its maximum (Michaelis-Menten), otherwise if using signal strength model, change differentiation rate to memory to decrease in ρ
    r_de_fit(ρ, r_de_max, sρ_de) = (!eff_to_mem && !bi_diff || div_diff) && no_stim ? 1/2*r_de_max : r_de_max*ρ/(sρ_de + ρ); # if removing stimulus dependence, set differentiation rate to that at the average stimulus concentration for a stimulus ratio of 1 (mass action) or half its maximum (Michaelis-Menten)
    r_te_fit(ρ, r_te_max, sρ_te) = r_te_max*ρ/(sρ_te + ρ); # rate of terminal exhaustion
    r_p_fit(r_p) = r_p; # proliferation rate

    f_n_m_fit(ρ, r_dm_max, sρ_dm) = !eff_to_mem && !bi_diff || div_diff ? r_dm_fit(ρ, r_dm_max, sρ_dm) : 0; # transition rate from naive to memory
    f_n_e_fit(ρ, r_de_max, sρ_de) = eff_to_mem || bi_diff || div_diff ? r_de_fit(ρ, r_de_max, sρ_de) : 0; # transition rate from naive to effector
    f_m_e_fit(ρ, r_de_max, sρ_de) = (!eff_to_mem || bi_diff) && !no_diff ? r_de_fit(ρ, r_de_max, sρ_de) : 0; # transition rate from memory to effector 
    f_e_m_fit(ρ, r_dm_max, sρ_dm) = (eff_to_mem || bi_diff) && !no_diff ? r_dm_fit(ρ, r_dm_max, sρ_dm) : 0; # transition rate from effector to memory
    f_n_d_fit(ρ, r_te_max, sρ_te) = naive_to_exh ? r_te_fit(ρ, r_te_max, sρ_te) : 0; # transition rate from naive to terminally exhausted
    f_m_d_fit(ρ, r_te_max, sρ_te) = mem_to_exh ? r_te_fit(ρ, r_te_max, sρ_te) : 0; # transition rate from memory to terminally exhausted 
    f_e_d_fit(ρ, r_te_max, sρ_te) = eff_to_exh ? r_te_fit(ρ, r_te_max, sρ_te) : 0; # transition rate from effector to terminally exhausted 
    f_d_e_fit(ρ, r_de_max, sρ_de) = back_diff ? r_de_fit(ρ, r_de_max, sρ_de) : 0; # transition rate from exhausted to effector 
    f_n_n_prolif_fit(r_p) = naive_prolif ? r_p_fit(r_p) : 0; # symmetric naive cell proliferation rate
    f_m_n_prolif_fit(r_p) = asymm_div ? r_p_fit(r_p) : 0; # asymmetric memory cell proliferation rate
    f_m_m_prolif_fit(r_p) = asymm_div ? 0 : r_p_fit(r_p); # symmetric memory cell proliferation rate
    f_e_n_prolif_fit(r_p) = asymm_div && eff_prolif ? r_p_fit(r_p) : 0; # asymmetric effector cell proliferation rate
    f_e_e_prolif_fit(r_p) = !asymm_div && eff_prolif ? r_p_fit(r_p) : 0; # symmetric effector cell proliferation rate
    f_e_death_fit = eff_death ? r_d : 0; # effector cell death rate
    f_d_death_fit = exh_death ? r_d : 0; # terminally exhausted cell death rate

    R_nn_fit(ρ, r_dm_max, r_de_max, r_te_max, sρ_dm, sρ_de, sρ_te, r_p) = -f_n_m_fit(ρ, r_dm_max, sρ_dm) - f_n_e_fit(ρ, r_de_max, sρ_de) - f_n_d_fit(ρ, r_te_max, sρ_te) + f_n_n_prolif_fit(r_p); # stimulus-dependent reaction term driven by naive cells for n (naive cells)
    R_nm_fit(r_p) = f_m_n_prolif_fit(r_p); # reaction term driven by memory cells for n
    R_ne_fit(r_p) = f_e_n_prolif_fit(r_p); # reaction term driven by effector cells for n
    R_mn_fit(ρ, r_dm_max, sρ_dm) = f_n_m_fit(ρ, r_dm_max, sρ_dm); # stimulus-dependent reaction term driven by naive cells for m (memory cells)
    R_mm_fit(ρ, r_de_max, r_te_max, sρ_de, sρ_te, r_p) = -f_m_e_fit(ρ, r_de_max, sρ_de) - f_m_d_fit(ρ, r_te_max, sρ_te) + f_m_m_prolif_fit(r_p); # stimulus-dependent reaction term driven by memory cells for m
    R_me_fit(ρ, r_dm_max, sρ_dm) = f_e_m_fit(ρ, r_dm_max, sρ_dm); # stimulus-dependent reaction term driven by effector cells for m
    R_en_fit(ρ, r_de_max, sρ_de) = f_n_e_fit(ρ, r_de_max, sρ_de); # stimulus-dependent reaction term driven by naive cells for e (effector cells)
    R_em_fit(ρ, r_de_max, sρ_de) = f_m_e_fit(ρ, r_de_max, sρ_de); # stimulus-dependent reaction term driven by memory cells for e
    R_ee_fit(ρ, r_dm_max, r_te_max, sρ_dm, sρ_te, r_p) = -f_e_m_fit(ρ, r_dm_max, sρ_dm) - f_e_d_fit(ρ, r_te_max, sρ_te) + f_e_e_prolif_fit(r_p) - f_e_death_fit; # stimulus-dependent reaction term driven by effector cells for e 
    R_ed_fit(ρ, r_de_max, sρ_de) = f_d_e_fit(ρ, r_de_max, sρ_de); # stimulus-dependent reaction term driven by exhausted cells for e
    R_dn_fit(ρ, r_te_max, sρ_te) = f_n_d_fit(ρ, r_te_max, sρ_te); # stimulus-dependent reaction term driven by naive cells for d (terminally exhausted/dysfunction-associated cells)
    R_dm_fit(ρ, r_te_max, sρ_te) = f_m_d_fit(ρ, r_te_max, sρ_te); # stimulus-dependent reaction term driven by memory cells for d
    R_de_fit(ρ, r_te_max, sρ_te) = f_e_d_fit(ρ, r_te_max, sρ_te); # stimulus-dependent reaction term driven by effector cells for d
    R_dd_fit(ρ, r_de_max, sρ_de) = -f_d_e_fit(ρ, r_de_max, sρ_de) - f_d_death_fit; # reaction term driven by exhausted cells for d

    ODE_funcs = [R_nn_fit, R_nm_fit, R_ne_fit, R_mn_fit, R_mm_fit, R_me_fit, R_en_fit, R_em_fit, R_ee_fit, R_ed_fit, R_dn_fit, R_dm_fit, R_de_fit, R_dd_fit]; # initialise functions for ODE system

    for donor_n = 1:2 # for each donor
        r_dm_max, r_de_max, r_te_max, sρ_dm, sρ_de, sρ_te, r_p_max = param_fits_both[donor_n][model_ind];

        A_s(s) = [ODE_funcs[1](ρ_avg_base*s,r_dm_max,r_de_max,r_te_max,sρ_dm,sρ_de,sρ_te,r_p_max) ODE_funcs[2](r_p_max) ODE_funcs[3](r_p_max) 0; 
            ODE_funcs[4](ρ_avg_base*s,r_dm_max,sρ_dm) ODE_funcs[5](ρ_avg_base*s,r_de_max,r_te_max,sρ_de,sρ_te,r_p_max) ODE_funcs[6](ρ_avg_base*s,r_dm_max,sρ_dm) 0; 
            ODE_funcs[7](ρ_avg_base*s,r_de_max,sρ_de) ODE_funcs[8](ρ_avg_base*s,r_de_max,sρ_de) ODE_funcs[9](ρ_avg_base*s,r_dm_max,r_te_max,sρ_dm,sρ_te,r_p_max) ODE_funcs[10](ρ_avg_base*s,r_de_max,sρ_de); 
            ODE_funcs[11](ρ_avg_base*s,r_te_max,sρ_te) ODE_funcs[12](ρ_avg_base*s,r_te_max,sρ_te) ODE_funcs[13](ρ_avg_base*s,r_te_max,sρ_te) ODE_funcs[14](ρ_avg_base*s,r_de_max,sρ_de)]; # coefficients in ODE system as a function of stimulus ration s

        for s_n = eachindex(s_range) # for each stimulus ratio
            model_outputs[model_n,s_n,:,donor_n] = exp(A_s(s_range[s_n])*T)*ICs[donor_n]; # solution at time T for the donor_nth donor's initial condition
        end # end for each stimulus ratio
    end # end for each donor
end

# plot results for cell type cell_n and donor donor_n 
plot_x = [stim for i in 1:10 for stim in stim_vals]; # arrange data to plot on x axis of boxplots
cell_names = ["Naive", "Memory-like", "Effector-like", "Dysfunction-associated"]; # cell names for plotting
cell_colours = [col_nai, col_mem, col_eff, col_exh]; # colours for cell types
donor_names = ["Healthy", "Patient"]; # donor names for plotting

for cell_n = 1:4 # for each cell type;   1: naive,  2: memory,  3: effector,  4: exhausted 
    for donor_n = 1:2 # for each donor type;   1: healthy,  2: patient 
        nonNaN_inds = []; # indices of non-NaN datapoints for the current cell type
        for i = eachindex(plot_x) # for each datapoint
            if !isnan(final_num[1,:,donor_n,:][i])
                append!(nonNaN_inds, i); # if index has a non-NaN value, add to the list of indices
            end
        end

        if cell_n == 4 && donor_n == 2
            plt = CairoMakie.Figure(fontsize=25, size=(590,300)); 
        else
            plt = CairoMakie.Figure(fontsize=25, size=(600,300)); 
        end
        ax = CairoMakie.Axis(plt[1, 1], xlabel="Stimulus density (ng/μg)", ylabel="$(cell_names[cell_n])\ncells (millions)", limits = (0, maximum(s_range), -0.15, nothing))#, title="$(donor_names[donor_n]) cell models");
        if T==8*24*60; CairoMakie.boxplot!(plot_x[nonNaN_inds], final_num[[1,3,4,6][cell_n],:,donor_n,:][nonNaN_inds]/1e6; color=cell_colours[cell_n]/0.8, whiskerwidth = 1, width = 1); end
        CairoMakie.lines!(s_range[start_ind:end], model_outputs[1,start_ind:end,cell_n,donor_n]/1e6, color=cell_colours[cell_n]*0.8, linewidth=5, label=solid_name);
        CairoMakie.lines!(s_range[start_ind:end], model_outputs[2,start_ind:end,cell_n,donor_n]/1e6, color=cell_colours[cell_n]*0.8, linewidth=5, linestyle=:dash, label=dash_name);
        if cell_n == 2 && donor_n == 2
            CairoMakie.axislegend(labelsize=20, patchsize = (40, 10), position = :lt);
        end
        display(plt)
        if fig_name == "Fig_2" || fig_name == "Fig_3"
            save("figs_and_videos\\$(fig_name)_$(cell_names[cell_n])_$(donor_names[donor_n]).png", plt)
        end
    end
end



# code to investigate pathways proposed in literature that do not include exhausted cell compartments (not shown in paper)

# specify model features relating to naive, memory and effector cells only,  below approximates the pathway found by "Disparate Individual Fates Compose Robust CD8+ T Cell Immunity", Buchholz (2013)
eff_to_mem = false; # (1a) switch linear differentiation from effector to memory
div_diff = false; # (1b) divergent differentiation from naive to memory and effector
bi_diff = false; # (1c) naive -> effector then cyclic from effector <-> memory
no_diff = false; # (1d) no differentiation between memory and effector (makes sense for a divergent differentiation model)
sig_stren = false; # (1e) alter differentiation rates such that memory cells are produced more by less antigen stimulation
no_stim = false; # (1f) turn off antigen-dependent differentiation after the first differentiation event
asymm_div = false; # (3a) proliferating cells produce a new naive cell
eff_prolif = false; # (3b) effector cells proliferate
naive_prolif = false; # (3c) naive cells proliferate
eff_death = false; # (3d) effector cells die

#pathway_plt = plot_pathway(0, eff_to_mem, div_diff, bi_diff, no_diff, sig_stren, no_stim, false, false, false, false, asymm_div, eff_prolif, naive_prolif, eff_death, false); # plot current model pathway

back_diff_options = [false, true]; # (1g) back-differentiation from exhausted to effector
naive_to_exh_options = [false, true]; # (2a) allow naive cells to become terminally exhausted 
mem_to_exh_options = [false, true]; # (2b) allow memory cells to become terminally exhausted 
eff_to_exh_options = [false, true]; # (2c) effector cells become terminally exhausted
exh_death_options = [false, true]; # (3e) terminally exhausted cells die

best_model = [true for i in 1:15]; # initialise best model out of those involving dysfunction-associated cells (this initialised model is not a simulated model, due to exclusions)
best_placement = Inf; # initialise the best model placement
for back_diff = back_diff_options
    for naive_to_exh = naive_to_exh_options
        for mem_to_exh = mem_to_exh_options
            for eff_to_exh = eff_to_exh_options
                for exh_death = exh_death_options # for each option involving dysfunction-associated cells 
                    model_curr = [eff_to_mem, div_diff, bi_diff, no_diff, sig_stren, no_stim, back_diff, naive_to_exh, mem_to_exh, eff_to_exh, asymm_div, eff_prolif, naive_prolif, eff_death, exh_death]; # store features for the current model in array
                    if model_curr in models
                        model_curr_ind = findfirst(x->x==model_curr,models); # index of current model in models
                        if findfirst(x->x==model_curr_ind,err_order_avg) < best_placement # if the current model places higher than the previous best
                            best_placement = findfirst(x->x==model_curr_ind,err_order_avg); # replace best placement 
                            best_model = model_curr; # replace the best model
                        end
                    end
                end
            end
        end
    end
end

best_model_ind = findfirst(x->x==best_model,models); # index for best model

if best_model != [true for i in 1:15] # if the best model is not still at its initialised value
    eff_to_mem, div_diff, bi_diff, no_diff, sig_stren, no_stim, back_diff, naive_to_exh, mem_to_exh, eff_to_exh, asymm_div, eff_prolif, naive_prolif, eff_death, exh_death = best_model; # extract features of best model
    #plot_pathway(1, eff_to_mem, div_diff, bi_diff, no_diff, sig_stren, no_stim, back_diff, naive_to_exh, mem_to_exh, eff_to_exh, asymm_div, eff_prolif, naive_prolif, eff_death, exh_death); # plot best model pathway
    println("Model is #$(findfirst(x->x==best_model_ind,err_order_h)) for healthy models, #$(findfirst(x->x==best_model_ind,err_order_p)) for patient models, and #$(findfirst(x->x==best_model_ind,err_order_avg)) on average.")
else
    println("No models with this arrangement of naive, memory and effector cells were simulated.")
end








## code to generate SI Figs (need to have run the code to generate main text figures first)

# S1 Fig: simulating the best models past the stimulus densities used in experiments
T = 8*24*60; # final time for expansion (min)
s_range = range(0,40,2001); # values for stimulus ratio s
start_ind = 20; # starting index to avoid asymptote at the start

param_fits_both = [param_fits_h, param_fits_p]; # parameter fits for both healthy and patient datasets
healthy_IC = mean_init_num[[1,3,4,6],1];#[0, 0.1, 0.4, 0.5]*sum(mean_init_num[[1,3,4,6],1]); # initial numbers of cells for healthy donor data 
patient_IC = mean_init_num[[1,3,4,6],2];#[0, 0.1, 0.4, 0.5]*sum(mean_init_num[[1,3,4,6],1]); # initial numbers of cells for patient data
ICs = [healthy_IC, patient_IC]; # initial conditions for both healthy and patient-derived samples

model_outputs = zeros(top_model_num, length(s_range), 4, 2); # [model, s, cell, donor]:  cell count predictions for each cell in each model for each donor
@showprogress 1 "Simulating best models..." for top_model_n = 1:top_model_num # for each of the top models
    model_n = err_order_avg[top_model_n]; # get the model number in original arrays

    eff_to_mem, div_diff, cyc_diff, no_diff, sig_stren, no_stim, back_diff, naive_to_exh, mem_to_exh, eff_to_exh, asymm_div, eff_prolif, naive_prolif, eff_death, exh_death = models[model_n]; 

    r_dm_fit(ρ, r_dm_max, sρ_dm) = (eff_to_mem || cyc_diff || div_diff) && no_stim ? 1/2*r_dm_max : sig_stren ? r_dm_max*(1-ρ/(sρ_dm + ρ)) : r_dm_max*ρ/(sρ_dm + ρ); # if removing stimulus dependence, set differentiation rate to that at the average stimulus concentration for a stimulus ratio of 1 (mass action) or half its maximum (Michaelis-Menten), otherwise if using signal strength model, change differentiation rate to memory to decrease in ρ
    r_de_fit(ρ, r_de_max, sρ_de) = (!eff_to_mem && !cyc_diff || div_diff) && no_stim ? 1/2*r_de_max : r_de_max*ρ/(sρ_de + ρ); # if removing stimulus dependence, set differentiation rate to that at the average stimulus concentration for a stimulus ratio of 1 (mass action) or half its maximum (Michaelis-Menten)
    r_te_fit(ρ, r_te_max, sρ_te) = r_te_max*ρ/(sρ_te + ρ); # rate of terminal exhaustion
    r_p_fit(r_p) = r_p; # proliferation rate

    f_n_m_fit(ρ, r_dm_max, sρ_dm) = !eff_to_mem && !cyc_diff || div_diff ? r_dm_fit(ρ, r_dm_max, sρ_dm) : 0; # transition rate from naive to memory
    f_n_e_fit(ρ, r_de_max, sρ_de) = eff_to_mem || cyc_diff || div_diff ? r_de_fit(ρ, r_de_max, sρ_de) : 0; # transition rate from naive to effector
    f_m_e_fit(ρ, r_de_max, sρ_de) = (!eff_to_mem || cyc_diff) && !no_diff ? r_de_fit(ρ, r_de_max, sρ_de) : 0; # transition rate from memory to effector 
    f_e_m_fit(ρ, r_dm_max, sρ_dm) = (eff_to_mem || cyc_diff) && !no_diff ? r_dm_fit(ρ, r_dm_max, sρ_dm) : 0; # transition rate from effector to memory
    f_n_d_fit(ρ, r_te_max, sρ_te) = naive_to_exh ? r_te_fit(ρ, r_te_max, sρ_te) : 0; # transition rate from naive to terminally exhausted
    f_m_d_fit(ρ, r_te_max, sρ_te) = mem_to_exh ? r_te_fit(ρ, r_te_max, sρ_te) : 0; # transition rate from memory to terminally exhausted 
    f_e_d_fit(ρ, r_te_max, sρ_te) = eff_to_exh ? r_te_fit(ρ, r_te_max, sρ_te) : 0; # transition rate from effector to terminally exhausted 
    #=NEW=#f_d_e_fit(ρ, r_de_max, sρ_de) = back_diff ? r_de_fit(ρ, r_de_max, sρ_de) : 0; # transition rate from exhausted to effector 
    f_n_n_prolif_fit(r_p) = naive_prolif ? r_p_fit(r_p) : 0; # symmetric naive cell proliferation rate
    f_m_n_prolif_fit(r_p) = asymm_div ? r_p_fit(r_p) : 0; # asymmetric memory cell proliferation rate
    f_m_m_prolif_fit(r_p) = asymm_div ? 0 : r_p_fit(r_p); # symmetric memory cell proliferation rate
    f_e_n_prolif_fit(r_p) = asymm_div && eff_prolif ? r_p_fit(r_p) : 0; # asymmetric effector cell proliferation rate
    f_e_e_prolif_fit(r_p) = !asymm_div && eff_prolif ? r_p_fit(r_p) : 0; # symmetric effector cell proliferation rate
    f_e_death_fit = eff_death ? r_d : 0; # effector cell death rate
    f_d_death_fit = exh_death ? r_d : 0; # terminally exhausted cell death rate

    R_nn_fit(ρ, r_dm_max, r_de_max, r_te_max, sρ_dm, sρ_de, sρ_te, r_p) = -f_n_m_fit(ρ, r_dm_max, sρ_dm) - f_n_e_fit(ρ, r_de_max, sρ_de) - f_n_d_fit(ρ, r_te_max, sρ_te) + f_n_n_prolif_fit(r_p); # stimulus-dependent (and potentially cytokine-dependent) reaction term driven by naive cells for n (naive cells)
    R_nm_fit(r_p) = f_m_n_prolif_fit(r_p); # cytokine-dependent reaction term driven by memory cells for n
    R_ne_fit(r_p) = f_e_n_prolif_fit(r_p); # cytokine-dependent reaction term driven by effector cells for n
    R_mn_fit(ρ, r_dm_max, sρ_dm) = f_n_m_fit(ρ, r_dm_max, sρ_dm); # stimulus-dependent reaction term driven by naive cells for m (memory cells)
    R_mm_fit(ρ, r_de_max, r_te_max, sρ_de, sρ_te, r_p) = -f_m_e_fit(ρ, r_de_max, sρ_de) - f_m_d_fit(ρ, r_te_max, sρ_te) + f_m_m_prolif_fit(r_p); # stimulus- and cytokine-dependent reaction term driven by memory cells for m
    R_me_fit(ρ, r_dm_max, sρ_dm) = f_e_m_fit(ρ, r_dm_max, sρ_dm); # stimulus-dependent reaction term driven by effector cells for m
    R_en_fit(ρ, r_de_max, sρ_de) = f_n_e_fit(ρ, r_de_max, sρ_de); # stimulus-dependent reaction term driven by naive cells for e (effector cells)
    R_em_fit(ρ, r_de_max, sρ_de) = f_m_e_fit(ρ, r_de_max, sρ_de); # stimulus-dependent reaction term driven by memory cells for e
    R_ee_fit(ρ, r_dm_max, r_te_max, sρ_dm, sρ_te, r_p) = -f_e_m_fit(ρ, r_dm_max, sρ_dm) - f_e_d_fit(ρ, r_te_max, sρ_te) + f_e_e_prolif_fit(r_p) - f_e_death_fit; # stimulus- and cytokine-dependent reaction term driven by effector cells for e 
    #=NEW=#R_ed_fit(ρ, r_de_max, sρ_de) = f_d_e_fit(ρ, r_de_max, sρ_de); # stimulus-dependent reaction term driven by exhausted cells for e
    R_dn_fit(ρ, r_te_max, sρ_te) = f_n_d_fit(ρ, r_te_max, sρ_te); # stimulus- and cytokine-dependent reaction term driven by naive cells for d (terminally exhausted/dysfunction-associated cells)
    R_dm_fit(ρ, r_te_max, sρ_te) = f_m_d_fit(ρ, r_te_max, sρ_te); # stimulus- and cytokine-dependent reaction term driven by memory cells for d
    R_de_fit(ρ, r_te_max, sρ_te) = f_e_d_fit(ρ, r_te_max, sρ_te); # stimulus- and cytokine-dependent reaction term driven by effector cells for d
    #=UPDATED, WAS CONSTANT=#R_dd_fit(ρ, r_de_max, sρ_de) = -f_d_e_fit(ρ, r_de_max, sρ_de) - f_d_death_fit; # reaction term driven by exhausted cells for d

    ODE_funcs = [R_nn_fit, R_nm_fit, R_ne_fit, R_mn_fit, R_mm_fit, R_me_fit, R_en_fit, R_em_fit, R_ee_fit, R_ed_fit, R_dn_fit, R_dm_fit, R_de_fit, R_dd_fit]; # initialise functions for ODE system

    for donor_n = 1:2 # for each donor
        r_dm_max, r_de_max, r_te_max, sρ_dm, sρ_de, sρ_te, r_p_max = param_fits_both[donor_n][model_n];

        A_s(s) = [ODE_funcs[1](ρ_avg_base*s,r_dm_max,r_de_max,r_te_max,sρ_dm,sρ_de,sρ_te,r_p_max) ODE_funcs[2](r_p_max) ODE_funcs[3](r_p_max) 0; 
            ODE_funcs[4](ρ_avg_base*s,r_dm_max,sρ_dm) ODE_funcs[5](ρ_avg_base*s,r_de_max,r_te_max,sρ_de,sρ_te,r_p_max) ODE_funcs[6](ρ_avg_base*s,r_dm_max,sρ_dm) 0; 
            ODE_funcs[7](ρ_avg_base*s,r_de_max,sρ_de) ODE_funcs[8](ρ_avg_base*s,r_de_max,sρ_de) ODE_funcs[9](ρ_avg_base*s,r_dm_max,r_te_max,sρ_dm,sρ_te,r_p_max) ODE_funcs[10](ρ_avg_base*s,r_de_max,sρ_de); 
            ODE_funcs[11](ρ_avg_base*s,r_te_max,sρ_te) ODE_funcs[12](ρ_avg_base*s,r_te_max,sρ_te) ODE_funcs[13](ρ_avg_base*s,r_te_max,sρ_te) ODE_funcs[14](ρ_avg_base*s,r_de_max,sρ_de)]; # coefficients in ODE system as a function of stimulus ration s

        for s_n = eachindex(s_range) # for each stimulus ratio
            model_outputs[top_model_n,s_n,:,donor_n] = exp(A_s(s_range[s_n])*T)*ICs[donor_n]; # solution at time T for the donor_nth donor's initial condition
        end # end for each stimulus ratio
    end # end for each donor
end # end for each model

# plot results for cell type cell_n and donor donor_n 
plot_x = [stim for i in 1:10 for stim in stim_vals]; # arrange data to plot on x axis of boxplots

tot_cells_data = sum(final_num[[1,3,4,6],:,:,:], dims=1)[1,:,:,:]; # total cells in data 
tot_cells_sims = sum(model_outputs, dims=3)[:,:,1,:]; # total cells in simulation outputs

for cell_n = 1:4 # for each cell type;   1: naive,  2: memory,  3: effector,  4: exhausted 
    for donor_n = 1:2 # for each donor type;   1: healthy,  2: patient 
        nonNaN_inds = []; # indices of non-NaN datapoints for the current cell type
        for i = eachindex(plot_x) # for each datapoint
            if !isnan(final_num[1,:,donor_n,:][i])
                append!(nonNaN_inds, i); # if index has a non-NaN value, add to the list of indices
            end
        end

        plt = CairoMakie.Figure(fontsize=32); 
        ax = CairoMakie.Axis(plt[1, 1], xlabel="Stimulus density (ng/μg)", ylabel="$(cell_names[cell_n]) cells (millions)", limits = (0, maximum(s_range), 0, nothing))#, title="$(donor_names[donor_n]) cell models");
        for model_n = 1:top_model_num 
            CairoMakie.lines!(s_range[start_ind:end], model_outputs[model_n,start_ind:end,cell_n,donor_n]/1e6, color=(cell_colours[cell_n]*0.8,0.1));
        end
        if T==8*24*60; CairoMakie.boxplot!(plot_x[nonNaN_inds], final_num[[1,3,4,6][cell_n],:,donor_n,:][nonNaN_inds]/1e6; color=cell_colours[cell_n]/0.8, whiskerwidth = 1, width = 1); end
        display(plt)
        save("figs_and_videos\\Fig_S1_$(cell_names[cell_n])_$(donor_names[donor_n]).png", plt)
    end
end




# S2 Fig: Yule's Q for each combination of model choices
# healthy models (left)
plt = CairoMakie.Figure(fontsize=32, size=(800,600)); 
ax = CairoMakie.Axis(plt[1, 1], xlabel="Model feature selection", ylabel="Model feature selection", xticks = (1:length(choice_names_plot), choice_names_plot), yticks = (1:length(choice_names_plot), choice_names_plot), xticklabelrotation=pi/2);
hmap = CairoMakie.heatmap!(1:length(choice_names_plot), 1:length(choice_names_plot), Q_h; colormap=cgrad(:vik, rev=true), colorrange=(-1,1))
CairoMakie.Colorbar(plt[1, 2], hmap; label="Yule's Q", width=15, ticksize=5, tickalign=1);
CairoMakie.colsize!(plt.layout, 1, Aspect(1, 1.0));
CairoMakie.colgap!(plt.layout, 10)
display(plt)
save("figs_and_videos\\Fig_S2_left.png", plt)

# patient models (right)
plt = CairoMakie.Figure(fontsize=32, size=(800,600)); 
ax = CairoMakie.Axis(plt[1, 1], xlabel="Model feature selection", ylabel="Model feature selection", xticks = (1:length(choice_names_plot), choice_names_plot), yticks = (1:length(choice_names_plot), choice_names_plot), xticklabelrotation=pi/2);
hmap = CairoMakie.heatmap!(1:length(choice_names_plot), 1:length(choice_names_plot), Q_p; colormap=cgrad(:vik, rev=true), colorrange=(-1,1))
CairoMakie.Colorbar(plt[1, 2], hmap; label="Yule's Q", width=15, ticksize=5, tickalign=1);
CairoMakie.colsize!(plt.layout, 1, Aspect(1, 1.0));
CairoMakie.colgap!(plt.layout, 10)
display(plt)
save("figs_and_videos\\Fig_S2_right.png", plt)



# S3 Fig: Yule's Q for models with increasing model loss, all features shown
tot_models = 2000;#top_model_num # total number of models (ordered from best to worst)
sweep_num = 500; # number of models to "sweep" over when checking proportion of models that contain a model choice 

choice_freq_sweep = zeros(num_choices+1, tot_models-sweep_num); # frequency of model choices as we sweep through the top tot_models models
choice_on = [false for i in 1:num_choices+1]; # true if the feature is on for any models in the top tot_models models
avg_avg_loss = zeros(tot_models-sweep_num); # average of the average normalised model loss for all models in each sweep
for start_ind = 1:tot_models-sweep_num # for each starting index (start of the current sweep)
    for model_n = err_order_avg[start_ind:start_ind+sweep_num-1] # for each model index in the current sweep
        choice_freq_sweep[2:end, start_ind] += models[model_n]; # increment the number of times each model feature was used for this sweep
        if !models[model_n][1] && !(model_n in ignore_1a)# if 1a is truly off (not excluding this model)
            choice_freq_sweep[1, start_ind] += 1; # increment number of 1a off
            choice_on[1] = true; # 1a is off in at least one model
        end
        for choice_n = 1:num_choices # for each model feature 
            if models[model_n][choice_n] # if it is on in the current model 
                choice_on[choice_n+1] = true; # feature is on in at least one model
            end
        end
        avg_avg_loss[start_ind] += model_err_avg[model_n]; # add error from the current model
    end
end
avg_avg_loss = avg_avg_loss/sweep_num; # divide by number of entries

PrX = pair_num_tot_1aoff[1:17:16^2]/num_models; # proportion of all models that have feature X turned on   Pr(X)
PrXA_sweep = choice_freq_sweep/sweep_num; # proportion of the accepted models that have feature X turned on   Pr(X|A)
PrA_sweep = sweep_num/num_models; # proportion of models that were accepted out of all models   Pr(A)

loss_x_axis = true; # true to make the x axis the model loss, otherwise the model ID (of the central model in each window)
plt = CairoMakie.Figure(fontsize=28, size=(750,500)); 
if loss_x_axis
    ax = CairoMakie.Axis(plt[1, 1], xlabel="Average normalised model loss", ylabel="Yule's Q");
else
    ax = CairoMakie.Axis(plt[1, 1], xlabel="Models ordered by increasing loss", ylabel="Yule's Q");
end
plot_cols = cgrad(:jet, num_choices+1, categorical = true);
for choice_n = 1:num_choices+1 # for each feature that we are not ignoring
    if choice_on[choice_n] # if at least one model has the current feature
        Q_sweep = (PrXA_sweep[choice_n,:] .- PrX[choice_n])./(PrXA_sweep[choice_n,:] .- 2*PrXA_sweep[choice_n,:].*(PrX[choice_n] .+ PrA_sweep .- PrA_sweep.*PrXA_sweep[choice_n,:]) .+ PrX[choice_n]); # Yule's Q (coefficient of association) for features in models
        if loss_x_axis
            CairoMakie.lines!(avg_avg_loss, Q_sweep, label="$(choice_names_plot[choice_n])", color=plot_cols[choice_n], linewidth=3);#model_err_avg[err_order_avg[(1:tot_models-sweep_num).+sweep_num÷2]] # plot proportion of models that contain the current model feature for all sweeps
        else
            CairoMakie.lines!((1:tot_models-sweep_num).+sweep_num÷2, Q_sweep, label="$(choice_names_plot[choice_n])", color=plot_cols[choice_n], linewidth=3); # plot proportion of models that contain the current model feature for all sweeps
        end
    end
end
leg = CairoMakie.axislegend(labelsize = 23); plt[1, 2] = leg; # make legend and place it next to plot
display(plt)
save("figs_and_videos\\Fig_S3.png", plt)




# S4 Fig: simulating the worst 1000 models past the stimulus densities used in experiments
T = 8*24*60; # final time for expansion (min)
s_range = range(0,40,2001); # values for stimulus ratio s
start_ind = 20; # starting index to avoid asymptote at the start

param_fits_both = [param_fits_h, param_fits_p]; # parameter fits for both healthy and patient datasets
healthy_IC = mean_init_num[[1,3,4,6],1];#[0, 0.1, 0.4, 0.5]*sum(mean_init_num[[1,3,4,6],1]); # initial numbers of cells for healthy donor data 
patient_IC = mean_init_num[[1,3,4,6],2];#[0, 0.1, 0.4, 0.5]*sum(mean_init_num[[1,3,4,6],1]); # initial numbers of cells for patient data
ICs = [healthy_IC, patient_IC]; # initial conditions for both healthy and patient-derived samples

worst_model_num = 1000; # number of bad models to check
worst_model_inds = err_order_avg[num_models-1000+1:num_models]; # indices of the worst models in the original array

model_outputs = zeros(1000, length(s_range), 4, 2); # [model, s, cell, donor]:  cell count predictions for each cell in each model for each donor
@showprogress 1 "Simulating worst models..." for worst_model_n = 1:worst_model_num # for each of the worst models
    model_n = worst_model_inds[worst_model_n]; # get the model number in original array

    eff_to_mem, div_diff, cyc_diff, no_diff, sig_stren, no_stim, back_diff, naive_to_exh, mem_to_exh, eff_to_exh, asymm_div, eff_prolif, naive_prolif, eff_death, exh_death = models[model_n]; 

    r_dm_fit(ρ, r_dm_max, sρ_dm) = (eff_to_mem || cyc_diff || div_diff) && no_stim ? 1/2*r_dm_max : sig_stren ? r_dm_max*(1-ρ/(sρ_dm + ρ)) : r_dm_max*ρ/(sρ_dm + ρ); # if removing stimulus dependence, set differentiation rate to that at the average stimulus concentration for a stimulus ratio of 1 (mass action) or half its maximum (Michaelis-Menten), otherwise if using signal strength model, change differentiation rate to memory to decrease in ρ
    r_de_fit(ρ, r_de_max, sρ_de) = (!eff_to_mem && !cyc_diff || div_diff) && no_stim ? 1/2*r_de_max : r_de_max*ρ/(sρ_de + ρ); # if removing stimulus dependence, set differentiation rate to that at the average stimulus concentration for a stimulus ratio of 1 (mass action) or half its maximum (Michaelis-Menten)
    r_te_fit(ρ, r_te_max, sρ_te) = r_te_max*ρ/(sρ_te + ρ); # rate of terminal exhaustion
    r_p_fit(r_p) = r_p; # proliferation rate

    f_n_m_fit(ρ, r_dm_max, sρ_dm) = !eff_to_mem && !cyc_diff || div_diff ? r_dm_fit(ρ, r_dm_max, sρ_dm) : 0; # transition rate from naive to memory
    f_n_e_fit(ρ, r_de_max, sρ_de) = eff_to_mem || cyc_diff || div_diff ? r_de_fit(ρ, r_de_max, sρ_de) : 0; # transition rate from naive to effector
    f_m_e_fit(ρ, r_de_max, sρ_de) = (!eff_to_mem || cyc_diff) && !no_diff ? r_de_fit(ρ, r_de_max, sρ_de) : 0; # transition rate from memory to effector 
    f_e_m_fit(ρ, r_dm_max, sρ_dm) = (eff_to_mem || cyc_diff) && !no_diff ? r_dm_fit(ρ, r_dm_max, sρ_dm) : 0; # transition rate from effector to memory
    f_n_d_fit(ρ, r_te_max, sρ_te) = naive_to_exh ? r_te_fit(ρ, r_te_max, sρ_te) : 0; # transition rate from naive to terminally exhausted
    f_m_d_fit(ρ, r_te_max, sρ_te) = mem_to_exh ? r_te_fit(ρ, r_te_max, sρ_te) : 0; # transition rate from memory to terminally exhausted 
    f_e_d_fit(ρ, r_te_max, sρ_te) = eff_to_exh ? r_te_fit(ρ, r_te_max, sρ_te) : 0; # transition rate from effector to terminally exhausted 
    f_d_e_fit(ρ, r_de_max, sρ_de) = back_diff ? r_de_fit(ρ, r_de_max, sρ_de) : 0; # transition rate from exhausted to effector 
    f_n_n_prolif_fit(r_p) = naive_prolif ? r_p_fit(r_p) : 0; # symmetric naive cell proliferation rate
    f_m_n_prolif_fit(r_p) = asymm_div ? r_p_fit(r_p) : 0; # asymmetric memory cell proliferation rate
    f_m_m_prolif_fit(r_p) = asymm_div ? 0 : r_p_fit(r_p); # symmetric memory cell proliferation rate
    f_e_n_prolif_fit(r_p) = asymm_div && eff_prolif ? r_p_fit(r_p) : 0; # asymmetric effector cell proliferation rate
    f_e_e_prolif_fit(r_p) = !asymm_div && eff_prolif ? r_p_fit(r_p) : 0; # symmetric effector cell proliferation rate
    f_e_death_fit = eff_death ? r_d : 0; # effector cell death rate
    f_d_death_fit = exh_death ? r_d : 0; # terminally exhausted cell death rate

    R_nn_fit(ρ, r_dm_max, r_de_max, r_te_max, sρ_dm, sρ_de, sρ_te, r_p) = -f_n_m_fit(ρ, r_dm_max, sρ_dm) - f_n_e_fit(ρ, r_de_max, sρ_de) - f_n_d_fit(ρ, r_te_max, sρ_te) + f_n_n_prolif_fit(r_p); # stimulus-dependent reaction term driven by naive cells for n (naive cells)
    R_nm_fit(r_p) = f_m_n_prolif_fit(r_p); # reaction term driven by memory cells for n
    R_ne_fit(r_p) = f_e_n_prolif_fit(r_p); # reaction term driven by effector cells for n
    R_mn_fit(ρ, r_dm_max, sρ_dm) = f_n_m_fit(ρ, r_dm_max, sρ_dm); # stimulus-dependent reaction term driven by naive cells for m (memory cells)
    R_mm_fit(ρ, r_de_max, r_te_max, sρ_de, sρ_te, r_p) = -f_m_e_fit(ρ, r_de_max, sρ_de) - f_m_d_fit(ρ, r_te_max, sρ_te) + f_m_m_prolif_fit(r_p); # stimulus- and cytokine-dependent reaction term driven by memory cells for m
    R_me_fit(ρ, r_dm_max, sρ_dm) = f_e_m_fit(ρ, r_dm_max, sρ_dm); # stimulus-dependent reaction term driven by effector cells for m
    R_en_fit(ρ, r_de_max, sρ_de) = f_n_e_fit(ρ, r_de_max, sρ_de); # stimulus-dependent reaction term driven by naive cells for e (effector cells)
    R_em_fit(ρ, r_de_max, sρ_de) = f_m_e_fit(ρ, r_de_max, sρ_de); # stimulus-dependent reaction term driven by memory cells for e
    R_ee_fit(ρ, r_dm_max, r_te_max, sρ_dm, sρ_te, r_p) = -f_e_m_fit(ρ, r_dm_max, sρ_dm) - f_e_d_fit(ρ, r_te_max, sρ_te) + f_e_e_prolif_fit(r_p) - f_e_death_fit; # stimulus- and cytokine-dependent reaction term driven by effector cells for e 
    R_ed_fit(ρ, r_de_max, sρ_de) = f_d_e_fit(ρ, r_de_max, sρ_de); # stimulus-dependent reaction term driven by exhausted cells for e
    R_dn_fit(ρ, r_te_max, sρ_te) = f_n_d_fit(ρ, r_te_max, sρ_te); # stimulus-dependent reaction term driven by naive cells for d (terminally exhausted/dysfunction-associated cells)
    R_dm_fit(ρ, r_te_max, sρ_te) = f_m_d_fit(ρ, r_te_max, sρ_te); # stimulus-dependent reaction term driven by memory cells for d
    R_de_fit(ρ, r_te_max, sρ_te) = f_e_d_fit(ρ, r_te_max, sρ_te); # stimulus-dependent reaction term driven by effector cells for d
    R_dd_fit(ρ, r_de_max, sρ_de) = -f_d_e_fit(ρ, r_de_max, sρ_de) - f_d_death_fit; # reaction term driven by exhausted cells for d

    ODE_funcs = [R_nn_fit, R_nm_fit, R_ne_fit, R_mn_fit, R_mm_fit, R_me_fit, R_en_fit, R_em_fit, R_ee_fit, R_ed_fit, R_dn_fit, R_dm_fit, R_de_fit, R_dd_fit]; # initialise functions for ODE system

    for donor_n = 1:2 # for each donor
        r_dm_max, r_de_max, r_te_max, sρ_dm, sρ_de, sρ_te, r_p_max = param_fits_both[donor_n][model_n];

        A_s(s) = [ODE_funcs[1](ρ_avg_base*s,r_dm_max,r_de_max,r_te_max,sρ_dm,sρ_de,sρ_te,r_p_max) ODE_funcs[2](r_p_max) ODE_funcs[3](r_p_max) 0; 
            ODE_funcs[4](ρ_avg_base*s,r_dm_max,sρ_dm) ODE_funcs[5](ρ_avg_base*s,r_de_max,r_te_max,sρ_de,sρ_te,r_p_max) ODE_funcs[6](ρ_avg_base*s,r_dm_max,sρ_dm) 0; 
            ODE_funcs[7](ρ_avg_base*s,r_de_max,sρ_de) ODE_funcs[8](ρ_avg_base*s,r_de_max,sρ_de) ODE_funcs[9](ρ_avg_base*s,r_dm_max,r_te_max,sρ_dm,sρ_te,r_p_max) ODE_funcs[10](ρ_avg_base*s,r_de_max,sρ_de); 
            ODE_funcs[11](ρ_avg_base*s,r_te_max,sρ_te) ODE_funcs[12](ρ_avg_base*s,r_te_max,sρ_te) ODE_funcs[13](ρ_avg_base*s,r_te_max,sρ_te) ODE_funcs[14](ρ_avg_base*s,r_de_max,sρ_de)]; # coefficients in ODE system as a function of stimulus ration s

        for s_n = eachindex(s_range) # for each stimulus ratio
            model_outputs[worst_model_n,s_n,:,donor_n] = exp(A_s(s_range[s_n])*T)*ICs[donor_n]; # solution at time T for the donor_nth donor's initial condition
        end # end for each stimulus ratio
    end # end for each donor
end # end for each model

# plot results for cell type cell_n and donor donor_n 
plot_x = [stim for i in 1:10 for stim in stim_vals]; # arrange data to plot on x axis of boxplots

tot_cells_data = sum(final_num[[1,3,4,6],:,:,:], dims=1)[1,:,:,:]; # total cells in data 
tot_cells_sims = sum(model_outputs, dims=3)[:,:,1,:]; # total cells in simulation outputs

for cell_n = 1:4 # for each cell type;   1: naive,  2: memory,  3: effector,  4: exhausted 
    for donor_n = 1:2 # for each donor type;   1: healthy,  2: patient 
        nonNaN_inds = []; # indices of non-NaN datapoints for the current cell type
        for i = eachindex(plot_x) # for each datapoint
            if !isnan(final_num[1,:,donor_n,:][i])
                append!(nonNaN_inds, i); # if index has a non-NaN value, add to the list of indices
            end
        end

        plt = CairoMakie.Figure(fontsize=32); 
        ax = CairoMakie.Axis(plt[1, 1], xlabel="Stimulus density (ng/μg)", ylabel="$(cell_names[cell_n]) cells (millions)", limits = (0, maximum(s_range), 0, nothing))#, title="$(donor_names[donor_n]) cell models");
        for model_n = 1:1000 
            CairoMakie.lines!(s_range[start_ind:end], model_outputs[model_n,start_ind:end,cell_n,donor_n]/1e6, color=(cell_colours[cell_n]*0.8,0.1));
        end
        if T==8*24*60; CairoMakie.boxplot!(plot_x[nonNaN_inds], final_num[[1,3,4,6][cell_n],:,donor_n,:][nonNaN_inds]/1e6; color=cell_colours[cell_n]/0.8, whiskerwidth = 1, width = 1); end
        display(plt)
        save("figs_and_videos\\Fig_S4_$(cell_names[cell_n])_$(donor_names[donor_n]).png", plt)
    end
end



# S5 Fig: results from the worst 1000 models 
# Fig S5a: model loss, with the worst 1000 indicated
colours = ifelse.(inbounds_h[err_order_avg] .&& inbounds_p[err_order_avg], :blue, :red); # define colours to differentiate between models that lie within data bounds (blue) and those that don't (red)
plt = CairoMakie.Figure(fontsize=28); 
ax = CairoMakie.Axis(plt[1, 1], xlabel="Models ordered by increasing average loss", ylabel="Average normalised model loss");
CairoMakie.barplot!(eachindex(err_order_avg), model_err_avg[err_order_avg], color=colours)
CairoMakie.vlines!(num_models-worst_model_num+1, linestyle=:dash, color=:black, linewidth=3)
CairoMakie.barplot!([0],[0], color=:blue, label="In bounds") # dummy plots for making the legend
CairoMakie.barplot!([0],[0], color=:red, label="Out of bounds")
CairoMakie.axislegend(position = :lt)
display(plt)
save("figs_and_videos\\Fig_S5a.png", plt)


# Fig S5b: Yule's Q for pairs of features in the worst models 
pair_num_worst = zeros(num_choices,num_choices); # number of times a pair of model choices appears in the worst models according to average loss
pair_num_worst_1aoff = zeros(num_choices+1,num_choices+1); # number of times a pair of model choices appears in the worst models according to average loss (also considering 1a off as a feature selection)
for model_n = worst_model_inds # for each of the worst models
    curr_model = models[model_n]; # current model
    for choice_n = findall(x->x==1,curr_model) # for each individual model choice in the current model 
        pair_num_worst[choice_n,choice_n] += 1; # add to the number of times this model choice was used
    end
    for comb = combinations(findall(x->x==1,curr_model),2) # for each combination of model choices in the current model
        pair_num_worst[comb[1],comb[2]] += 1; # add to the number of times this pair of model choices is used 
    end
    if !curr_model[1] && !(model_n in ignore_1a) # if 1a is truly off (not overwritten)
        pair_num_worst_1aoff[1,1] += 1; # increment the number of times feature 1a is off
        for feat_n = 2:num_choices # for each other feature
            if curr_model[feat_n] # if that feature is on (alongside 1a off)
                pair_num_worst_1aoff[1,feat_n+1] += 1; # increment that pair
            end
        end
    end
end
pair_num_worst_1aoff[2:num_choices+1,2:num_choices+1] = pair_num_worst;
for choice_n = 1:num_choices 
    pair_num_worst[choice_n, 1:choice_n-1] .= NaN; # remove the lower triangular portion for plotting
    pair_num_worst_1aoff[choice_n, 1:choice_n-1] .= NaN; # remove the lower triangular portion for plotting
end
pair_num_worst_1aoff[num_choices+1, 1:num_choices] .= NaN; # remove the lower triangular portion for plotting

choice_names_plot = ["1a off","1a on","1b on","1c on","1d on","1e on","1f on","1g on","2a on","2b on","2c on","3a on","3b on","3c on","3d on","3e on"];

# for events X: "feature X is turned on" and A: "model is accepted as a good model (by healthy, patient or average model loss)"
PrX_worst = pair_num_tot_1aoff/num_models; # proportion of all models that have feature X turned on   Pr(X)
PrXA_worst = pair_num_worst_1aoff/worst_model_num; # Pr(X|A) for worst models
PrA_worst = worst_model_num/num_models; # Pr(A) for worst models

Q_worst = (PrXA_worst - PrX_worst)./(PrXA_worst - 2*PrXA_worst.*(PrX_worst .+ PrA_worst - PrA_worst.*PrXA_worst) + PrX_worst)#(OR_worst.-1)./(OR_worst.+1); # Yule's Q for features in models ordered by average loss

plt = CairoMakie.Figure(fontsize=32, size=(800,600)); 
ax = CairoMakie.Axis(plt[1, 1], xlabel="Model feature selection", ylabel="Model feature selection", xticks = (1:length(choice_names_plot), choice_names_plot), yticks = (1:length(choice_names_plot), choice_names_plot), xticklabelrotation=pi/2);
hmap = CairoMakie.heatmap!(1:length(choice_names_plot), 1:length(choice_names_plot), Q_worst; colormap=cgrad(:vik, rev=true), colorrange=(-1,1))
CairoMakie.Colorbar(plt[1, 2], hmap; label="Yule's Q", width=15, ticksize=5, tickalign=1);
CairoMakie.colsize!(plt.layout, 1, Aspect(1, 1.0));
CairoMakie.colgap!(plt.layout, 10)
display(plt)
save("figs_and_videos\\Fig_S5b.png", plt)


# Fig S5c: Yule's Q compared to literature for the worst 1000 models 
plt = CairoMakie.Figure(fontsize=26, size=(600,400)); ax = CairoMakie.Axis(plt[1, 1], xticks = (1:16, choice_names_plot), xticklabelrotation=pi/2, xlabel="Model feature selection", ylabel="Yule's Q / Weighted average", limits=(0,17,-1.1,1.1));
plot_cols = [(:gray,0.9), (col_nai,0.75)]; # colours for barplot
CairoMakie.barplot!([i for i in 1:num_choices+1 for j in 1:2], [lit_support_ratio Q_worst[1:length(choice_names_plot)+1:length(choice_names_plot)^2]]'[:], color=plot_cols[[j for i in 1:num_choices+1 for j in 1:2]])
CairoMakie.barplot!([0], color=plot_cols[1], label="Literature")
CairoMakie.barplot!([0], color=plot_cols[2], label="Worst models")
CairoMakie.axislegend(halign=:left, valign=:top);
display(plt)
save("figs_and_videos\\Fig_S5c.png", plt)


# Fig S5d: parameter values in worst 1000 Models
# plot for parameters with units 1/day
plt = CairoMakie.Figure(fontsize=28, size=(350,500)); ax = CairoMakie.Axis(plt[1, 1], xticks = (1.5:3:3*4-0.5, param_names_short[[1,2,3,7]]), xticklabelrotation=pi/4, xlabel="Parameter", ylabel="Parameter value (day⁻¹)", limits=(0,nothing,0,nothing));
for param_n_plot = 1:4
    param_n = [1,2,3,7][param_n_plot];
    CairoMakie.boxplot!((3*param_n_plot-2)*ones(worst_model_num), [param_fits_h[i][param_n] for i in worst_model_inds]*param_scales[param_n]; whiskerwidth = 1, width = 1, color=donor_cols[1]);#, show_outliers=false
    CairoMakie.boxplot!((3*param_n_plot-1)*ones(worst_model_num), [param_fits_p[i][param_n] for i in worst_model_inds]*param_scales[param_n]; whiskerwidth = 1, width = 1, color=donor_cols[2]);
end
CairoMakie.boxplot!([-1], [0]; whiskerwidth = 1, width = 1, color=donor_cols[1], label="Healthy");
CairoMakie.boxplot!([-1], [0]; whiskerwidth = 1, width = 1, color=donor_cols[2], label="Patient");
CairoMakie.axislegend(position = :rt)
display(plt)
save("figs_and_videos\\Fig_S5d_left.png", plt)

# plot for parameters with units ng/μg
plt = CairoMakie.Figure(fontsize=28, size=(280,500)); ax = CairoMakie.Axis(plt[1, 1], xticks = (1.5:3:3*3-0.5, param_names_short[[4,5,6]]), xticklabelrotation=pi/4, xlabel="Parameter", ylabel="Parameter value (ng/μg)", limits=(0,nothing,0,nothing));
for param_n_plot = 1:3
    param_n = [4,5,6][param_n_plot];
    CairoMakie.boxplot!((3*param_n_plot-2)*ones(worst_model_num), [param_fits_h[i][param_n] for i in worst_model_inds]*param_scales[param_n]; whiskerwidth = 1, width = 1, color=donor_cols[1]);#, show_outliers=false
    CairoMakie.boxplot!((3*param_n_plot-1)*ones(worst_model_num), [param_fits_p[i][param_n] for i in worst_model_inds]*param_scales[param_n]; whiskerwidth = 1, width = 1, color=donor_cols[2]);
end
CairoMakie.boxplot!([-1], [0]; whiskerwidth = 1, width = 1, color=donor_cols[1], label="Healthy");
CairoMakie.boxplot!([-1], [0]; whiskerwidth = 1, width = 1, color=donor_cols[2], label="Patient");
display(plt)
save("figs_and_videos\\Fig_S5d_right.png", plt)


# Fig S5e: difference in parameter values in worst 1000 models 
# plot for parameters with units 1/day
plt = CairoMakie.Figure(fontsize=28, size=(300,500)); ax = CairoMakie.Axis(plt[1, 1], xticks = (1:4, param_names_short[[1,2,3,7]]), xticklabelrotation=pi/4, xlabel="Parameter", ylabel="Parameter difference, h-p (day⁻¹)", limits=(nothing, nothing, -10, 10));
for param_n_plot = 1:4
    param_n = [1,2,3,7][param_n_plot];
    CairoMakie.boxplot!(param_n_plot*ones(worst_model_num), ([param_fits_h[i][param_n] for i in worst_model_inds]-[param_fits_p[i][param_n] for i in worst_model_inds])*param_scales[param_n]; whiskerwidth = 1, width = 1, color=col_nai);#, show_outliers=false
end
display(plt)
save("figs_and_videos\\Fig_S5e_left.png", plt)

# plot for parameters with units ng/μg
plt = CairoMakie.Figure(fontsize=28, size=(250,500)); ax = CairoMakie.Axis(plt[1, 1], xticks = (1:3, param_names_short[[4,5,6]]), xticklabelrotation=pi/4, xlabel="Parameter", ylabel="Parameter difference, h-p (ng/μg)", limits=(nothing, nothing, -20, 20));
for param_n_plot = 1:3
    param_n = [4,5,6][param_n_plot];
    CairoMakie.boxplot!(param_n_plot*ones(worst_model_num), ([param_fits_h[i][param_n] for i in worst_model_inds]-[param_fits_p[i][param_n] for i in worst_model_inds])*param_scales[param_n]; whiskerwidth = 1, width = 1, color=col_nai);#, show_outliers=false
end
display(plt)
save("figs_and_videos\\Fig_S5e_right.png", plt)







## code to generate SI Videos (these will save in the current directory)

order_by_loss = false; # true to order animation by model loss, choose order below. otherwise, order by array model_inds_animate
order_name = "average"; # name of dataset/fitted models for title:  choose "healthy" to reproduce S1 Video (saves as "S1_Video"), "patient" to reproduce S2 Video (saves as "S2_Video"), and "average" to reproduce S3 Video (saves as "S3_Video")
num_models_animate = 300; # number of models to animate (either ordered by loss or otherwise)
model_inds_animate = 1:num_models_animate; # if not ordering by model loss, define order here

show_choices = false; # true to show model choices in animation
framerate = 10; # frames per second in animation

if order_by_loss
    if order_name == "healthy"
        model_inds_animate = err_order_h[1:num_models_animate]; # indices of models to animate (select models from the list ordered by loss)
        anim_name = "S1_Video"; # animation file name
    elseif order_name == "patient" 
        model_inds_animate = err_order_p[1:num_models_animate]; 
        anim_name = "S2_Video"; # animation file name
    elseif order_name == "average"
        model_inds_animate = err_order_avg[1:num_models_animate]; 
        anim_name = "S3_Video"; # animation file name
    else
        error("Pick order_name as 'healthy', 'patient', or 'average'")
    end
    title_text = " (ordered by $order_name model loss)"; # text to go in animation title
else
    anim_name = "all_models";
    title_text = "";
end


col_diff = RGB(0/255,0/255,0/255); # colour for state transitions/differentiation and exhaustion arrows
col_prolif = RGB(33/255,69/255,166/255); # colour for proliferation event arrows
col_death = RGB(145/255,33/255,33/255); # colour for death event arrows
arrow_width = 5; # linewidth for arrows

cell_pos = [-1 0; 0 0; 1 0; 0.5 -0.5]; # positions of each cell in the graph [cell type, x/y]  (origin on memory cell)

scale_down = 1; # factor to scale the plot down by
state_graph = CairoMakie.Figure(fontsize=25/scale_down, size=(600,350)./scale_down); # initialise graph
ax = CairoMakie.Axis(state_graph[1, 1], limits=(minimum(cell_pos[:,1])-0.25,maximum(cell_pos[:,1])+0.26,minimum(cell_pos[:,2])-0.26,maximum(cell_pos[:,2])+0.4), aspect=DataAspect()); # define axis with title
hidespines!(ax); # remove axis grid and other background elements
hidedecorations!(ax);

# create cells in figure
CairoMakie.scatter!(cell_pos[:,1], cell_pos[:,2]; markersize=100/scale_down, color=[col_nai,col_mem,col_eff,col_exh], strokecolor=:black, strokewidth=3/scale_down); # draw cells 
CairoMakie.scatter!(cell_pos[:,1]+(rand(4).-0.5)/20, cell_pos[:,2]+(rand(4).-0.5)/20; markersize=60/scale_down, color=[col_nai,col_mem,col_eff,col_exh]/1.5); # draw cell nuclei

# draw all lines and arrows into figure first
# all state transitions: differentiation and exhaustion 
n2m = CairoMakie.lines!([cell_pos[1,1]+0.2,cell_pos[2,1]-0.22], [cell_pos[1,2],cell_pos[2,2]], color=col_diff, linewidth=arrow_width/scale_down); # line and arrowhead for naive -> memory
n2e = CairoMakie.lines!(cell_pos[1,1]+0.14:0.01:cell_pos[3,1]-0.15, -0.11*(cell_pos[1,1]+0.14:0.01:cell_pos[3,1]-0.15).^2 .+ 0.22, color=col_diff, linewidth=arrow_width/scale_down); # naive -> effector
n2d = CairoMakie.lines!(cell_pos[1,1]+0.16:0.01:cell_pos[4,1]-0.22, 0.29*((cell_pos[1,1]+0.16:0.01:cell_pos[4,1]-0.22).-0.3).^2 .- 0.5, color=col_diff, linewidth=arrow_width/scale_down); # naive -> exhausted
m2e = CairoMakie.lines!([cell_pos[2,1]+0.2,cell_pos[3,1]-0.22], [cell_pos[2,2]+0.03,cell_pos[3,2]+0.03], color=col_diff, linewidth=arrow_width/scale_down); # memory -> effector
m2d = CairoMakie.lines!([cell_pos[2,1]+0.14,cell_pos[4,1]-0.15], [cell_pos[2,2]-0.14,cell_pos[4,2]+0.15], color=col_diff, linewidth=arrow_width/scale_down); # memory -> exhausted
e2m = CairoMakie.lines!([cell_pos[2,1]+0.22,cell_pos[3,1]-0.2], [cell_pos[2,2]-0.03,cell_pos[3,2]-0.03], color=col_diff, linewidth=arrow_width/scale_down); # effector -> memory
e2d = CairoMakie.lines!([cell_pos[3,1]-0.125,cell_pos[4,1]+0.165], [cell_pos[3,2]-0.165,cell_pos[4,2]+0.125], color=col_diff, linewidth=arrow_width/scale_down); # effector -> exhausted
d2e = CairoMakie.lines!([cell_pos[3,1]-0.165,cell_pos[4,1]+0.125], [cell_pos[3,2]-0.125,cell_pos[4,2]+0.165], color=col_diff, linewidth=arrow_width/scale_down); # exhausted -> effector

# all proliferation events 
snp = CairoMakie.arc!(cell_pos[1,:]+[-0.1,0.2], 0.12, 1.3*pi, -0.06*pi, color=col_prolif, linewidth=arrow_width/scale_down)#, linestyle=(:dash,1.5)); # symmetric naive proliferation
smp = CairoMakie.arc!(cell_pos[2,:]+[-0.05,-0.22], 0.12, -1.25*pi, 0.05*pi, color=col_prolif, linewidth=arrow_width/scale_down)#, linestyle=(:dash,1.5)); # symmetric memory proliferation
amp = CairoMakie.lines!(cell_pos[1,1]+0.2:0.01:cell_pos[2,1]-0.18, 0.3*((cell_pos[1,1]+0.2:0.01:cell_pos[2,1]-0.18).+0.5).^2 .- 0.09, color=col_prolif, linewidth=arrow_width/scale_down)#, linestyle=(:dash,1.5)); # asymmetric memory proliferation
sep = CairoMakie.arc!(cell_pos[3,:]+[0.13,0.2], 0.12, 1.0*pi, -0.4*pi, color=col_prolif, linewidth=arrow_width/scale_down)#, linestyle=(:dash,1.5)); # symmetric effector proliferation
aep = CairoMakie.lines!(cell_pos[1,1]+0.1:0.01:cell_pos[3,1]-0.08, -0.17*(cell_pos[1,1]+0.1:0.01:cell_pos[3,1]-0.08).^2 .+ 0.33, color=col_prolif, linewidth=arrow_width/scale_down)#, linestyle=(:dash,1.5)); # asymmetric effector proliferation

# all death events 
ed = CairoMakie.lines!([cell_pos[3,1]+0.13,cell_pos[3,1]+0.25], [cell_pos[3,2]-0.13,cell_pos[3,2]-0.25], color=col_death, linewidth=arrow_width/scale_down)#, linestyle=(:dash,1.5)); # effector death
dd = CairoMakie.lines!([cell_pos[4,1]+0.13,cell_pos[4,1]+0.25], [cell_pos[4,2]-0.13,cell_pos[4,2]-0.25], color=col_death, linewidth=arrow_width/scale_down)#, linestyle=(:dash,1.5)); # exhausted death

if use_arrows2d # use arrows2d function
    n2ma = CairoMakie.arrows2d!([cell_pos[2,1]-0.26], [cell_pos[2,2]], [0.1], [0]; tipwidth=arrow_width/scale_down*4, tiplength=arrow_width/scale_down*3, shaftwidth=0); # arrow for naive -> memory
    n2ea = CairoMakie.arrows2d!([cell_pos[3,1]-0.19], [-0.11*(cell_pos[3,1]-0.19)^2+0.225], [0.1], [-0.04]; tipwidth=arrow_width/scale_down*4, tiplength=arrow_width/scale_down*3, shaftwidth=0); # etc...
    n2da = CairoMakie.arrows2d!([cell_pos[4,1]-0.26], [0.29*(cell_pos[4,1]-0.54)^2-0.5], [0.1], [0]; tipwidth=arrow_width/scale_down*4, tiplength=arrow_width/scale_down*3, shaftwidth=0);
    m2ea = CairoMakie.arrows2d!([cell_pos[3,1]-0.26], [cell_pos[3,2]+0.03], [0.1], [0]; tipwidth=arrow_width/scale_down*4, tiplength=arrow_width/scale_down*3, shaftwidth=0);
    m2da = CairoMakie.arrows2d!([cell_pos[4,1]-0.18], [cell_pos[4,2]+0.18], [0.08], [-0.08]; tipwidth=arrow_width/scale_down*4, tiplength=arrow_width/scale_down*3, shaftwidth=0); 
    e2ma = CairoMakie.arrows2d!([cell_pos[2,1]+0.26], [cell_pos[2,2]-0.03], [-0.1], [0]; tipwidth=arrow_width/scale_down*4, tiplength=arrow_width/scale_down*3, shaftwidth=0);
    e2da = CairoMakie.arrows2d!([cell_pos[4,1]+0.2], [cell_pos[4,2]+0.16], [-0.08], [-0.08]; tipwidth=arrow_width/scale_down*4, tiplength=arrow_width/scale_down*3, shaftwidth=0);
    d2ea = CairoMakie.arrows2d!([cell_pos[3,1]-0.2], [cell_pos[3,2]-0.16], [0.08], [0.08]; tipwidth=arrow_width/scale_down*4, tiplength=arrow_width/scale_down*3, shaftwidth=0);
    snpa = CairoMakie.arrows2d!([cell_pos[1,1]-0.2], [cell_pos[1,2]+0.13], [0.1], [-0.08]; tipwidth=arrow_width/scale_down*4, tiplength=arrow_width/scale_down*3, shaftwidth=0, color=col_prolif); # proliferation arrows...
    smpa = CairoMakie.arrows2d!([cell_pos[2,1]+0.065], [cell_pos[2,2]-0.23], [-0.02], [0.1]; tipwidth=arrow_width/scale_down*4, tiplength=arrow_width/scale_down*3, shaftwidth=0, color=col_prolif);
    ampa = CairoMakie.arrows2d!([cell_pos[1,1]+0.24], [0.3*(cell_pos[1,1]+0.74)^2-0.09], [-0.1], [0.02]; tipwidth=arrow_width/scale_down*4, tiplength=arrow_width/scale_down*3, shaftwidth=0, color=col_prolif);
    sepa = CairoMakie.arrows2d!([cell_pos[3,1]+0.01], [cell_pos[3,2]+0.24], [0], [-0.1]; tipwidth=arrow_width/scale_down*4, tiplength=arrow_width/scale_down*3, shaftwidth=0, color=col_prolif);
    aepa = CairoMakie.arrows2d!([cell_pos[1,1]+0.13], [-0.17*(cell_pos[1,1]+0.13)^2+0.33], [-0.1], [-0.07]; tipwidth=arrow_width/scale_down*4, tiplength=arrow_width/scale_down*3, shaftwidth=0, color=col_prolif);
    eda = CairoMakie.arrows2d!([cell_pos[3,1]+0.215], [cell_pos[3,2]-0.215], [0.08], [-0.08]; tipwidth=arrow_width/scale_down*4, tiplength=arrow_width/scale_down*3, shaftwidth=0, color=col_death); # death arrows...
    dda = CairoMakie.arrows2d!([cell_pos[4,1]+0.215], [cell_pos[4,2]-0.215], [0.08], [-0.08]; tipwidth=arrow_width/scale_down*4, tiplength=arrow_width/scale_down*3, shaftwidth=0, color=col_death); 
else # use arrows function
    n2ma = CairoMakie.arrows!([cell_pos[2,1]-0.24], [cell_pos[2,2]], [0.01], [0]; arrowsize = arrow_width/scale_down*4); # arrow for naive -> memory
    n2ea = CairoMakie.arrows!([cell_pos[3,1]-0.175], [-0.11*(cell_pos[3,1]-0.175)^2+0.22], [0.01], [-0.004]; arrowsize = arrow_width/scale_down*4); # etc...
    n2da = CairoMakie.arrows!([cell_pos[4,1]-0.24], [0.29*(cell_pos[4,1]-0.52)^2-0.5], [0.01], [0]; arrowsize = arrow_width/scale_down*4);
    m2ea = CairoMakie.arrows!([cell_pos[3,1]-0.24], [cell_pos[3,2]+0.03], [0.01], [0]; arrowsize = arrow_width/scale_down*4);
    m2da = CairoMakie.arrows!([cell_pos[4,1]-0.17], [cell_pos[4,2]+0.17], [0.01], [-0.01]; arrowsize = arrow_width/scale_down*4);
    e2ma = CairoMakie.arrows!([cell_pos[2,1]+0.24], [cell_pos[2,2]-0.03], [-0.01], [0]; arrowsize = arrow_width/scale_down*4); 
    e2da = CairoMakie.arrows!([cell_pos[4,1]+0.19], [cell_pos[4,2]+0.15], [-0.01], [-0.01]; arrowsize = arrow_width/scale_down*4); 
    d2ea = CairoMakie.arrows!([cell_pos[3,1]-0.19], [cell_pos[3,2]-0.15], [0.01], [0.01]; arrowsize = arrow_width/scale_down*4); 
    snpa = CairoMakie.arrows!([cell_pos[1,1]-0.19], [cell_pos[1,2]+0.12], [0.01], [-0.008]; arrowsize = arrow_width/scale_down*4, color=col_prolif); # proliferation arrows...
    smpa = CairoMakie.arrows!([cell_pos[2,1]+0.065], [cell_pos[2,2]-0.21], [-0.0002], [0.001]; arrowsize = arrow_width/scale_down*4, color=col_prolif); 
    ampa = CairoMakie.arrows!([cell_pos[1,1]+0.23], [0.3*(cell_pos[1,1]+0.73)^2-0.09], [-0.01], [0.002]; arrowsize = arrow_width/scale_down*4, color=col_prolif); 
    sepa = CairoMakie.arrows!([cell_pos[3,1]+0.01], [cell_pos[3,2]+0.22], [0.0], [-0.001]; arrowsize = arrow_width/scale_down*4, color=col_prolif); 
    aepa = CairoMakie.arrows!([cell_pos[1,1]+0.12], [-0.17*(cell_pos[1,1]+0.12)^2+0.33], [-0.01], [-0.007]; arrowsize = arrow_width/scale_down*4, color=col_prolif); 
    eda = CairoMakie.arrows!([cell_pos[3,1]+0.225], [cell_pos[3,2]-0.225], [0.01], [-0.01]; arrowsize = arrow_width/scale_down*4, color=col_death); # death arrows... 
    dda = CairoMakie.arrows!([cell_pos[4,1]+0.225], [cell_pos[4,2]-0.225], [0.01], [-0.01]; arrowsize = arrow_width/scale_down*4, color=col_death);
end

# legend 
CairoMakie.scatter!([-2], [-2]; markersize=25/scale_down, color=col_nai, strokecolor=:black, strokewidth=2/scale_down, label="Naive"); # dummy plots to use for legend
CairoMakie.scatter!([-2], [-2]; markersize=25/scale_down, color=col_mem, strokecolor=:black, strokewidth=2/scale_down, label="Memory-like");
CairoMakie.scatter!([-2], [-2]; markersize=25/scale_down, color=col_eff, strokecolor=:black, strokewidth=2/scale_down, label="Effector-like");
CairoMakie.scatter!([-2], [-2]; markersize=25/scale_down, color=col_exh, strokecolor=:black, strokewidth=2/scale_down, label="Dysfunction-associated");
CairoMakie.lines!([-2], [-2], color=col_diff, linewidth=arrow_width/scale_down*0.75, label="Transition");
CairoMakie.lines!([-2], [-2], color=col_prolif, linewidth=arrow_width/scale_down*0.75, label="Proliferation"); # , linestyle=(:dot,1.2)
CairoMakie.lines!([-2], [-2], color=col_death, linewidth=arrow_width/scale_down*0.75, label="Death"); # , linestyle=(:dot,1.2)
CairoMakie.axislegend(position=(-0.03*scale_down,-0.2*scale_down), framevisible=false, nbanks=2, labelsize=18/scale_down, colgap=15/scale_down, patchsize = (15/scale_down, 10)); # construct legend

#display(state_graph); # view full graph with all arrows

# animate the figure where each frame turns off various arrows in the full pathway
nframes = num_models; # number of frames in animation 
record(state_graph, "figs_and_videos\\$anim_name.mp4", eachindex(model_inds_animate); framerate=framerate) do model_n_ordered # for each model with ordered index model_n_ordered, gif saved as "state_graph"
    
    model_n = model_inds_animate[model_n_ordered]; # extract model index in models array

    if show_choices # if showing model choices on the graphs
        if model_n in ignore_1a # if 1a is not considered for the current model, ignore 1a in the title
            if sum(models[model_n]) <= num_choices÷2 # display the on model choices
                ax.title = "Model $model_n_ordered$title_text:\n$(join(["$(choice_names[i+1])" for i in findall(x->x==true,models[model_n][2:end])],", ")) on, the rest off (1a overwritten)"; # set current title for on choices 
            else # display the off model choices
                ax.title = "Model $model_n_ordered$title_text:\n$(join(["$(choice_names[i+1])" for i in findall(x->x==false,models[model_n][2:end])],", ")) off, the rest on (1a overwritten)"; # set current title for off choices
            end
        else
            if sum(models[model_n]) <= num_choices÷2 # display the on model choices
                ax.title = "Model $model_n_ordered$title_text:\n$(join(["$(choice_names[i])" for i in findall(x->x==true,models[model_n])],", ")) on, the rest off"; # set current title for on choices 
            else # display the off model choices
                ax.title = "Model $model_n_ordered$title_text:\n$(join(["$(choice_names[i])" for i in findall(x->x==false,models[model_n])],", ")) off, the rest on"; # set current title for off choices
            end
        end
    else # not showing model choices
        ax.title = "Model $model_n_ordered$title_text"; # title is model number
    end

    c1a, c1b, c1c, c1d, c1e, c1f, c1g, c2a, c2b, c2c, c3a, c3b, c3c, c3d, c3e = models[model_n]; # extract all model choices

    if (!c1a && !c1c) || c1b; n2m.visible = true; n2ma.visible = true; else; n2m.visible = false; n2ma.visible = false; end
    if c1a || c1b || c1c; n2e.visible = true; n2ea.visible = true; else; n2e.visible = false; n2ea.visible = false; end
    if c2a; n2d.visible = true; n2da.visible = true; else; n2d.visible = false; n2da.visible = false; end
    if (!c1a || c1c) && !c1d; m2e.visible = true; m2ea.visible = true; else; m2e.visible = false; m2ea.visible = false; end
    if c2b; m2d.visible = true; m2da.visible = true; else; m2d.visible = false; m2da.visible = false; end 
    if (c1a || c1c) && !c1d; e2m.visible = true; e2ma.visible = true; else; e2m.visible = false; e2ma.visible = false; end
    if c2c; e2d.visible = true; e2da.visible = true; else; e2d.visible = false; e2da.visible = false; end
    if c1g; d2e.visible = true; d2ea.visible = true; else; d2e.visible = false; d2ea.visible = false; end

    if c3c; snp.visible = true; snpa.visible = true; else; snp.visible = false; snpa.visible = false; end
    if !c3a; smp.visible = true; smpa.visible = true; else; smp.visible = false; smpa.visible = false; end
    if c3a; amp.visible = true; ampa.visible = true; else; amp.visible = false; ampa.visible = false; end
    if !c3a && c3b; sep.visible = true; sepa.visible = true; else; sep.visible = false; sepa.visible = false; end
    if c3a && c3b; aep.visible = true; aepa.visible = true; else; aep.visible = false; aepa.visible = false; end

    if c3d; ed.visible = true; eda.visible = true; else; ed.visible = false; eda.visible = false; end
    if c3e; dd.visible = true; dda.visible = true; else; dd.visible = false; dda.visible = false; end

    println("For $order_name animation: $model_n_ordered out of $(length(model_inds_animate))")
end
