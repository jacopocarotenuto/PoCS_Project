using Graphs, Plots, StatsBase, LsqFit, Distributions, Pkg

function SetUpEnviroment()
    current_dir = pwd()
    if current_dir[end-11:end] == "PoCS_Project"
        Pkg.activate("./")
        Pkg.instantiate()
    elseif current_dir[end-7:end] == "task_15"
        Pkg.activate("../../")
        Pkg.instantiate()
    else
        print("Please run script from the main directory or the script parent directory")
    end

    # Get the current directory
    current_dir = pwd()
    # Check if the current directory is the one with the script
    if current_dir[end-11:end] == "PoCS_Project"
        # If not change to the directory with the script
        cd("code/task_15")
    end
end

SetUpEnviroment()
gr()

file_to_write_results = "./FitResults.txt"
AllStruct = []

# Data Structure to store all relevant information of a Sandpile Simulation
struct SandPile
    graph::SimpleGraph
    graph_type::String
    shedding_probability::Float64
    total_grains::Int64
    max_iter::Int64
    NumberOfCascades::Int64
    CascadeSizes::Array{Int64}
    CascadeDurations::Array{Int64}
    Sizes::Array{Int64}
    SizesProbabilities::Array{Float64}
    TauSizes::Float64
    CharacSize::Float64
    Times::Array{Int64}
    TimesProbabilities::Array{Float64}
    TauTimes::Float64
    CharacTime::Float64
    TauSizesTimes::Float64
end

# Function so simulate the Sandpile Dynamics
function SandPile_Dynamic(graph, shedding_probability::Float64, total_grains::Int; max_iter::Int=1000)
    # Initialize variables from graph
    N = nv(graph)
    ids = 1:N
    thresholds = degree(graph) # Thresholds are set equal to the degree of the nodes
    heights = zeros(N)
    neighbours_list = [neighbors(graph, i) for i in ids] # Precompute neighbours list

    # Initialize statistics collection variables
    NumberOfCascades = 0
    CascadeSizes = Array{Int64}(undef, 0)
    CascadeDurations = Array{Int64}(undef, 0)


    # Main Dynamics Loop
    for step in 1:total_grains
        heights[rand(1:N)] += 1 # 1) Add one grain to a random node

        overloaded = ids[heights .> thresholds] # 2) Check for overloaded nodes
        if length(overloaded) > 0 # 3) If there are overloaded nodes
            
            NumberOfCascades += 1 # Start a new cascade
            append!(CascadeSizes,1)
            append!(CascadeDurations,0)
            iter = 0
            while (length(overloaded) > 0) & (iter < max_iter) # 4) While there are overloaded nodes
                CascadeDurations[end] += 1
                # Move grains to neighbours
                for i in overloaded # For each overloaded node
                    neighbours = neighbours_list[i]
                    CascadeSizes[end] += length(neighbours)
                    for j in neighbours
                        heights[j] += 1 # Move one grain to each neighbour
                    end
                    heights[i] = heights[i] - length(neighbours) - ceil(heights[i]*shedding_probability) # Update height of the overloaded node and loose some grains according to the shedding probability
                end
                overloaded = ids[heights .> thresholds]
                iter += 1
            end
        end
    end

    return heights, NumberOfCascades, CascadeSizes, CascadeDurations
end

# Function to analyze Cascade Sizes Distribution
function AnalyzeCascadeSizes(CascadeSizes)
    h = fit(Histogram, CascadeSizes, nbins=floor(maximum(CascadeSizes)/5)) # Fit a histogram to the cascade sizes
    Probabilities = h.weights[h.weights .> 0] # Eliminate empty bins
    Probabilities = Probabilities / sum(Probabilities) # Normalize the probabilities
    Sizes = collect(h.edges[1])[2:end][h.weights .> 0] 

    Sizes = Sizes[Probabilities .> minimum(Probabilities)*5] # Eliminate small probabilities
    Probabilities = Probabilities[Probabilities .> minimum(Probabilities)*5] # Eliminate small probabilities

    
    @. model(x, p) = -p[1]*x - exp(x)/p[2]
    fitt = curve_fit(model, log.(Sizes), log.(Probabilities), [1.5,10000]) # Fit a model to the data
    tau, sc = fitt.param
    return Sizes, Probabilities, tau, sc
end

# Function to analyze Cascade Times Distribution
function AnalyzeCascadeTimes(CascadeTimes)
    h = fit(Histogram, CascadeTimes, nbins=floor(maximum(CascadeTimes))) # Fit a histogram to the cascade times
    Probabilities = h.weights[h.weights .> 0] # Eliminate empty bins
    Probabilities = Probabilities / sum(Probabilities) # Normalize the probabilities
    Times = collect(h.edges[1])[2:end][h.weights .> 0]

    Times = Times[Probabilities .> minimum(Probabilities)*5] # Eliminate small probabilities
    Probabilities = Probabilities[Probabilities .> minimum(Probabilities)*5] # Eliminate small probabilities

    @. model(x, p) = -p[1]*x - exp(x)/p[2]
    fitt = curve_fit(model, log.(Times), log.(Probabilities), [1.5,100])
    tau = fitt.param[1]
    sc = fitt.param[2]

    return Times, Probabilities, tau, sc
end

# Function to analyze the relation between Cascade Sizes and Times
function AnalyzeTimes_vs_Sizes(CascadeSizes, CascadeTimes)
    @. model(x, p) = p[1]*x
    fitt = curve_fit(model, log.(CascadeTimes), log.(CascadeSizes), [1.5]) # Fit a model to the data
    tau = fitt.param[1]
    return tau
end

# Function to run the complete pipeline of the Sandpile Simulation
function CompleteSandpilePipeline(graph, graph_type::String; shedding_probability::Float64=0.0001, total_grains::Int=100000, max_iter::Int=1000)
    z, NumberOfCascades, CascadeSizes, CascadeDurations = SandPile_Dynamic(graph, shedding_probability, total_grains, max_iter=max_iter) # Run the sandpile dynamics
    Sizes, SizesProbabilities, TauSizes, CharacSize = AnalyzeCascadeSizes(CascadeSizes) # Analyze the cascade sizes
    Times, TimesProbabilities, TauTimes, CharacSizeTimes = AnalyzeCascadeTimes(CascadeDurations) # Analyze the cascade times
    tauSizesTimes = AnalyzeTimes_vs_Sizes(CascadeSizes, CascadeDurations) # Analyze the relation between cascade sizes and times
    return SandPile(g, graph_type, shedding_probability, total_grains, max_iter, NumberOfCascades, CascadeSizes, CascadeDurations, Sizes, SizesProbabilities, TauSizes, CharacSize, Times, TimesProbabilities, TauTimes, CharacSizeTimes, tauSizesTimes)
end

# Function to plot the results of the Sandpile Simulation
function PlotResults(file_name, SandPile)
    @. ExponentialThresholdModel(x, p) = -p[1]*x - exp(x)/p[2]
    @. SimpleExponentialModel(x, p) = -p[1]*x
    
    
    ## Plot Cascade Size Distribution
    CascadeDistributionPlot = plot(SandPile.Sizes[3:end], SandPile.SizesProbabilities[3:end], legend=false, xaxis=:log, yaxis=:log,linewidth=6)
    CascadeDistributionPlot = plot!(SandPile.Sizes[3:end], exp.(ExponentialThresholdModel(log.(SandPile.Sizes), [SandPile.TauSizes, SandPile.CharacSize]))[3:end], xaxis=:log, yaxis=:log, legend=false, xlabel="Cascade Sizes", ylabel="Probability", title="Cascade Size Distribution", yguidefontsize=25,xguidefontsize=25, titlefontsize=25,xtickfontsize=18,ytickfontsize=18, linewidth=6)
    
    ## Plot Cascade Times Distribution
    CascadeTimesDistributionPlot = plot(SandPile.Times[3:end], SandPile.TimesProbabilities[3:end], legend=false, xaxis=:log, yaxis=:log,linewidth=6)
    CascadeTimesDistributionPlot = plot!(SandPile.Times[3:end], exp.(ExponentialThresholdModel(log.(SandPile.Times), [SandPile.TauTimes, SandPile.CharacTime]))[3:end], xaxis=:log, yaxis=:log, legend=false, xlabel="Cascade Durations", ylabel="Probability", title="Cascade Durations Distribution", yguidefontsize=25,xguidefontsize=25, titlefontsize=25,xtickfontsize=18,ytickfontsize=18,linewidth=6)
    
    ## Plot Cascade Times vs Cascade Sizes
    CascadeTimesVsSizesPlot = scatter(SandPile.CascadeDurations, SandPile.CascadeSizes, xlabel="Cascade Sizes", ylabel="Cascade Durations", legend=false, xaxis=:log,yaxis=:log, title = "Cascade Durations vs Cascade Sizes", markershape=:x, markerstrokewidth=0.1, yguidefontsize=25,xguidefontsize=25, titlefontsize=25,xtickfontsize=18,ytickfontsize=18)
    CascadeTimesVsSizesPlot = plot!(range(minimum(SandPile.CascadeDurations),maximum(SandPile.CascadeDurations),length=100), range(minimum(SandPile.CascadeDurations),maximum(SandPile.CascadeDurations),length=100).^(SandPile.TauSizesTimes), xaxis=:log, yaxis=:log ,legend=false, linewidth = 3)

    ## Layout Plots
    plot(CascadeDistributionPlot, CascadeTimesDistributionPlot, CascadeTimesVsSizesPlot, layout=(1,3), size=(2880,1080),margin = 1.8Plots.cm, dpi=400)

    ## Save Plots
    savefig(file_name)
    return nothing
end

# Function to save the results of the fits to a file
function SaveResultsToFile(file_name, SandPile; io_mode = "a")
    results = 
    """
    ###########################
    GRAPH TYPE: $(SandPile.graph_type)
    NUMBER OF CASCADES: $(SandPile.NumberOfCascades)
    TAU SIZE: $(SandPile.TauSizes)
    CHARACTERISTIC SIZE: $(SandPile.CharacSize)
    TAU TIMES: $(SandPile.TauTimes)
    TAU SIZES vs TIMES: $(SandPile.TauSizesTimes)

    """
    
    open(file_name, io_mode) do f
        write(f, results)
    end
end

# Function to plot all the results in a single plot
function TotalPlot(AllSandpiles)
    plot(AllSandpiles[1].Sizes[3:end], AllSandpiles[1].SizesProbabilities[3:end], xaxis=:log, yaxis=:log, title="Cascade Size Distribution", xlabel="Cascade Sizes", ylabel="Probability", label=AllSandpiles[1].graph_type, legendfontsize=6)
    for sp in AllSandpiles[2:end]
        plot!(sp.Sizes[3:end], sp.SizesProbabilities[3:end], xaxis=:log, yaxis=:log, label=sp.graph_type)
    end
    savefig("./Plots/TotalPlot.pdf")
    return nothing
end


###### WARNING: All this file will take a while to run.


### Gaussian Degree Distribution
print("\rSimulating Gaussian Degree Distribution Graph                                                               ")

N = 100000
degree_sequence = abs.(Int.(floor.(rand(Normal(20,6), N))))
if sum(degree_sequence) % 2 != 0
    degree_sequence[rand(1:N)] += 1
end
g = random_configuration_model(N, degree_sequence)
GaussianDegreeDistribution = CompleteSandpilePipeline(g, "Gaussian Degree Distribution"; total_grains=N*20)
PlotResults("./Plots/Gaussian.png", GaussianDegreeDistribution)
SaveResultsToFile(file_to_write_results, GaussianDegreeDistribution; io_mode="w")

append!(AllStruct, [GaussianDegreeDistribution])

#### Uniform Degree Distribution
print("\rSimulating Uniform Degree Distribution Graph                                                               ")
N = 100000
degree_sequence = abs.(Int.(floor.(rand(1:40, N))))
if sum(degree_sequence) % 2 != 0
    degree_sequence[rand(1:N)] += 1
end
g = random_configuration_model(N, degree_sequence)
UniformDegreeDistribution = CompleteSandpilePipeline(g, "Uniform Degree Distribution"; total_grains=N*20)
PlotResults("./Plots/Uniform.png", UniformDegreeDistribution)
SaveResultsToFile(file_to_write_results, UniformDegreeDistribution)

append!(AllStruct, [UniformDegreeDistribution])

#### Power Degree Distribution Gamma 2
print("\rSimulating Scale Free Graph, γ = 2                                                               ")
N = 100000
g = static_scale_free(N, N*2, 2)
PowerDegreeDistribution2 = CompleteSandpilePipeline(g, "Power Degree Distribution Gamma 2"; total_grains=N*20)
PlotResults("./Plots/ScaleFree2.png", PowerDegreeDistribution2)
SaveResultsToFile(file_to_write_results, PowerDegreeDistribution2)

append!(AllStruct, [PowerDegreeDistribution2])

#### Power Degree Distribution Gamma 2.5
print("\rSimulating Scale Free Graph, γ = 2.5                                                               ")
N = 100000
g = static_scale_free(N, N*2, 2.5)
PowerDegreeDistribution2_5 = CompleteSandpilePipeline(g, "Power Degree Distribution Gamma 2.5"; total_grains=N*20)
PlotResults("./Plots/ScaleFree25.png", PowerDegreeDistribution2_5)
SaveResultsToFile(file_to_write_results, PowerDegreeDistribution2_5)

append!(AllStruct, [PowerDegreeDistribution2_5])

#### Power Degree Distribution Gamma 3
print("\rSimulating Scale Free Graph, γ = 3                                                               ")
N = 100000
g = static_scale_free(N, N*2, 3)
PowerDegreeDistribution3 = CompleteSandpilePipeline(g, "Power Degree Distribution Gamma 3"; total_grains=N*20)
PlotResults("./Plots/ScaleFree3.png", PowerDegreeDistribution3)
SaveResultsToFile(file_to_write_results, PowerDegreeDistribution3)

append!(AllStruct, [PowerDegreeDistribution3])

#### Power Degree Distribution Gamma 4
print("\rSimulating Scale Free Graph, γ = 4                                                               ")
N = 100000
g = static_scale_free(N, N*2, 4)
PowerDegreeDistribution4 = CompleteSandpilePipeline(g, "Power Degree Distribution Gamma 4"; total_grains=N*20)
PlotResults("./Plots/ScaleFree4.png", PowerDegreeDistribution4)
SaveResultsToFile(file_to_write_results, PowerDegreeDistribution4)

append!(AllStruct, [PowerDegreeDistribution4])

#### Erdos Renyi Network
print("\rSimulating Erdos Renyi Graph                                                               ")
N = 100000
g = erdos_renyi(N, 0.0002)
ErdosRenyi = CompleteSandpilePipeline(g, "Erdos Renyi Network"; total_grains=N*20)
PlotResults("./Plots/ErdosRenyi.png", ErdosRenyi)
SaveResultsToFile(file_to_write_results, ErdosRenyi)

append!(AllStruct, [ErdosRenyi])

#### Barabasi Albert Network
print("\rSimulating Barabasi Albert Graph                                                               ")
N = 100000
g = barabasi_albert(N, 10)
BarabasiAlbert = CompleteSandpilePipeline(g, "Barabasi Albert Network"; total_grains=N*20)
PlotResults("./Plots/BarabasiAlbert.png", BarabasiAlbert)
SaveResultsToFile(file_to_write_results, BarabasiAlbert)

append!(AllStruct, [BarabasiAlbert])

TotalPlot(AllStruct)

print("\rTask 15 Completed                                                                                 ")
