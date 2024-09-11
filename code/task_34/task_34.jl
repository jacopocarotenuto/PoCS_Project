using GraphPlot, Graphs, Plots, Distributions, Pkg

function SetUpEnviroment()
    current_dir = pwd()
    if current_dir[end-11:end] == "PoCS_Project"
        Pkg.activate("./")
        Pkg.instantiate()
    elseif current_dir[end-7:end] == "task_34"
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
        cd("code/task_34")
    end
end


SetUpEnviroment()

gr()

# Rounds to record for each simulation
rounds_to_record_NS = [1,100,1000,10000,20000]
rounds_to_record_SP =  [1,100,1000,10000,100000]


# Data structure to store the results of the simulations
struct UG_Stats
    graph_type::String
    TotalSimulations::Int
    N::Int
    p_A_NS::Array{Float64, 2}
    q_A_NS::Array{Float64, 2}
    p_B_NS::Array{Float64, 2}
    q_B_NS::Array{Float64, 2}
    p_C_NS::Array{Float64, 2}
    q_C_NS::Array{Float64, 2}
    p_A_SP::Array{Float64, 2}
    q_A_SP::Array{Float64, 2}
    p_B_SP::Array{Float64, 2}
    q_B_SP::Array{Float64, 2}
    p_C_SP::Array{Float64, 2}
    q_C_SP::Array{Float64, 2}
end

# Function to perform one round of the Ultimatum Game between two players
function OneRoundUG!(Player_i, Player_j, strategies, payoffs)
    ## Player_i is the proposer
    offer = strategies[Player_i, 1]
    response = strategies[Player_j, 2]
    if offer >= response
        payoffs[Player_i] += 1 - offer
        payoffs[Player_j] += offer
    end

    ## Player_i is the respondent
    offer = strategies[Player_j, 1]
    response = strategies[Player_i, 2]
    if offer >= response
        payoffs[Player_j] += 1 - offer
        payoffs[Player_i] += offer
    end

    return nothing
end

# Function to perform the Ultimatum Game on the entire graph
function UG_on_entire_graph!(list_of_neighbors, strategies::Array{Float64, 2}, payoffs)
    for i in 1:N::Int
        for j in list_of_neighbors[i]::Vector{Int}
            OneRoundUG!(i, j, strategies, payoffs)
        end
    end
    return nothing
end

# Function to update the strategies of the players after each cycle
function UpdateStrategies!(strategies::Array{Float64, 2}, list_of_neighbors::Array{Array{Int64, 1}, 1}, payoffs::Array{Float64, 1}; modality = "Natural Selection")
    if modality == "Natural Selection"
        # Choose a random neighbor to copy the strategy from with a probability proportional to the difference in payoffs, vectorized version
        chosen_neighbors = rand.(list_of_neighbors)
        chosen_neighbors_payoffs = payoffs[chosen_neighbors]
        p = (chosen_neighbors_payoffs .- payoffs) ./ (2*maximum([length.(list_of_neighbors) length.(list_of_neighbors[chosen_neighbors])], dims = 2))
        strategies .= ifelse.(rand(length(list_of_neighbors)) .< p, strategies[chosen_neighbors,:], strategies)
    end

    if modality == "Social Penalty"
        # The player with the lowest payoff and its neighbors die and are replaced by random strategies (of the same overall type)

        lowest_id = argmin(payoffs)
        if strategies[lowest_id, 1] == strategies[lowest_id, 2] # Checking the type of player like this is a slight abuse, but probability of being wrong is very low
            temp = rand(0:0.05:1)
            strategies[lowest_id, :] = [temp, temp]
            for i in list_of_neighbors[lowest_id]::Vector{Int64}
                temp = rand(0:0.05:1)
                strategies[i, :] = [temp, temp]
            end
        elseif strategies[lowest_id, 1] == (1 - strategies[lowest_id, 2])
            temp = rand(0:0.05:1)
            strategies[lowest_id, :] = [temp, 1-temp]
            for i in list_of_neighbors[lowest_id]::Vector{Int64}
                temp = rand(0:0.05:1)
                strategies[i, :] = [temp, 1-temp]
            end
        else
            strategies[lowest_id, :] = rand(0:0.05:1, 2)
            for i in list_of_neighbors[lowest_id]::Vector{Int64}
                strategies[i, :] = rand(0:0.05:1, 2)
            end
        end
    end
    return nothing

end

# Function to play the Ultimatum Game on a graph for a given number of rounds
function PlayUG(g; rounds_to_play = 1000, rounds_to_record = [1], modality = "Natural Selection", strategy_choice = "Type C")
    if strategy_choice == "Type C" # P and Q are indipendent
        strategies = rand(0:0.05:1, N, 2)
    end
    if strategy_choice == "Type B" # P is random and Q = 1 - P
        p = rand(0:0.05:1, N)
        strategies = hcat(p, 1 .- p)
    end
    if strategy_choice == "Type A" # P and Q are the same
        p = rand(0:0.05:1, N)
        strategies = hcat(p, p)
    end

    list_of_neighbors = [neighbors(g, i) for i in 1:N] # Precompute the neighbors of each node


    ## Add self loops for nodes with no neighbors, to avoid problems with random picking of neighbors.
    findall(length.(list_of_neighbors) .==0)
    for i in findall(length.(list_of_neighbors) .==0)
        list_of_neighbors[i] = [i]
    end

    p_data = zeros(length(rounds_to_record), 21)
    q_data = zeros(length(rounds_to_record), 21)
    payoffs = zeros(N)
    counter = 1
    for round in 1:rounds_to_play
        if round in rounds_to_record

            h_p = fit(Histogram, strategies[:,1], -0.025:0.05:1.025)
            h_q = fit(Histogram, strategies[:,2], -0.025:0.05:1.025)
            p_data[counter, :] = h_p.weights / sum(h_p.weights)
            q_data[counter, :] = h_q.weights / sum(h_q.weights)
            counter += 1
        end

        payoffs .= 0
        UG_on_entire_graph!(list_of_neighbors, strategies, payoffs)
        UpdateStrategies!(strategies, list_of_neighbors, payoffs, modality = modality)

    end

    return p_data, q_data
end

# Function to plot the results of a simulation
function PlotData(data,rounds_to_record, title="Histogram of", xlabel="YOU FORGOT";file_name="", legend=:topright)
    scatter(0:0.05:1, data[1,:], label = "Round 1", xlabel = xlabel, ylabel = "Probability", title = title, legend = legend; palette = :Dark2_5, legendfontsize = 20, yguidefontsize=25,xguidefontsize=25, titlefontsize=25,xtickfontsize=18,ytickfontsize=18, markersize = 7)
    for i in eachindex(rounds_to_record)[2:end]
        scatter!(0:0.05:1, data[i,:], label = "Round $(rounds_to_record[i])", markersize = 7)
    end

    for i in eachindex(rounds_to_record)
        plot!(0:0.05:1, data[i,:], label = "", linewidth=4)
    end
    plot1 = plot!(size=(2880,1080), margin = 2Plots.cm, dpi = 600, markerstrokewidth = 0.6)

    if file_name != ""
        savefig(file_name)
    end

    return plot1
end

# Function to plot P and Q for a given simulation
function PlotPQ(p_data, q_data, rounds_to_record, file_name)
    
    plot1 = PlotData(p_data,rounds_to_record, "Histogram of p", "p")
    plot2 = PlotData(q_data,rounds_to_record, "Histogram of q", "q")

    plot(plot1, plot2, layout = (1,2), size=(2880,1080), margin = 2Plots.cm, dpi = 600, markerstrokewidth = 0.6)
    savefig(file_name)
end

# Function to perform multiple simulations of the Ultimatum Game withe every type of player and both update rules
function MultipleUG(total_simulations, graph_type, N)

    rounds_to_record_NS = [1,100,1000,10000,20000]
    rounds_to_record_SP =  [1,100,1000,10000,100000]

    p_data_A_NS = zeros(length(rounds_to_record_NS), 21)
    q_data_A_NS = zeros(length(rounds_to_record_NS), 21)
    p_data_B_NS = zeros(length(rounds_to_record_NS), 21)
    q_data_B_NS = zeros(length(rounds_to_record_NS), 21)
    p_data_C_NS = zeros(length(rounds_to_record_NS), 21)
    q_data_C_NS = zeros(length(rounds_to_record_NS), 21)

    # First we do "Natural Selection"
    for i in 1:total_simulations
        print("\rSimulation $i out of $total_simulations of Natural Selection")
        g = GenerateGraph(graph_type, N)
        p_data_i_A, q_data_i_A = PlayUG(g, rounds_to_play = 2e4, rounds_to_record = rounds_to_record_NS, modality = "Natural Selection", strategy_choice = "Type A")
        p_data_i_B, q_data_i_B = PlayUG(g, rounds_to_play = 2e4, rounds_to_record = rounds_to_record_NS, modality = "Natural Selection", strategy_choice = "Type B")
        p_data_i_C, q_data_i_C = PlayUG(g, rounds_to_play = 2e4, rounds_to_record = rounds_to_record_NS, modality = "Natural Selection", strategy_choice = "Type C")
        p_data_A_NS += p_data_i_A
        q_data_A_NS += q_data_i_A
        p_data_B_NS += p_data_i_B
        q_data_B_NS += q_data_i_B
        p_data_C_NS += p_data_i_C
        q_data_C_NS += q_data_i_C
    end
    p_data_A_NS /= total_simulations
    q_data_A_NS /= total_simulations
    p_data_B_NS /= total_simulations
    q_data_B_NS /= total_simulations
    p_data_C_NS /= total_simulations
    q_data_C_NS /= total_simulations

    ## Then we do "Social Penalty"
    p_data_A_SP = zeros(length(rounds_to_record_SP), 21)
    q_data_A_SP = zeros(length(rounds_to_record_SP), 21)
    p_data_B_SP = zeros(length(rounds_to_record_SP), 21)
    q_data_B_SP = zeros(length(rounds_to_record_SP), 21)
    p_data_C_SP = zeros(length(rounds_to_record_SP), 21)
    q_data_C_SP = zeros(length(rounds_to_record_SP), 21)

    for i in 1:total_simulations
        print("\rSimulation $i out of $total_simulations of Social Penalty                ")
        g = GenerateGraph(graph_type, N)
        p_data_i_A, q_data_i_A = PlayUG(g, rounds_to_play = 1e5, rounds_to_record = rounds_to_record_SP, modality = "Social Penalty", strategy_choice = "Type A")
        p_data_i_B, q_data_i_B = PlayUG(g, rounds_to_play = 1e5, rounds_to_record = rounds_to_record_SP, modality = "Social Penalty", strategy_choice = "Type B")
        p_data_i_C, q_data_i_C = PlayUG(g, rounds_to_play = 1e5, rounds_to_record = rounds_to_record_SP, modality = "Social Penalty", strategy_choice = "Type C")
        p_data_A_SP += p_data_i_A
        q_data_A_SP += q_data_i_A
        p_data_B_SP += p_data_i_B
        q_data_B_SP += q_data_i_B
        p_data_C_SP += p_data_i_C
        q_data_C_SP += q_data_i_C
    end
    p_data_A_SP /= total_simulations
    q_data_A_SP /= total_simulations
    p_data_B_SP /= total_simulations
    q_data_B_SP /= total_simulations
    p_data_C_SP /= total_simulations
    q_data_C_SP /= total_simulations

    return UG_Stats(graph_type, total_simulations, N, p_data_A_NS, q_data_A_NS, p_data_B_NS, q_data_B_NS, p_data_C_NS, q_data_C_NS, p_data_A_SP, q_data_A_SP, p_data_B_SP, q_data_B_SP, p_data_C_SP, q_data_C_SP)

end

# Function to generate the graph
function GenerateGraph(graph_type, N)
    if graph_type == "Erdos Renyi"
        return erdos_renyi(N, 4/N)
    elseif graph_type == "Barabasi Albert"
        return barabasi_albert(N, 4)
    elseif graph_type == "Scale Free"
        return static_scale_free(N, N*4, 3)
    end
end

# Function to plot all the results of a simulation
function PlotAllPQ(UG)
    PlotPQ(UG.p_A_NS, UG.q_A_NS, rounds_to_record_NS, "Plots/"*UG.graph_type*"_TypeA_NS.png")
    PlotPQ(UG.p_B_NS, UG.q_B_NS, rounds_to_record_NS, "Plots/"*UG.graph_type*"_TypeB_NS.png")
    PlotPQ(UG.p_C_NS, UG.q_C_NS, rounds_to_record_NS, "Plots/"*UG.graph_type*"_TypeC_NS.png")

    PlotPQ(UG.p_A_SP, UG.q_A_SP, rounds_to_record_SP, "Plots/"*UG.graph_type*"_TypeA_SP.png")
    PlotPQ(UG.p_B_SP, UG.q_B_SP, rounds_to_record_SP, "Plots/"*UG.graph_type*"_TypeB_SP.png")
    PlotPQ(UG.p_C_SP, UG.q_C_SP, rounds_to_record_SP, "Plots/"*UG.graph_type*"_TypeC_SP.png")
end

# Global Variables
N = 10000

### WARNING: This will take a long time to run

## Erdos Renyi
ErdosRenyi = MultipleUG(10, "Erdos Renyi", N)

## Barabasi Albert
BarabasiAlbert = MultipleUG(10, "Barabasi Albert", N)

## Scale Free
ScaleFree = MultipleUG(10, "Scale Free", N)


####### Relevant Plots #######

# More plots can be done with the PlotData function, it's very versatile

# This plot could be done for every type of player
ER_p_A_NS = PlotData(ErdosRenyi.p_A_NS,rounds_to_record_NS, "ER (Natural Selection)", "p")
BA_p_A_NS = PlotData(BarabasiAlbert.p_A_NS,rounds_to_record_NS, "BA (Natural Selection)", "p")
SF_p_A_NS = PlotData(ScaleFree.p_A_NS,rounds_to_record_NS, "SF (Natural Selection)", "p")
ER_p_A_SP = PlotData(ErdosRenyi.p_A_SP,rounds_to_record_SP, "ER (Social Penalty)", "p"; legend=:bottom)
BA_p_A_SP = PlotData(BarabasiAlbert.p_A_SP,rounds_to_record_SP, "BA (Social Penalty)", "p"; legend=:bottom)
SF_p_A_SP = PlotData(ScaleFree.p_A_SP,rounds_to_record_SP, "SF (Social Penalty)", "p"; legend=:bottom)
plot(ER_p_A_NS, BA_p_A_NS, SF_p_A_NS, ER_p_A_SP, BA_p_A_SP, SF_p_A_SP, layout = (2,3), size=(4320,2160), margin = 3Plots.cm, dpi = 600, markerstrokewidth = 0.6)
savefig("Plots/TypeA_AllGraphs_AllRules.png")


# This plot could be done for every type of graph
SF_p_A_NS = PlotData(ScaleFree.p_A_NS,rounds_to_record_NS, "Only Type A", "p")
SF_p_B_NS = PlotData(ScaleFree.p_B_NS,rounds_to_record_NS, "Only Type B", "p")
SF_p_C_NS = PlotData(ScaleFree.p_C_NS,rounds_to_record_NS, "Only Type C", "p")
SF_p_A_SP = PlotData(ScaleFree.p_A_SP,rounds_to_record_SP, "Only Type A", "p"; legend=:bottom)
SF_p_B_SP = PlotData(ScaleFree.p_B_SP,rounds_to_record_SP, "Only Type B", "p"; legend=:bottom)
SF_p_C_SP = PlotData(ScaleFree.p_C_SP,rounds_to_record_SP, "Only Type C", "p"; legend=:bottom)
plot(SF_p_A_NS, SF_p_B_NS, SF_p_C_NS, SF_p_A_SP, SF_p_B_SP, SF_p_C_SP, layout = (2,3), size=(4320,2160), margin = 3Plots.cm, dpi = 600, markerstrokewidth = 0.6)
savefig("Plots/AllTypes_SF_AllRules.png")

