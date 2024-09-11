using GraphPlot, SimpleWeightedGraphs, Graphs, Plots, DataFrames, CSV, GraphDataFrameBridge, Statistics, Glob, Pkg

function SetUpEnviroment()
    current_dir = pwd()
    if current_dir[end-11:end] == "PoCS_Project"
        Pkg.activate("./")
        Pkg.instantiate()
    elseif current_dir[end-7:end] == "task_44"
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
        cd("code/task_44")
    end
end


SetUpEnviroment()
gr()

# Relevant Functions

# Function to extract the graph for a specific country from edge_list
function extract_country_graph(edge_list)
    country_edges = CSV.read(edge_list, DataFrame)[:, [1,2,5]]
    country_graph = MetaDiGraph(country_edges, :nodeID_from, :nodeID_to, weight=:scaled_sci)
    return country_graph, country_edges
end

# Function to analyze the graph for a list of countries, computing relevant metrics
function analyze_country_graph(countries_edges_list)
    results = DataFrame(
        :ISO => String[],
        :total_SCI => Float64[],
        :total_edges => Int64[], 
        :total_nodes => Int64[], 
        :average_indegree => Float64[],
        :average_outdegree => Float64[],
        :average_SCI => Float64[],
        :std_SCI => Float64[],
        :global_clustering => Float64[],
        :average_betweenness => Float64[],
        :max_betweenness => Float64[], 
        :average_closeness => Float64[],
        :max_closeness => Float64[])
    for edge_list in countries_edges_list
        print("\rAnalyzing $(edge_list[end-6:end-4]), $(length(countries_edges_list) - findfirst(isequal(edge_list), countries_edges_list)) countries left              ")
        ISO = edge_list[end-6:end-4]
        country_graph, country_edges = extract_country_graph(edge_list)
        total_SCI = sum(country_edges.scaled_sci)
        total_edges = ne(country_graph)
        total_nodes = nv(country_graph)
        average_indegree = mean(indegree(country_graph))
        average_outdegree = mean(outdegree(country_graph))
        average_SCI = mean(country_edges.scaled_sci)
        std_SCI = std(country_edges.scaled_sci)
        global_clustering = global_clustering_coefficient(country_graph)
        beetwenness = betweenness_centrality(country_graph)
        average_betweenness = mean(beetwenness)
        max_betweenness = maximum(beetwenness)
        closeness = closeness_centrality(country_graph)
        average_closeness = mean(closeness)
        max_closeness = maximum(closeness)
        push!(results, [ISO, total_SCI, total_edges, total_nodes, average_indegree, average_outdegree, average_SCI, std_SCI, global_clustering, average_betweenness, max_betweenness, average_closeness, max_closeness])
    end
    return results
end

# Some plot options
graph_plot_options = Dict(:EDGELINEWIDTH => 0.1, :arrowlengthfrac => 0.02, :plot_size => (12cm,12cm));
plot_options = Dict(:legendfontsize => 20, :yguidefontsize=>25,:xguidefontsize=>25, :titlefontsize=>25,:xtickfontsize=>18,:ytickfontsize=>18, :markersize => 7, :margin=>1Plots.cm)

# Load node_list
node_list = CSV.read("Output/node_list.csv", DataFrame)
# Extract all the country codes
countries_edges_list = glob("./Output/edge_list_*.csv")
print("There are ", length(countries_edges_list), " countries in the dataset.") # There are 200 countries in the dataset

# Extracting the graph for every single country
graph_data = analyze_country_graph(countries_edges_list) # This is really really long, provided there are the pre-computed results

graph_data = CSV.read("Data/graph_data.csv", DataFrame) # Load the pre-computed results
graph_data = graph_data[graph_data.total_nodes .> 1, :]; # Remove countries with only one node




# Total Nodes vs Total Edges
scatter(graph_data.total_nodes, graph_data.total_edges, yscale=:log10, xlabel="Total Nodes", ylabel="Total Edges", title="Total Nodes vs Total Edges", size=(1000,800), label = "Data"; plot_options...)
plot!(1:700, (1:700).^2, color=:red, label="Max Edges", linewidth=2, legend=:bottomright; plot_options...)
savefig("Plots/TotalNodesVsTotalEdges.png")

# Total Nodes vs Average Betweenness
scatter(graph_data.total_nodes, graph_data.average_betweenness, xlabel="Total Nodes", xscale=:log10,ylabel="Average Betweenness", title="Total Nodes vs Average Betweenness", size=(1000,800), legend=false; plot_options...)
savefig("Plots/TotalNodesVsAverageBetweenness.png")

# Total Nodes vs Max Betweenness
scatter(graph_data.total_nodes, graph_data.max_betweenness, xlabel="Total Nodes", xscale=:log10,ylabel="Max Betweenness", title="Total Nodes vs Max Betweenness", size=(1000,800), legend=false; plot_options...)
savefig("Plots/TotalNodesVsMaxBetweenness.png")

# Total Nodes vs Average Closeness
scatter(graph_data.total_nodes, graph_data.average_closeness, xlabel="Total Nodes", xscale=:log10,yscale=:log10,ylabel="Average Closeness", title="Total Nodes vs Average Closeness", size=(1000,800), legend=false; plot_options...)
savefig("Plots/TotalNodesVsAverageCloseness.png")

# Total Nodes vs Max Closeness
scatter(graph_data.total_nodes, graph_data.total_SCI, xlabel="Total Nodes", ylabel="Total SCI",yscale=:log10,xscale=:log10, title="Total Nodes vs Total SCI", size=(1000,800), legend=false; plot_options...)
savefig("Plots/TotalNodesVsTotalSCI.png")