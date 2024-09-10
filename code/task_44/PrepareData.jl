using DataFrames, CSV, Countries, GeoIO, GeoStats

function FixWorkingDirectory()
    # Get the current directory
    current_dir = pwd()
    # Check if the current directory is the one with the script
    if current_dir[end-7:end] != "task_44"
        # If not change to the directory with the script
        cd("code/task_44")
    end
end


# Function to convert code to ISO3 format based on the type
function from_code_to_ISO3(code, type)
    if type == "country"
        return get_country(code).alpha3
    elseif type == "nuts3"
        if code[1:2] == "EL"
            return "GR"
        elseif code[1:2] == "UK"
            return "GB"
        end
        return get_country(code[1:2]).alpha3
    elseif type == "nuts2"
        if code[1:2] == "EL"
            return "GR"
        elseif code[1:2] == "UK"
            return "GB"
        end
        return get_country(code[1:2]).alpha3
    elseif (type == "gadm1") | (type == "gadm2")
        if code[1:3] == "XKO"
            return code[1:3]
        end
        return get_country(code[1:3]).alpha3
    elseif type == "county"
        return "USA"
    end
end

# Function to get the country name based on the code
function get_country_name(code)
    if code == "XKO"
        return "Kosovo"
    else
        return get_country(code).name
    end
end

# Function to convert ISO code to Facebook standard used
function CodeToFacebookStandard(ISO, ID1 = " ", ID2 = " ")
    code = ISO
    if ID1 != " "
        code *= ID1
    end
    if ID2 != " "
        code *= "_" * ID2
    end
    return code
end

# Function to save all countries edge_list on different files
function SaveAllCountriesEdgeList(edge_list)
    # DO NOT SAVE USA
    # Extract all the country codes
    all_country_codes = unique(edge_list.country_ISO3)
    for code in all_country_codes
        if code == "USA"
            continue
        end
        country_edges = edge_list[edge_list.country_ISO3 .== code, :]
        CSV.write("Output/edge_list_$code.csv", country_edges)
    end
end

FixWorkingDirectory()
################################### CREATION OF EDGE LIST #######################################
# Load the data
data = CSV.read("Data/gadm1_nuts3_counties-gadm1_nuts3_counties - FB Social Connectedness Index - October 2021.tsv", DataFrame)

# Create the dictionaries for the IDs
unique_codes = unique(data[:, 1])
GADM3_to_ID = Dict(unique_codes[i] => i for i in eachindex(unique_codes))

# Load the data for Facebook levels
FB_LV = CSV.read("Data/gadm1_nuts3_counties_levels.csv", DataFrame)
FB_LV = transform(FB_LV, [:key, :level] => ByRow((x,y) -> from_code_to_ISO3(x, y)) => :country) # Transform the data adding the country name based on code

# Create the dictionary for code to ISO3 mapping
code_to_ISO3 = Dict(FB_LV[:, 1] .=> FB_LV[:, 3])

# Get all the countries
all_countries_codes = DataFrame(all_countries())
all_countries_codes = select(all_countries_codes, 1:3) # Get the names of the countries

# Add the country names to the dictionary
for i in 1:size(all_countries_codes, 1)
    code_to_ISO3[all_countries_codes[i, 2]] = all_countries_codes[i, 2]
end

# Add Kosovo and Mauritius to the dictionary
code_to_ISO3["XKO"] = "XKO"
code_to_ISO3["MUS1"] = "MUS"
# The above code is necessary as the provided data has some inconsistencies in the codes, and not all code are present in the levels file.


# Transform the data by adding the country (from) name based on code
EdgeList = transform(data, :user_loc => ByRow(x -> code_to_ISO3[x]) => :country_ISO3)
# Transform the data by adding the country (to) name based on code
EdgeList = transform(EdgeList, :fr_loc => ByRow(x -> code_to_ISO3[x]) => :country_ISO3_to)

################################### IF WE DONT NEED INTER COUNTRY CONNECTIONS #######################################
# Filter the data to keep only the connections within the same country
EdgeList = filter(row -> row[:country_ISO3] == row[:country_ISO3_to], EdgeList)
#####################################################################################################################

# Transform the data by adding the country name based on code
EdgeList = transform(EdgeList, :country_ISO3 => ByRow(x -> get_country_name(x)) => :country_name)
# Transform the data by adding the node ID (from) based on code
EdgeList = transform(EdgeList, :user_loc => ByRow(x -> GADM3_to_ID[x]) => :nodeID_from)
# Transform the data by adding the node ID (to) based on code
EdgeList = transform(EdgeList, :fr_loc => ByRow(x -> GADM3_to_ID[x]) => :nodeID_to)
# Select the columns we need
EdgeList = select(EdgeList, [:nodeID_from, :nodeID_to, :country_name, :country_ISO3, :scaled_sci])

SaveAllCountriesEdgeList(EdgeList)

################################### CREATION OF NODE LIST #######################################

# Create a dictionary to store code to coordinates mapping
CODE_to_COORD = Dict()

#### NUTS DATA ####
NUTS = DataFrame(GeoIO.load("Data/NUTS_RG_60M_2016_4326.geojson")) # Load the data
NUTS3 = filter(row -> row[:LEVL_CODE] == 3, NUTS) # Filter the data to keep only NUTS3 regions
NUTS3_Coord = [centroid(i) for i in NUTS3.geometry]

# Add NUTS3 code to coordinates mapping to the dictionary
for i in eachindex(NUTS3.NUTS_ID)
    CODE_to_COORD[NUTS3[i,:].NUTS_ID] = NUTS3_Coord[i]
end

#### GADM DATA ####
GADM0_DATA = DataFrame(GeoIO.load("Data/GADM_DATA/gadm28_adm0.shp"))
GADM0_DATA = select(GADM0_DATA, [:ISO, :geometry])
GADM0_DATA = transform(GADM0_DATA, [:ISO] => ByRow((ISO) -> CodeToFacebookStandard(ISO)) => :CODE)
GADM1_DATA = DataFrame(GeoIO.load("Data/GADM_DATA/gadm28_adm1.shp"))
GADM1_DATA = select(GADM1_DATA, [:ISO, :ID_1, :geometry])
GADM1_DATA = transform(GADM1_DATA, [:ISO, :ID_1] => ByRow((ISO, ID1) -> CodeToFacebookStandard(ISO, string(ID1))) => :CODE)
GADM2_DATA = DataFrame(GeoIO.load("Data/GADM_DATA/gadm28_adm2.shp"))
GADM2_DATA = select(GADM2_DATA, [:ISO, :ID_1, :ID_2, :geometry])
GADM2_DATA = transform(GADM2_DATA, [:ISO, :ID_1, :ID_2] => ByRow((ISO, ID1, ID2) -> CodeToFacebookStandard(ISO, string(ID1), string(ID2))) => :CODE)

# Filter GADM data based on Facebook levels
FB_GADM0 = filter(row -> row[:level] == "country", FB_LV)
FB_GADM1 = filter(row -> row[:level] == "gadm1", FB_LV)
FB_GADM2 = filter(row -> row[:level] == "gadm2", FB_LV)
GADM0_DATA = filter(row -> row.CODE in FB_GADM0.key, GADM0_DATA)
GADM1_DATA = filter(row -> row.CODE in FB_GADM1.key, GADM1_DATA)
GADM2_DATA = filter(row -> row.CODE in FB_GADM2.key, GADM2_DATA)

# Add GADM0 code to coordinates mapping to the dictionary
for i in eachindex(GADM0_DATA.CODE)
    CODE_to_COORD[GADM0_DATA[i,:].CODE] = centroid(first(rings(GADM0_DATA[i,:].geometry)))
end

# Add GADM1 code to coordinates mapping to the dictionary
for i in eachindex(GADM1_DATA.CODE)
    CODE_to_COORD[GADM1_DATA[i,:].CODE] = centroid(first(rings(GADM1_DATA[i,:].geometry)))
end

# Add GADM2 code to coordinates mapping to the dictionary
for i in eachindex(GADM2_DATA.CODE)
    CODE_to_COORD[GADM2_DATA[i,:].CODE] = centroid(first(rings(GADM2_DATA[i,:].geometry)))
end

# Get the IDs from the previously created dictionary
NodeList = DataFrame("label" => collect(keys(GADM3_to_ID)), "ID" => collect(values(GADM3_to_ID)))
LAT = zeros(length(NodeList.label))
LON = zeros(length(NodeList.label))

# Add latitude and longitude coordinates to the NodeList
for i in eachindex(NodeList.label)
    if haskey(CODE_to_COORD, NodeList.label[i])
        LAT[i] = GeoStats.coords(CODE_to_COORD[NodeList.label[i]]).lat
        LON[i] = GeoStats.coords(CODE_to_COORD[NodeList.label[i]]).lon
    else
        LAT[i] = 0
        LON[i] = 0
    end
end

NodeList = hcat(NodeList, DataFrame("LAT" => LAT, "LON" => LON))

# Save to CSV
CSV.write("Output/node_list.csv", NodeList)