using Pkg
current_dir = pwd()
if current_dir[end-11:end] == "PoCS_Project"
    Pkg.activate("./")
    Pkg.instantiate()
elseif current_dir[end-3:end] == "code"
    Pkg.activate("../")
    Pkg.instantiate()
else
    print("Please run script from the main directory or the script parent directory")
end
