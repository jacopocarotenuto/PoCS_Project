# PoCS_Project
This repository contains the code for the completion of all the assigned task for the final project of the course "Physics of Complex Networks: Structure and Dynamics" at the University of Padova. The project was solo developed by me, [Jacopo Carotenuto](https://github.com/jacopocarotenuto) in the academic year 2023/2024 in [Julia](https://julialang.org/).

## Repository Structure
The repository is divided into 3 main folders:
- `code`: Contains all the code for the tasks.
- `data`: Contains all the resulting data from the "Data" task (Task 44).
- `latex`: Contains the LaTeX files for the final report.

Additionaly there are some files in the root directory:
- `README.md`: This file.
- `final_report.pdf`: The final report of the project.
- `LICENSE`: The license file.
- `.gitignore`: The gitignore file.
- `Project.toml`: The Julia environment files.

To run the code it is necessary to have Julia 1.10 installed on the machine. Before attempting to run the task code, to setup the environment run the following command from the root directory of the repository: `julia EnvironmentSetup.jl`.

## The Tasks
Inside the `code` folder there is a folder for each task and a `README.md` file with the instruction to run the code.

### Task 15
In this task the "Sandpile Model" on different types of networks is studied. This model is a type of dynamic process that can be simulated on a network that exhibits self-organized criticality (SOC). The model has been studied in different types of networks and some interesting relations have emerged, especially about the distribution of the avalanche size and duration. In this task the model will be simulated on different networks and calculate the corresponding distributions.

### Task 34

In this task the evolution of strategies according to different "survival" rules will be studied. The "Ultimatum Game" on different networks with different types of players and strategies will be considered.

### Task 44

In this task the goal was to extract a network for each country present in the Facebook Social Connectedness Index data. The data was collected by Facebook and it is based on the number of friendships between people in different countries.
> The Social Connectedness Index uses an anonymized snapshot of active Facebook users and their friendship networks to measure the intensity of social connectedness between locations. Users are assigned to locations based on their information and activity on Facebook, including the stated city on their Facebook profile, and device and connection information.


