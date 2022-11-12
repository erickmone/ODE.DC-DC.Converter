# This design is based on the book by Perez-Ramirez and Beristain-Jimenez 2016
# If this code helps you please consider on cite the original work:

# Pérez Ramírez, J., & Beristain Jimenez, J. A. (2016).
# Electrónica de potencia: modelado y control de convertidores cd-cd.

# Copyright 2022
# Author Erick Moreno Negrete
# Visit my Github page
# https://github.com/erickmone

# ----------------------------------------------------------------
#                     Sepic Converter:
# ----------------------------------------------------------------

using OrdinaryDiffEq 
using ControlSystems
using OrdinaryDiffEq
using LinearAlgebra
using LaTeXStrings

using Plots
theme(:dao)

# Parameters of the converter

E  = 5;
R  = 5;
L1 = 40e-6;
L2 = 40e-6;
C1 = 10e-6;
C2 = 100e-6;
D  = 0.8;
Dₚ = 1-D;
x₀ = [0; 0; 0; 0];

# Dynamical Equations

function dxdt(dx,x,p,t)
    global L1,L2,C1,C2,E,R,Dₚ;
    dx[1] = (1/L1)*(E - Dₚ*x[3] - Dₚ*x[4]);
    dx[2] = (1/L2)*(D*x[3] - Dₚ*x[4]);
    dx[3] = (1/C1)*(Dₚ*x[1] - D*x[2]);
    dx[4] = (1/C2)*(Dₚ*x[1] + Dₚ*x[2] - (1/R)*x[4]);
end 

# Simulate using ODE Solver

t = (0,1)
prob = ODEProblem(dxdt,x₀,t)
sol = solve(prob, Tsit5())
 
# Plotting results

p = plot(sol,
         lw=2,
         tspan=(0,0.004),
         layout=grid(4,1), 
         size=(600, 800),
         xlabel=["" "" "" "Time \$t\$ [s]"],
         ylabel=["Current [A]" "Current [A]" "Voltage [V]" "Voltage [V]"],
         label = ["\$i_{L1}(t)\$" "\$i_{L2}(t)\$" "\$v_{C1}(t)\$" "\$v_{C2}(t)\$"], 
         title = ["Response of the Sepic Converter" "" "" ""],
         color = [:red :red :blue :blue]
         )