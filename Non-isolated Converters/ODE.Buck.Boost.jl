# This design is based on the book by Perez-Ramirez and Beristain-Jimenez 2016
# If this code helps you please consider on cite the original work:

# Pérez Ramírez, J., & Beristain Jimenez, J. A. (2016).
# Electrónica de potencia: modelado y control de convertidores cd-cd.

# Copyright 2022
# Author Erick Moreno Negrete
# Visit my Github page
# https://github.com/erickmone

# ----------------------------------------------------------------
#                     Buck-Boost Converter:
# ----------------------------------------------------------------

using OrdinaryDiffEq 
using ControlSystems
using OrdinaryDiffEq
using LinearAlgebra
using LaTeXStrings

using Plots
theme(:dao)

# Parameters of the converter

E  = 48;
R  = 10;
L  = 100e-6;
C  = 100e-6;
D  = 0.5;
Dₚ = 1-D;
x₀ = [0; 0];

# Dynamical Equations

function dxdt(dx,x,p,t)
    global L,C,E,R,Dₚ;
    dx[1] = (1/L)*(D*E - Dₚ*x[2]);
    dx[2] = (1/C)*(Dₚ*x[1] - (1/R)*x[2]);
end 

# Simulate using ODE Solver

t = (0,1)
prob = ODEProblem(dxdt,x₀,t)
sol = solve(prob, Tsit5())
 
# Plotting results

p = plot(sol,
         lw=2,
         tspan=(0,0.008),
         layout=(2,1), 
         xlabel=["" "Time \$t\$ [s]"],
         ylabel=["Current [A]" "Voltage [V]"],
         label = ["\$i_L(t)\$" "\$v_C(t)\$"], 
         title = ["Response of the Buck-Boost Converter" ""],
         color = [:red :blue]
         )
