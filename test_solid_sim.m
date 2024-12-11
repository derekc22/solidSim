clear;clc
close all

rho = 2000;
T0 = 2500;
k = 0.002;
n = 0.5;
gamma = 1.2;
Mbar = 10;
ri = 0.1;
rf = 0.2;
L = 2;
Astar = pi*ri^2;



[tList, P0List] = solidSim(rho, T0, k, n, gamma, Mbar, ri, rf, L, Astar)