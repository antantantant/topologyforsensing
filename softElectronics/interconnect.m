%%% Main code for interconnect design %%%

%% set parameters
E = 1; % Young's modulus
mu = 0.3; % Poisson's ratio
nx = 10; % horizontal number of elements
ny = 10; % vertical number of elements

%% create topology
x = @(p)par2top(p, nx, ny);

%% K assembly
K = @(x)kassemble(x, E, mu);

%% set boundary conditions
% Need Xu's input

%% simulation


