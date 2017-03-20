%%%%%%%%% Calculate input observability for various structures %%%%%%%%%
clear;
close all;
% load structure.mat; % load xPhys

T = 1e1; % maximum time
nsteps = 1e1; % time interval

%% set up finite element method
% Set parameters
nelx = 10; % horizontal number of elements (left to right)
nely = 1; % vertical number of elements (top to down)
penal = 3; % polynomial order to define density-young's modulus relationship
E0 = 1; % young's modulus at density=1
Emin = 1e-9; % young's modulus at density=0, keep this small
nu = 0.3; % poisson's ratio

% element-wise stiffness matrix for a quadrilateral element (square in shape)
A11 = [12  3 -6 -3;  3 12  3  0; -6  3 12 -3; -3  0 -3 12]; %stiffness matrix
A12 = [-6 -3  0  3; -3 -6 -3 -6;  0 -3 -6  3;  3 -6  3 -6];
B11 = [-4  3 -2  9;  3 -4 -9  4; -2 -9 -4 -3;  9  4 -3 -4];
B12 = [ 2 -3  4 -9; -3  2  9 -2;  4  9  2  3; -9 -2  3  2];
KE = 1/(1-nu^2)/24*([A11 A12;A12' A11]+nu*[B11 B12;B12' B11]);

% element-wise consistent mass matrix for a quadrilateral element (square in shape)
ME = (4*eye(8) + [zeros(4),eye(4);eye(4),zeros(4)] + repmat(kron([0 1;1 0],2*eye(2)),2,2))/9;

% define material density
% this is now loaded from an existing file structure.mat
xPhys = ones(nely,nelx);
[nely,nelx] = size(xPhys);

% element-to-global assembly
nodenrs = reshape(1:(1+nelx)*(1+nely),1+nely,1+nelx);
edofVec = reshape(2*nodenrs(1:end-1,1:end-1)+1,nelx*nely,1);
edofMat = repmat(edofVec,1,8)+repmat([0 1 2*nely+[2 3 0 1] -2 -1],nelx*nely,1);
iK = reshape(kron(edofMat,ones(8,1))',64*nelx*nely,1);
jK = reshape(kron(edofMat,ones(1,8))',64*nelx*nely,1);
sK = reshape(KE(:)*(Emin+xPhys(:)'.^penal*(E0-Emin)),64*nelx*nely,1);
K = sparse(iK,jK,sK); K = (K+K')/2;
sM = reshape(ME(:)*xPhys(:)',64*nelx*nely,1);
M = sparse(iK,jK,sM); M = (M+M')/2;

% define loads and dofs
F = sparse(2*(nelx+1)*(nely+1),1,-1,2*(nely+1)*(nelx+1),1);
fixeddofs = [1:2*(nely+1)];
alldofs = [1:2*(nely+1)*(nelx+1)];
freedofs = setdiff(alldofs,fixeddofs);

% dynamic solution
p = size(freedofs,2); % DOF of the system, number of state is 2p


% prepare matrices for input (force) estimation
dt = T/nsteps;
Kb = K(freedofs,freedofs);
Mb = M(freedofs,freedofs);
Fb = F(freedofs);
A = [zeros(p), -Mb\Kb; eye(p), zeros(p)];
B = [inv(Mb); zeros(p)];

% set observer
% S = zeros(p,1); %Sa in the document
% S(2:2:p) = 1;
% S = diag(S);
num_observer = 5;
S = zeros(num_observer,p); 
S(:,randperm(p,num_observer))=eye(num_observer);

% convert continuous time to discrete time
AA = expm(A*dt);
BB = inv(A)*(AA-eye(2*p))*B;
% check AA and BB calculation, OK this is correct
y = zeros(nsteps,2*p);
for i = 2:size(y,1)
    y(i,:) = (AA*y(i-1,:)' + BB*Fb)';
end
   
% calculate matrices for observability
D = S*inv(Mb);
C1 = -S*inv(Mb)*Kb;
B1 = inv(Mb);
A1 = -inv(Mb)*Kb;
C = [zeros(num_observer,p), C1];

% M2 follows corollary 5 eq 35
% M1 follows eq 13b

% M2 = [];
% rankM2 = zeros(floor((nsteps-1)/2)+2,1);
% for i = 0:floor((nsteps-1)/2)+1
%     M2 = [M2; S*A1^i*B1];
%     rankM2(i+1) = rank(M2);
% end
M2 = D;
M1 = [zeros(num_observer*nsteps,p);M2];
rankM2 = zeros(nsteps,1);
rankM1 = zeros(nsteps,1);
for i = 0:nsteps-1
    M2 = [M2; C*A^i*B];
    M1 = [[zeros(num_observer*(nsteps-1-i),p);M2],M1];
    rankM2(i+1) = rank(M2);
    rankM1(i+1) = rank(M1);
end

plot(1:nsteps,rankM1);




