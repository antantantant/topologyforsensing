%% dynamics using finite element method
% Set parameters
nelx = 1; % horizontal number of elements (left to right)
nely = 1; % vertical number of elements (top to down)
penal = 2; % polynomial order to define density-young's modulus relationship
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
xPhys = ones(nely,nelx);

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

% static solution
U = zeros(2*(nely+1)*(nelx+1),1);
U(freedofs) = K(freedofs,freedofs)\F(freedofs);

% dynamic solution
p = size(freedofs,2);
T = 1e2; % maximum time
nsteps = 1e3;
tspan = linspace(0,T,nsteps); % t0 = 0; tf = 10, nsteps = 100; 
y = ode4(@(t,y)testode(t,y,...
    K(freedofs,freedofs),M(freedofs,freedofs),F(freedofs)),tspan,zeros(2*p,1));

% t0 = 0:T/(length(t)-1):T;
% y0 = interp1(t0,t,y);

figure; hold on;
plot(tspan,y(:,end)); plot([0,tspan(end)],[U(end),U(end)],'-.'); 

% figure; hold on;
% plot(t0,y0(:,end)); plot([0,t0(end)],[U(end),U(end)],'-.'); 

