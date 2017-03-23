% time invariant force, linear time invariant system

clear;
close all;
% load structure.mat; % load xPhys

draw_deformation = 0; % draw deformation of the structure
T = 2; % maximum time
nsteps = 2; % time interval
num_observer = 2; % number of observers

%% define the structure
% Set parameters
nelx = 4; % horizontal number of elements (left to right)
nely = 4; % vertical number of elements (top to down)
penal = 3; % polynomial order to define density-young's modulus relationship
E0 = 1; % young's modulus at density=1
Emin = 1e-9; % young's modulus at density=0, keep this small
M0 = 1;
Mmin = 1e-9;
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
% xPhys = rand(nely,nelx);
xPhys = ones(nely,nelx); 
% xPhys(2:end-1,2:end-1)=0;
% xPhys = [1 1 0 0; 1 1 1 0; 0 1 1 1; 0 0 1 1];
[nely,nelx] = size(xPhys);

% element-to-global assembly
nodenrs = reshape(1:(1+nelx)*(1+nely),1+nely,1+nelx);
edofVec = reshape(2*nodenrs(1:end-1,1:end-1)+1,nelx*nely,1);
edofMat = repmat(edofVec,1,8)+repmat([0 1 2*nely+[2 3 0 1] -2 -1],nelx*nely,1);
iK = reshape(kron(edofMat,ones(8,1))',64*nelx*nely,1);
jK = reshape(kron(edofMat,ones(1,8))',64*nelx*nely,1);
sK = reshape(KE(:)*(Emin+xPhys(:)'.^penal*(E0-Emin)),64*nelx*nely,1);
K = sparse(iK,jK,sK); K = (K+K')/2;
sM = reshape(ME(:)*(Mmin+xPhys(:)'*(M0-Mmin)),64*nelx*nely,1);
M = sparse(iK,jK,sM); M = (M+M')/2;


%% dynamic simulation
% define loads and dofs
% F = sparse(2*(nelx+1)*(nely+1),1,-1e-2,2*(nely+1)*(nelx+1),1);
fixeddofs = [1:2*(nely+1)];
alldofs = [1:2*(nely+1)*(nelx+1)];
freedofs = setdiff(alldofs,fixeddofs);

% prepare matrices for input (force) estimation
p = size(freedofs,2);
dt = T/nsteps;
Kb = K(freedofs,freedofs);
Mb = M(freedofs,freedofs);
Fb = [0.1;0.1];
Sp = zeros(p,2); % Sp specifies the loading location
Sp(end,1)=1; Sp(end-1,2)=1;
nf = size(Sp,2); % number of forces

A = [zeros(p), -Mb\Kb; eye(p), zeros(p)];
B = [inv(Mb)*Sp; zeros(p,nf)];

% calculate some big matrices
%The following uses a theoretical solution
AA = expm(A*dt);
BB = inv(A)*(AA-eye(2*p))*B;

% time variant input Fb*i
Fb_set = kron((1:nsteps-1)',eye(nf))*Fb;
y = zeros(nsteps,2*p);
for i = 2:nsteps
    y(i,:) = (AA*y(i-1,:)' + BB*Fb_set(nf*(i-2)+(1:nf)))';
end
% calculate the acceleration (since we will use acc as the observation)
acc = (inv(Mb)*Sp*reshape(Fb_set,nf,nsteps-1) -...
    inv(Mb)*Kb*y(1:end-1,p+1:end)')';

% draw the structure
if draw_deformation
    for tt = 1:nsteps
        if mod(tt,10)==1
            figure; hold on; axis equal;
            sx = repmat(0:nelx,nely+1,1);
            sy = repmat((nely:-1:0)',1,nelx+1);
            ally = zeros((nely+1)*(nelx+1)*2,1);
            ally(freedofs) = y(tt,p+1:end)';
            dx = reshape(ally(1:2:end),nely+1,nelx+1);
            dy = reshape(ally(2:2:end),nely+1,nelx+1);
            dx_normal = dx/max(abs(y(:))+eps);
            dy_normal = dy/max(abs(y(:))+eps);
            sxx = sx+dx_normal;
            syy = sy+dy_normal;
            for ii=1:nelx+1
                plot(sxx(:,ii),syy(:,ii));
            end
            for ii=1:nely+1
                plot(sxx(ii,:),syy(ii,:));
            end
        end
    end
end

%% model reduction
% [Phi,Omega] = eig(Mb\Kb); % Phi are the eigenvectors and Omega are the eigenvalues
% omega = diag(Omega);
% plot(omega);
% nm = 5; % num_modal
% Phi_r = Phi(:,1:nm);
% omega_r = omega(1:nm);
% Omega_r = diag(omega_r);
% A_r = [zeros(nm), -Omega_r^2; eye(nm), zeros(nm)];
% B_r = [Phi_r'*Sp; zeros(nm,nf)]; 
% AA_r = expm(A_r*dt);
% BB_r = inv(A_r)*(AA_r-eye(2*nm))*B_r;

%% input estimation
% set observer
S = zeros(num_observer,p); 
% S(:,randperm(p,num_observer))=eye(num_observer);
S(1,end) = 1; S(2,end-1) = 1;

% calculate matrices for observability
D = S*inv(Mb)*Sp;
C1 = -S*inv(Mb)*Kb;
C = [zeros(num_observer,p), C1];
% M2 follows corollary 5 eq 35
% M1 follows eq 13b
M2 = D;
M1 = M2;

rankM2 = zeros(nsteps,1);
rankM1 = zeros(nsteps,1);
num_test = 100;
error = zeros(num_test, nsteps-1,nf);
fbar = zeros(nsteps-1,nf);

% add noise to state
noise_std = (max(abs(acc(:)))-min(abs(acc(:))))*0.1;
theoretical_estimation_std = noise_std*sqrt(diag(inv(D'*D)))

for test_id = 1:num_test
    acc_noise = acc + normrnd(0,noise_std,size(acc,1),size(acc,2));
    for i = 1:nsteps-1
        Ytemp = S*acc_noise(1:i,:)';
        Ytemp = Ytemp(:);
        temp = ((M1'*M1)\(M1'*Ytemp))';
        fbar(i,:) = temp((i-1)*nf+(1:nf));
        error(test_id, i,:) = (fbar(i,:)'-Fb_set((i-1)*nf+(1:nf)))';
    %     error(i) = norm(fbar-Fb_set(1:i*nf))/i;
    %     rankM2(i) = rank(M2);
    %     rankM1(i) = rank(M1);
%         M2 = [M2; C*AA^(i-1)*BB];
%         M1 = [M2,[zeros(num_observer,nf*i);M1]];
    end
    % rankM2(i) = rank(M2);
    % rankM1(i) = rank(M1);
end

% figure;
% plot(1:nsteps-1,error);
std(error,1)
% log(det(D'*D))
% eig(D'*D)