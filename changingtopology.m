% This code tests the observability of various topologies under the same
% sensor setting.
% For now, all displacements are observable.

clear;
close all;
load structure.mat; % load a set of xPhys 
number_cases = size(xPhys_set,1);
deflection_set = zeros(number_cases,1);
estimation_error_set = zeros(number_cases,1);
% observability_set = zeros(number_cases,1);

nelx = 40;
nely = 10;

for case_id = 1:number_cases
    xPhys = reshape(xPhys_set(case_id,:), nely, nelx);
    
    %% dynamics using finite element method
    % Set parameters
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
    T = 1e1; % maximum time
    nsteps = 1e1;
    tspan = linspace(0,T,nsteps); % t0 = 0; tf = 10, nsteps = 100; 
    y = ode4(@(t,y)testode(t,y,...
        K(freedofs,freedofs),M(freedofs,freedofs),F(freedofs)),tspan,zeros(2*p,1));

    % store max deflection
    deflection = y(:,p+1:end);
    deflection_set(case_id) = max(abs(deflection(:)));
    
    % plot solution
%     figure; hold on;
%     plot(tspan,y(:,end)); plot([0,tspan(end)],[U(end),U(end)],'-.'); 

    % prepare matrices for input (force) estimation
    dt = T/nsteps;
    Kb = K(freedofs,freedofs);
    Mb = M(freedofs,freedofs);
    Fb = F(freedofs);
    A = [zeros(p), -inv(Mb)*Kb; eye(p), zeros(p)];
    B = [inv(Mb); zeros(p)];

    % set observer
    C = zeros(2*p);
    C((p+1):4:end,(p+1):4:end) = eye(p/4); 
    C((p+2):4:end,(p+2):4:end) = eye(p/4); % assume odd nodes' displacements are observable
    
    % calculate some big matrices
    Ab = eye(2*p)+dt/2*A;
    AA = eye(2*p)+(dt/6*(eye(2*p) + ...
        2*Ab + 2*(eye(2*p)+dt/2*A*Ab) + (eye(2*p)+dt*A*(eye(2*p)+dt/2*A*Ab))))*A;
    BB = (dt/6*(eye(2*p) + ...
        2*Ab + 2*(eye(2*p)+dt/2*A*Ab) + (eye(2*p)+dt*A*(eye(2*p)+dt/2*A*Ab))))*B;

    % calculate more matrices, see README 
    %TODO: write README
    D = zeros(2*p*(nsteps-1),p);
    for i=1:(nsteps-1)
        D((i-1)*2*p+(1:2*p),:) = C*AA^(i-1)*BB;
    end
    D = kron(sparse(triu(ones(nsteps-1))'),eye(2*p))*D;

    % input estimation 
    Y = y(2:end,:)'; Y = Y(:);
    fbar = (D'*D)\(D'*Y);
    
    estimation_error_set(case_id) = norm(fbar-Fb);
%     observability_set(case_id) = det(D'*D);
end
plot(estimation_error_set, deflection_set,'.','MarkerSize',20);
[pareto,pareto_ID] = prtp([estimation_error_set,deflection_set]);

figure;colormap(gray); imagesc(1-reshape(xPhys_set(pareto_ID(3),:),nely,nelx)); caxis([0 1]); axis equal; axis off; drawnow;
