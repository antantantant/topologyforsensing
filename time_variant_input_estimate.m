% time invariant force, linear time invariant system

% clear;
% close all;
% load structure.mat; % load xPhys

draw_deformation = 0; % draw deformation of the structure
test_variance = 1; % test theoretical variance, TESTED

optimize = 0; % topology optimization for observability
T = 1e2; % maximum time
nsteps = 2; % time interval

volfrac = 0.3; % volume fraction

% Set parameters
nelx = 10; % horizontal number of elements (left to right)
nely = 4; % vertical number of elements (top to down)
penal = 3; % polynomial order to define density-young's modulus relationship
E0 = 1; % young's modulus at density=1
Emin = 1e-9; % young's modulus at density=0, keep this small
M0 = 1;
Mmin = 1e-9;
nu = 0.3; % poisson's ratio

nf = nelx; % number of forces
num_observer = nf; % number of observers

% Case 1
% Left size of the beam is fixed to the ground
fixeddofs = [1:2*(nely+1)];
alldofs = [1:2*(nely+1)*(nelx+1)];
freedofs = setdiff(alldofs,fixeddofs);
p = size(freedofs,2);

% set observer
S = zeros(num_observer,p); 
% S(:,randperm(p,num_observer))=eye(num_observer);
S(:,2*(1:nelx)*(nely+1)) = eye(nelx); % put y-axis sensors on top of the beam

Sp = zeros(p,nf); % Sp specifies the loading location
Sp(2*(1:nelx)*(nely+1),:) = eye(nelx); % put loads at the bottom of the beam
Fb = -1*ones(nelx,1);

% % Case 2
% % Left size of the beam is fixed to the ground
% % fixeddofs = [1:2*(nely+1)];
% fixeddofs = [1:2*(nely+1),(2*(nelx-1)*(nely+1)+1):2*nelx*(nely+1)];
% alldofs = [1:2*(nely+1)*(nelx+1)];
% freedofs = setdiff(alldofs,fixeddofs);
% p = size(freedofs,2);
% 
% % set observer
% S = zeros(num_observer,p); 
% % S(:,2*(1:nelx)*(nely+1)) = eye(nelx); % put y-axis sensors on top of the beam
% S(:,[nelx*(nely+1)-1,nelx*(nely+1)]) = eye(num_observer); % put y-axis sensors on top of the beam
% 
% % set loads
% Sp = zeros(p,nf); % Sp specifies the loading location
% Sp([nelx*(nely+1)-1,nelx*(nely+1)],:) = eye(nf); % put loads at the bottom of the beam
% Fb = -1*ones(nf,1);



% define material density
xPhys = ones(nely,nelx);
% xPhys = reshape(x_soln(:,end),nely*nelx,1);
% xPhys = rand(nely,nelx);
% xPhys = xPhys/sum(xPhys(:))*volfrac*nelx*nely;
% xPhys = volfrac*ones(nely,nelx); 
% xPhys(2:end-1,2:end-1)=0;
% xPhys = [1 1 0 0; 1 1 1 0; 0 1 1 1; 0 0 1 1];

%% define the structure
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
sM = reshape(ME(:)*(Mmin+xPhys(:)'*(M0-Mmin)),64*nelx*nely,1);
M = sparse(iK,jK,sM); M = (M+M')/2;


%% dynamic simulation
% prepare matrices for input (force) estimation
dt = T/nsteps;
Kb = K(freedofs,freedofs);
Mb = M(freedofs,freedofs);

A = [zeros(p), -Mb\Kb; eye(p), zeros(p)];
B = [Mb\Sp; zeros(p,nf)];

% calculate some big matrices
%The following uses a theoretical solution
AA = expm(A*dt);
BB = A\(AA-eye(2*p))*B;

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
%         if mod(tt,10)==1
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
%         end
    end
    figure;
    xPhys_all = zeros(nely+1,nelx+1);xPhys_all(2:end,2:end)=reshape(xPhys,nely,nelx);
    y_all = zeros(2*(nely+1)*(nelx+1),1); y_all(2*(nely+1)+1:end) = y(end,p+1:end)';
    dy = reshape(y_all(2:2:end),nely+1,nelx+1); 
    dy_normal = dy/max(abs(y(:))+eps);
    final_deformation = (xPhys_all>0).*-dy_normal;
    imagesc(final_deformation);
    axis equal; axis off; drawnow;
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
% calculate matrices for observability
D = S*inv(Mb)*Sp;
C1 = -S*inv(Mb)*Kb;
C = [zeros(num_observer,p), C1];
% M2 follows corollary 5 eq 35
% M1 follows eq 13b


%% check theoretical variance with simulation results
if test_variance
    % rankM2 = zeros(nsteps,1);
    % rankM1 = zeros(nsteps,1);
    num_test = 1000;
    error = zeros(num_test, nsteps-1,nf);
    fbar = zeros(nsteps-1,nf);

    for test_id = 1:num_test
        % add noise to state
        noise_std = 0.1;
        acc_noise = acc + normrnd(0,noise_std,size(acc,1),size(acc,2));
        M2 = D;
        M1 = M2;

        for i = 1:nsteps-1
            Ytemp = S*acc_noise(1:i,:)';
            Ytemp = Ytemp(:);
            temp = ((M1'*M1)\(M1'*Ytemp))';
            fbar(i,:) = temp((i-1)*nf+(1:nf));
            error(test_id, i,:) = (fbar(i,:)'-Fb_set((i-1)*nf+(1:nf)))';
        %     error(i) = norm(fbar-Fb_set(1:i*nf))/i;
        %     rankM2(i) = rank(M2);
        %     rankM1(i) = rank(M1);
            M2 = [M2; C*AA^(i-1)*BB];
            M1 = [M2,[zeros(num_observer,nf*i);M1]];
        end
        % rankM2(i) = rank(M2);
        % rankM1(i) = rank(M1);
    end
    theoretical_estimation_std = noise_std*sqrt(diag(inv(M1'*M1)));


    % figure;
    % plot(1:nsteps-1,error);
    std(error,1)
    % log(det(D'*D))
    % eig(D'*D)
end
figure;shadedErrorBar(1:80,mean(error,1),std(error,1),'g');
hold on; plot(theoretical_estimation_std(1:80));
hold on; plot(-theoretical_estimation_std(1:80));
%% calculate theoretical variance gradient
if optimize
    %%%%%%%%%%%%%%%% just to double check: only need to consider D
    % M2 = D;
    % M1 = M2;
    % logdetM1 = zeros(nsteps,1);
    % logdetM1(1) = log(det(M1'*M1));
    % for i = 1:nsteps-1
    %     M2 = [M2; C*AA^(i-1)*BB];
    %     M1 = [M2,[zeros(num_observer,nf*i);M1]];
    %     logdetM1(i+1) = log(det(M1'*M1));
    % end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % This is how mass is calculated
    % iK = reshape(kron(edofMat,ones(8,1))',64*nelx*nely,1);
    % jK = reshape(kron(edofMat,ones(1,8))',64*nelx*nely,1);
    % sK = reshape(KE(:)*(Emin+xPhys(:)'.^penal*(E0-Emin)),64*nelx*nely,1);
    % K = sparse(iK,jK,sK); K = (K+K')/2;
    % sM = reshape(ME(:)*(Mmin+xPhys(:)'*(M0-Mmin)),64*nelx*nely,1);
    % M = sparse(iK,jK,sM);

    % All dofs:
    % fixeddofs = [1:2*(nely+1)];
    % alldofs = [1:2*(nely+1)*(nelx+1)];
    % freedofs = setdiff(alldofs,fixeddofs);
    Sp_all = sparse(repmat(freedofs,1,nf)',kron(1:nf,ones(1,length(freedofs)))',...
        Sp(:),2*(nely+1)*(nelx+1),nf);
    S_all = sparse(repmat(1:num_observer,1,length(freedofs))',kron(freedofs,ones(1,nf))',...
        S(:),num_observer,2*(nely+1)*(nelx+1));

    % gradient descent for observability
    % ***zero density of the element where input (or observer???) is at maximizes the trace
    % but this violates the stress/compliance requirement
    
    
    % use multistart to deal with local solutions
    ntrial = 1000;
    x_soln = zeros(nelx*nely,ntrial);
    f_soln = zeros(1,ntrial);
    for trial = 1:ntrial
        x0 = rand(nely,nelx);
        x0 = x0/sum(x0(:))*volfrac*nelx*nely;
        x_old = x0(:);
        x = x_old;
        dx = ones(nelx*nely,1);
        free_var = 1:(nelx*nely-1); % fix the last element where observers/inputs are currently at

        loop = 0;
        change = 1;
        detD = [];

        %% START ITERATION
        figure;
        while change > 1e-6
            loop = loop + 1;

    %         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %         % calculation of dtr(D)/dx = dtr(SM^-1Sp)/d(M)*d(M)/dx
    %         dtr_dx = zeros(nelx*nely,1);
    %         dtr_dM = -(inv(Mb)*Sp*S*inv(Mb))';
    %         count = 0;
    %         while count<nelx*nely
    %             dM_dx = sparse(iK(count*64+(1:64)), jK(count*64+(1:64)), ...
    %                 ME(:)*(M0-Mmin), 2*(nely+1)*(nelx+1), 2*(nely+1)*(nelx+1));
    %             dtr_dx(count+1) = trace(dtr_dM*dM_dx(freedofs,freedofs));
    %             count = count + 1;
    %         end
    %         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %         
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % calculation of ddet(D)/dx = det(SM^-1Sp)M'*d(M^-1)/dx
            ddet_dx = zeros(nelx*nely,1);
            ddet_dinvM = det(S*inv(Mb)*Sp)*Mb';
            count = 0;
            while count<nelx*nely
                dM_dx = sparse(iK(count*64+(1:64)), jK(count*64+(1:64)), ...
                    ME(:)*(M0-Mmin), 2*(nely+1)*(nelx+1), 2*(nely+1)*(nelx+1));
                ddet_dx(count+1) = trace(ddet_dinvM*...
                    (-inv(Mb)*dM_dx(freedofs,freedofs)*inv(Mb))');
                count = count + 1;
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            dc = ddet_dx(free_var);
            dv = ones(length(free_var),1);
            % iterate on step size
            alpha = 1e-1;
            vol_gap = 1;
            while vol_gap>1e-3
                alpha = alpha/2;
                dx = alpha*(dc - dv'*dc*dv/(dv'*dv));
                x(free_var) = max(min(x_old(free_var)+dx,1),0);
                vol_gap = abs(sum(x)-volfrac*nelx*nely);
            end
            change = max(abs(x-x_old));
            x_old = x;
            % update mass
            sM = reshape(ME(:)*(Mmin+x(:)'*(M0-Mmin)),64*nelx*nely,1);
            M = sparse(iK,jK,sM);  
            Mb = M(freedofs,freedofs);

            % update trace of D
            detD(loop) = log(det(S*inv(Mb)*Sp));

        end
        x_soln(:,trial) = x;
        f_soln(trial) = detD(end);
    end
    x_best = x_soln(:,f_soln==max(f_soln));
    colormap(gray); imagesc(1-reshape(x_best,nely,nelx)); 
    caxis([0 1]); axis equal; axis off; drawnow;
%     figure;
%     plot(detD)
    save('topopt_obs.mat','x_best','nelx','nely','S','Sp');
end

%% check difference between topopt_obs and default design
x_default = xPhys;
x_topopt = x_best;
sM = reshape(ME(:)*(Mmin+x_default(:)'*(M0-Mmin)),64*nelx*nely,1);
M_default = sparse(iK,jK,sM); 
sK = reshape(KE(:)*(Emin+x_topopt(:)'.^penal*(E0-Emin)),64*nelx*nely,1);
sM = reshape(ME(:)*(Mmin+x_topopt(:)'*(M0-Mmin)),64*nelx*nely,1);
K_topopt = sparse(iK,jK,sK); K_topopt = (K_topopt+K_topopt')/2;
M_topopt = sparse(iK,jK,sM);

num_test = 1000;
error_default = zeros(num_test, nsteps-1,nf);
error_topopt = zeros(num_test, nsteps-1,nf);

% test default first
Mb_default = M_default(freedofs,freedofs);
D_default = S*inv(Mb_default)*Sp;
C1_default = -S*inv(Mb_default)*Kb;
C_default = [zeros(num_observer,p), C1_default];
fbar_default = zeros(nsteps-1,nf);
acc_default = acc; 
for test_id = 1:num_test
    % add noise to state
    noise_std = 0.1;
    acc_noise_default = acc_default + normrnd(0,noise_std,...
        size(acc_default,1),size(acc_default,2));
    M2_default = D_default;
    M1_default = M2_default;

    for i = 1:nsteps-1
        Ytemp = S*acc_noise_default(1:i,:)';
        Ytemp = Ytemp(:);
        temp = ((M1_default'*M1_default)\(M1_default'*Ytemp))';
        fbar_default(i,:) = temp((i-1)*nf+(1:nf));
        error_default(test_id, i,:) = (fbar_default(i,:)'-Fb_set((i-1)*nf+(1:nf)))';
        M2_default = [M2_default; C_default*AA^(i-1)*BB];
        M1_default = [M2_default,[zeros(num_observer,nf*i);M1_default]];
    end
end
theoretical_estimation_std_default = ...
    noise_std*sqrt(diag(inv(M1_default'*M1_default)));
std(error_default,1)

% now test topopt
Mb_topopt = M_topopt(freedofs,freedofs);
Kb_topopt = K_topopt(freedofs,freedofs);
% rerun dynamic simulation
A_topopt = [zeros(p), -Mb_topopt\Kb_topopt; eye(p), zeros(p)];
B_topopt = [inv(Mb_topopt)*Sp; zeros(p,nf)];
AA_topopt = expm(A_topopt*dt);
BB_topopt = inv(A_topopt)*(AA_topopt-eye(2*p))*B_topopt;
y = zeros(nsteps,2*p);
for i = 2:nsteps
    y(i,:) = (AA_topopt*y(i-1,:)' + BB_topopt*Fb_set(nf*(i-2)+(1:nf)))';
end

% calculate the acceleration (since we will use acc as the observation)
acc_topopt = (inv(Mb_topopt)*Sp*reshape(Fb_set,nf,nsteps-1) -...
    inv(Mb_topopt)*Kb_topopt*y(1:end-1,p+1:end)')';

D_topopt = S*inv(Mb_topopt)*Sp;
C1_topopt = -S*inv(Mb_topopt)*Kb_topopt;
C_topopt = [zeros(num_observer,p), C1_topopt];
fbar_topopt = zeros(nsteps-1,nf);
for test_id = 1:num_test
    % add noise to state
    noise_std = 0.1;
    acc_noise_topopt = acc_topopt + ...
        normrnd(0,noise_std,size(acc_topopt,1),size(acc_topopt,2));
    M2_topopt = D_topopt;
    M1_topopt = M2_topopt;

    for i = 1:nsteps-1
        Ytemp = S*acc_noise_topopt(1:i,:)';
        Ytemp = Ytemp(:);
        temp = ((M1_topopt'*M1_topopt)\(M1_topopt'*Ytemp))';
        fbar_topopt(i,:) = temp((i-1)*nf+(1:nf));
        error_topopt(test_id, i,:) = ...
            (fbar_topopt(i,:)'-Fb_set((i-1)*nf+(1:nf)))';
        M2_topopt = [M2_topopt; C_topopt*AA_topopt^(i-1)*BB_topopt];
        M1_topopt = [M2_topopt,[zeros(num_observer,nf*i);M1_topopt]];
    end
end
theoretical_estimation_std_topopt = ...
    noise_std*sqrt(diag(inv(M1_topopt'*M1_topopt)));
std(error_topopt,1)


