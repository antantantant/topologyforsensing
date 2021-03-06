% draw pareto curve for compliance and observability
% time invariant force, linear time invariant system

% clear;
% close all;

T = 1e1; % maximum time
nsteps = 2; % time interval
volfrac = 0.3; % volume fraction

% Set parameters
nelx = 80; % horizontal number of elements (left to right)
nely = 20; % vertical number of elements (top to down)
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

%% all designs
% x = ones(nely*nelx,1)*volfrac;
x = x_soln(:,end);

%% create all data
xPhys = reshape(x,nely,nelx);
% xPhys(nely,2:end) = 0;
figure; colormap(gray); imagesc(1-reshape(xPhys,nely,nelx)); 

sK = reshape(KE(:)*(Emin+xPhys(:)'.^penal*(E0-Emin)),64*nelx*nely,1);
K = sparse(iK,jK,sK); K = (K+K')/2;
sM = reshape(ME(:)*(Mmin+xPhys(:)'*(M0-Mmin)),64*nelx*nely,1);
M = sparse(iK,jK,sM); M = (M+M')/2;

% calculate observability
Mb = M(freedofs,freedofs);
Kb = K(freedofs,freedofs);
% rerun dynamic simulation
D = S*inv(Mb)*Sp;
observability = log(det(D'*D));

% calculate compliance
Fb = sparse(find(sum(Sp,2)==1),1,-1,p,1);
u = Kb\Fb;
compliance = u'*Kb*u;
figure; hold on; axis equal;
sx = repmat(0:nelx,nely+1,1);
sy = repmat((nely:-1:0)',1,nelx+1);
ally = zeros((nely+1)*(nelx+1)*2,1);
ally(freedofs) = u;
dx = reshape(ally(1:2:end),nely+1,nelx+1);
dy = reshape(ally(2:2:end),nely+1,nelx+1);
dx_normal = dx/1e4;
dy_normal = dy/1e4;
sxx = sx+dx_normal;
syy = sy+dy_normal;

xx = reshape(x,nely,nelx);
for i = 1:nelx
    for j = 1:nely
        if xx(j,i)>0
            posx = [sxx(j,i),sxx(j,i+1),sxx(j+1,i+1),sxx(j+1,i)];
            posy = [syy(j,i),syy(j,i+1),syy(j+1,i+1),syy(j+1,i)];
            f = fill(posx,posy,ones(1,3)*(1-xx(j,i)));
%             set(f,'EdgeColor','none')
        end
    end
end
% for ii=1:nelx+1
%     plot(sxx(:,ii),syy(:,ii));
% end
% for ii=1:nely+1
%     plot(sxx(ii,:),syy(ii,:));
% end