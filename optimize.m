% clear;
% close all;
% 
% %% Set parameters
% nelx = 80; % horizontal number of elements (left to right)
% nely = 20; % vertical number of elements (top to down)
% penal = 3; % polynomial order to define density-young's modulus relationship
% E0 = 1; % young's modulus at density=1
% Emin = 1e-9; % young's modulus at density=0, keep this small
% M0 = 1;
% Mmin = 1e-9;
% nu = 0.3; % poisson's ratio
% rmin = 1.5; % filter size
% volfrac = 0.3; % volume fraction
% 
% nf = nelx; % number of forces
% % nf = 2;
% num_observer = nf; % number of observers
% 
% % Case 1
% % Left size of the beam is fixed to the ground
% fixeddofs = [1:2*(nely+1)];
% alldofs = [1:2*(nely+1)*(nelx+1)];
% freedofs = setdiff(alldofs,fixeddofs);
% p = size(freedofs,2);
% 
% % set observer
% S = zeros(num_observer,p); 
% % S(:,randperm(p,num_observer))=eye(num_observer);
% S(:,2*(1:nelx)*(nely+1)) = eye(nelx); % put y-axis sensors on top of the beam
% 
% Sp = zeros(p,nf); % Sp specifies the loading location
% Sp(2*(1:nelx)*(nely+1),:) = eye(nelx); % put loads at the bottom of the beam
% Fb = -1*ones(nelx,1);
% 
% % % Case 2
% % % Left size of the beam is fixed to the ground
% % % fixeddofs = [1:2*(nely+1)];
% % fixeddofs = [1:2*(nely+1),(2*(nelx-1)*(nely+1)+1):2*nelx*(nely+1)];
% % alldofs = [1:2*(nely+1)*(nelx+1)];
% % freedofs = setdiff(alldofs,fixeddofs);
% % p = size(freedofs,2);
% % 
% % % set observer
% % S = zeros(num_observer,p); 
% % % S(:,2*(1:nelx)*(nely+1)) = eye(nelx); % put y-axis sensors on top of the beam
% % S(:,[nelx*(nely+1)-1,nelx*(nely+1)]) = eye(num_observer); % put y-axis sensors on top of the beam
% % 
% % % set loads
% % Sp = zeros(p,nf); % Sp specifies the loading location
% % Sp([nelx*(nely+1)-1,nelx*(nely+1)],:) = eye(nf); % put loads at the bottom of the beam
% 
% %% define the structure
% % element-wise stiffness matrix for a quadrilateral element (square in shape)
% A11 = [12  3 -6 -3;  3 12  3  0; -6  3 12 -3; -3  0 -3 12]; %stiffness matrix
% A12 = [-6 -3  0  3; -3 -6 -3 -6;  0 -3 -6  3;  3 -6  3 -6];
% B11 = [-4  3 -2  9;  3 -4 -9  4; -2 -9 -4 -3;  9  4 -3 -4];
% B12 = [ 2 -3  4 -9; -3  2  9 -2;  4  9  2  3; -9 -2  3  2];
% KE = 1/(1-nu^2)/24*([A11 A12;A12' A11]+nu*[B11 B12;B12' B11]);
% 
% % element-wise consistent mass matrix for a quadrilateral element (square in shape)
% ME = (4*eye(8) + [zeros(4),eye(4);eye(4),zeros(4)] + repmat(kron([0 1;1 0],2*eye(2)),2,2))/9;
% 
% % element-to-global assembly
% nodenrs = reshape(1:(1+nelx)*(1+nely),1+nely,1+nelx);
% edofVec = reshape(2*nodenrs(1:end-1,1:end-1)+1,nelx*nely,1);
% edofMat = repmat(edofVec,1,8)+repmat([0 1 2*nely+[2 3 0 1] -2 -1],nelx*nely,1);
% iK = reshape(kron(edofMat,ones(8,1))',64*nelx*nely,1);
% jK = reshape(kron(edofMat,ones(1,8))',64*nelx*nely,1);
% 
% % static input
% Fb = sparse(find(sum(Sp,2)==1),1,-1,p,1);
% % Fb1 = sparse([nelx*(nely+1)-1,nelx*(nely+1)],1,[cos(0),sin(0)],p,1);
% % Fb2 = sparse([nelx*(nely+1)-1,nelx*(nely+1)],1,[cos(pi/6),sin(pi/6)],p,1);
% % Fb3 = sparse([nelx*(nely+1)-1,nelx*(nely+1)],1,[cos(pi/3),sin(pi/3)],p,1);
% % Fb4 = sparse([nelx*(nely+1)-1,nelx*(nely+1)],1,[cos(pi/2),sin(pi/2)],p,1);
% % Fb5 = sparse([nelx*(nely+1)-1,nelx*(nely+1)],1,[cos(4*pi/6),sin(4*pi/6)],p,1);
% % Fb6 = sparse([nelx*(nely+1)-1,nelx*(nely+1)],1,[cos(5*pi/6),sin(5*pi/6)],p,1);
% % Fb7 = sparse([nelx*(nely+1)-1,nelx*(nely+1)],1,[cos(6*pi/6),sin(6*pi/6)],p,1);
% U = zeros(2*(nely+1)*(nelx+1),1);
% 
% %% PREPARE FILTER
% iH = ones(nelx*nely*(2*(ceil(rmin)-1)+1)^2,1);
% jH = ones(size(iH));
% sH = zeros(size(iH));
% k = 0;
% for i1 = 1:nelx
%   for j1 = 1:nely
%     e1 = (i1-1)*nely+j1;
%     for i2 = max(i1-(ceil(rmin)-1),1):min(i1+(ceil(rmin)-1),nelx)
%       for j2 = max(j1-(ceil(rmin)-1),1):min(j1+(ceil(rmin)-1),nely)
%         e2 = (i2-1)*nely+j2;
%         k = k+1;
%         iH(k) = e1;
%         jH(k) = e2;
%         sH(k) = max(0,rmin-sqrt((i1-i2)^2+(j1-j2)^2));
%       end
%     end
%   end
% end
% H = sparse(iH,jH,sH);
% Hs = sum(H,2);
% 
% %% Optimization, use multistart to deal with local solutions
% max_iter = 200;
% 
% load dM_dX_case1.mat;
% % count = 0;
% % dM_dX = sparse(p*nelx*nely,p);
% % while count<nelx*nely
% %     fprintf('calculate d|D|/dx, count=%d \n',count);
% %     dM_dx = sparse(jK(count*64+(1:64)), iK(count*64+(1:64)),...
% %         ME(:)*(M0-Mmin), 2*(nely+1)*(nelx+1), 2*(nely+1)*(nelx+1));
% %     dM_dX(count*p+(1:p),:) = dM_dx(freedofs,freedofs);
% %     count = count + 1;
% % end
% 
% load x0_obs_case1.mat; 
% x0 = x;
% o0_set = exp(51:1:122);
% ntrial = length(o0_set);
% x_soln = zeros(nelx*nely,ntrial);
% f_soln = zeros(2,ntrial);
for trial = (length(o0_set)-1):-1:1
%     x0 = ones(nely,nelx)*volfrac;
%     x0 = x0/sum(x0(:))*volfrac*nelx*nely;
    x0 = x_soln(:,trial+1);
    x_old = x0(:);
    x = x_old;
    dx = ones(nelx*nely,1);

    o0 = o0_set(trial); % target minimal detD

    loop = 0;
    change = 1;
    o_set = [];
    c_set = [];
    
    w_old = 0;
    
    %% START ITERATION
    while change > 1e-6 && loop < max_iter
        loop = loop + 1;
        sK = reshape(KE(:)*(Emin+x(:)'.^penal*(E0-Emin)),64*nelx*nely,1);
        K = sparse(iK,jK,sK); K = (K+K')/2;
        sM = reshape(ME(:)*(Mmin+x(:)'*(M0-Mmin)),64*nelx*nely,1);
        M = sparse(iK,jK,sM); M = (M+M')/2;
        Kb = K(freedofs,freedofs);
        Mb = M(freedofs,freedofs);
        o = det(S*inv(Mb)*Sp); % current detD
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % calculation of ddet(D)/dx = det(SM^-1Sp)M'*d(M^-1)/dx
        % calculate this only when the observability constraint is violated
        if o<o0
            temp = (-dM_dX*inv(Mb)*o)';
            temp3 = reshape(temp,p*p,[]);
            ddet_dx = sum(temp3(1:(p+1):end,:))';

            dd = ddet_dx;
            A = -dd';
            b = -o0+o;
        else
            A = [];
            b = [];
            dd = zeros(nelx*nely,1);
        end
        dv = ones(nelx*nely,1);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % calculate gradient for compliance
        U(freedofs) = Kb\Fb;
%         U1(freedofs) = Kb\Fb1;
%         U2(freedofs) = Kb\Fb2;
%         U3(freedofs) = Kb\Fb3;
%         U4(freedofs) = Kb\Fb4;
%         U5(freedofs) = Kb\Fb5;
%         U6(freedofs) = Kb\Fb6;
%         U7(freedofs) = Kb\Fb7;
        ce = reshape(sum((U(edofMat)*KE).*U(edofMat),2),nely,nelx);
%         ce = reshape(sum((U1(edofMat)*KE).*U1(edofMat)+...
%             (U2(edofMat)*KE).*U2(edofMat)+...
%             (U3(edofMat)*KE).*U3(edofMat)+...
%             (U4(edofMat)*KE).*U4(edofMat)+...
%             (U5(edofMat)*KE).*U5(edofMat)+...
%             (U6(edofMat)*KE).*U6(edofMat)+...
%             (U7(edofMat)*KE).*U7(edofMat),2),nely,nelx);
        c = sum(sum((Emin+reshape(x,nely,nelx).^penal*(E0-Emin)).*ce));
        dc = -penal*(E0-Emin)*reshape(x,nely,nelx).^(penal-1).*ce;
        dc(:) = H*(x(:).*dc(:))./Hs./max(1e-3,x(:)); % use filter
        dc = dc(:);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % solve a linear programming problem for dx
        f = dc/1e3;
        Aeq = ones(1,nelx*nely);
        beq = 0;
        lb = zeros(nelx*nely,1)-x;
        ub = ones(nelx*nely,1)-x;
        [dx,fval,exitflag,output,lambda_all] = linprog(f,A,b,Aeq,beq,lb,ub,zeros(nelx*nely,1));
        
%         f = -ddet_dx;
%         f = f/max(abs(f));
%         Aeq = ones(1,nelx*nely);
%         beq = 0;
%         lb = zeros(nelx*nely,1)-x;
%         ub = ones(nelx*nely,1)-x;
%         [dx,fval,exitflag,output,lambda_all] = linprog(f,[],[],Aeq,beq,lb,ub,zeros(nelx*nely,1));
%         lambda = lambda_all.eqlin;
%         mu = lambda_all.ineqlin;
%         if isempty(mu) mu = 0;end
%         
%         w = max(mu, 0.5*(mu+w_old))*1e3;
%         w_old = w;
%         
        % do a line search
        alpha = 0.2;
%         L0 = c/1e3+w*(-o+o0)*(-o+o0>0); % Lagrangian at alpha = 0
% %         L0 = -o;
%         L = 1e16; 
%         dLdx = dc/1e3 -w*ddet_dx;
% %         dLdx = -ddet_dx;
        
        o_temp = -1e16;
        while (o_temp<o0) && alpha>1e-12 
%         while (L>L0+0.1*dLdx'*alpha*dx) && alpha>1e-12 
%         while c_temp>c
            alpha = alpha/2;
fprintf('take a step: alpha=%d \n',alpha)

            x_temp = x + alpha*dx;
            % update compliance
            sK = reshape(KE(:)*(Emin+x_temp(:)'.^penal*(E0-Emin)),64*nelx*nely,1);
            K = sparse(iK,jK,sK); K = (K+K')/2;Kb = K(freedofs,freedofs);
            U(freedofs) = Kb\Fb;
            ce = reshape(sum((U(edofMat)*KE).*U(edofMat),2),nely,nelx);
%             U1(freedofs) = Kb\Fb1;
%             U2(freedofs) = Kb\Fb2;
%             U3(freedofs) = Kb\Fb3;
%             U4(freedofs) = Kb\Fb4;
%             U5(freedofs) = Kb\Fb5;
%             U6(freedofs) = Kb\Fb6;
%             U7(freedofs) = Kb\Fb7;
%             ce = reshape(sum((U1(edofMat)*KE).*U1(edofMat)+...
%                 (U2(edofMat)*KE).*U2(edofMat)+...
%                 (U3(edofMat)*KE).*U3(edofMat)+...
%                 (U4(edofMat)*KE).*U4(edofMat)+...
%                 (U5(edofMat)*KE).*U5(edofMat)+...
%                 (U6(edofMat)*KE).*U6(edofMat)+...
%                 (U7(edofMat)*KE).*U7(edofMat),2),nely,nelx);
            c_temp = sum(sum((Emin+reshape(x_temp,nely,nelx).^penal*(E0-Emin)).*ce));
            
            % update observability
            sM = reshape(ME(:)*(Mmin+x_temp(:)'*(M0-Mmin)),64*nelx*nely,1);
            M = sparse(iK,jK,sM);  
            Mb = M(freedofs,freedofs);
            o_temp = det(S*inv(Mb)*Sp);
            
%             L = c_temp/1e3 + w*(-o_temp+o0)*(-o_temp+o0>0);
%             L = -o_temp;
            
            
            
        end
        
        % update x
        x_old = x;
        x = x_temp;
        change = max(abs(x-x_old))
        
        % store performance
        o_set(loop) = o_temp;
        c_set(loop) = c_temp;

    end
    x_soln(:,trial) = x;
    f_soln(:,trial) = [o_set(end);c_set(end)];
end
x_best = x_soln(:,f_soln(2,:)==max(f_soln(2,:)));
colormap(gray); imagesc(1-reshape(x_best,nely,nelx)); 
caxis([0 1]); axis equal; axis off; drawnow;
% colormap(gray); imagesc(1-reshape(x,nely,nelx)); 
% caxis([0 1]); axis equal; axis off; drawnow;
figure;
countx = 0;
county = 0;
margin = 2;
% image_all = zeros(ceil(length(o_set)/4)*nely+(ceil(length(o_set)/4)-1)*margin,4*nelx+3*margin);
image_all = zeros(3*nely+2*margin,4*nelx+3*margin);

for i = 1:6:ntrial
    image_all(county*(nely+margin)+(1:nely),countx*(nelx+margin)+(1:nelx)) = reshape(x_soln(:,i),nely,nelx);
    countx = countx+1;
    if countx==4
        county = county+1;
        countx = 0;
    end
end
colormap(gray); imagesc(1-image_all);
caxis([0 1]); axis equal; axis off; drawnow;

figure;plot(-log(f_soln(1,:)),log(f_soln(2,:)),'.','MarkerSize',20);