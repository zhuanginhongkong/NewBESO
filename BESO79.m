function BESO79(nelx,nely,vol_con,er,rmin)
scale = 1;
nelx = 160*scale; nely = 100*scale;
E0 = 1; nu = 0.3; vol_con = 0.5;
er = 0.01; rmin = 5*scale^1.5;
%% INITIALIZE PATTERN
x = ones(nely,nelx); vol = 1; iter = 0; change = 1; c = []; index = 0;
%% INDEXING NODES AND ELEMENT
nodenrs = reshape(1:(1+nelx)*(1+nely),1+nely,1+nelx);
edofVec = reshape(2*nodenrs(1:end-1,1:end-1)+1,nelx*nely,1);
edofMat = repmat(edofVec,1,8)+repmat([0 1 2*nely+[2 3 0 1] -2 -1],nelx*nely,1);
%% STIFFNESS MATRIX
k = [1/2-nu/6 1/8+nu/8 -1/4-nu/12 -1/8+3*nu/8 -1/4+nu/12 -1/8-nu/8 nu/6 1/8-3*nu/8];
KE = 1/(1-nu^2)*[k(1) k(2) k(3) k(4) k(5) k(6) k(7) k(8)
    k(2) k(1) k(8) k(7) k(6) k(5) k(4) k(3)
    k(3) k(8) k(1) k(6) k(7) k(4) k(5) k(2)
    k(4) k(7) k(6) k(1) k(8) k(3) k(2) k(5)
    k(5) k(6) k(7) k(8) k(1) k(2) k(3) k(4)
    k(6) k(5) k(4) k(3) k(2) k(1) k(8) k(7)
    k(7) k(4) k(5) k(2) k(3) k(8) k(1) k(6)
    k(8) k(3) k(2) k(5) k(4) k(7) k(6) k(1)];
iK = reshape(kron(edofMat,ones(8,1))',64*nelx*nely,1);
jK = reshape(kron(edofMat,ones(1,8))',64*nelx*nely,1);
%% DEFINE LOADS AND SUPPORTS (CANTILEVER)
F = sparse(2*(nelx+1)*(nely+1)-nely,1,-1,2*(nely+1)*(nelx+1),1);
U = zeros(2*(nely+1)*(nelx+1),1);
fixeddofs = 1:2*(nely+1);
alldofs = 1:2*(nely+1)*(nelx+1);
freedofs = setdiff(alldofs,fixeddofs);
%% PREPARE FILTER
iH = ones(nelx*nely*(2*(ceil(rmin)-1)+1)^2,1);
jH = ones(size(iH));
sH = zeros(size(iH));
for i1 = 1:nelx
    for j1 = 1:nely
        e1 = (i1-1)*nely+j1;
        [i2,j2] = ndgrid(max(i1-(ceil(rmin)-1),1):min(i1+(ceil(rmin)-1),nelx),max(j1-(ceil(rmin)-1),1):min(j1+(ceil(rmin)-1),nely));
        e2 = (i2(:)-1)*nely+j2(:);
        iH(index + (1:numel(e2))) = e1;
        jH(index + (1:numel(e2))) = e2;
        sH(index + (1:numel(e2))) = max(0,rmin-sqrt((i1-i2(:)).^2+(j1-j2(:)).^2));
        index = index + numel(e2);
    end
end
H = sparse(iH,jH,sH); Hs = sum(H,2);
%% MAIN OPTIMIZATION LOOP
while change > 0.001
    iter = iter + 1;
    vol = max(vol*(1-er),vol_con);
    if iter >1; olddc = dc; end
    % FINITE ELEMENT ANALYSIS
    sK = reshape(KE(:)*E0*x(:)',64*nelx*nely,1);
    K = sparse(iK,jK,sK);  K = (K+K')/2;
    U(freedofs) = decomposition(K(freedofs,freedofs),'chol','lower')\F(freedofs);  
    ce = reshape(sum((U(edofMat)*KE).*U(edofMat),2),nely,nelx);
    c = [c, 0.5.*sum(sum(x*E0.*ce))];
    % SENSITIVITY ANALYSIS
    dc = -x.*ce;
    dv = ones(nely,nelx);
    dc(:) = H*dc(:)./Hs;
    if iter > 1; dc = (dc+olddc)/2.; end  % STABILIZATION OF EVOLUTIONARY PROCESS
    % BESO DESIGN UPDATE
    l1 = min(min(-dc./dv)); l2 = max(max(-dc./dv));
    while ((l2-l1)/l2 > 1e-5)
        th = (l1+l2)/2;
        x = max(1e-9,sign(-dc./dv-th));
        if sum(sum(x))-vol*(nelx*nely) > 0
            l1 = th;
        else
            l2 = th;
        end
    end
    % PRINT RESULTS AND PLOT DENSITIES
    if iter > 10
        change = abs(sum(c(iter-9:iter-5))-sum(c(iter-4:iter)))/sum(c(iter-4:iter));
    end
    disp([' It.: ' sprintf('%2i',iter) ' Obj.: ' sprintf('%6.4f',c(iter)) ' Vol.: ' sprintf('%5.3f',sum(sum(x))/(nelx*nely)) ' ch.: ' sprintf('%5.3f',change)])
    clf; colormap(summer); imagesc(x); axis equal tight off; pause(1e-6);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This Matlab code was written by Zicheng Zhuang and Yi Min Xie            %
% The Centre for Innovative Structures and Materials (CISM),               %
% RMIT University, Melbourne, VIC, Australia                               %
% Please send your comments to: zhuanginhongkong@outlook.com               %
%                                                                          %
% This code is proposed for educational purposes.                          %
%                                                                          %
% Disclaimer:                                                              %
% The authors reserve all rights but do not guarantee that the code is     %
% free from errors. Furthermore, we shall not be liable in any event.      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%