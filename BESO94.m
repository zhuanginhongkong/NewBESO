function BESO94(nelx,nely,nelz,vol_con,er,rmin)
scale = 1;
nelx = 48*scale; nely = 24*scale; nelz = 24*scale;
E0 = 1; nu = 0.3; vol_con = 0.12;
er = 0.02;
rmin = 2;
%% INITIALIZE PATTERN
nElem = nely*nelx*nelz;  % TOTAL ELEMENT NUMBER
x = ones(nElem,1);
dv0 = ones(nElem,1);
vol = 1; iter = 0; change = 1; c = []; maxiter = 1000;
%% INDEXING NODES AND ELEMENT
nodenrs = reshape(1:(1+nelx)*(1+nely)*(1+nelz),1+nely,1+nelz,1+nelx);
edofVec = reshape(3*nodenrs(1:end-1,1:end-1,1:end-1)+1,nElem,1);
edofMat = repmat(edofVec,1,24)+repmat([0,1,2,3*(nely+1)*(nelz+1)+[0,1,2,-3,-2,-1],-3,-2,-1,3*(nely+...
    1)+[0,1,2],3*(nely+1)*(nelz+2)+[0,1,2,-3,-2,-1],3*(nely+1)+[-3,-2,-1]],nElem,1);
%% 3D STIFFNESS MATRIX
k = 1/(1+nu)/(2*nu-1)/144 *([-32;-6;-6;8;6;6;10;6;3;-4;-6;-3;-4;-3;-6;10;3;6;8;3;3;4;-3;-3; -32;-6;-6;-4;-3;6;10;3;6;8;6;-3;-4;-6;-3;4;-3;3;8;3;
    3;10;6;-32;-6;-3;-4;-3;-3;4;-3;-6;-4;6;6;8;6;3;10;3;3;8;3;6;10;-32;6;6;-4;6;3;10;-6;-3;10;-3;-6;-4;3;6;4;3;3;8;-3;-3;-32;-6;-6;8;6;-6;10;3;3;4;
    -3;3;-4;-6;-3;10;6;-3;8;3;-32;3;-6;-4;3;-3;4;-6;3;10;-6;6;8;-3;6;10;-3;3;8;-32;-6;6;8;6;-6;8;3;-3;4;-3;3;-4;-3;6;10;3;-6;-32;6;-6;-4;3;3;8;-3;
    3;10;-6;-3;-4;6;-3;4;3;-32;6;3;-4;-3;-3;8;-3;-6;10;-6;-6;8;-6;-3;10;-32;6;-6;4;3;-3;8;-3;3;10;-3;6;-4;3;-6;-32;6;-3;10;-6;-3;8;-3;3;4;3;3;-4;6;
    -32;3;-6;10;3;-3;8;6;-3;10;6;-6;8;-32;-6;6;8;6;-6;10;6;-3;-4;-6;3;-32;6;-6;-4;3;6;10;-3;6;8;-6;-32;6;3;-4;3;3;4;3;6;-4;-32;6;-6;-4;6;-3;10;-6;3;
    -32;6;-6;8;-6;-6;10;-3;-32;-3;6;-4;-3;3;4;-32;-6;-6;8;6;6;-32;-6;-6;-4;-3;-32;-6;-3;-4;-32;6;6;-32;-6;-32] ...
    +nu*[48;0;0;0;-24;-24;-12;0;-12;0;24;0;0;0;24;-12;-12;0;-12;0;0;-12;12;12;48;0;24;0;0;0;-12;-12;-24;0;-24;0;0;24;12;-12;12;0;-12;0;-12;-12;0;
    48;24;0;0;12;12;-12;0;24;0;-24;-24;0;0;-12;-12;0;0;-12;-12;0;-12;48;0;0;0;-24;0;-12;0;12;-12;12;0;0;0;-24;-12;-12;-12;-12;0;0;48;0;24;0;-24;0;
    -12;-12;-12;-12;12;0;0;24;12;-12;0;0;-12;0;48;0;24;0;-12;12;-12;0;-12;-12;24;-24;0;12;0;-12;0;0;-12;48;0;0;0;-24;24;-12;0;0;-12;12;-12;0;0;-24;
    -12;-12;0;48;0;24;0;0;0;-12;0;-12;-12;0;0;0;-24;12;-12;-12;48;-24;0;0;0;0;-12;12;0;-12;24;24;0;0;12;-12;48;0;0;-12;-12;12;-12;0;0;-12;12;0;0;0;
    24;48;0;12;-12;0;0;-12;0;-12;-12;-12;0;0;-24;48;-12;0;-12;0;0;-12;0;12;-12;-24;24;0;48;0;0;0;-24;24;-12;0;12;0;24;0;48;0;24;0;0;0;-12;12;-24;0;
    24;48;-24;0;0;-12;-12;-12;0;-24;0;48;0;0;0;-24;0;-12;0;-12;48;0;24;0;24;0;-12;12;48;0;-24;0;12;-12;-12;48;0;0;0;-24;-24;48;0;24;0;0;48;24;0;0;48;0;0;48;0;48]);
KE(tril(ones(24))==1) = k';
KE = reshape(KE,24,24);
KE = KE+KE'-diag(diag(KE));
index_ik = arrayfun(@(i) i:24,1:24,'UniformOutput',false);
index_ik = horzcat(index_ik{:});
index_jk = arrayfun(@(j) repmat(j,1,24-j+1),1:24,'UniformOutput',false);
index_jk = horzcat(index_jk{:});
iK = edofMat(:, index_ik)';
jK = edofMat(:, index_jk)';
index_k = sort([iK(:), jK(:)],2,'descend');  clear iK jK;
%% DEFINE LOADS AND SUPPORTS (CANTILEVER)
F = sparse(3*nodenrs(1:nely+1,1,nelx+1),1,-1,3*(nely+1)*(nelx+1)*(nelz+1),1);
fixeddofs = 1:3*(nely+1)*(nelz+1);
U = zeros(3*(nely+1)*(nelx+1)*(nelz+1),1);
alldofs = 1:3*(nely+1)*(nelx+1)*(nelz+1);
freedofs = setdiff(alldofs,fixeddofs);
%% PREPARE 3D FILTER
[dy,dz,dx]=meshgrid(-ceil(rmin)+1:ceil(rmin)-1,-ceil(rmin)+1:ceil(rmin)-1,-ceil(rmin)+1:ceil(rmin)-1);
h = max(0,rmin-sqrt(dx.^2+dy.^2+dz.^2));
Hs = imfilter(ones(nely,nelz,nelx),h,'symmetric');
% Hs = imfilter(ones(nely,nelz,nelx),h,0);
%% MAIN OPTIMIZATION LOOP
while change > 0.001 && iter < maxiter
    iter = iter + 1;
    vol = max(vol*(1-er),vol_con);
    % FINITE ELEMENT ANALYSIS
    sK = reshape(k(:)*E0*x(:)',length(k)*nElem,1);
    K = sparse(index_k(:,1),index_k(:,2),sK);
    LK = chol(K(freedofs,freedofs), 'lower');
    U(freedofs) = LK'\(LK\F(freedofs));
    ce = reshape(sum((U(edofMat)*KE).*U(edofMat),2),nElem,1);
    c = [c, 0.5.*sum(x*E0.*ce)];
    % SENSITIVITY ANALYSIS
    dc = -x.*ce*E0;
    dc = imfilter(reshape(dc,nely,nelz,nelx)./Hs,h,'symmetric');
    dv = imfilter(reshape(dv0,nely,nelz,nelx)./Hs,h,'symmetric');
    if iter > 1; dc = (dc+olddc)/2.; end  % STABILIZATION OF EVOLUTIONARY PROCESS
    % BESO DESIGN UPDATE
    l1 = min(min(min(-dc./dv))); l2 = max(max(max(-dc./dv)));
    while ((l2-l1)/l2 > 1e-5)
        th = (l1+l2)/2;
        x = max(1e-9,sign(-dc./dv-th));
        if sum(sum(sum(x)))-vol*nElem > 0
            l1 = th;
        else
            l2 = th;
        end
    end
    % PRINT RESULTS AND PLOT DENSITIES
    if iter > 10
        change = abs(sum(c(iter-9:iter-5))-sum(c(iter-4:iter)))/sum(c(iter-4:iter));
    end
    disp([' It.: ' sprintf('%2i',iter) ' Obj.: ' sprintf('%6.4f',c(iter)) ' Vol.: ' sprintf('%5.3f',sum(sum(sum(x)))/(nElem)) ' ch.: ' sprintf('%5.3f',change)])
    surf = shiftdim(reshape(x,nely,nelz,nelx),2);
    surf = smooth3(surf,'box',1);
    clf; patch(isosurface(surf,0.5),'FaceColor',[0 127 102]/255,'EdgeColor','none','FaceAlpha', 0.9);
    patch(isocaps(surf,0.5),'FaceColor',[255 255 102]/255,'EdgeColor','none','FaceAlpha', 0.7);
    light('Position', [1 1 1], 'Style', 'infinite');   % ADD LIGHTING
    light('Position', [-1 -1 1], 'Style', 'infinite');
    lighting phong; 
    material([0.5 0.6 0.4]);
    drawnow; view([110,20]); axis equal tight off; pause(1e-6);
    x = reshape(x,nElem,1);
    olddc = dc;
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