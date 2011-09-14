function W = test_FistaTree(Y,X, W0)

param.num_threads=-1; % all cores (-1 by default)
param.verbose=false;   % verbosity, false by default
param.it0=10;      % frequency for duality gap computations
param.max_it=500; % maximum number of iterations
param.L0=0.1;
param.tol=1e-5;
param.intercept=false;

param.lambda=1e-1; % regularization parameter
param.pos=true;

%% Example of tree structure
%% tree structured groups:
%% g1= {0 1 2 3 4 5 6}    root(g1) = {0};
%% g2= {1 2 3}            root(g2) = {1};
%% g3= {4 5 6}            root(g3) = {4};
%% g4= {2}                root(g4) = {2};
%% g5= {3}                root(g5) = {3};
%% g6= {5}                root(g6) = {5};
%% g7= {6}                root(g7) = {6};
%tree.own_variables=  int32([0 1 4 2 3 5 6]);   % pointer to the first variable of each group
%tree.N_own_variables=int32([1 1 1 1 1 1 1]); % number of "root" variables in each group
%%tree.eta_g=[1 2 2 3.5 3.5 3.5 3.5];
%tree.eta_g=[1 1 1 1 1 1 1];
%tree.groups=sparse([0 0 0 0 0 0 0; ...
%                    1 0 0 0 0 0 0; ...
%                    1 0 0 0 0 0 0; ...
%                    0 1 0 0 0 0 0; ...
%                    0 1 0 0 0 0 0; ...
%                    0 0 1 0 0 0 0; ...
%                    0 0 1 0 0 0 0]);  % first group should always be the root of the tree

% Example of tree structure
% tree structured groups:
% g1= {0 1 2 3 4}    root(g1) = {0};
% g2= {1 2}          root(g2) = {1};
% g3= {3 4}          root(g3) = {3};
% g4= {2}            root(g4) = {2};
% g5= {4}            root(g5) = {4};
tree.own_variables=  int32([0 1 3 2 4]);   % pointer to the first variable of each group
tree.N_own_variables=int32([1 1 1 1 1]); % number of "root" variables in each group
tree.eta_g=[1 2 2 3 3];
tree.groups=sparse([0 0 0 0 0; ...
                    1 0 0 0 0; ...
                    1 0 0 0 0; ...
                    0 1 0 0 0; ...
                    0 0 1 0 0]);  % first group should always be the root of the tree

%X=randn(100,10);
%X=X-repmat(mean(X),[size(X,1) 1]);
%X=mexNormalize(X);
%Y=randn(100,100);
%Y=Y-repmat(mean(Y),[size(Y,1) 1]);
%Y=mexNormalize(Y);
%W0=zeros(size(X,2),size(Y,2));

param.compute_gram=true;
param.loss='square';
param.regul='tree-l2';
[W optim_info]=mexFistaTree(Y,X,W0,tree,param);
