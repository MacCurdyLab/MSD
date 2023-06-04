%% Shell Runner

clear; clc; close all

A = stlread('fingerDemo.stl');

rho     = 1e-4;        %[kg/mm^3] Density
E       = 200;         %[MPa] Elastic Modulus
th      = 2;           %[mm] Thickness
p       = 1;           %[MPa] applied pressure
zeta    = 300;         %damping ratio
beta    = 1e-2;        %bending stiffness ratio
tend    = 6;          %end time

ampfunc = @(t) 0.5*(erf(2*(t-1))+1);

CM = A.ConnectivityList;
NC = A.Points;

%extract edge list
LI = edges(A);

%Wing Node List [facet1 wingnode1 facet2 wingnode2]
WN = cell2mat(edgeAttachments(A,LI));
W1 = arrayfun(@(i) setdiff(CM(WN(i,1),:),LI(i,:)), 1:size(LI,1));
W2 = arrayfun(@(i) setdiff(CM(WN(i,2),:),LI(i,:)), 1:size(LI,1));
WN = [WN(:,1) W1' WN(:,2) W2'];

%compute facenormals
[N,~]=patchNormal(CM,NC);

%compute initial spring lengths
L0 = rssq(NC(LI(:,1),:)-NC(LI(:,2),:),2)';

%compute extension spring stiffnesses
ke = L0'*E;
d = zeta*ones(size(ke));

%compute torsion spring stiffnesses
kt = beta*th*ones(size(LI,1),1)*E;

%compute initial facet areas
A0 = tri_area(CM,NC);

%compute masses of each node by accumulating volume of each element
m = accumarray(CM(:),reshape(repmat(rho*th*A0/4,1,3),[],1));

%get angle associated with each edge
e0 = computeEdgeAngles(WN,N);

%set up time vector
ti = linspace(0,tend,200);

%Initial Conditions
X0 = reshape([reshape(NC',[],1)'; zeros(1,numel(NC))],[],1);

%Find Fixed nodes
nF0 = find(NC(:,3)<7);

%perform integration
tic
[t1,y1] = odeVc(@(t,X) shellMSDFunc(t,X,CM,WN,L0,e0,LI,m,ke,kt,d,p,nF0,ampfunc,false),[0 tend],X0,1e-3);
X1 = permute(reshape(y1(:,1:2:end)',3,[],length(t1)),[2 1 3]);
fprintf('Explicit (Verlet) Integration complete. %i timesteps solved in %1.2f s.\n',length(t1),toc);
plotme3D(ti,t1,X1,CM,nF0,false)

% ode45 runs about 2000x slower than odeVc... I estimate it would take about
% 3.5 hours on my machine. 

% [t2,y2] = ode45(@(t,X) shellMSDFunc(t,X,CM,WN,L0,e0,LI,m,ke,kt,d,p,nF0,ampfunc,true),[0 tend],X0);
% X2 = permute(reshape(y2(:,1:2:end)',3,[],length(t2)),[2 1 3]);
% fprintf('Explicit (ode45) Integration complete. %i timesteps solved in %1.2f s.\n',length(t2),toc);
% plotme3D(ti,t2,X2,CM,nF0,false)

