%% Vectorized MSD Simulator on Arbitrary Tets in 110 lines
%Lawrence Smith | lasm4254@colorado.edu | 2 April 2021
clear; clc; close all
load small_bunny_mesh.mat

%Extract mesh connectivity
CM = T.ConnectivityList;  
NC = T.Points; 
LI = edges(T); 
ID = edgeAttachments(T,edges(T)); 
ID2 = cell2mat(ID')';

rho     = 10;           %Density
E       = 50;           %Elastic Modulus
nu      = 0.25;         %Poisson's Rato
kappa   = E*(4*nu-1)/(2*(1+nu)*(1-nu)); %Eqn 5.36 of Golec

%compute squares of initial spring lengths
L0 = arrayfun(@(i) sumsqr(diff(NC(LI(i,:),:))), 1:size(LI,1));

%compute initial tetrahedral volumes
V0 = arrayfun(@(i) abs(det([NC(CM(i,:)',:) ones(4,1)]))/6,1:size(CM,1));

%Compute effective lengths of each unit cell
L0eff = (V0*12/sqrt(2)).^(1/3);

%Compute and assemble spring stiffnesses using effective lengths
I = cell2mat(arrayfun(@(i) repmat(i,1,numel(ID{i})),1:size(LI,1),'UniformOutput',false));
Ke = 8*sqrt(2)/105*L0eff'*E; 
K0 = accumarray(I',Ke(ID2)); %EQ 14 of Harders

%compute masses of each node by accumulating volume of each element
m = accumarray(CM(:),reshape(repmat(rho*V0/4,1,4),[],1));

%find indices of points which are to be fixed in space, and loaded
fixed =  find(NC(:,3)<0);
forced = find(NC(:,3)>35);

%Prepare the initial conditions
X0 = reshape([reshape(NC',[],1)'; zeros(1,numel(NC))],[],1);

tic; 
opts = odeset('RelTol',1e-3,'AbsTol',1e-4);
[t,y] = ode23(@(t,X) integrate(t,X,L0,V0,K0,LI,CM,kappa,m,fixed,forced),linspace(0,10,100),X0,opts);
fprintf('\n Simulation complete. Clock Time Elapsed: %f s.\n',toc);
animateTet(t,y,CM,fixed,forced)

function [Xdot] = integrate(t,X,L0,V0,K0,LI,CM,kappa,m,fixed,forced)
fprintf('Sim time = %f s\n',t); 
g = 0; 
d = 0.25;

%extract positions, velocities of nodes [x1 y1 z1; x2 y2 z2; ... ]
x = reshape(X(1:2:end),3,[])';
v = reshape(X(2:2:end),3,[])';

%compute the length of all links
L = arrayfun(@(i) sumsqr(diff(x(LI(i,:),:))), 1:size(LI,1));

%compute the direction of each link
D = x(LI(:,1),:)-x(LI(:,2),:);
D = D'./vecnorm(D');

%assemble spring force vector
f_s = -K0'.*D.*(L - L0); 

%assemble damping force vecor
dv = v(LI(:,1),:)-v(LI(:,2),:); 
f_d = -d*K0'.*D.*dot(dv',D);

%compute volumes of unit cells
V = arrayfun(@(i) abs(det([x(CM(i,:)',:) ones(4,1)]))/6,1:size(CM,1));

%create a 3D matrix with a row for each unit cell
co = cat(3,CM(:,[2 3 4]),CM(:,[1 3 4]),CM(:,[1 2 4]),CM(:,[1 2 3]));

%compute the derivative of each volume wrt position of each node
dVdx = sum(cross(reshape(x(co,3),[],3,4),reshape(x(co,2),[],3,4),2),2);
dVdy = sum(cross(reshape(x(co,1),[],3,4),reshape(x(co,3),[],3,4),2),2);
dVdz = sum(cross(reshape(x(co,2),[],3,4),reshape(x(co,1),[],3,4),2),2);
dVdX = permute(cat(2,dVdx,dVdy,dVdz),[2 3 1]);

%compute change in volume from reference, the volumetric forces
dV = reshape(repmat((V-V0)./V0,1,12),[],1); 
f_v = -kappa*dV.*dVdX(:);

%append the inverted forces at the oppostie end of each link
f = [f_s(:); -f_s(:); f_d(:); -f_d(:); f_v];

%apply forces to model
index = reshape(1:numel(x),3,[])';
F = accumarray([reshape(repmat(index(LI,:),2,1)',[],1); reshape(index(CM,:)',[],1)],f);
Xdot = zeros(size(X)); Xdot(2:2:end) = F;

if t<2 %apply traction to ears of bunny
    Xdot(6*forced) = Xdot(6*forced) + 6e2;
end

%store incoming velocities as outgoing position derivatives
Xdot(1:2:end) = X(2:2:end);

%apply gravity loading
Xdot(6:6:end) = Xdot(6:6:end) - g; 

%divide by masses in order to produce accelerations
Xdot(2:2:end) = Xdot(2:2:end)./reshape(repmat(m,1,3),[],1);

%remove accelerations from nodes which are fixed in space
Xdot([6*fixed-4; 6*fixed-2; 6*fixed]) = 0;
end