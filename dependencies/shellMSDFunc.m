function [Xdot] = shellMSDFunc(t,X,CM,WN,L0,e0,LI,m,ke,kt,d,p,nF0,afn,printtime)

%REQUIRED:
%t - time
%X - state vector containing x1, vx1, y1, vy1, etc.
%CM- connectivity matrix
%WN- wing node matrix
%L0- initial spring lengths
%A0- initial facet areas
%e0- initial edge angles
%LI- list of springs, indicated [startnode endnode]
%m - mass vector, giving the mass of each particle
%p - applied pressure

%OPTIONAL:
%Ff- load amplitude vs. time function
%nF0 - nodes which are fixed

%the loaded nodes
% nFe = (1:size(X,1)/6)';

if printtime
fprintf('Sim Time is %f\n',t)
end
%initialize derivative of state vector with zeros
Xdot = zeros(size(X)); 

%extract positions, velocities of nodes [x1 y1 z1; x2 y2 z2; ... ]
x = reshape(X(1:2:end),3,[])';
v = reshape(X(2:2:end),3,[])';

%compute the length of all links
L = rssq(x(LI(:,1),:)-x(LI(:,2),:),2)';

%compute the direction of each link
D = x(LI(:,1),:)-x(LI(:,2),:);
D = D'./vecnorm(D');

%assemble spring force vector
f_s = -ke'.*D.*(L - L0); 

%assemble damping force vecor
dv = v(LI(:,1),:)-v(LI(:,2),:); 
f_d = -d'.*D.*dot(dv',D);

%determine the normals of each facet
[N,~]=patchNormal(CM,x);

%compute facet areas
A = tri_area(CM,x);

%compute the pressure vector on each face
pt = p*afn(t);
Fp = A.*pt.*N/3;
Fp = reshape(repmat(Fp',3,1),3,[])';

%compute edge angles
e = computeEdgeAngles(WN,N);

%compute torque at edge
Te = kt.*(e-e0);

%compute wingtip force2
Fw1 = N(WN(:,1),:).*Te./A(WN(:,1));
Fw2 = N(WN(:,3),:).*Te./A(WN(:,3));
Fw = [Fw1 Fw2];

%compute spine forces to preserve equilibrium
Fsp = -(Fw1+Fw2)/2;
Fsp = [Fsp Fsp];

%add forces at the oppostie end of each link, and append Fp and Fb
f = [f_s(:); -f_s(:); ...
     f_d(:); -f_d(:); ...
     Fp(:);...
     Fw(:);...
     Fsp(:)];

%Create a Matrix which maps degrees of freedom to nodes
ID = reshape(1:numel(x),3,[])';

idsd = repmat(ID(LI,:),2,1)';   %ID's for spring and damper forces
idfp = ID(reshape(CM',[],1),:); %ID's for pressure forces
idfwn = [ID(WN(:,2),:) ID(WN(:,4),:)];      %ID's for Wing Node forces
idfsp = [ID(LI(:,1),:) ID(LI(:,2),:)]; %ID's for Spine Nodes

id = [idsd(:); idfp(:); idfwn(:); idfsp(:)];

F = accumarray(id,f);

%divide by masses in order to produce accelerations
%we duplicate mass vector to match dimension of nDOF
Xdot(2:2:end) = F./reshape(repmat(m',3,1),[],1);

%zero out accelerations for fixed nodes
Xdot(reshape(2*ID(nF0,:),[],1))=0;

%store incoming velocities as outgoing position derivatives
Xdot(1:2:end) = X(2:2:end);

% fprintf('sim time = %f\n',t)

% cla
% gpatch(CM,x,'gw','k',0.5);hold on
% patchNormPlot(CM,x);
% view(3)
% drawnow

% %plot the endpoints of an edge in question
% plot3(x(18,1),x(18,2),x(18,3),'b.','markersize',40)
% plot3(x(32,1),x(32,2),x(32,3),'b.','markersize',40)
% i = 65;
% plot3(x(WN(i,2),1),x(WN(i,2),2),x(WN(i,2),3),'r.','markersize',40)
% plot3(x(WN(i,4),1),x(WN(i,4),2),x(WN(i,4),3),'r.','markersize',40)
% quiver3(x(WN(i,2),1),x(WN(i,2),2),x(WN(i,2),3),Fw1(i,1),Fw1(i,2),Fw1(i,3),'r','filled','LineWidth',2)
% quiver3(x(WN(i,4),1),x(WN(i,4),2),x(WN(i,4),3),Fw2(i,1),Fw2(i,2),Fw2(i,3),'r','filled','LineWidth',2)
% quiver3(x(18,1),x(18,2),x(18,3),Fsp(i,1),Fsp(i,2),Fsp(i,3),'b','filled','LineWidth',2)
% quiver3(x(32,1),x(32,2),x(32,3),Fsp(i,1),Fsp(i,2),Fsp(i,3),'b','filled','LineWidth',2)

end