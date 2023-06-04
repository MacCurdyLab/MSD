function [] = plotme3D(t,t1,X1,CM,nF0,writeMov)

X = interpDisp(t1,X1,t);

figure
set(gcf,'Position',[385 65 1041 904])
M = max(cat(1,X(:,:,1),X(:,:,end)))*1.1;
m = min(cat(1,X(:,:,1),X(:,:,end)))*1.1;

% Fb = freeBoundary(triangulation(CM,X(:,:,1)));

clor = cool(10);

for i = 1:length(t)
    cla
%     patch('faces',CM,'Vertices',X(:,:,i),'facecolor',clor(1,:),'Marker','.','markersize',10,'linewidth',1,'edgecolor','k'); 
    D = rssq(X(:,:,i)-X(:,:,1),2);
    patch('faces',CM,'Vertices',X(:,:,i),'facecolor','interp','facevertexcdata',D,'Marker','.','markersize',10,'linewidth',1,'edgecolor','k');
    hold on
    plot3(X(nF0,1,1),X(nF0,2,1),X(nF0,3,1),'r.','markersize',20)
    colormap(clor)
    colorbar
    axis equal
    axis([m(1) M(1) m(2) M(2) m(3) M(3)])
    title(sprintf('Sim. Time = %1.1f',t(i)))
    campos([-509.3266 -141.2313  160.0748])
    camlight headlight
    material dull
    gdrawnow()
    F(i) = getframe();
end

if writeMov
v = VideoWriter('shellFing.avi');
open(v)
writeVideo(v,F)
close(v)
end

end