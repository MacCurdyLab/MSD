function []  = animateTet(t,y,CM,fixed,forced)
close all
vlog = [];
figure
set(gcf,'position',[400 50 800 800])
    for i = 1:size(t,1)
        X = reshape(y(i,1:2:end),3,[])';
        V = arrayfun(@(i) abs(det([X(CM(i,:)',:) ones(4,1)]))/6,1:size(CM,1));
        vlog(i) = sum(V);
        [F,P] = freeBoundary(triangulation(CM,X));
        plot3(X(fixed,1),X(fixed,2),X(fixed,3),'r.','markersize',20)
        hold on
        if t(i)<2
            plot3(X(forced,1),X(forced,2),X(forced,3),'b.','markersize',20)
        end
        h = patch('faces',F,'vertices',P,'facecolor','w'); axis equal;
        view([0 -1 0])
        axis([-30 30 -30 30 -5 70])
        hold off
        pause(0.01); delete(h)
    end
    h = patch('faces',F,'vertices',P,'facecolor','w'); axis equal;
    axis([-30 30 -30 30 -5 70])
%     trisurf(F,P(:,1),P(:,2),P(:,3),'FaceColor','white'); axis equal
    
%     figure
%     plot(t,vlog,'linewidth',1.5); xlabel('time [s]'); ylabel('volume []');
end