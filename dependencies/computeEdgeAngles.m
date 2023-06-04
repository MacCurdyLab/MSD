function [eA] = computeEdgeAngles(WN,N)

% eA = acos(dot(N(LI(:,1),:),N(LI(:,2),:),2));

eA = real(acos(dot(N(WN(:,1),:), N(WN(:,3),:), 2)));

end