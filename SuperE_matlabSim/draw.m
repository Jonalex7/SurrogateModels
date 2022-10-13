function[M] = draw(propertiesInit,dispTot,dofs,solveOption,plasticity,time)

maxDimensions = max([max(propertiesInit.nodes(:,1))-min(propertiesInit.nodes(:,1)) max(propertiesInit.nodes(:,2))-min(propertiesInit.nodes(:,2)) max(propertiesInit.nodes(:,3))-min(propertiesInit.nodes(:,3))]);

nodesDisp = zeros(size(propertiesInit.nodes,1),3);

% Coordinates of displaced nodes

for i = 1:size(nodesDisp,1)
    for j = 1:3
        nodesDisp(i,j) = propertiesInit.nodes(i,j) + solveOption.amplification*dispTot(propertiesInit.nodes(i,j+3));
    end
end

% Coordinates of elements

originDisp      = zeros(dofs.nbElements,3); % Coordinates of origin of displaced elements
endDisp         = zeros(dofs.nbElements,3); % Coordinates of end of displaced element
for i = 1:dofs.nbElements
    originDisp(i,:)      = nodesDisp(propertiesInit.nodein(i),1:3);
    endDisp(i,:)         = nodesDisp(propertiesInit.nodeout(i),1:3);
end
LDisp     = zeros(dofs.nbElements,1); % Length
angleDisp = zeros(3,3,dofs.nbElements); % vector of local axis
for i = 1:dofs.nbElements
    Deltax = endDisp(i,1) - originDisp(i,1);
    Deltay = endDisp(i,2) - originDisp(i,2);
    Deltaz = endDisp(i,3) - originDisp(i,3);
    LDisp(i) = sqrt(Deltax^2 + Deltay^2 + Deltaz^2);
    
    angleDisp(1,:,i) = [Deltax Deltay Deltaz];
    angleDisp(1,:,i) = angleDisp(1,:,i) / norm(angleDisp(1,:,i));
end

% Draw

[xsphere,ysphere,zsphere] = sphere;

figure;
hold on
for i = 1:dofs.nbElements
    % Initial structure
%     plot3([propertiesInit.origin(i,1) propertiesInit.end(i,1)],[propertiesInit.origin(i,2) propertiesInit.end(i,2)],[propertiesInit.origin(i,3) propertiesInit.end(i,3)]);
    
    % Deformed structure
    if plasticity(1,i) == -1
        lineColor = 'g';
        lineStyle = '-';    
    elseif plasticity(1,i) == -2
        lineColor = 'g';
        lineStyle = '--';
    elseif plasticity(1,i) == 1;
        lineColor = 'r';
        lineStyle = '-';
    else
        lineColor = 'r';
        lineStyle = '--';
    end
    plot3([originDisp(i,1) endDisp(i,1)],[originDisp(i,2) endDisp(i,2)],[originDisp(i,3) endDisp(i,3)],'LineWidth',2,'color',lineColor,'LineStyle',lineStyle);
    
    % Plastic hinges
    radius = min(maxDimensions/100,LDisp(i)/3);
    if plasticity(2,i) == 1
        xCenterSphere = originDisp(i,1) + radius*angleDisp(1,1,i);
        yCenterSphere = originDisp(i,2) + radius*angleDisp(1,2,i);
        zCenterSphere = originDisp(i,3) + radius*angleDisp(1,3,i);
        x = xsphere*radius + xCenterSphere;
        y = ysphere*radius + yCenterSphere;
        z = zsphere*radius + zCenterSphere;
        colormap([1 0 0])
        surf(x,y,z,'EdgeColor','none');
    end
    if plasticity(3,i) == 1
        xCenterSphere = endDisp(i,1) - radius*angleDisp(1,1,i);
        yCenterSphere = endDisp(i,2) - radius*angleDisp(1,2,i);
        zCenterSphere = endDisp(i,3) - radius*angleDisp(1,3,i);
        x = xsphere*radius + xCenterSphere;
        y = ysphere*radius + yCenterSphere;
        z = zsphere*radius + zCenterSphere;
        colormap([1 0 0])
        surf(x,y,z,'EdgeColor','none');
    end
end
view(solveOption.az,solveOption.el)
axis equal
axis([min(propertiesInit.nodes(:,1))-1 max(propertiesInit.nodes(:,1))+1 min(propertiesInit.nodes(:,2))-1 max(propertiesInit.nodes(:,2))+1 min(propertiesInit.nodes(:,3))-1 max(propertiesInit.nodes(:,3))+1])
title(['\fontsize{24}After ',num2str(time),' s'])
% hold off

M = getframe(gcf);