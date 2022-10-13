function[punchProp,punchPropLeft,punchPropRight,resultPunch] = punchingScenarios2(properties,data,force,forceElem,punchProp,punchPropLeft,punchPropRight)

%% Localise punching at each node

% Left
for i = 1:size(punchPropLeft.nodeProp,1)
    if punchPropLeft.leg(i,3) == 1
        if forceElem(1,punchPropLeft.leg(i,2)) > punchPropLeft.nodeProp(i,13)
%             disp(['Punching at node ',num2str(punchPropLeft.leg(i,1)),' from brace ',num2str(punchPropLeft.leg(i,2)),' at left side'])
            if i <= punchProp.nodesLevels
                punchPropLeft.behaviourNodes(i,2) = 1;
            else
                for j = 1:punchProp.nodesLevels
                    if punchPropLeft.behaviourNodes(j,4) == punchPropLeft.leg(i,1)
                        punchPropLeft.behaviourNodes(j,5) = 1;
                    end
                    if punchPropLeft.behaviourNodes(j,7) == punchPropLeft.leg(i,1)
                        punchPropLeft.behaviourNodes(j,8) = 1;
                    end
                end
            end
        end
    end
    if punchPropLeft.leg(i,5) == 1
        if forceElem(1,punchPropLeft.leg(i,4)) > punchPropLeft.nodeProp(i,13)
%             disp(['Punching at node ',num2str(punchPropLeft.leg(i,1)),' from brace ',num2str(punchPropLeft.leg(i,4)),' at left side'])
            if i <= punchProp.nodesLevels
                punchPropLeft.behaviourNodes(i,2) = 1;
            else
                for j = 1:punchProp.nodesLevels
                    if punchPropLeft.behaviourNodes(j,4) == punchPropLeft.leg(i,1)
                        punchPropLeft.behaviourNodes(j,5) = 1;
                    end
                    if punchPropLeft.behaviourNodes(j,7) == punchPropLeft.leg(i,1)
                        punchPropLeft.behaviourNodes(j,8) = 1;
                    end
                end
            end
        end
    end
end

% Right
for i = 1:size(punchPropRight.nodeProp,1)
    if punchPropRight.leg(i,3) == 1
        if forceElem(1,punchPropRight.leg(i,2)) > punchPropRight.nodeProp(i,13)
%             disp(['Punching at node ',num2str(punchPropRight.leg(i,1)),' from brace ',num2str(punchPropRight.leg(i,2)),' at right side'])
            if i <= punchProp.nodesLevels
                punchPropRight.behaviourNodes(i,2) = 1;
            else
                for j = 1:punchProp.nodesLevels
                    if punchPropRight.behaviourNodes(j,4) == punchPropRight.leg(i,1)
                        punchPropRight.behaviourNodes(j,5) = 1;
                    end
                    if punchPropRight.behaviourNodes(j,7) == punchPropRight.leg(i,1)
                        punchPropRight.behaviourNodes(j,8) = 1;
                    end
                end
            end
        end
    end
    
    if punchPropRight.leg(i,5) == 1
        if forceElem(1,punchPropRight.leg(i,4)) > punchPropRight.nodeProp(i,13)
%             disp(['Punching at node ',num2str(punchPropRight.leg(i,1)),' from brace ',num2str(punchPropRight.leg(i,4)),' at right side'])
            if i <= punchProp.nodesLevels
                punchPropRight.behaviourNodes(i,2) = 1;
            else
                for j = 1:punchProp.nodesLevels
                    if punchPropRight.behaviourNodes(j,4) == punchPropRight.leg(i,1)
                        punchPropRight.behaviourNodes(j,5) = 1;
                    end
                    if punchPropRight.behaviourNodes(j,7) == punchPropRight.leg(i,1)
                        punchPropRight.behaviourNodes(j,8) = 1;
                    end
                end
            end
        end
    end
end

% No punching on impacted node

for i = 1:punchProp.nodesLevels
    if data.impactedZone == 0 && punchPropLeft.behaviourNodes(i,1) == data.impactedNode
        punchPropLeft.behaviourNodes(i,2) = -1;
        punchPropRight.behaviourNodes(i,2) = -1;
    end
end

%% Localise punching for each level

% 1st column : number of node on impacted leg
% 2nd column : displacement allowed on impacted leg ?
% 3rd column : displacement allowed on back leg ?

for i = 1:punchProp.nodesLevels
    % Left
    
    if punchPropLeft.behaviourNodes(i,2) == 1 % Impacted leg
        punchPropLeft.behaviourLevels(i,2) = 1;
    end
    
    if i ~= 1 && i ~= punchProp.nodesLevels % Rear leg
        if punchPropLeft.behaviourNodes(i,4) ~= 0 && punchPropLeft.behaviourNodes(i,7) ~= 0
            if punchPropLeft.behaviourNodes(i,5) == 1 && punchPropLeft.behaviourNodes(i,8) == 1
                punchPropLeft.behaviourLevels(i,3) = 1;
            end
        elseif punchPropLeft.behaviourNodes(i,4) ~= 0 && punchPropLeft.behaviourNodes(i,7) == 0
            if punchPropLeft.behaviourNodes(i,5) == 1
                punchPropLeft.behaviourLevels(i,3) = 1;
            end
        elseif punchPropLeft.behaviourNodes(i,4) == 0 && punchPropLeft.behaviourNodes(i,7) ~= 0
            if punchPropLeft.behaviourNodes(i,8) == 1
                punchPropLeft.behaviourLevels(i,3) = 1;
            end
        end
    end
    punchPropLeft.behaviourLevels(i,4) = max([punchPropLeft.behaviourLevels(i,2) punchPropLeft.behaviourLevels(i,3)]);
    if max(punchPropLeft.behaviourLevels(:,4)) == 0
        punchPropLeft.activated = 0;
    else
        punchPropLeft.activated = 1;
    end
    
    % Right
    
    if punchPropRight.behaviourNodes(i,2) == 1 % Impacted leg
        punchPropRight.behaviourLevels(i,2) = 1;
    end
    
    if i ~= 1 && i ~= punchProp.nodesLevels % Rear leg
        if punchPropRight.behaviourNodes(i,4) ~= 0 && punchPropRight.behaviourNodes(i,7) ~= 0
            if punchPropRight.behaviourNodes(i,5) == 1 && punchPropRight.behaviourNodes(i,8) == 1
                punchPropRight.behaviourLevels(i,3) = 1;
            end
        elseif punchPropRight.behaviourNodes(i,4) ~= 0 && punchPropRight.behaviourNodes(i,7) == 0
            if punchPropRight.behaviourNodes(i,5) == 1
                punchPropRight.behaviourLevels(i,3) = 1;
            end
        elseif punchPropRight.behaviourNodes(i,4) == 0 && punchPropRight.behaviourNodes(i,7) ~= 0
            if punchPropRight.behaviourNodes(i,8) == 1
                punchPropRight.behaviourLevels(i,3) = 1;
            end
        end
    end
    punchPropRight.behaviourLevels(i,4) = max([punchPropRight.behaviourLevels(i,2) punchPropRight.behaviourLevels(i,3)]);
    if max(punchPropRight.behaviourLevels(:,4)) == 0
        punchPropRight.activated = 0;
    else
        punchPropRight.activated = 1;
    end
end

%% Punching distance at each punched node

ddeltaShip = norm(force.valueDisp);

% ----------
% Left
% ----------

punchedLeft = 0;
belowBasisIndexLeft = 0;
aboveBasisIndexLeft = punchProp.nodesLevels;
for i = 1:punchProp.nodesLevels
    if punchedLeft == 0 && punchPropLeft.behaviourLevels(i,4) == 1
        belowBasisIndexLeft = i-1;
        punchedLeft = 1;
    end
    if punchedLeft == 1 && punchPropLeft.behaviourLevels(i,4) == 0
        aboveBasisIndexLeft = i;
        punchedLeft = 0;
    end
end
if belowBasisIndexLeft == 0
    belowBasisHeightLeft = 0;
else
    belowBasisHeightLeft = properties.nodes(punchPropLeft.behaviourLevels(belowBasisIndexLeft,1),3);
end
aboveBasisHeightLeft = properties.nodes(punchPropLeft.behaviourLevels(aboveBasisIndexLeft,1),3);

ddeltaLeft = ddeltaShip * (cos(punchPropLeft.angle)*cos(data.shipTrajectory*pi/180) + sin(punchPropLeft.angle)*sin(data.shipTrajectory*pi/180));

ddeltaLeftVector = zeros(punchProp.nodesLevels,1);
if data.impactedZone == 0 % on node
    
    impactedNodeHeight = properties.nodes(punchPropLeft.behaviourLevels(punchProp.impactedNodesIndex,1),3);
    for i = 1:punchProp.nodesLevels
        if punchPropLeft.behaviourLevels(i,4) == 1
            nodeHeight = properties.nodes(punchPropLeft.behaviourLevels(i,1),3);
            if i == punchProp.impactedNodesIndex
                ddeltaLeftVector(i) = ddeltaLeft;
            elseif i < punchProp.impactedNodesIndex
                ddeltaLeftVector(i) = ddeltaLeft * (nodeHeight-belowBasisHeightLeft)/(impactedNodeHeight-belowBasisHeightLeft);
            else
                ddeltaLeftVector(i) = ddeltaLeft * (aboveBasisHeightLeft-nodeHeight)/(aboveBasisHeightLeft-impactedNodeHeight);
            end
            if isinf(ddeltaLeftVector(i)) || isnan(ddeltaLeftVector(i))
                ddeltaLeftVector(i) = 0;
            end
            if ddeltaLeftVector(i) < 0
                ddeltaLeftVector(i) = 0;
            end
        end
    end
    
else % on element
    
    if punchPropLeft.behaviourLevels(punchProp.impactedNodesIndex(1),4) == 1 && punchPropLeft.behaviourLevels(punchProp.impactedNodesIndex(2),4) == 1
        dispBelowMag = 1;
        dispAboveMag = 1;
    else
        dispBelowMag = (punchProp.lprime+punchProp.lprimeprime)/punchProp.lprime;
        dispAboveMag = (punchProp.lprime+punchProp.lprimeprime)/punchProp.lprimeprime;
    end
    
    impactedNodeBelowHeight = properties.nodes(punchPropLeft.behaviourLevels(punchProp.impactedNodesIndex(1),1),3);
    impactedNodeAboveHeight = properties.nodes(punchPropLeft.behaviourLevels(punchProp.impactedNodesIndex(2),1),3);
    for i = 1:punchProp.nodesLevels
        if punchPropLeft.behaviourLevels(i,4) == 1
            nodeHeight = properties.nodes(punchPropLeft.behaviourLevels(i,1),3);
            if i <= punchProp.impactedNodesIndex(1)
                ddeltaLeftVector(i) = ddeltaLeft * dispBelowMag * (nodeHeight-belowBasisHeightLeft)/(impactedNodeBelowHeight-belowBasisHeightLeft);
            elseif i >= punchProp.impactedNodesIndex(2)
                ddeltaLeftVector(i) = ddeltaLeft * dispAboveMag * (aboveBasisHeightLeft-nodeHeight)/(aboveBasisHeightLeft-impactedNodeAboveHeight);
            end
            if isinf(ddeltaLeftVector(i)) || isnan(ddeltaLeftVector(i))
                ddeltaLeftVector(i) = 0;
            end
        end
    end
    
end

% ----------
% Right
% ----------

punchedRight = 0;
belowBasisIndexRight = 0;
aboveBasisIndexRight = punchProp.nodesLevels;
for i = 1:punchProp.nodesLevels
    if punchedRight == 0 && punchPropRight.behaviourLevels(i,4) == 1
        belowBasisIndexRight = i-1;
        punchedRight = 1;
    end
    if punchedRight == 1 && punchPropRight.behaviourLevels(i,4) == 0
        aboveBasisIndexRight = i;
        punchedRight = 0;
    end
end
if belowBasisIndexRight == 0
    belowBasisHeightRight = 0;
else
    belowBasisHeightRight = properties.nodes(punchPropRight.behaviourLevels(belowBasisIndexRight,1),3);
end
aboveBasisHeightRight = properties.nodes(punchPropRight.behaviourLevels(aboveBasisIndexRight,1),3);

ddeltaRight = ddeltaShip * (cos(punchPropRight.angle)*cos(data.shipTrajectory*pi/180) + sin(punchPropRight.angle)*sin(data.shipTrajectory*pi/180));

ddeltaRightVector = zeros(punchProp.nodesLevels,1);
if data.impactedZone == 0 % on node
    
    impactedNodeHeight = properties.nodes(punchPropRight.behaviourLevels(punchProp.impactedNodesIndex,1),3);
    for i = 1:punchProp.nodesLevels
        if punchPropRight.behaviourLevels(i,4) == 1
            nodeHeight = properties.nodes(punchPropRight.behaviourLevels(i,1),3);
            if i == punchProp.impactedNodesIndex
                ddeltaRightVector(i) = ddeltaRight;
            elseif i < punchProp.impactedNodesIndex
                ddeltaRightVector(i) = ddeltaRight * (nodeHeight-belowBasisHeightRight)/(impactedNodeHeight-belowBasisHeightRight);
            else
                ddeltaRightVector(i) = ddeltaRight * (aboveBasisHeightRight-nodeHeight)/(aboveBasisHeightRight-impactedNodeHeight);
            end
            if isinf(ddeltaRightVector(i)) || isnan(ddeltaRightVector(i))
                ddeltaRightVector(i) = 0;
            end
            if ddeltaRightVector(i) < 0
                ddeltaRightVector(i) = 0;
            end
        end
    end
    
else % on element
    
    if punchPropRight.behaviourLevels(punchProp.impactedNodesIndex(1),4) == 1 && punchPropRight.behaviourLevels(punchProp.impactedNodesIndex(2),4) == 1
        dispBelowMag = 1;
        dispAboveMag = 1;
    else
        dispBelowMag = (punchProp.lprime+punchProp.lprimeprime)/punchProp.lprime;
        dispAboveMag = (punchProp.lprime+punchProp.lprimeprime)/punchProp.lprimeprime;
    end
    
    impactedNodeBelowHeight = properties.nodes(punchPropRight.behaviourLevels(punchProp.impactedNodesIndex(1),1),3);
    impactedNodeAboveHeight = properties.nodes(punchPropRight.behaviourLevels(punchProp.impactedNodesIndex(2),1),3);
    for i = 1:punchProp.nodesLevels
        if punchPropRight.behaviourLevels(i,4) == 1
            nodeHeight = properties.nodes(punchPropRight.behaviourLevels(i,1),3);
            if i <= punchProp.impactedNodesIndex(1)
                ddeltaRightVector(i) = ddeltaRight * dispBelowMag * (nodeHeight-belowBasisHeightRight)/(impactedNodeBelowHeight-belowBasisHeightRight);
            elseif i >= punchProp.impactedNodesIndex(2)
                ddeltaRightVector(i) = ddeltaRight * dispAboveMag * (aboveBasisHeightRight-nodeHeight)/(aboveBasisHeightRight-impactedNodeAboveHeight);
            end
            if isinf(ddeltaRightVector(i)) || isnan(ddeltaRightVector(i))
                ddeltaRightVector(i) = 0;
            end
        end
    end
    
end

%% Compute punching force

% ----------
% Left
% ----------

FLeft = 0;
for i = 1:punchProp.nodesLevels
    if punchPropLeft.behaviourLevels(i,4) == 1
        ddeltaPunch = ddeltaLeftVector(i);
        punchImpactActivated = 0;
        punchRearActivated = 0;
        
        % punching on impacted leg
        if punchPropLeft.behaviourLevels(i,2) == 1
            punchImpactActivated = 1;
            deltaPunchImpact = punchPropLeft.behaviourNodes(i,3) + ddeltaPunch;
            if deltaPunchImpact > punchProp.maxDe*data.De(1)
                deltaPunchImpact = punchProp.maxDe*data.De(1);
            end
            [resultPunchImpactLeft,punchModeNodesImpactLeft] = solvePunch2(punchPropLeft,i,deltaPunchImpact,ddeltaPunch);
            FnodeImpact = resultPunchImpactLeft.PtotPunch;
        end
        
        % punching on rear leg
        if punchPropLeft.behaviourLevels(i,3) == 1
            punchRearActivated = 1;
            if punchPropLeft.behaviourNodes(i,5) == 1
                deltaPunchRearBottom = punchPropLeft.behaviourNodes(i,6) + ddeltaPunch; % bottom
                if deltaPunchRearBottom > punchProp.maxDe*data.De(1)
                    deltaPunchRearBottom = punchProp.maxDe*data.De(1);
                end
                [resultPunchRearBottomLeft,punchModeNodesRearBottomLeft] = solvePunch2(punchPropLeft,i+punchProp.nodesLevels-1,deltaPunchRearBottom,ddeltaPunch);
            else
                deltaPunchRearBottom = 0;
                resultPunchRearBottomLeft.PtotPunch = 0;
                punchModeNodesRearBottomLeft = ones(1,4);
            end
            
            if punchPropLeft.behaviourNodes(i,8) == 1
                deltaPunchRearTop = punchPropLeft.behaviourNodes(i,9) + ddeltaPunch; % top
                if deltaPunchRearTop > punchProp.maxDe*data.De(1)
                    deltaPunchRearTop = punchProp.maxDe*data.De(1);
                end
                [resultPunchRearTopLeft,punchModeNodesRearTopLeft] = solvePunch2(punchPropLeft,i+punchProp.nodesLevels+1,deltaPunchRearTop,ddeltaPunch);
            else
                deltaPunchRearTop = 0;
                resultPunchRearTopLeft.PtotPunch = 0;
                punchModeNodesRearTopLeft = ones(1,4);
            end
            
            FnodeRear = resultPunchRearBottomLeft.PtotPunch + resultPunchRearTopLeft.PtotPunch;
        end
        
        % comparison of forces
        if punchImpactActivated == 1 && punchRearActivated == 1
            if FnodeImpact < FnodeRear
                Fnode = FnodeImpact;
                punchPropLeft.behaviourNodes(i,3) = deltaPunchImpact;
                punchPropLeft.punchModeNodes(i,:) = punchModeNodesImpactLeft(1,:);
            else
                Fnode = FnodeRear;
                punchPropLeft.behaviourNodes(i,6) = deltaPunchRearBottom;
                punchPropLeft.behaviourNodes(i,9) = deltaPunchRearTop;
                punchPropLeft.punchModeNodes(i+punchProp.nodesLevels-1,:) = punchModeNodesRearBottomLeft(1,:);
                punchPropLeft.punchModeNodes(i+punchProp.nodesLevels+1,:) = punchModeNodesRearTopLeft(1,:);
            end
            
        elseif punchImpactActivated == 1 && punchRearActivated == 0
            Fnode = FnodeImpact;
            punchPropLeft.behaviourNodes(i,3) = deltaPunchImpact;
            punchPropLeft.punchModeNodes(i,:) = punchModeNodesImpactLeft(1,:);
            
        elseif punchImpactActivated == 0 && punchRearActivated == 1
            Fnode = FnodeRear;
            punchPropLeft.behaviourNodes(i,6) = deltaPunchRearBottom;
            punchPropLeft.behaviourNodes(i,9) = deltaPunchRearTop;
            punchPropLeft.punchModeNodes(i+punchProp.nodesLevels-1,:) = punchModeNodesRearBottomLeft(1,:);
            punchPropLeft.punchModeNodes(i+punchProp.nodesLevels+1,:) = punchModeNodesRearTopLeft(1,:);
            
        end
        
        
        nodeHeight = properties.nodes(punchPropLeft.behaviourLevels(i,1),3);
        if nodeHeight < data.shipHeight
            ampl = (nodeHeight-belowBasisHeightLeft) / (data.shipHeight-belowBasisHeightLeft);
        else
            ampl = (aboveBasisHeightLeft-nodeHeight) / (aboveBasisHeightLeft-data.shipHeight);
        end
        FLeft = FLeft + Fnode*ampl;
    end
    
end

% ----------
% Right
% ----------

FRight = 0;
for i = 1:punchProp.nodesLevels
    if punchPropRight.behaviourLevels(i,4) == 1
        ddeltaPunch = ddeltaRightVector(i);
        punchImpactActivated = 0;
        punchRearActivated = 0;
        
        % punching on impacted leg
        if punchPropRight.behaviourLevels(i,2) == 1
            punchImpactActivated = 1;
            deltaPunchImpact = punchPropRight.behaviourNodes(i,3) + ddeltaPunch;
            if deltaPunchImpact > punchProp.maxDe*data.De(1)
                deltaPunchImpact = punchProp.maxDe*data.De(1);
            end
            [resultPunchImpactRight,punchModeNodesImpactRight] = solvePunch2(punchPropRight,i,deltaPunchImpact,ddeltaPunch);
            FnodeImpact = resultPunchImpactRight.PtotPunch;
        end
        
        % punching on rear leg
        if punchPropRight.behaviourLevels(i,3) == 1
            punchRearActivated = 1;
            if punchPropRight.behaviourNodes(i,5) == 1
                deltaPunchRearBottom = punchPropRight.behaviourNodes(i,6) + ddeltaPunch; % bottom
                if deltaPunchRearBottom > punchProp.maxDe*data.De(1)
                    deltaPunchRearBottom = punchProp.maxDe*data.De(1);
                end
                [resultPunchRearBottomRight,punchModeNodesRearBottomRight] = solvePunch2(punchPropRight,i+punchProp.nodesLevels-1,deltaPunchRearBottom,ddeltaPunch);
            else
                deltaPunchRearBottom = 0;
                resultPunchRearBottomRight.PtotPunch = 0;
                punchModeNodesRearBottomRight = ones(1,4);
            end
            
            if punchPropRight.behaviourNodes(i,8) == 1
                deltaPunchRearTop = punchPropRight.behaviourNodes(i,9) + ddeltaPunch; % top
                if deltaPunchRearTop > punchProp.maxDe*data.De(1)
                    deltaPunchRearTop = punchProp.maxDe*data.De(1);
                end
                [resultPunchRearTopRight,punchModeNodesRearTopRight] = solvePunch2(punchPropRight,i+punchProp.nodesLevels+1,deltaPunchRearTop,ddeltaPunch);
            else
                deltaPunchRearTop = 0;
                resultPunchRearTopRight.PtotPunch = 0;
                punchModeNodesRearTopRight = ones(1,4);
            end
            
            FnodeRear = resultPunchRearBottomRight.PtotPunch + resultPunchRearTopRight.PtotPunch;
        end
        
        % comparison of forces
        if punchImpactActivated == 1 && punchRearActivated == 1
            if FnodeImpact < FnodeRear
                Fnode = FnodeImpact;
                punchPropRight.behaviourNodes(i,3) = deltaPunchImpact;
                punchPropRight.punchModeNodes(i,:) = punchModeNodesImpactRight(1,:);
            else
                Fnode = FnodeRear;
                punchPropRight.behaviourNodes(i,6) = deltaPunchRearBottom;
                punchPropRight.behaviourNodes(i,9) = deltaPunchRearTop;
                punchPropRight.punchModeNodes(i+punchProp.nodesLevels-1,:) = punchModeNodesRearBottomRight(1,:);
                punchPropRight.punchModeNodes(i+punchProp.nodesLevels+1,:) = punchModeNodesRearTopRight(1,:);
            end
            
        elseif punchImpactActivated == 1 && punchRearActivated == 0
            Fnode = FnodeImpact;
            punchPropRight.behaviourNodes(i,3) = deltaPunchImpact;
            punchPropRight.punchModeNodes(i,:) = punchModeNodesImpactRight(1,:);
            
        elseif punchImpactActivated == 0 && punchRearActivated == 1
            Fnode = FnodeRear;
            punchPropRight.behaviourNodes(i,6) = deltaPunchRearBottom;
            punchPropRight.behaviourNodes(i,9) = deltaPunchRearTop;
            punchPropRight.punchModeNodes(i+punchProp.nodesLevels-1,:) = punchModeNodesRearBottomRight(1,:);
            punchPropRight.punchModeNodes(i+punchProp.nodesLevels+1,:) = punchModeNodesRearTopRight(1,:);
            
        end
        
        
        nodeHeight = properties.nodes(punchPropRight.behaviourLevels(i,1),3);
        if nodeHeight < data.shipHeight
            ampl = (nodeHeight-belowBasisHeightRight) / (data.shipHeight-belowBasisHeightRight);
        else
            ampl = (aboveBasisHeightRight-nodeHeight) / (aboveBasisHeightRight-data.shipHeight);
        end
        FRight = FRight + Fnode*ampl;
    end
    
end

% if FRight == 0
%     FRight = 10*punchPropRight.nodeProp(2,13)*10^(-6);
% end

%% Reduction of section

for i = 1:2*punchProp.nodesLevels
    if punchPropLeft.angle == 0
        % Plastic moment reduction
        punchProp.redSect(punchPropLeft.leg(i,1),1) = punchPropLeft.punchModeNodes(i,3);
        punchProp.redSect(punchPropRight.leg(i,1),2) = punchPropRight.punchModeNodes(i,3);
        % Inertia reduction
        punchProp.redSect(punchPropLeft.leg(i,1),3) = punchPropLeft.punchModeNodes(i,4);
        punchProp.redSect(punchPropRight.leg(i,1),4) = punchPropRight.punchModeNodes(i,4);
    else
        % Plastic moment reduction
        punchProp.redSect(punchPropRight.leg(i,1),1) = punchPropRight.punchModeNodes(i,3);
        punchProp.redSect(punchPropLeft.leg(i,1),2) = punchPropLeft.punchModeNodes(i,3);
        % Inertia reduction
        punchProp.redSect(punchPropRight.leg(i,1),3) = punchPropRight.punchModeNodes(i,4);
        punchProp.redSect(punchPropLeft.leg(i,1),4) = punchPropLeft.punchModeNodes(i,4);
    end
end

%% Results

resultPunch.FLeft = FLeft;
resultPunch.FRight = FRight;
FpunchX = FLeft*cos(punchPropLeft.angle) + FRight*cos(punchPropRight.angle);
FpunchY = FLeft*sin(punchPropLeft.angle) + FRight*sin(punchPropRight.angle);
resultPunch.FPunch = [FpunchX FpunchY];