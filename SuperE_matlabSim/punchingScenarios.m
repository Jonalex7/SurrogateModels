function[punchProp,punchPropLeft,punchPropRight,resultPunch] = punchingScenarios(properties,data,force,forceElem,punchProp,punchPropLeft,punchPropRight)

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

%% Localise punching for each level

% 1st column : number of node on impacted leg
% 2nd column : displacement allowed on impacted leg ?
% 3rd column : displacement allowed on back leg ?

for i = 1:punchProp.nodesLevels
    % Left
    
    if punchPropLeft.behaviourNodes(i,2) == 1
        punchPropLeft.behaviourLevels(i,2) = 1;
    end
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
    punchPropLeft.behaviourLevels(i,4) = max([punchPropLeft.behaviourLevels(i,2) punchPropLeft.behaviourLevels(i,3)]);
    if max(punchPropLeft.behaviourLevels(:,4)) == 0
        punchPropLeft.activated = 0;
    else
        punchPropLeft.activated = 1;
    end
    
    % Right
    
    if punchPropRight.behaviourNodes(i,2) == 1
        punchPropRight.behaviourLevels(i,2) = 1;
    end
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
    punchPropRight.behaviourLevels(i,4) = max([punchPropRight.behaviourLevels(i,2) punchPropRight.behaviourLevels(i,3)]);
    if max(punchPropRight.behaviourLevels(:,4)) == 0
        punchPropRight.activated = 0;
    else
        punchPropRight.activated = 1;
    end
end

%% Compute punching

if data.impactedZone == 0 %On node
    
    
    
elseif data.impactedZone == 1 %On element
    
    % --------------------------------------------------------------------
    % --------------------------------------------------------------------
    % --------------------------------------------------------------------
    % --------------------------------------------------------------------
    % ----
    % Left
    % ----
    % --------------------------------------------------------------------
    % --------------------------------------------------------------------
    % --------------------------------------------------------------------
    % --------------------------------------------------------------------
    
    % Below the impact
    indexBelow = punchProp.impactedNodesIndex(1);
    numberNodesBelow = 0;
    while indexBelow-numberNodesBelow > 0 && punchPropLeft.behaviourLevels(indexBelow-numberNodesBelow,4) == 1
        numberNodesBelow = numberNodesBelow + 1;
    end
    punchPropLeft.numberNodesBelow = numberNodesBelow;
    punchPropLeft.dispNodesBelow = zeros(numberNodesBelow,2);
    punchPropLeft.dispNodesBelow(:,1) = (indexBelow-numberNodesBelow+1:indexBelow)';
    if indexBelow-numberNodesBelow == 0
        belowBasis = 0;
    else
        belowBasis = properties.nodes(punchPropLeft.behaviourLevels(indexBelow-numberNodesBelow,1),3);
    end
    belowHeight = properties.nodes(punchPropLeft.behaviourLevels(indexBelow,1),3);
    belowNodeHeight = zeros(numberNodesBelow,1);
    for i = 1:numberNodesBelow
        belowNodeHeight(i) = properties.nodes(punchPropLeft.behaviourLevels(punchPropLeft.dispNodesBelow(i,1),1),3);
        punchPropLeft.dispNodesBelow(i,2) = (belowNodeHeight(i)-belowBasis)/(belowHeight-belowBasis);
    end
    
    % Above the impact
    indexAbove = punchProp.impactedNodesIndex(2);
    numberNodesAbove = 0;
    while indexAbove+numberNodesAbove <= punchProp.nodesLevels && punchPropLeft.behaviourLevels(indexAbove+numberNodesAbove,4) == 1
        numberNodesAbove = numberNodesAbove + 1;
    end
    punchPropLeft.numberNodesAbove = numberNodesAbove;
    punchPropLeft.dispNodesAbove = zeros(numberNodesAbove,2);
    punchPropLeft.dispNodesAbove(:,1) = (indexAbove:indexAbove+numberNodesAbove-1)';
    if indexAbove+numberNodesAbove >= punchProp.nodesLevels
        aboveBasis = properties.nodes(punchPropLeft.behaviourLevels(punchProp.nodesLevels,1),3);
    else
        aboveBasis = properties.nodes(punchPropLeft.behaviourLevels(indexAbove+numberNodesAbove,1),3);
    end
    aboveHeight = properties.nodes(punchPropLeft.behaviourLevels(indexAbove,1),3);
    aboveNodeHeight = zeros(numberNodesAbove,1);
    for i = 1:numberNodesAbove
        aboveNodeHeight(i) = properties.nodes(punchPropLeft.behaviourLevels(punchPropLeft.dispNodesAbove(i,1),1),3);
        punchPropLeft.dispNodesAbove(i,2) = (aboveBasis-aboveNodeHeight(i))/(aboveBasis-aboveHeight);
    end
    
    % Scenarios
    
    if numberNodesBelow >= 1 && numberNodesAbove >= 1
        dispBelow = norm(force.valueDisp);
        dispAbove = norm(force.valueDisp);
    elseif numberNodesBelow >= 1 && numberNodesAbove == 0
        dispBelow = norm(force.valueDisp) * (punchProp.lprime+punchProp.lprimeprime)/punchProp.lprime;
        dispAbove = 0;
    elseif numberNodesBelow == 0 && numberNodesAbove >= 1
        dispBelow = 0;
        dispAbove = norm(force.valueDisp) * (punchProp.lprime+punchProp.lprimeprime)/punchProp.lprimeprime;
    else
        dispBelow = 0;
        dispAbove = 0;
    end
    
    dispBelow = dispBelow * (cos(punchPropLeft.angle)*cos(data.shipTrajectory*pi/180) + sin(punchPropLeft.angle)*sin(data.shipTrajectory*pi/180));
    dispAbove = dispAbove * (cos(punchPropLeft.angle)*cos(data.shipTrajectory*pi/180) + sin(punchPropLeft.angle)*sin(data.shipTrajectory*pi/180));
    
    % ------------------------------
    % Compute resitant force - below
    % ------------------------------
    FnodeBelow = 0;
    if numberNodesBelow > 0
        
        % for each level on which there is punching
        for i = 1:numberNodesBelow
            
            ddeltaPunch = dispBelow * punchPropLeft.dispNodesBelow(i,2);
            
            % Punching of impacted leg
            if punchPropLeft.behaviourLevels(punchPropLeft.dispNodesBelow(i,1),2) == 1
                punchImpactLeg = 1;
                deltaPunchImpact = punchPropLeft.behaviourNodes(punchPropLeft.dispNodesBelow(i,1),3) + ddeltaPunch;
                if deltaPunchImpact > punchProp.maxDe*data.De(1)
                    deltaPunchImpact = punchProp.maxDe*data.De(1);
                end
                [resultPunch,punchPropLeft] = solvePunch(punchPropLeft,punchPropLeft.dispNodesBelow(i,1),deltaPunchImpact,ddeltaPunch);
                FnodeImpactLeg = resultPunch.Ptotpunch * (belowNodeHeight(i)-belowBasis)/(data.shipHeight-belowBasis);
            else
                punchImpactLeg = 0;
            end
            
            % Punching on rear leg
            if punchPropLeft.behaviourLevels(punchPropLeft.dispNodesBelow(i,1),3) == 1 && punchPropLeft.dispNodesBelow(i,1) ~= 1
                punchRearLeg = 1;
                
                if punchPropLeft.behaviourNodes(punchPropLeft.dispNodesBelow(i,1),4) ~= 0
                    deltaPunchRearDown = punchPropLeft.behaviourNodes(punchPropLeft.dispNodesBelow(i,1),6) + ddeltaPunch;
                    if deltaPunchRearDown > punchProp.maxDe*data.De(1)
                        deltaPunchRearDown = punchProp.maxDe*data.De(1);
                    end
                    [resultPunch,punchPropLeft] = solvePunch(punchPropLeft,punchProp.nodesLevels+punchPropLeft.dispNodesBelow(i,1)-1,deltaPunchRearDown,ddeltaPunch);
                    FnodeRearLegDown = resultPunch.Ptotpunch;
                else
                    deltaPunchRearDown = punchPropRight.behaviourNodes(punchPropLeft.dispNodesBelow(i,1),6);
                    FnodeRearLegDown = 0;
                end
                
                if punchPropLeft.behaviourNodes(punchPropLeft.dispNodesBelow(i,1),7) ~= 0
                    deltaPunchRearUp = punchPropLeft.behaviourNodes(punchPropLeft.dispNodesBelow(i,1),9) + ddeltaPunch;
                    if deltaPunchRearUp > punchProp.maxDe*data.De(1)
                        deltaPunchRearUp = punchProp.maxDe*data.De(1);
                    end
                    [resultPunch,punchPropLeft] = solvePunch(punchPropLeft,punchProp.nodesLevels+punchPropLeft.dispNodesBelow(i,1)+1,deltaPunchRearUp,ddeltaPunch);
                    FnodeRearLegUp = resultPunch.Ptotpunch;
                else
                    deltaPunchRearUp = punchPropRight.behaviourNodes(punchPropRight.dispNodesBelow(i,1),9);
                    FnodeRearLegUp = 0;
                end
                FnodeRearLeg = (FnodeRearLegDown + FnodeRearLegUp) * (belowNodeHeight(i)-belowBasis)/(data.shipHeight-belowBasis);
            else
                punchRearLeg = 0;
            end
            
            % Comparison of both forces
            if punchImpactLeg == 1 && punchRearLeg == 1
                Fnode = min(FnodeImpactLeg,FnodeRearLeg);
                if FnodeImpactLeg < FnodeRearLeg
                    punchPropLeft.behaviourNodes(punchPropLeft.dispNodesBelow(i,1),3) = deltaPunchImpact;
                else
                    punchPropLeft.behaviourNodes(punchPropLeft.dispNodesBelow(i,1),6) = deltaPunchRearDown;
                    punchPropLeft.behaviourNodes(punchPropLeft.dispNodesBelow(i,1),9) = deltaPunchRearUp;
                end
            elseif punchImpactLeg == 1 && punchRearLeg == 0
                Fnode = FnodeImpactLeg;
                punchPropLeft.behaviourNodes(punchPropLeft.dispNodesBelow(i,1),3) = deltaPunchImpact;
            elseif punchImpactLeg == 0 && punchRearLeg == 1
                Fnode = FnodeRearLeg;
                punchPropLeft.behaviourNodes(punchPropLeft.dispNodesBelow(i,1),6) = deltaPunchRearDown;
                punchPropLeft.behaviourNodes(punchPropLeft.dispNodesBelow(i,1),9) = deltaPunchRearUp;
            else
                Fnode = 0;
            end
            
            % For each node
            FnodeBelow = FnodeBelow + Fnode;
            
        end
        
    end
    
    % -------------------------------
    % Compute resistant force - above
    % -------------------------------
    FnodeAbove = 0;
    if numberNodesAbove > 0
        
        % for each level on which there is punching
        for i = 1:numberNodesAbove
            
            ddeltaPunch = dispAbove * punchPropLeft.dispNodesAbove(i,2);
            
            % Punching of impacted leg
            if punchPropLeft.behaviourLevels(punchPropLeft.dispNodesAbove(i,1),2) == 1
                punchImpactLeg = 1;
                deltaPunchImpact = punchPropLeft.behaviourNodes(punchPropLeft.dispNodesAbove(i,1),3) + ddeltaPunch;
                if deltaPunchImpact > punchProp.maxDe*data.De(1)
                    deltaPunchImpact = punchProp.maxDe*data.De(1);
                end
                [resultPunch,punchPropLeft] = solvePunch(punchPropLeft,punchPropLeft.dispNodesAbove(i,1),deltaPunchImpact,ddeltaPunch);
                FnodeImpactLeg = resultPunch.Ptotpunch * (aboveBasis-aboveNodeHeight(i))/(aboveBasis-data.shipHeight);
            else
                punchImpactLeg = 0;
            end
            
            % Punching on rear leg
            if punchPropLeft.behaviourLevels(punchPropLeft.dispNodesAbove(i,1),3) == 1 && punchPropLeft.dispNodesAbove(i,1) ~= punchProp.nodesLevels
                punchRearLeg = 1;
                
                if punchPropLeft.behaviourNodes(punchPropLeft.dispNodesAbove(i,1),4) ~= 0
                    deltaPunchRearDown = punchPropLeft.behaviourNodes(punchPropLeft.dispNodesAbove(i,1),6) + ddeltaPunch;
                    if deltaPunchRearDown > punchProp.maxDe*data.De(1)
                        deltaPunchRearDown = punchProp.maxDe*data.De(1);
                    end
                    [resultPunch,punchPropLeft] = solvePunch(punchPropLeft,punchProp.nodesLevels+punchPropLeft.dispNodesAbove(i,1)-1,deltaPunchRearDown,ddeltaPunch);
                    FnodeRearLegDown = resultPunch.Ptotpunch;
                else
                    deltaPunchRearDown = punchPropRight.behaviourNodes(punchPropLeft.dispNodesAbove(i,1),6);
                    FnodeRearLegDown = 0;
                end
                
                if punchPropLeft.behaviourNodes(punchPropLeft.dispNodesAbove(i,1),7) ~= 0
                    deltaPunchRearUp = punchPropLeft.behaviourNodes(punchPropLeft.dispNodesAbove(i,1),9) + ddeltaPunch;
                    if deltaPunchRearUp > punchProp.maxDe*data.De(1)
                        deltaPunchRearUp = punchProp.maxDe*data.De(1);
                    end
                    [resultPunch,punchPropLeft] = solvePunch(punchPropLeft,punchProp.nodesLevels+punchPropLeft.dispNodesAbove(i,1)+1,deltaPunchRearUp,ddeltaPunch);
                    FnodeRearLegUp = resultPunch.Ptotpunch;
                else
                    deltaPunchRearUp = punchPropLeft.behaviourNodes(punchPropLeft.dispNodesAbove(i,1),9);
                    FnodeRearLegUp = 0;
                end
                FnodeRearLeg = (FnodeRearLegDown + FnodeRearLegUp) * (aboveBasis-aboveNodeHeight(i))/(aboveBasis-data.shipHeight);
            else
                punchRearLeg = 0;
            end
            
            % Comparison of both forces
            if punchImpactLeg == 1 && punchRearLeg == 1
                Fnode = min(FnodeImpactLeg,FnodeRearLeg);
                if FnodeImpactLeg < FnodeRearLeg
                    punchPropLeft.behaviourNodes(punchPropLeft.dispNodesAbove(i,1),3) = deltaPunchImpact;
                else
                    punchPropLeft.behaviourNodes(punchPropLeft.dispNodesAbove(i,1),6) = deltaPunchRearDown;
                    punchPropLeft.behaviourNodes(punchPropLeft.dispNodesAbove(i,1),9) = deltaPunchRearUp;
                end
            elseif punchImpactLeg == 1 && punchRearLeg == 0
                Fnode = FnodeImpactLeg;
                punchPropLeft.behaviourNodes(punchPropLeft.dispNodesAbove(i,1),3) = deltaPunchImpact;
            elseif punchImpactLeg == 0 && punchRearLeg == 1
                Fnode = FnodeRearLeg;
                punchPropLeft.behaviourNodes(punchPropLeft.dispNodesAbove(i,1),6) = deltaPunchRearDown;
                punchPropLeft.behaviourNodes(punchPropLeft.dispNodesAbove(i,1),9) = deltaPunchRearUp;
            else
                Fnode = 0;
            end
            
            % For each node
            FnodeAbove = FnodeAbove + Fnode;
            
        end
        
    end
    
%     FnodeBelow
%     FnodeAbove
    
    FLeft = FnodeBelow + FnodeAbove;
%     if FLeft == 0
%         FLeft = 5*punchPropRight.nodeProp(2,13)*10^(-6);
%     end
    
    % --------------------------------------------------------------------
    % --------------------------------------------------------------------
    % --------------------------------------------------------------------
    % --------------------------------------------------------------------
    % -----
    % Right
    % -----
    % --------------------------------------------------------------------
    % --------------------------------------------------------------------
    % --------------------------------------------------------------------
    % --------------------------------------------------------------------
    
    % Below the impact
    indexBelow = punchProp.impactedNodesIndex(1);
    numberNodesBelow = 0;
    while indexBelow-numberNodesBelow > 0 && punchPropRight.behaviourLevels(indexBelow-numberNodesBelow,4) == 1
        numberNodesBelow = numberNodesBelow + 1;
    end
    punchPropRight.numberNodesBelow = numberNodesBelow;
    punchPropRight.dispNodesBelow = zeros(numberNodesBelow,2);
    punchPropRight.dispNodesBelow(:,1) = (indexBelow-numberNodesBelow+1:indexBelow)';
    if indexBelow-numberNodesBelow == 0
        belowBasis = 0;
    else
        belowBasis = properties.nodes(punchPropRight.behaviourLevels(indexBelow-numberNodesBelow,1),3);
    end
    belowHeight = properties.nodes(punchPropRight.behaviourLevels(indexBelow,1),3);
    belowNodeHeight = zeros(numberNodesBelow,1);
    for i = 1:numberNodesBelow
        belowNodeHeight(i) = properties.nodes(punchPropRight.behaviourLevels(punchPropRight.dispNodesBelow(i,1),1),3);
        punchPropRight.dispNodesBelow(i,2) = (belowNodeHeight(i)-belowBasis)/(belowHeight-belowBasis);
    end
    
    % Above the impact
    indexAbove = punchProp.impactedNodesIndex(2);
    numberNodesAbove = 0;
    while indexAbove+numberNodesAbove <= punchProp.nodesLevels && punchPropRight.behaviourLevels(indexAbove+numberNodesAbove,4) == 1
        numberNodesAbove = numberNodesAbove + 1;
    end
    punchPropRight.numberNodesAbove = numberNodesAbove;
    punchPropRight.dispNodesAbove = zeros(numberNodesAbove,2);
    punchPropRight.dispNodesAbove(:,1) = (indexAbove:indexAbove+numberNodesAbove-1)';
    if indexAbove+numberNodesAbove >= punchProp.nodesLevels
        aboveBasis = properties.nodes(punchPropRight.behaviourLevels(punchProp.nodesLevels,1),3);
    else
        aboveBasis = properties.nodes(punchPropRight.behaviourLevels(indexAbove+numberNodesAbove,1),3);
    end
    aboveHeight = properties.nodes(punchPropRight.behaviourLevels(indexAbove,1),3);
    aboveNodeHeight = zeros(numberNodesAbove,1);
    for i = 1:numberNodesAbove
        aboveNodeHeight(i) = properties.nodes(punchPropRight.behaviourLevels(punchPropRight.dispNodesAbove(i,1),1),3);
        punchPropRight.dispNodesAbove(i,2) = (aboveBasis-aboveNodeHeight(i))/(aboveBasis-aboveHeight);
    end
    
    % Scenarios
    
    if numberNodesBelow >= 1 && numberNodesAbove >= 1
        dispBelow = norm(force.valueDisp);
        dispAbove = norm(force.valueDisp);
    elseif numberNodesBelow >= 1 && numberNodesAbove == 0
        dispBelow = norm(force.valueDisp) * (punchProp.lprime+punchProp.lprimeprime)/punchProp.lprime;
        dispAbove = 0;
    elseif numberNodesBelow == 0 && numberNodesAbove >= 1
        dispBelow = 0;
        dispAbove = norm(force.valueDisp) * (punchProp.lprime+punchProp.lprimeprime)/punchProp.lprimeprime;
    else
        dispBelow = 0;
        dispAbove = 0;
    end
    
    dispBelow = dispBelow * (cos(punchPropRight.angle)*cos(data.shipTrajectory*pi/180) + sin(punchPropRight.angle)*sin(data.shipTrajectory*pi/180));
    dispAbove = dispAbove * (cos(punchPropRight.angle)*cos(data.shipTrajectory*pi/180) + sin(punchPropRight.angle)*sin(data.shipTrajectory*pi/180));
    
    % ------------------------------
    % Compute resitant force - below
    % ------------------------------
    FnodeBelow = 0;
    if numberNodesBelow > 0
        
        % for each level on which there is punching
        for i = 1:numberNodesBelow
            
            ddeltaPunch = dispBelow * punchPropRight.dispNodesBelow(i,2);
            
            % Punching of impacted leg
            if punchPropRight.behaviourLevels(punchPropRight.dispNodesBelow(i,1),2) == 1
                punchImpactLeg = 1;
                deltaPunchImpact = punchPropRight.behaviourNodes(punchPropRight.dispNodesBelow(i,1),3) + ddeltaPunch;
                if deltaPunchImpact > punchProp.maxDe*data.De(1)
                    deltaPunchImpact = punchProp.maxDe*data.De(1);
                end
                [resultPunch,punchPropRight] = solvePunch(punchPropRight,punchPropRight.dispNodesBelow(i,1),deltaPunchImpact,ddeltaPunch);
                FnodeImpactLeg = resultPunch.Ptotpunch * (belowNodeHeight(i)-belowBasis)/(data.shipHeight-belowBasis);
            else
                punchImpactLeg = 0;
            end
            
            % Punching on rear leg
            if punchPropRight.behaviourLevels(punchPropRight.dispNodesBelow(i,1),3) == 1 && punchPropRight.dispNodesBelow(i,1) ~= 1
                punchRearLeg = 1;
                
                if punchPropRight.behaviourNodes(punchPropRight.dispNodesBelow(i,1),4) ~= 0
                    deltaPunchRearDown = punchPropRight.behaviourNodes(punchPropRight.dispNodesBelow(i,1),6) + ddeltaPunch;
                    if deltaPunchRearDown > punchProp.maxDe*data.De(1)
                        deltaPunchRearDown = punchProp.maxDe*data.De(1);
                    end
                    [resultPunch,punchPropRight] = solvePunch(punchPropRight,punchProp.nodesLevels+punchPropRight.dispNodesBelow(i,1)-1,deltaPunchRearDown,ddeltaPunch);
                    FnodeRearLegDown = resultPunch.Ptotpunch;
                else
                    deltaPunchRearDown = punchPropRight.behaviourNodes(punchPropRight.dispNodesBelow(i,1),6);
                    FnodeRearLegDown = 0;
                end
                
                if punchPropRight.behaviourNodes(punchPropRight.dispNodesBelow(i,1),7) ~= 0
                    deltaPunchRearUp = punchPropRight.behaviourNodes(punchPropRight.dispNodesBelow(i,1),9) + ddeltaPunch;
                    if deltaPunchRearUp > punchProp.maxDe*data.De(1)
                        deltaPunchRearUp = punchProp.maxDe*data.De(1);
                    end
                    [resultPunch,punchPropRight] = solvePunch(punchPropRight,punchProp.nodesLevels+punchPropRight.dispNodesBelow(i,1)+1,deltaPunchRearUp,ddeltaPunch);
                    FnodeRearLegUp = resultPunch.Ptotpunch;
                else
                    deltaPunchRearUp = punchPropRight.behaviourNodes(punchPropRight.dispNodesBelow(i,1),9);
                    FnodeRearLegUp = 0;
                end
                FnodeRearLeg = (FnodeRearLegDown + FnodeRearLegUp) * (belowNodeHeight(i)-belowBasis)/(data.shipHeight-belowBasis);
            else
                punchRearLeg = 0;
            end
            
            % Comparison of both forces
            if punchImpactLeg == 1 && punchRearLeg == 1
                Fnode = min(FnodeImpactLeg,FnodeRearLeg);
                if FnodeImpactLeg < FnodeRearLeg
                    punchPropRight.behaviourNodes(punchPropRight.dispNodesBelow(i,1),3) = deltaPunchImpact;
                else
                    punchPropRight.behaviourNodes(punchPropRight.dispNodesBelow(i,1),6) = deltaPunchRearDown;
                    punchPropRight.behaviourNodes(punchPropRight.dispNodesBelow(i,1),9) = deltaPunchRearUp;
                end
            elseif punchImpactLeg == 1 && punchRearLeg == 0
                Fnode = FnodeImpactLeg;
                punchPropRight.behaviourNodes(punchPropRight.dispNodesBelow(i,1),3) = deltaPunchImpact;
            elseif punchImpactLeg == 0 && punchRearLeg == 1
                Fnode = FnodeRearLeg;
                punchPropRight.behaviourNodes(punchPropRight.dispNodesBelow(i,1),6) = deltaPunchRearDown;
                punchPropRight.behaviourNodes(punchPropRight.dispNodesBelow(i,1),9) = deltaPunchRearUp;
            else
                Fnode = 0;
            end
            
            % For each node
            FnodeBelow = FnodeBelow + Fnode;
            
        end
        
    end
    
    % -------------------------------
    % Compute resistant force - above
    % -------------------------------
    FnodeAbove = 0;
    if numberNodesAbove > 0
        
        % for each level on which there is punching
        for i = 1:numberNodesAbove
            
            ddeltaPunch = dispAbove * punchPropRight.dispNodesAbove(i,2);
            
            % Punching of impacted leg
            if punchPropRight.behaviourLevels(punchPropRight.dispNodesAbove(i,1),2) == 1
                punchImpactLeg = 1;
                deltaPunchImpact = punchPropRight.behaviourNodes(punchPropRight.dispNodesAbove(i,1),3) + ddeltaPunch;
                if deltaPunchImpact > punchProp.maxDe*data.De(1)
                    deltaPunchImpact = punchProp.maxDe*data.De(1);
                end
                [resultPunch,punchPropRight] = solvePunch(punchPropRight,punchPropRight.dispNodesAbove(i,1),deltaPunchImpact,ddeltaPunch);
                FnodeImpactLeg = resultPunch.Ptotpunch * (aboveBasis-aboveNodeHeight(i))/(aboveBasis-data.shipHeight);
            else
                punchImpactLeg = 0;
            end
            
            % Punching on rear leg
            if punchPropRight.behaviourLevels(punchPropRight.dispNodesAbove(i,1),3) == 1 && punchPropRight.dispNodesAbove(i,1) ~= punchProp.nodesLevels
                punchRearLeg = 1;
                
                if punchPropRight.behaviourNodes(punchPropRight.dispNodesAbove(i,1),4) ~= 0
                    deltaPunchRearDown = punchPropRight.behaviourNodes(punchPropRight.dispNodesAbove(i,1),6) + ddeltaPunch;
                    if deltaPunchRearDown > punchProp.maxDe*data.De(1)
                        deltaPunchRearDown = punchProp.maxDe*data.De(1);
                    end
                    [resultPunch,punchPropRight] = solvePunch(punchPropRight,punchProp.nodesLevels+punchPropRight.dispNodesAbove(i,1)-1,deltaPunchRearDown,ddeltaPunch);
                    FnodeRearLegDown = resultPunch.Ptotpunch;
                else
                    deltaPunchRearDown = punchPropRight.behaviourNodes(punchPropRight.dispNodesAbove(i,1),6);
                    FnodeRearLegDown = 0;
                end
                
                if punchPropRight.behaviourNodes(punchPropRight.dispNodesAbove(i,1),7) ~= 0
                    deltaPunchRearUp = punchPropRight.behaviourNodes(punchPropRight.dispNodesAbove(i,1),9) + ddeltaPunch;
                    if deltaPunchRearUp > punchProp.maxDe*data.De(1)
                        deltaPunchRearUp = punchProp.maxDe*data.De(1);
                    end
                    [resultPunch,punchPropRight] = solvePunch(punchPropRight,punchProp.nodesLevels+punchPropRight.dispNodesAbove(i,1)+1,deltaPunchRearUp,ddeltaPunch);
                    FnodeRearLegUp = resultPunch.Ptotpunch;
                else
                    deltaPunchRearUp = punchPropRight.behaviourNodes(punchPropRight.dispNodesAbove(i,1),9);
                    FnodeRearLegUp = 0;
                end
                FnodeRearLeg = (FnodeRearLegDown + FnodeRearLegUp) * (aboveBasis-aboveNodeHeight(i))/(aboveBasis-data.shipHeight);
            else
                punchRearLeg = 0;
            end
            
            % Comparison of both forces
            if punchImpactLeg == 1 && punchRearLeg == 1
                Fnode = min(FnodeImpactLeg,FnodeRearLeg);
                if FnodeImpactLeg < FnodeRearLeg
                    punchPropRight.behaviourNodes(punchPropRight.dispNodesAbove(i,1),3) = deltaPunchImpact;
                else
                    punchPropRight.behaviourNodes(punchPropRight.dispNodesAbove(i,1),6) = deltaPunchRearDown;
                    punchPropRight.behaviourNodes(punchPropRight.dispNodesAbove(i,1),9) = deltaPunchRearUp;
                end
            elseif punchImpactLeg == 1 && punchRearLeg == 0
                Fnode = FnodeImpactLeg;
                punchPropRight.behaviourNodes(punchPropRight.dispNodesAbove(i,1),3) = deltaPunchImpact;
            elseif punchImpactLeg == 0 && punchRearLeg == 1
                Fnode = FnodeRearLeg;
                punchPropRight.behaviourNodes(punchPropRight.dispNodesAbove(i,1),6) = deltaPunchRearDown;
                punchPropRight.behaviourNodes(punchPropRight.dispNodesAbove(i,1),9) = deltaPunchRearUp;
            else
                Fnode = 0;
            end
            
            % For each node
            FnodeAbove = FnodeAbove + Fnode;
            
        end
        
    end
    
%     FnodeBelow
%     FnodeAbove
    
    FRight = FnodeBelow + FnodeAbove;
%     if FRight == 0
%         FRight = 5*punchPropRight.nodeProp(2,13)*10^(-6);
%     end
    
end

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