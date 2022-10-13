function[punchPropLeft,punchPropRight] = punchingNode2(properties,data,forceElem,punchProp,punchPropLeft,punchPropRight)

%% Traction or compression in all braces

for j = 1:size(punchPropLeft.leg,1)
    punchPropLeft.leg(j,3) = sign(forceElem(1,punchPropLeft.leg(j,2)));
    punchPropLeft.leg(j,5) = sign(forceElem(1,punchPropLeft.leg(j,4)));
end
for j = 1:size(punchPropRight.leg,1)
    punchPropRight.leg(j,3) = sign(forceElem(1,punchPropRight.leg(j,2)));
    punchPropRight.leg(j,5) = sign(forceElem(1,punchPropRight.leg(j,4)));
end

%% Properties of nodes for punching

nodesLevels = length(data.h) - 1;

dirAllLegs = [properties.angle(1,:,2);
    properties.angle(1,:,nodesLevels+2);
    properties.angle(1,:,2*nodesLevels+2);
    properties.angle(1,:,3*nodesLevels+2)];

% Left
% 12 columns : Rleg, L, tp, fy, alpha, xi, Db, L1, L2, NB, Gap, multiplcation with regard to delta, F0
punchPropLeft.nodeProp = zeros(size(punchPropLeft.leg,1),13);

for i = 1:size(punchPropLeft.leg,1)
    punchPropLeft.nodeProp(i,1) = (data.De(1) - data.t(1)) / 2;
    punchPropLeft.nodeProp(i,3) = data.t(1);
    punchPropLeft.nodeProp(i,4) = data.fy(1);
    
    % Look for the compressed brace
    if punchPropLeft.leg(i,3) > 0
        dirBrace = properties.angle(1,:,punchPropLeft.leg(i,2));
        dirLeg = dirAllLegs(punchPropLeft.leg(i,6),:);
        alpha = acos(dirBrace*dirLeg');
        if alpha > pi/2 && alpha <= pi
            alpha = pi - alpha;
        elseif alpha > pi && alpha <= 3*pi/2
            alpha = alpha - pi;
        elseif alpha > 3*pi/2 && alpha <= 2*pi
            alpha = 2*pi - alpha;
        end
    elseif punchPropLeft.leg(i,5) > 0
        dirBrace = properties.angle(1,:,punchPropLeft.leg(i,4));
        dirLeg = dirAllLegs(punchPropLeft.leg(i,6),:);
        alpha = acos(dirBrace*dirLeg');
        if alpha > pi/2 && alpha <= pi
            alpha = pi - alpha;
        elseif alpha > pi && alpha <= 3*pi/2
            alpha = alpha - pi;
        elseif alpha > 3*pi/2 && alpha <= 2*pi
            alpha = 2*pi - alpha;
        end
    else
        punchPropLeft.nodeProp(i,5) = -100;
    end
    
    punchPropLeft.nodeProp(i,5) = pi/2 - alpha;
    
    punchPropLeft.nodeProp(i,6) = 1;
    punchPropLeft.nodeProp(i,7) = data.De(4);
    if i <= punchProp.nodesLevels
        punchPropLeft.nodeProp(i,8) = punchProp.LLegs(i);
    else
        punchPropLeft.nodeProp(i,8) = punchProp.LLegs(i-punchProp.nodesLevels);
    end
    if i < punchProp.nodesLevels
        punchPropLeft.nodeProp(i,9) = punchProp.LLegs(i+1);
    elseif i > punchProp.nodesLevels && i ~= 2*punchProp.nodesLevels
        punchPropLeft.nodeProp(i,9) = punchProp.LLegs(i+1-punchProp.nodesLevels);
    else
        punchPropLeft.nodeProp(i,9) = punchProp.LLegs(punchProp.nodesLevels);
    end
    punchPropLeft.nodeProp(i,2) = punchPropLeft.nodeProp(i,8) + punchPropLeft.nodeProp(i,9);
    punchPropLeft.nodeProp(i,10) = 1;
    punchPropLeft.nodeProp(i,11) = 0;
    
    if punchPropLeft.leg(i,1) == 33
        punchPropLeft.nodeProp(i,12) = (punchProp.lprime + punchProp.lprimeprime)/punchProp.lprime;
    else
        punchPropLeft.nodeProp(i,12) = (punchProp.lprime + punchProp.lprimeprime)/punchProp.lprimeprime;
    end
    
    [resultPunchLeft,punchPropLeft] = solvePunch(punchPropLeft,i,0,0);
    
    punchPropLeft.nodeProp(i,13) = resultPunchLeft.Ptotpunch * 10^6;
end

% Right
% 12 columns : Rleg, L, tp, fy, alpha, xi, Db, L1, L2, NB, Gap, multiplcation with regard to delta, F0
punchPropRight.nodeProp = zeros(size(punchPropRight.leg,1),13);

for i = 1:size(punchPropRight.leg,1)
    punchPropRight.nodeProp(i,1) = (data.De(1) - data.t(1)) / 2;
    punchPropRight.nodeProp(i,3) = data.t(1);
    punchPropRight.nodeProp(i,4) = data.fy(1);
    
    % Look for the compressed brace
    if punchPropRight.leg(i,3) > 0
        dirBrace = properties.angle(1,:,punchPropRight.leg(i,2));
        dirLeg = dirAllLegs(punchPropRight.leg(i,6),:);
        alpha = acos(dirBrace*dirLeg');
        if alpha > pi/2 && alpha <= pi
            alpha = pi - alpha;
        elseif alpha > pi && alpha <= 3*pi/2
            alpha = alpha - pi;
        elseif alpha > 3*pi/2 && alpha <= 2*pi
            alpha = 2*pi - alpha;
        end
    elseif punchPropRight.leg(i,5) > 0
        dirBrace = properties.angle(1,:,punchPropRight.leg(i,4));
        dirLeg = dirAllLegs(punchPropRight.leg(i,6),:);
        alpha = acos(dirBrace*dirLeg');
        if alpha > pi/2 && alpha <= pi
            alpha = pi - alpha;
        elseif alpha > pi && alpha <= 3*pi/2
            alpha = alpha - pi;
        elseif alpha > 3*pi/2 && alpha <= 2*pi
            alpha = 2*pi - alpha;
        end
    else
        punchPropRight.nodeProp(i,5) = -100;
    end
    
    punchPropRight.nodeProp(i,5) = pi/2 - alpha;
    
    punchPropRight.nodeProp(i,6) = 1;
    punchPropRight.nodeProp(i,7) = data.De(4);
    if i <= punchProp.nodesLevels
        punchPropRight.nodeProp(i,8) = punchProp.LLegs(i);
    else
        punchPropRight.nodeProp(i,8) = punchProp.LLegs(i-punchProp.nodesLevels);
    end
    if i < punchProp.nodesLevels
        punchPropRight.nodeProp(i,9) = punchProp.LLegs(i+1);
    elseif i > punchProp.nodesLevels && i ~= 2*punchProp.nodesLevels
        punchPropRight.nodeProp(i,9) = punchProp.LLegs(i+1-punchProp.nodesLevels);
    else
        punchPropRight.nodeProp(i,9) = punchProp.LLegs(punchProp.nodesLevels);
    end
    punchPropRight.nodeProp(i,2) = punchPropRight.nodeProp(i,8) + punchPropRight.nodeProp(i,9);
    punchPropRight.nodeProp(i,10) = 1;
    punchProp.nodePropRight(i,11) = 0;
    
    if punchPropRight.leg(i,1) == 33
        punchPropRight.nodeProp(i,12) = (punchProp.lprime + punchProp.lprimeprime)/punchProp.lprime;
    else
        punchPropRight.nodeProp(i,12) = (punchProp.lprime + punchProp.lprimeprime)/punchProp.lprimeprime;
    end
    
    [resultPunchRight,punchPropRight] = solvePunch(punchPropRight,i,0,0);
    
    punchPropRight.nodeProp(i,13) = resultPunchRight.Ptotpunch * 10^6;
end

%% Behaviour of the node
% Column 1 : number of node on the impacted leg
% Column 2 : activated or not ? (yes if F_brace > F_Rd,punch)
% Column 3 : delta of punching on that node
% Columns (1,2,3)+i : the same for opposite nodes (remains 0 if brace in
% traction)

punchPropLeft.behaviourNodes = zeros(nodesLevels,9);
punchPropRight.behaviourNodes = zeros(nodesLevels,9);

punchPropLeft.behaviourNodes(:,1) = punchPropLeft.leg(1:nodesLevels,1);
punchPropRight.behaviourNodes(:,1) = punchPropRight.leg(1:nodesLevels,1);

for i = 1:nodesLevels
    % Left
    if punchPropLeft.leg(i,3) == 1
        if i == 1
            punchPropLeft.behaviourNodes(i,4) = punchProp.posNodes(i,punchPropLeft.leg(nodesLevels+i,6));
        else
            punchPropLeft.behaviourNodes(i,4) = punchProp.posNodes(i-1,punchPropLeft.leg(nodesLevels+i,6));
        end
    end
    if punchPropLeft.leg(i,5) == 1
        if i == nodesLevels
            punchPropLeft.behaviourNodes(i,7) = punchProp.posNodes(i,punchPropLeft.leg(nodesLevels+i,6));
        else
            punchPropLeft.behaviourNodes(i,7) = punchProp.posNodes(i+1,punchPropLeft.leg(nodesLevels+i,6));
        end
    end
    
    % Right
    if punchPropRight.leg(i,3) == 1
        if i == 1
            punchPropRight.behaviourNodes(i,4) = punchProp.posNodes(i,punchPropRight.leg(nodesLevels+i,6));
        else
            punchPropRight.behaviourNodes(i,4) = punchProp.posNodes(i-1,punchPropRight.leg(nodesLevels+i,6));
        end
    end
    if punchPropRight.leg(i,5) == 1
        if i == nodesLevels
            punchPropRight.behaviourNodes(i,7) = punchProp.posNodes(i,punchPropRight.leg(nodesLevels+i,6));
        else
            punchPropRight.behaviourNodes(i,7) = punchProp.posNodes(i+1,punchPropRight.leg(nodesLevels+i,6));
        end
    end
end

%% Behaviour of levels

punchPropLeft.behaviourLevels = zeros(nodesLevels,4);
punchPropRight.behaviourLevels = zeros(nodesLevels,4);
punchPropLeft.behaviourLevels(:,1) = punchPropLeft.behaviourNodes(:,1);
punchPropRight.behaviourLevels(:,1) = punchPropRight.behaviourNodes(:,1);