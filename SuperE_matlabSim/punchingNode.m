function[punchProp,punchPropLeft,punchPropRight] = punchingNode(properties,data)

%% Identification of nodes

nodesLevels = length(data.h) - 1;

punchProp.posNodes = zeros(nodesLevels,4);
for i = 1:nodesLevels
    punchProp.posNodes(i,1) = 8*i+1;
    punchProp.posNodes(i,2) = 8*i+2;
    punchProp.posNodes(i,3) = 8*i+3;
    punchProp.posNodes(i,4) = 8*i+4;
end

if data.impactedZone == 1
    pos = properties.origin(data.impactedElement,:);
else
    pos = properties.nodes(data.impactedNode,1:3);
end

if pos(1) < 0 && pos(2) < 0 % leg 1
    startLeg = 1;
    startLeft = 12;
    startImpact = 9;
    startRight = 10;
    
    horizLeft = data.nbLegs + 4;
    horizRight = data.nbLegs + 1;
    startBraceLeft = data.nbLegs + data.nbBottomHoriz + data.nbTopHoriz + (length(data.h)-2)*4*3 + 1;
    startBraceImpact1 = startBraceLeft + 1;
    startBraceImpact2 = data.nbLegs + data.nbBottomHoriz + data.nbTopHoriz + 1;
    startBraceRight = startBraceImpact2 + 1;
elseif pos(1) > 0 && pos(2) < 0 % leg 2
    startLeg = 2;
    startLeft = 9;
    startImpact = 10;
    startRight = 11;
    
    horizLeft = data.nbLegs + 1;
    horizRight = data.nbLegs + 2;
    startBraceLeft = data.nbLegs + data.nbBottomHoriz + data.nbTopHoriz + 1;
    startBraceImpact1 = startBraceLeft + 1;
    startBraceImpact2 = data.nbLegs + data.nbBottomHoriz + data.nbTopHoriz + (length(data.h)-2)*4*1 + 1;
    startBraceRight = startBraceImpact2 + 1;
elseif pos(1) > 0 && pos(2) > 0 % leg 3
    startLeg = 3;
    startLeft = 10;
    startImpact = 11;
    startRight = 12;
    
    horizLeft = data.nbLegs + 2;
    horizRight = data.nbLegs + 3;
    startBraceLeft = data.nbLegs + data.nbBottomHoriz + data.nbTopHoriz + (length(data.h)-2)*4*1 + 1;
    startBraceImpact1 = startBraceLeft + 1;
    startBraceImpact2 = data.nbLegs + data.nbBottomHoriz + data.nbTopHoriz + (length(data.h)-2)*4*2 + 1;
    startBraceRight = startBraceImpact2 + 1;
else % leg 4
    startLeg = 4;
    startLeft = 11;
    startImpact = 12;
    startRight = 9;
    
    horizLeft = data.nbLegs + 3;
    horizRight = data.nbLegs + 4;
    startBraceLeft = data.nbLegs + data.nbBottomHoriz + data.nbTopHoriz + (length(data.h)-2)*4*2 + 1;
    startBraceImpact1 = startBraceLeft + 1;
    startBraceImpact2 = data.nbLegs + data.nbBottomHoriz + data.nbTopHoriz + (length(data.h)-2)*4*3 + 1;
    startBraceRight = startBraceImpact2 + 1;
end

% 1 matrix for each leg :
% - leg at left of the impacted one
% - impacted leg, one side
% - impacted leg, other side
% - leg at right of the impacted one

% 1st column : number of node
% 2nd and 4th : numbers of connected braces
% 3rd and 5th : status ofthat braces (compression, tension,...
% 6th for the orientation of leg
legLeft = zeros(nodesLevels,6);
legImpact1 = zeros(nodesLevels,6);
legImpact2 = zeros(nodesLevels,6);
legRight = zeros(nodesLevels,6);

for i = 1:nodesLevels
    legLeft(i,1) = startLeft + (i-1)*8;
    legImpact1(i,1) = startImpact + (i-1)*8;
    legImpact2(i,1) = startImpact + (i-1)*8;
    legRight(i,1) = startRight + (i-1)*8;
    
    legLeft(i,2) = startBraceLeft - 2 + (i-1)*4;
    legLeft(i,4) = startBraceLeft + (i-1)*4;
    
    legImpact1(i,2) = startBraceImpact1 - 2 + (i-1)*4;
    legImpact1(i,4) = startBraceImpact1 + (i-1)*4;
    
    legImpact2(i,2) = startBraceImpact2 - 2 + (i-1)*4;
    legImpact2(i,4) = startBraceImpact2 + (i-1)*4;
    
    legRight(i,2) = startBraceRight - 2 + (i-1)*4;
    legRight(i,4) = startBraceRight + (i-1)*4;
    
    if startLeg == 1
        legLeft(i,6) = 4;
        legRight(i,6) = 2;
    elseif startLeg == 4
        legLeft(i,6) = 3;
        legRight(i,6) = 1;
    else
        legLeft(i,6) = startLeg - 1;
        legRight(i,6) = startLeg + 1;
    end
    legImpact1(i,6) = startLeg;
    legImpact2(i,6) = startLeg;
end

legLeft(1,2) = horizLeft;
legLeft(nodesLevels,4) = horizLeft + 4;

legImpact1(1,2) = horizLeft;
legImpact1(nodesLevels,4) = horizLeft + 4;

legImpact2(1,2) = horizRight;
legImpact2(nodesLevels,4) = horizRight + 4;

legRight(1,2) = horizRight;
legRight(nodesLevels,4) = horizRight + 4;

punchPropLeft.leg = [legImpact1 ; legLeft];
punchPropRight.leg = [legImpact2 ; legRight];

dirLeg = properties.angle(1,:,2);
alphaAxisShip = acos(dirLeg*[1 ; 0 ; 0]);% * 180/pi;

if data.impactedZone == 0
    punchProp.lprime = 0;
    punchProp.lprimeprime = 0;
else
    H = properties.end(data.impactedElement,3) - properties.origin(data.impactedElement,3);
    d = properties.end(data.impactedElement,3) - data.shipHeight;
    
    punchProp.lprime = (H-d)/(cos(pi/2-alphaAxisShip));
    punchProp.lprimeprime = d/(cos(pi/2-alphaAxisShip));
end

%% Impacted nodes

if data.impactedZone == 0
    punchProp.impactedNodes = data.impactedNode;
    i = 1;
    while data.h(i) < properties.nodes(data.impactedNode,3)
        i = i+1;
    end
    punchProp.impactedNodesIndex = i-1;
elseif data.impactedZone == 1
    punchProp.impactedNodes = [properties.nodein(data.impactedElement) properties.nodeout(data.impactedElement)];
    i = 1;
    while data.h(i) < properties.origin(data.impactedElement,3)
        i = i+1;
    end
    punchProp.impactedNodesIndex = [i-1 i];
end

punchProp.nodesLevels = nodesLevels;

%% Length of legs

punchProp.LLegs = properties.L(1:nodesLevels);

%% Orientation

if startLeg == 1 || startLeg == 3
    punchPropLeft.angle = pi/2;
    punchPropRight.angle = 0;
else
    punchPropLeft.angle = 0;
    punchPropRight.angle = pi/2;
end