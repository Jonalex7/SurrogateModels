%% Nodes

deltaH = 0;
% data.h = [0 2.5-deltaH 22-deltaH 36-deltaH 46.5-deltaH 55.5-deltaH];
data.h = [0 2.5 22 36 46.5 55.5];
data.b = [25 6.4];
data.inclinedFeet = 1; % 0->No ; 1->Yes

properties.nodes = zeros((2*length(data.h)-1)*4,9);

if data.inclinedFeet == 0
    data.alpha = atan(((data.b(1)-data.b(2))/2)/(data.h(end)-data.h(2)));
else
    data.alpha = atan(((data.b(1)-data.b(2))/2)/(data.h(end)-data.h(1)));
end

data.posLeg    = zeros(length(data.h),1);
data.posBrace  = zeros(length(data.h)-2,1);
data.posBraceZ = zeros(length(data.h)-2,1);

for i = 1:length(data.posLeg)
    if data.inclinedFeet == 1
        data.posLeg(i) = data.b(1)/2 - data.h(i)*tan(data.alpha);
    else
        if i == 1
            data.posLeg(i) = data.b(1)/2;
        else
            data.posLeg(i) = data.b(1)/2 - (data.h(i)-data.h(2))*tan(data.alpha);
        end
    end
end
for i = 1:length(data.posBrace)
    deltah = data.h(i+2)-data.h(i+1);
    beta = atan(deltah / (2*data.posLeg(i+1)-deltah*tan(data.alpha)));
    data.posBraceZ(i) = data.h(i+1) + data.posLeg(i+1)*tan(beta);
    if data.inclinedFeet == 1
        data.posBrace(i) = data.b(1)/2 - data.posBraceZ(i)*tan(data.alpha);
    else
        data.posBrace(i) = data.b(1)/2 - (data.posBraceZ(i)-data.h(2))*tan(data.alpha);
    end
end

level = 1;
legLevel = 1;
braceLevel = 1;
for i = 1:size(properties.nodes,1)
    if level == 1 %fictive elements
        if mod(i,4) == 1
            properties.nodes(i,1) = -data.posLeg(legLevel);
            properties.nodes(i,2) = -data.posLeg(legLevel);
            properties.nodes(i,3) = data.h(1) - data.h(2);
        elseif mod(i,4) == 2
            properties.nodes(i,1) = data.posLeg(legLevel);
            properties.nodes(i,2) = -data.posLeg(legLevel);
            properties.nodes(i,3) = data.h(1) - data.h(2);
        elseif mod(i,4) == 3
            properties.nodes(i,1) = data.posLeg(legLevel);
            properties.nodes(i,2) = data.posLeg(legLevel);
            properties.nodes(i,3) = data.h(1) - data.h(2);
        elseif mod(i,4) == 0
            properties.nodes(i,1) = -data.posLeg(legLevel);
            properties.nodes(i,2) = data.posLeg(legLevel);
            properties.nodes(i,3) = data.h(1) - data.h(2);
            
            level = level + 1;
        end
        
    elseif level == 2 || mod(level,2) == 1 %leg
        if mod(i,4) == 1
            properties.nodes(i,1) = -data.posLeg(legLevel);
            properties.nodes(i,2) = -data.posLeg(legLevel);
            properties.nodes(i,3) = data.h(legLevel);
        elseif mod(i,4) == 2
            properties.nodes(i,1) = data.posLeg(legLevel);
            properties.nodes(i,2) = -data.posLeg(legLevel);
            properties.nodes(i,3) = data.h(legLevel);
        elseif mod(i,4) == 3
            properties.nodes(i,1) = data.posLeg(legLevel);
            properties.nodes(i,2) = data.posLeg(legLevel);
            properties.nodes(i,3) = data.h(legLevel);
        elseif mod(i,4) == 0
            properties.nodes(i,1) = -data.posLeg(legLevel);
            properties.nodes(i,2) = data.posLeg(legLevel);
            properties.nodes(i,3) = data.h(legLevel);
            
            level = level + 1;
            legLevel = legLevel + 1;
        end
    else %brace
        if mod(i,4) == 1
            properties.nodes(i,1) = 0;
            properties.nodes(i,2) = -data.posBrace(braceLevel);
            properties.nodes(i,3) = data.posBraceZ(braceLevel);
        elseif mod(i,4) == 2
            properties.nodes(i,1) = data.posBrace(braceLevel);
            properties.nodes(i,2) = 0;
            properties.nodes(i,3) = data.posBraceZ(braceLevel);
        elseif mod(i,4) == 3
            properties.nodes(i,1) = 0;
            properties.nodes(i,2) = data.posBrace(braceLevel);
            properties.nodes(i,3) = data.posBraceZ(braceLevel);
        elseif mod(i,4) == 0
            properties.nodes(i,1) = -data.posBrace(braceLevel);
            properties.nodes(i,2) = 0;
            properties.nodes(i,3) = data.posBraceZ(braceLevel);
            
            level = level + 1;
            braceLevel = braceLevel + 1;
        end
    end
end

for i = 1:size(properties.nodes,1)
    for j = 4:9
        properties.nodes(i,j) = 6*i+j-9;
    end
end


%% Elements - Data to be provided

% Each line correspond to an element
% Element 1 : nodes 1-2
% Element 2 : nodes 2-3
% Element 3 : nodes 3-4
% ...

% Geometric properties

data.nbLegs        = (length(data.h)-1)*4;
data.nbBottomHoriz = 4;
data.nbTopHoriz    = 4;
data.nbBraces      = (length(data.h)-2)*16;
data.nFictiveElements = 4;
data.nbElements = data.nbLegs + data.nbBottomHoriz + data.nbTopHoriz + data.nbBraces + data.nFictiveElements;

% Type : 1->leg ; 2->bottom horiz ; 3->top horiz ; 4->brace

properties.type = [ones(data.nbLegs,1) ; ones(data.nbBottomHoriz,1)*2 ; ones(data.nbTopHoriz,1)*3 ; ones(data.nbBraces,1)*4 ; ones(data.nFictiveElements,1)*5];
%                     1    2    3   4    5    6   7    8    9  10   11   12
properties.nodein = zeros(data.nbElements,1);
properties.nodeout = zeros(data.nbElements,1);

levelBrace = 1;
shift = 0;
brace = 1;

for i = 1:data.nbElements
    
    % Braces
    
    if i <= data.nbLegs
        if floor((i-1)/(data.nbLegs/4)) == 0 %1st leg
            if mod(i,data.nbLegs/4) == 1
                properties.nodein(i) = 5;
                properties.nodeout(i) = properties.nodein(i)+4;
            else
                properties.nodein(i) = properties.nodeout(i-1);
                properties.nodeout(i) = properties.nodeout(i-1)+8;
            end
        elseif floor((i-1)/(data.nbLegs/4)) == 1
            if mod(i,data.nbLegs/4) == 1
                properties.nodein(i) = 6;
                properties.nodeout(i) = properties.nodein(i)+4;
            else
                properties.nodein(i) = properties.nodeout(i-1);
                properties.nodeout(i) = properties.nodeout(i-1)+8;
            end
        elseif floor((i-1)/(data.nbLegs/4)) == 2
            if mod(i,data.nbLegs/4) == 1
                properties.nodein(i) = 7;
                properties.nodeout(i) = properties.nodein(i)+4;
            else
                properties.nodein(i) = properties.nodeout(i-1);
                properties.nodeout(i) = properties.nodeout(i-1)+8;
            end
        elseif floor((i-1)/(data.nbLegs/4)) == 3
            if mod(i,data.nbLegs/4) == 1
                properties.nodein(i) = 8;
                properties.nodeout(i) = properties.nodein(i)+4;
            else
                properties.nodein(i) = properties.nodeout(i-1);
                properties.nodeout(i) = properties.nodeout(i-1)+8;
            end
        end
    
    % Bottom Horiz
    
    elseif i <= data.nbLegs + data.nbBottomHoriz
        if mod(i,4) == 1
            properties.nodein(i) = 9;
            properties.nodeout(i) = 10;
        elseif mod(i,4) == 2
            properties.nodein(i) = 10;
            properties.nodeout(i) = 11;
        elseif mod(i,4) == 3
            properties.nodein(i) = 11;
            properties.nodeout(i) = 12;
        elseif mod(i,4) == 0
            properties.nodein(i) = 12;
            properties.nodeout(i) = 9;
        end
    
    % Top Horiz
    
    elseif i <= data.nbLegs + data.nbBottomHoriz + data.nbTopHoriz
        if mod(i,4) == 1
            properties.nodein(i) = size(properties.nodes,1)-3;
            properties.nodeout(i) = size(properties.nodes,1)-2;
        elseif mod(i,4) == 2
            properties.nodein(i) = size(properties.nodes,1)-2;
            properties.nodeout(i) = size(properties.nodes,1)-1;
        elseif mod(i,4) == 3
            properties.nodein(i) = size(properties.nodes,1)-1;
            properties.nodeout(i) = size(properties.nodes,1);
        elseif mod(i,4) == 0
            properties.nodein(i) = size(properties.nodes,1);
            properties.nodeout(i) = size(properties.nodes,1)-3;
        end
        
    % Braces
    
    elseif i <= data.nbLegs + data.nbBottomHoriz + data.nbTopHoriz + data.nbBraces
        if brace == 1
            properties.nodein(i) = 13+(levelBrace-1)*8+shift;
            properties.nodeout(i) = 9+(levelBrace-1)*8+shift;
            brace = 2;
        elseif brace == 2
            properties.nodein(i) = 13+(levelBrace-1)*8+shift;
            if shift == 3
                properties.nodeout(i) = 10+(levelBrace-1)*8-1;
            else
                properties.nodeout(i) = 10+(levelBrace-1)*8+shift;
            end
            brace = 3;
        elseif brace == 3
            properties.nodein(i) = 13+(levelBrace-1)*8+shift;
            properties.nodeout(i) = 17+(levelBrace-1)*8+shift;
            brace = 4;
        elseif brace == 4
            properties.nodein(i) = 13+(levelBrace-1)*8+shift;
            if shift == 3
                properties.nodeout(i) = 18+(levelBrace-1)*8-1;
            else
                properties.nodeout(i) = 18+(levelBrace-1)*8+shift;
            end
            brace = 1;
            levelBrace = levelBrace + 1;
        end
        if levelBrace == length(data.h)-1
            levelBrace = 1;
            shift = shift + 1;
        end
        
    % Fictive elements
    
    else
        if mod(i,4) == 1
            properties.nodein(i) = 1;
            properties.nodeout(i) = 5;
        elseif mod(i,4) == 2
            properties.nodein(i) = 2;
            properties.nodeout(i) = 6;
        elseif mod(i,4) == 3
            properties.nodein(i) = 3;
            properties.nodeout(i) = 7;
        elseif mod(i,4) == 0
            properties.nodein(i) = 4;
            properties.nodeout(i) = 8;
        end
    
    end
    
end

% Cross section properties

data.De = [ 1.3 ;  1.3 ;  1.3 ; 0.65 ;  1.3];
data.t  = [0.05 ; 0.05 ; 0.05 ; 0.05 ; 0.05];

% Material properties

data.fy = [ones(4,1)*317*10^6 ; 317*10^9];
data.E  = [ones(4,1)*210000*10^6 ; 210000*10^9];

% Boundary conditions

data.rho = [1 1 ;
            1 1 ; 
            1 1 ;
            1 1 ;
            1 1];
% % properties.rho = ones(length(properties.nodein),2); % Fixity factors
% % properties.rho(9:16,:) = 0.6;

%% Boundary Conditions

dofs.blokedDofs = 1:24; % Numbers of bloked Dofs

%% Ship

data.impactedZone = 1; % 0->node ; 1->element
if data.impactedZone == 0
    data.impactedNode = 25;%17
    if ismember(properties.nodes(data.impactedNode,3),data.h(2:end-1)) % Node on a leg
%         properties.impactedNodeModel = 1;
        properties.defLoc = 1;
    else % Crossing of braces
        properties.defLoc = 0;
    end
    data.shipHeight = properties.nodes(data.impactedNode,3);
else
    data.shipHeight = 43;%43;%44.5-deltaH;%43
    data.impactedElement = 4;%86 - 4 - 5
    properties.impactedElemModel = [1 2];
    properties.defLoc = 1;
end

data.VInit = 1; % 0->cstt speed ; 1->initial speed
data.dispmax = 1.0; %%%%% 0.3 Originally in Timothee's Code
%data.shipTrajectory = 0;%0; % angle (in degrees) with regard to the X axis
%data.pointTrajectory = [0 0];%[0 -7];%[0 -2.5];
% data.shipSpeed = 0.3/3; %m/s
%data.shipSpeed = 5; %m/s
data.shipWeight = 6*10^6; %kg

% Ship

ship.bulbous = 0; %0->no ; 1-> yes

ship.p = 5; % m %5    related with psib
ship.q = 12; % m     related with phib
ship.phib = 75;% 75 degrees     % related with q
ship.psib = 78.7; % 78.7 degrees    % related with p
ship.hb = 8; % m

ship.hBulb = 0.1;%3
ship.bBulb = 0.1;%2
ship.pBulb = 0.1;%4
ship.shiftCentreBulb = 5;

% ship.p = 11.5; % m
% ship.q = 32; % m
% ship.phib = 59.7;%80; % degrees
% ship.psib = 82.9; % degrees
% ship.hb = 25.9; % m = hship
% 
% ship.hBulb = 9.1;%3
% ship.bBulb = 1.9;%2 = pbulb
% ship.pBulb = 6;%4 = qbulb
% ship.shiftCentreBulb = 22.5;