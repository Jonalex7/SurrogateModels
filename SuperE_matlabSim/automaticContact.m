function[contactElement,contactNode] = automaticContact(data,properties,ship)

nbCheckPoints = 10; %10

bowSlices = 10; %5
bulbSlices = 50;

% % % % disp('Computation of potential impacted elements or nodes')

%% Non-bulbous ship


if ship.bulbous == 0
    
    contactElement.listBow = zeros(size(properties.origin,1),1);
    contactElement.listBulb = 0;
    
    zmax = data.shipHeight;
    zmin = data.shipHeight - ship.hb;
    
    slopeTraj = tan(data.shipTrajectory*pi/180);
    shiftTraj = data.pointTrajectory(2) - slopeTraj*data.pointTrajectory(1);
    pb = ship.p - (ship.hb/tan(ship.psib*pi/180));  %Lower semi-beam of the ship
%     shiftTrajUpBow = shiftTraj + ship.p/(cos(data.shipTrajectory*pi/180));
%     shiftTrajDownBow = shiftTraj - ship.p/(cos(data.shipTrajectory*pi/180));
    
    numberImpactedElementsBow = 1;
    
    for i = 1:size(properties.origin,1)
        
        impactElemBow = 0;
        
        coordOrigin = properties.origin(i,:);
        coordEnd = properties.end(i,:);
        
        deltax = (coordEnd(1) - coordOrigin(1))/nbCheckPoints;
        deltay = (coordEnd(2) - coordOrigin(2))/nbCheckPoints;
        deltaz = (coordEnd(3) - coordOrigin(3))/nbCheckPoints;
        
        for checkPoint = 1:nbCheckPoints+1 % if one extremity in the box
            
            xPoint = coordOrigin(1) + (checkPoint-1)*deltax;
            yPoint = coordOrigin(2) + (checkPoint-1)*deltay;
            zPoint = coordOrigin(3) + (checkPoint-1)*deltaz;
            
            
            % Check according to Z axis
            
            if zPoint > zmin && zPoint < zmax
                
                % Check according to in-plane displacement

%                 ymin = slopeTraj*xPoint + shiftTrajDownBow;
%                 ymax = slopeTraj*xPoint + shiftTrajUpBow;
                zpby=zPoint-zmin;
                pby = (((ship.p-pb)/ship.hb)*zpby)+pb;  %beam at zPoint
                shiftTrajUpBow = shiftTraj + pby/(cos(data.shipTrajectory*pi/180));
                shiftTrajDownBow = shiftTraj - pby/(cos(data.shipTrajectory*pi/180));

                ymin = (slopeTraj*xPoint + shiftTrajDownBow);
                ymax = (slopeTraj*xPoint + shiftTrajUpBow);
                
                if yPoint > ymin && yPoint < ymax
                    impactElemBow = 1;
                end
                
            end
            
        end
        
        if impactElemBow == 1
            contactElement.listBow(numberImpactedElementsBow) = i;
            numberImpactedElementsBow = numberImpactedElementsBow + 1;
        end
        
    end
    
    
    
    contactNode.listBow = zeros(size(properties.nodes,1),1);
    contactNode.listBulb = 0;
    
    numberImpactedNodesBow = 1;
    
    for i = 1:size(properties.nodes,1)
        
        impactNodeBow = 0;
        
        xPoint = properties.nodes(i,1);
        yPoint = properties.nodes(i,2);
        zPoint = properties.nodes(i,3);
        
        % Check according to Z axis
        
        if zPoint > zmin && zPoint < zmax
            
            % Check according to in-plane displacement
            zpby=zPoint-zmin;
            pby = (((ship.p-pb)/ship.hb)*zpby)+pb;  %beam at zPoint
            shiftTrajUpBow = shiftTraj + pby/(cos(data.shipTrajectory*pi/180));
            shiftTrajDownBow = shiftTraj - pby/(cos(data.shipTrajectory*pi/180));            
            ymin = slopeTraj*xPoint + shiftTrajDownBow;
            ymax = slopeTraj*xPoint + shiftTrajUpBow;
            
            if yPoint > ymin && yPoint < ymax
                impactNodeBow = 1;
            end
            
        end
        
        if impactNodeBow == 1
            contactNode.listBow(numberImpactedNodesBow) = i;
            numberImpactedNodesBow = numberImpactedNodesBow + 1;
        end
        
    end
    
    
    if numberImpactedElementsBow == 1
        contactElement.listBow = contactElement.listBow(1:numberImpactedElementsBow);
    else
        contactElement.listBow = contactElement.listBow(1:numberImpactedElementsBow-1);
    end
    
    if numberImpactedNodesBow == 1
        contactNode.listBow = contactNode.listBow(1:numberImpactedNodesBow);
    else
        contactNode.listBow = contactNode.listBow(1:numberImpactedNodesBow-1);
    end
    
    
%% Bulbous ship
    
else
    
    contactElement.listBow = zeros(size(properties.origin,1),1);
    contactElement.listBulb = zeros(size(properties.origin,1),1);
    
    zmax = data.shipHeight;
    zmin = data.shipHeight - ship.hb;
    zTopBulb = zmin + ship.hBulb;
    
    slopeTraj = tan(data.shipTrajectory*pi/180);
    shiftTraj = data.pointTrajectory(2) - slopeTraj*data.pointTrajectory(1);
    shiftTrajUpBow = shiftTraj + ship.p/(cos(data.shipTrajectory*pi/180));
    shiftTrajDownBow = shiftTraj - ship.p/(cos(data.shipTrajectory*pi/180));
    shiftTrajUpBulb = shiftTraj + ship.bBulb/2/(cos(data.shipTrajectory*pi/180));
    shiftTrajDownBulb = shiftTraj - ship.bBulb/2/(cos(data.shipTrajectory*pi/180));
    
    numberImpactedElementsBow = 1;
    numberImpactedElementsBulb = 1;
    
    for i = 1:size(properties.origin,1)
        
        impactElemBow = 0;
        impactElemBulb = 0;
        
        coordOrigin = properties.origin(i,:);
        coordEnd = properties.end(i,:);
        
        deltax = (coordEnd(1) - coordOrigin(1))/nbCheckPoints;
        deltay = (coordEnd(2) - coordOrigin(2))/nbCheckPoints;
        deltaz = (coordEnd(3) - coordOrigin(3))/nbCheckPoints;
        
        for checkPoint = 1:nbCheckPoints+1 % if one extremity in the box
            
            xPoint = coordOrigin(1) + (checkPoint-1)*deltax;
            yPoint = coordOrigin(2) + (checkPoint-1)*deltay;
            zPoint = coordOrigin(3) + (checkPoint-1)*deltaz;
            
            % Check according to Z axis - Bow
            
            if zPoint > zTopBulb && zPoint < zmax
                
                % Check according to in-plane displacement
                
                ymin = slopeTraj*xPoint + shiftTrajDownBow;
                ymax = slopeTraj*xPoint + shiftTrajUpBow;
                
                if yPoint > ymin && yPoint < ymax
                    impactElemBow = 1;
                end
                
            end
            
            % Check according to Z axis - Bulb
            
            if zPoint > zmin && zPoint < zTopBulb
                
                % Check according to in-plane displacement
                
                ymin = slopeTraj*xPoint + shiftTrajDownBulb;
                ymax = slopeTraj*xPoint + shiftTrajUpBulb;
                
                if yPoint > ymin && yPoint < ymax
                    impactElemBulb = 1;
                end
                
            end
            
        end
        
        if impactElemBow == 1
            contactElement.listBow(numberImpactedElementsBow) = i;
            numberImpactedElementsBow = numberImpactedElementsBow + 1;
        end
        if impactElemBulb == 1
            contactElement.listBulb(numberImpactedElementsBulb) = i;
            numberImpactedElementsBulb = numberImpactedElementsBulb + 1;
        end
        
    end
    
    
    contactNode.listBow = zeros(size(properties.origin,1),1);
    contactNode.listBulb = zeros(size(properties.origin,1),1);
    
    numberImpactedNodesBow = 1;
    numberImpactedNodesBulb = 1;
    
    for i = 1:size(properties.nodes,1)
        
        impactNodeBow = 0;
        impactNodeBulb = 0;
        
        xPoint = properties.nodes(i,1);
        yPoint = properties.nodes(i,2);
        zPoint = properties.nodes(i,3);
        
        % Check according to Z axis - Bow
        
        if zPoint > zTopBulb && zPoint < zmax
            
            % Check according to in-plane displacement
            
            ymin = slopeTraj*xPoint + shiftTrajDownBow;
            ymax = slopeTraj*xPoint + shiftTrajUpBow;
            
            if yPoint > ymin && yPoint < ymax
                impactNodeBow = 1;
            end
            
        end
        
        % Check according to Z axis - Bulb
        
        if zPoint > zmin && zPoint < zTopBulb
            
            % Check according to in-plane displacement
            
            ymin = slopeTraj*xPoint + shiftTrajDownBulb;
            ymax = slopeTraj*xPoint + shiftTrajUpBulb;
            
            if yPoint > ymin && yPoint < ymax
                impactNodeBulb = 1;
            end
            
        end
        
        if impactNodeBow == 1
            contactNode.listBow(numberImpactedNodesBow) = i;
            numberImpactedNodesBow = numberImpactedNodesBow + 1;
        end
        if impactNodeBulb == 1
            contactNode.listBulb(numberImpactedNodesBulb) = i;
            numberImpactedNodesBulb = numberImpactedNodesBulb + 1;
        end
        
    end
    
    
    if numberImpactedElementsBow == 1
        contactElement.listBow = contactElement.listBow(1:numberImpactedElementsBow);
    else
        contactElement.listBow = contactElement.listBow(1:numberImpactedElementsBow-1);
    end
    if numberImpactedElementsBulb == 1
        contactElement.listBulb = contactElement.listBulb(1:numberImpactedElementsBulb);
    else
        contactElement.listBulb = contactElement.listBulb(1:numberImpactedElementsBulb-1);
    end
    
    if numberImpactedNodesBow == 1
        contactNode.listBow = contactNode.listBow(1:numberImpactedNodesBow);
    else
        contactNode.listBow = contactNode.listBow(1:numberImpactedNodesBow-1);
    end
    if numberImpactedNodesBulb == 1
        contactNode.listBulb = contactNode.listBulb(1:numberImpactedNodesBulb);
    else
        contactNode.listBulb = contactNode.listBulb(1:numberImpactedNodesBulb-1);
    end
    
    
    
end

listBowElem = contactElement.listBow;
listBulbElem = contactElement.listBulb;
listBowNode = contactNode.listBow;
listBulbNode = contactNode.listBulb;

%% Initial position of the ship

% % % % disp('Computation of distance')

dataCyl.a2 = ship.q;
dataCyl.b2 = ship.p;
dataCyl.alphaDeg2 = data.shipTrajectory;

xi = atan((data.h(end)-data.h(1)) / sqrt(2*((data.b(2)-data.b(1))/2)^2));
Rleg = data.De(1)/2;
a1 = Rleg/sin(xi);
heightReport = data.shipHeight / data.h(end);
posLegHeight = -(data.b(1)*(1-heightReport) + data.b(2)*heightReport);

if abs(dataCyl.alphaDeg2) < 45
    xS = posLegHeight - a1 - dataCyl.a2;
    yS = slopeTraj*xS + shiftTraj;    
else
    yS = posLegHeight - a1 - dataCyl.a2;
    xS = (yS-shiftTraj)/slopeTraj;
end

shiftx2Init = xS;
shifty2Init = yS;

distMax = sqrt(2)*data.b(1) + ship.q;


%% Bow-jacket distances

contactElement.distBow = zeros(length(contactElement.listBow),5);

if ship.bulbous == 0
    deltaSlice = ship.hb / bowSlices;
else
    deltaSlice = (ship.hb - ship.hBulb) / bowSlices;
end

for elem = 1:length(contactElement.listBow)
    
    nbElem = contactElement.listBow(elem);
    
    zmin = min([properties.origin(nbElem,3) properties.end(nbElem,3)]);
    zmax = max([properties.origin(nbElem,3) properties.end(nbElem,3)]);
    
    deltaX = properties.end(nbElem,1) - properties.origin(nbElem,1);
    deltaY = properties.end(nbElem,2) - properties.origin(nbElem,2);
    deltaZ = properties.end(nbElem,3) - properties.origin(nbElem,3);
    
    xi = atan(abs(deltaZ) / sqrt(deltaX^2 + deltaY^2));
    dataCyl.alphaDeg1 = atan(deltaY/deltaX)*180/pi;
    Relem = data.De(properties.type(nbElem))/2;
    dataCyl.a1 = Relem/(sin(xi));
    dataCyl.b1 = Relem;
    
    % Compute distance
    
    dist = 0;
    contactEllipses = 0;
    slicesIn = zeros(bowSlices,1);
    
    while contactEllipses == 0 && dist <= distMax
        dist = dist + Relem;
        
        contactVecBow = zeros(bowSlices+1,1);
        
        for nbSlice = 1:bowSlices+1
            zSlice = data.shipHeight - (nbSlice-1)*deltaSlice;
            
            if zSlice > zmin && zSlice < zmax
                
                heightReport = (zSlice - properties.origin(nbElem,3)) / deltaZ;
                dataCyl.shiftx1 = properties.origin(nbElem,1)*(1-heightReport) + properties.end(nbElem,1)*heightReport;
                dataCyl.shifty1 = properties.origin(nbElem,2)*(1-heightReport) + properties.end(nbElem,2)*heightReport;
                
                dataCyl.shiftx2 = shiftx2Init + dist*cos(dataCyl.alphaDeg2*pi/180);
                dataCyl.shifty2 = shifty2Init + dist*sin(dataCyl.alphaDeg2*pi/180);
                
                dataCyl.a2 = ship.q - (nbSlice-1)*deltaSlice*cos(ship.phib*pi/180);
                dataCyl.b2 = ship.p - (nbSlice-1)*deltaSlice*cos(ship.psib*pi/180);
                
                [contactEllipses,solveX,solveY] = Ellipses(dataCyl);
                
                contactVecBow(nbSlice) = contactEllipses;
                
            else
                
                contactVecBow(nbSlice) = -1;
                
            end
            
        end
        
        [contactEllipses,labelZ] = max(contactVecBow);
        
        if contactEllipses == 1
            labelZOK = labelZ;
        end
        
        countSlicesIn = 1;
        for i = 1:bowSlices+1
            if contactVecBow(i) == 1
                slicesIn(countSlicesIn) = i;
                countSlicesIn = countSlicesIn + 1;
            end
        end
        if countSlicesIn == 1
            slicesIn1 = 0;
        else
            slicesIn1 = slicesIn(1:countSlicesIn-1);
        end
        
    end
    
    while contactEllipses == 1 && dist <= distMax
        dist = dist - 0.1;
        
        for nbSlice = slicesIn1(1):slicesIn1(end)
            zSlice = data.shipHeight - (nbSlice-1)*deltaSlice;
            
            heightReport = (zSlice - properties.origin(nbElem,3)) / deltaZ;
            dataCyl.shiftx1 = properties.origin(nbElem,1)*(1-heightReport) + properties.end(nbElem,1)*heightReport;
            dataCyl.shifty1 = properties.origin(nbElem,2)*(1-heightReport) + properties.end(nbElem,2)*heightReport;
            
            dataCyl.shiftx2 = shiftx2Init + dist*cos(dataCyl.alphaDeg2*pi/180);
            dataCyl.shifty2 = shifty2Init + dist*sin(dataCyl.alphaDeg2*pi/180);
            
            dataCyl.a2 = ship.q - (nbSlice-1)*deltaSlice*cos(ship.phib*pi/180);
            dataCyl.b2 = ship.p - (nbSlice-1)*deltaSlice*cos(ship.psib*pi/180);
            
            [contactEllipses,solveX,solveY] = Ellipses(dataCyl);
            
            contactVecBow(nbSlice) = contactEllipses;
            
        end
        
        slicesIn2 = slicesIn1;
        
        [contactEllipses,labelZ] = max(contactVecBow);
        
        if contactEllipses == 1
            labelZOK = labelZ;
            
            countSlicesIn = 1;
            for i = slicesIn1(1):slicesIn1(end)
                if contactVecBow(i) == 1
                    slicesIn(countSlicesIn) = i;
                    countSlicesIn = countSlicesIn + 1;
                end
            end
            if countSlicesIn == 1
                slicesIn2 = 0;
            else
                slicesIn2 = slicesIn(1:countSlicesIn-1);
            end
        end
        
    end
    
    while contactEllipses == 0 && dist <= distMax
        dist = dist + 0.01;
        
        for nbSlice = slicesIn2(1):slicesIn2(end)
            zSlice = data.shipHeight - (nbSlice-1)*deltaSlice;
            
            heightReport = (zSlice - properties.origin(nbElem,3)) / deltaZ;
            dataCyl.shiftx1 = properties.origin(nbElem,1)*(1-heightReport) + properties.end(nbElem,1)*heightReport;
            dataCyl.shifty1 = properties.origin(nbElem,2)*(1-heightReport) + properties.end(nbElem,2)*heightReport;
            
            dataCyl.shiftx2 = shiftx2Init + dist*cos(dataCyl.alphaDeg2*pi/180);
            dataCyl.shifty2 = shifty2Init + dist*sin(dataCyl.alphaDeg2*pi/180);
            
            dataCyl.a2 = ship.q - (nbSlice-1)*deltaSlice*cos(ship.phib*pi/180);
            dataCyl.b2 = ship.p - (nbSlice-1)*deltaSlice*cos(ship.psib*pi/180);
            
            [contactEllipses,solveX,solveY] = Ellipses(dataCyl);
            
            contactVecBow(nbSlice) = contactEllipses;
            
        end
        
        [contactEllipses,labelZ] = max(contactVecBow);
        
        if contactEllipses == 1
            labelZOK = labelZ;
        end
        
        slicesIn3 = slicesIn2;
        
        countSlicesIn = 1;
        for i = slicesIn2(1):slicesIn2(end)
            if contactVecBow(i) == 1
                slicesIn(countSlicesIn) = i;
                countSlicesIn = countSlicesIn + 1;
            end
        end
        if countSlicesIn == 1
            slicesIn3 = 0;
        else
            slicesIn3 = slicesIn(1:countSlicesIn-1);
        end
        
    end
    
    solveResult = zeros(bowSlices+1,2);
    
    while contactEllipses == 1 && dist <= distMax
        dist = dist - 0.001;
        
        for nbSlice = slicesIn3(1):slicesIn3(end)
            zSlice = data.shipHeight - (nbSlice-1)*deltaSlice;
            
            heightReport = (zSlice - properties.origin(nbElem,3)) / deltaZ;
            dataCyl.shiftx1 = properties.origin(nbElem,1)*(1-heightReport) + properties.end(nbElem,1)*heightReport;
            dataCyl.shifty1 = properties.origin(nbElem,2)*(1-heightReport) + properties.end(nbElem,2)*heightReport;
            
            dataCyl.shiftx2 = shiftx2Init + dist*cos(dataCyl.alphaDeg2*pi/180);
            dataCyl.shifty2 = shifty2Init + dist*sin(dataCyl.alphaDeg2*pi/180);
            
            dataCyl.a2 = ship.q - (nbSlice-1)*deltaSlice*cos(ship.phib*pi/180);
            dataCyl.b2 = ship.p - (nbSlice-1)*deltaSlice*cos(ship.psib*pi/180);
            
            [contactEllipses,solveX,solveY] = Ellipses(dataCyl);
            
            contactVecBow(nbSlice) = contactEllipses;
            if contactEllipses == 1
                solveResult(nbSlice,1) = solveX;
                solveResult(nbSlice,2) = solveY;
            end
            
        end
        
        [contactEllipses,labelZ] = max(contactVecBow);
        
        if contactEllipses == 1
            labelZOK = labelZ;
        end
        
    end
    
    contactElement.distBow(elem,1) = dist;
    contactElement.distBow(elem,2) = data.shipHeight - (labelZOK-1)*deltaSlice;
    
    solveX = solveResult(labelZOK,1);
    solveY = solveResult(labelZOK,2);
    
    contactElement.distBow(elem,3) = (solveX - dataCyl.shiftx2)*cos(pi+data.shipTrajectory*pi/180) + (solveY - dataCyl.shifty2)*sin(pi+data.shipTrajectory*pi/180);
    contactElement.distBow(elem,4) = -(solveX - dataCyl.shiftx2)*sin(pi+data.shipTrajectory*pi/180) + (solveY - dataCyl.shifty2)*cos(pi+data.shipTrajectory*pi/180);
    contactElement.distBow(elem,5) =  - (labelZOK-1)*deltaSlice;
    
    if dist > distMax
        contactElement.distBow(elem,1) = 1000;
    end
    
end

%% Bow-node distance

contactNode.distBow = zeros(length(contactNode.listBow),5);

if contactNode.listBow(1) ~= 0
    
    for node = 1:length(contactNode.listBow)
        
        nbNode = contactNode.listBow(node);
        
        xPoint = properties.nodes(nbNode,1);
        yPoint = properties.nodes(nbNode,2);
        zPoint = properties.nodes(nbNode,3);
        
        shiftz = data.shipHeight - zPoint;
        
        dataCyl.a2 = ship.q - shiftz*sin(ship.phib*pi/180);
        dataCyl.b2 = ship.p - shiftz*sin(ship.psib*pi/180);
        dataCyl.alphaDeg2 = data.shipTrajectory;
        
        if xPoint ~= 0 && yPoint ~= 0 %Leg
            
            xi = atan((data.h(end)-data.h(1)) / sqrt(2*((data.b(2)-data.b(1))/2)^2));
            Rleg = data.De(1)/2;
            dataCyl.a1 = Rleg/(sin(xi));
            dataCyl.b1 = Rleg;
            dataCyl.alphaDeg1 = 45 * sign(xPoint) * sign(yPoint);
            
        else %Crossing of braces
            
            Rbrace = data.De(4)/2;
            level = 1;
            while data.h(level) < zPoint
                level = level + 1;
            end
            deltaZ = data.h(level) - data.h(level-1);
            bottomReport = data.h(level-1)/data.h(end);
            bottomWidth = (data.b(1)*(1-bottomReport) + data.b(end)*bottomReport)/2;
            topReport = data.h(level)/data.h(end);
            topWidth = (data.b(1)*(1-topReport) + data.b(end)*topReport)/2;
            deltaX = bottomWidth + topWidth;
            
            xi = atan(deltaZ/deltaX);
            
            dataCyl.a1 = 2*Rbrace*cos(xi)/(cos(pi/2-2*xi));
            dataCyl.b1 = Rbrace;
            
            if xPoint == 0
                dataCyl.alphaDeg1 = 0;
            else
                dataCyl.alphaDeg1 = 90;
            end
            
        end
        
        dataCyl.shiftx1 = xPoint;
        dataCyl.shifty1 = yPoint;
        
        
        % Compute distance
        
        dist = 0;
        contactEllipses = 0;
        
        while contactEllipses == 0 && dist <= distMax
            dist = dist + Rleg;
            
            dataCyl.shiftx2 = shiftx2Init + dist*cos(dataCyl.alphaDeg2*pi/180);
            dataCyl.shifty2 = shifty2Init + dist*sin(dataCyl.alphaDeg2*pi/180);
            [contactEllipses,solveX,solveY] = Ellipses(dataCyl);
            
        end
        
        while contactEllipses == 1 && dist <= distMax
            dist = dist - 0.1;
            
            dataCyl.shiftx2 = shiftx2Init + dist*cos(dataCyl.alphaDeg2*pi/180);
            dataCyl.shifty2 = shifty2Init + dist*sin(dataCyl.alphaDeg2*pi/180);
            [contactEllipses,solveX,solveY] = Ellipses(dataCyl);
            
        end
        
        while contactEllipses == 0 && dist <= distMax
            dist = dist + 0.01;
            
            dataCyl.shiftx2 = shiftx2Init + dist*cos(dataCyl.alphaDeg2*pi/180);
            dataCyl.shifty2 = shifty2Init + dist*sin(dataCyl.alphaDeg2*pi/180);
            [contactEllipses,solveX,solveY] = Ellipses(dataCyl);
            
        end
        
        while contactEllipses == 1 && dist <= distMax
            dist = dist - 0.001;
            
            dataCyl.shiftx2 = shiftx2Init + dist*cos(dataCyl.alphaDeg2*pi/180);
            dataCyl.shifty2 = shifty2Init + dist*sin(dataCyl.alphaDeg2*pi/180);
            [contactEllipses,solveX,solveY] = Ellipses(dataCyl);
            
        end
        
        contactNode.distBow(node,1) = dist;
        contactNode.distBow(node,2) = zPoint;
        
        contactNode.distBow(node,3) = (solveX - dataCyl.shiftx2)*cos(pi+data.shipTrajectory*pi/180) + (solveY - dataCyl.shifty2)*sin(pi+data.shipTrajectory*pi/180);
        contactNode.distBow(node,4) = -(solveX - dataCyl.shiftx2)*sin(pi+data.shipTrajectory*pi/180) + (solveY - dataCyl.shifty2)*cos(pi+data.shipTrajectory*pi/180);
        contactNode.distBow(node,5) =  zPoint - data.shipHeight;
        
        if dist > distMax
            contactNode.distBow(node,1) = 1000;
        end
        
    end
    
end

%% Bulb-jacket distances

contactElement.distBulb = zeros(length(contactElement.listBulb),5);

deltaSlice = ship.hBulb / bulbSlices;

hEq = ship.hBulb/2;
bEq = ship.bBulb/2;
pEq = ship.pBulb;

if contactElement.listBulb(1) ~= 0
    
    for elem = 1:length(contactElement.listBulb)
        
        nbElem = contactElement.listBulb(elem);
        
        zmin = min([properties.origin(nbElem,3) properties.end(nbElem,3)]);
        zmax = max([properties.origin(nbElem,3) properties.end(nbElem,3)]);
        
        deltaX = properties.end(nbElem,1) - properties.origin(nbElem,1);
        deltaY = properties.end(nbElem,2) - properties.origin(nbElem,2);
        deltaZ = properties.end(nbElem,3) - properties.origin(nbElem,3);
        
        xi = atan(abs(deltaZ) / sqrt(deltaX^2 + deltaY^2));
        dataCyl.alphaDeg1 = atan(deltaY/deltaX)*180/pi;
        Relem = data.De(properties.type(nbElem))/2;
        dataCyl.a1 = Relem/(sin(xi));
        dataCyl.b1 = Relem;
        
        % Compute distance
        
        dist = 0;
        contactEllipses = 0;
        slicesIn = zeros(bulbSlices,1);
        
        while contactEllipses == 0 && dist <= distMax
            dist = dist + Relem;
            
            contactVecBulb = zeros(bulbSlices-1,1);
            
            for nbSlice = 1:bulbSlices-1
                zSlice = data.shipHeight - ship.hb + ship.hBulb - nbSlice*deltaSlice;
                
                if zSlice > zmin && zSlice < zmax
                    
                    heightReport = (zSlice - properties.origin(nbElem,3)) / deltaZ;
                    dataCyl.shiftx1 = properties.origin(nbElem,1)*(1-heightReport) + properties.end(nbElem,1)*heightReport;
                    dataCyl.shifty1 = properties.origin(nbElem,2)*(1-heightReport) + properties.end(nbElem,2)*heightReport;
                    
                    dataCyl.shiftx2 = shiftx2Init + (ship.shiftCentreBulb+dist)*cos(dataCyl.alphaDeg2*pi/180);
                    dataCyl.shifty2 = shifty2Init + (ship.shiftCentreBulb+dist)*sin(dataCyl.alphaDeg2*pi/180);
                    
                    zEq = zSlice - (data.shipHeight - ship.hb + hEq);
                    dataCyl.a2 = pEq * sqrt(1-(zEq/hEq)^2);
                    dataCyl.b2 = bEq * sqrt(1-(zEq/hEq)^2);
                    
                    [contactEllipses,solveX,solveY] = Ellipses(dataCyl);
                    
                    contactVecBulb(nbSlice) = contactEllipses;
                    
                else
                    
                    contactVecBulb(nbSlice) = -1;
                    
                end
                
            end
            
%             efg = contactVecBulb
            
            [contactEllipses,labelZ] = max(contactVecBulb);
            
            if contactEllipses == 1
                labelZOK = labelZ;
            end
            
            countSlicesIn = 1;
            for i = 1:bulbSlices-1
                if contactVecBulb(i) == 1
                    slicesIn(countSlicesIn) = i;
                    countSlicesIn = countSlicesIn + 1;
                end
            end
            if countSlicesIn == 1
                slicesIn1 = 0;
            else
                slicesIn1 = slicesIn(1:countSlicesIn-1);
            end
            
            %         contactVecBulb
            %         labelZ
            
        end
        
        while contactEllipses == 1 && dist <= distMax
            dist = dist - 0.1;
            
            contactVecBulb = zeros(bulbSlices-1,1);
            
            for nbSlice = slicesIn1(1):slicesIn1(end)
                zSlice = data.shipHeight - ship.hb + ship.hBulb - nbSlice*deltaSlice;
                
                if zSlice > zmin && zSlice < zmax
                    
                    heightReport = (zSlice - properties.origin(nbElem,3)) / deltaZ;
                    dataCyl.shiftx1 = properties.origin(nbElem,1)*(1-heightReport) + properties.end(nbElem,1)*heightReport;
                    dataCyl.shifty1 = properties.origin(nbElem,2)*(1-heightReport) + properties.end(nbElem,2)*heightReport;
                    
                    dataCyl.shiftx2 = shiftx2Init + (ship.shiftCentreBulb+dist)*cos(dataCyl.alphaDeg2*pi/180);
                    dataCyl.shifty2 = shifty2Init + (ship.shiftCentreBulb+dist)*sin(dataCyl.alphaDeg2*pi/180);
                    
                    zEq = zSlice - (data.shipHeight - ship.hb + hEq);
                    dataCyl.a2 = pEq * sqrt(1-(zEq/hEq)^2);
                    dataCyl.b2 = bEq * sqrt(1-(zEq/hEq)^2);
                    
                    [contactEllipses,solveX,solveY] = Ellipses(dataCyl);
                    
                    contactVecBulb(nbSlice) = contactEllipses;
                    
                else
                    
                    contactVecBulb(nbSlice) = -1;
                    
                end
                
            end
            
            [contactEllipses,labelZ] = max(contactVecBulb);
            
            if contactEllipses == 1
                labelZOK = labelZ;
                
                countSlicesIn = 1;
                for i = slicesIn1(1):slicesIn1(end)
                    if contactVecBulb(i) == 1
                        slicesIn(countSlicesIn) = i;
                        countSlicesIn = countSlicesIn + 1;
                    end
                end
                if countSlicesIn == 1
                    slicesIn2 = 0;
                else
                    slicesIn2 = slicesIn(1:countSlicesIn-1);
                end
            end
            
            %         contactVecBulb
            %         labelZ
            
        end
        
        while contactEllipses == 0 && dist <= distMax
            dist = dist + 0.01;
            
            contactVecBulb = zeros(bulbSlices-1,1);
            
            for nbSlice = slicesIn2(1):slicesIn2(end)
                zSlice = data.shipHeight - ship.hb + ship.hBulb - nbSlice*deltaSlice;
                
                if zSlice > zmin && zSlice < zmax
                    
                    heightReport = (zSlice - properties.origin(nbElem,3)) / deltaZ;
                    dataCyl.shiftx1 = properties.origin(nbElem,1)*(1-heightReport) + properties.end(nbElem,1)*heightReport;
                    dataCyl.shifty1 = properties.origin(nbElem,2)*(1-heightReport) + properties.end(nbElem,2)*heightReport;
                    
                    dataCyl.shiftx2 = shiftx2Init + (ship.shiftCentreBulb+dist)*cos(dataCyl.alphaDeg2*pi/180);
                    dataCyl.shifty2 = shifty2Init + (ship.shiftCentreBulb+dist)*sin(dataCyl.alphaDeg2*pi/180);
                    
                    zEq = zSlice - (data.shipHeight - ship.hb + hEq);
                    dataCyl.a2 = pEq * sqrt(1-(zEq/hEq)^2);
                    dataCyl.b2 = bEq * sqrt(1-(zEq/hEq)^2);
                    
                    [contactEllipses,solveX,solveY] = Ellipses(dataCyl);
                    
                    contactVecBulb(nbSlice) = contactEllipses;
                    
                else
                    
                    contactVecBulb(nbSlice) = -1;
                    
                end
                
            end
            
            [contactEllipses,labelZ] = max(contactVecBulb);
            
            if contactEllipses == 1
                labelZOK = labelZ;
            end
            
            countSlicesIn = 1;
            for i = slicesIn2(1):slicesIn2(end)
                if contactVecBulb(i) == 1
                    slicesIn(countSlicesIn) = i;
                    countSlicesIn = countSlicesIn + 1;
                end
            end
            if countSlicesIn == 1
                slicesIn3 = 0;
            else
                slicesIn3 = slicesIn(1:countSlicesIn-1);
            end
            
            %         contactVecBulb
            %         labelZ
            
        end
        
        while contactEllipses == 1 && dist <= distMax
            dist = dist - 0.001;
            
            contactVecBulb = zeros(bulbSlices-1,1);
            
            for nbSlice = slicesIn3(1):slicesIn3(end)
                zSlice = data.shipHeight - ship.hb + ship.hBulb - nbSlice*deltaSlice;
                
                if zSlice > zmin && zSlice < zmax
                    
                    heightReport = (zSlice - properties.origin(nbElem,3)) / deltaZ;
                    dataCyl.shiftx1 = properties.origin(nbElem,1)*(1-heightReport) + properties.end(nbElem,1)*heightReport;
                    dataCyl.shifty1 = properties.origin(nbElem,2)*(1-heightReport) + properties.end(nbElem,2)*heightReport;
                    
                    dataCyl.shiftx2 = shiftx2Init + (ship.shiftCentreBulb+dist)*cos(dataCyl.alphaDeg2*pi/180);
                    dataCyl.shifty2 = shifty2Init + (ship.shiftCentreBulb+dist)*sin(dataCyl.alphaDeg2*pi/180);
                    
                    zEq = zSlice - (data.shipHeight - ship.hb + hEq);
                    dataCyl.a2 = pEq * sqrt(1-(zEq/hEq)^2);
                    dataCyl.b2 = bEq * sqrt(1-(zEq/hEq)^2);
                    
                    [contactEllipses,solveX,solveY] = Ellipses(dataCyl);
                    
                    contactVecBulb(nbSlice) = contactEllipses;
                    
                else
                    
                    contactVecBulb(nbSlice) = -1;
                    
                end
                
            end
            
            [contactEllipses,labelZ] = max(contactVecBulb);
            
            if contactEllipses == 1
                labelZOK = labelZ;
                solveResultX = solveX;
                solveResultY = solveY;
            end
            
            %         contactVecBulb
            %         labelZOK
            
        end
        
        contactElement.distBulb(elem,1) = dist;
        contactElement.distBulb(elem,2) = data.shipHeight - ship.hb + ship.hBulb - labelZOK*deltaSlice;
        
        solveX = solveResultX;
        solveY = solveResultY;
        
        contactElement.distBulb(elem,3) = (solveX - dataCyl.shiftx2)*cos(pi+data.shipTrajectory*pi/180) + (solveY - dataCyl.shifty2)*sin(pi+data.shipTrajectory*pi/180);
        contactElement.distBulb(elem,4) = -(solveX - dataCyl.shiftx2)*sin(pi+data.shipTrajectory*pi/180) + (solveY - dataCyl.shifty2)*cos(pi+data.shipTrajectory*pi/180);
        contactElement.distBulb(elem,5) =  - (labelZOK-1)*deltaSlice + ship.hBulb/2;
        
        if dist > distMax
            contactElement.distBulb(elem,1) = 1000;
        end
        
    end
    
end

%% Bulb-node distances

contactNode.distBulb = zeros(length(contactNode.listBulb),5);

if contactNode.listBulb(1) ~= 0
    
    for node = 1:length(contactNode.listBulb)
        
        nbNode = contactNode.listBulb(node);
        
        xPoint = properties.nodes(nbNode,1);
        yPoint = properties.nodes(nbNode,2);
        zPoint = properties.nodes(nbNode,3);
        
        zEq = zPoint - (data.shipHeight - ship.hb + hEq);
        
        dataCyl.a2 = pEq * sqrt(1-(zEq/hEq)^2);
        dataCyl.b2 = bEq * sqrt(1-(zEq/hEq)^2);
        dataCyl.alphaDeg2 = data.shipTrajectory;
        
        if xPoint ~= 0 && yPoint ~= 0 %Leg
            
            xi = atan((data.h(end)-data.h(1)) / sqrt(2*((data.b(2)-data.b(1))/2)^2));
            Rleg = data.De(1)/2;
            dataCyl.a1 = Rleg/(sin(xi));
            dataCyl.b1 = Rleg;
            dataCyl.alphaDeg1 = 45 * sign(xPoint) * sign(yPoint);
            
        else %Crossing of braces
            
            Rbrace = data.De(4)/2;
            level = 1;
            while data.h(level) < zPoint
                level = level + 1;
            end
            deltaZ = data.h(level) - data.h(level-1);
            bottomReport = data.h(level-1)/data.h(end);
            bottomWidth = (data.b(1)*(1-bottomReport) + data.b(end)*bottomReport)/2;
            topReport = data.h(level)/data.h(end);
            topWidth = (data.b(1)*(1-topReport) + data.b(end)*topReport)/2;
            deltaX = bottomWidth + topWidth;
            
            xi = atan(deltaZ/deltaX);
            
            dataCyl.a1 = 2*Rbrace*cos(xi)/(cos(pi/2-2*xi));
            dataCyl.b1 = Rbrace;
            
            if xPoint == 0
                dataCyl.alphaDeg1 = 0;
            else
                dataCyl.alphaDeg1 = 90;
            end
            
        end
        
        dataCyl.shiftx1 = xPoint;
        dataCyl.shifty1 = yPoint;
        
        
        % Compute distance
        
        dist = 0;
        contactEllipses = 0;
        
        while contactEllipses == 0 && dist <= distMax
            dist = dist + Rleg;
            
            dataCyl.shiftx2 = shiftx2Init + (ship.shiftCentreBulb+dist)*cos(dataCyl.alphaDeg2*pi/180);
            dataCyl.shifty2 = shifty2Init + (ship.shiftCentreBulb+dist)*sin(dataCyl.alphaDeg2*pi/180);
            [contactEllipses,solveX,solveY] = Ellipses(dataCyl);
            
        end
        
        while contactEllipses == 1 && dist <= distMax
            dist = dist - 0.1;
            
            dataCyl.shiftx2 = shiftx2Init + (ship.shiftCentreBulb+dist)*cos(dataCyl.alphaDeg2*pi/180);
            dataCyl.shifty2 = shifty2Init + (ship.shiftCentreBulb+dist)*sin(dataCyl.alphaDeg2*pi/180);
            [contactEllipses,solveX,solveY] = Ellipses(dataCyl);
            
        end
        
        while contactEllipses == 0 && dist <= distMax
            dist = dist + 0.01;
            
            dataCyl.shiftx2 = shiftx2Init + (ship.shiftCentreBulb+dist)*cos(dataCyl.alphaDeg2*pi/180);
            dataCyl.shifty2 = shifty2Init + (ship.shiftCentreBulb+dist)*sin(dataCyl.alphaDeg2*pi/180);
            [contactEllipses,solveX,solveY] = Ellipses(dataCyl);
            
        end
        
        while contactEllipses == 1 && dist <= distMax
            dist = dist - 0.001;
            
            dataCyl.shiftx2 = shiftx2Init + (ship.shiftCentreBulb+dist)*cos(dataCyl.alphaDeg2*pi/180);
            dataCyl.shifty2 = shifty2Init + (ship.shiftCentreBulb+dist)*sin(dataCyl.alphaDeg2*pi/180);
            [contactEllipses,solveX,solveY] = Ellipses(dataCyl);
            
        end
        
        contactNode.distBulb(node,1) = dist;
        contactNode.distBulb(node,2) = zPoint;
        
        contactNode.distBulb(node,3) = (solveX - dataCyl.shiftx2)*cos(pi+data.shipTrajectory*pi/180) + (solveY - dataCyl.shifty2)*sin(pi+data.shipTrajectory*pi/180);
        contactNode.distBulb(node,4) = -(solveX - dataCyl.shiftx2)*sin(pi+data.shipTrajectory*pi/180) + (solveY - dataCyl.shifty2)*cos(pi+data.shipTrajectory*pi/180);
        contactNode.distBulb(node,5) =  zPoint - (data.shipHeight - ship.hb + ship.hBulb/2);
        
        if dist > distMax
            contactNode.distBulb(node) = 1000;
        end
        
    end
    
end

%% Minimum displacement before first contact set to zero

% minDistBowElem = min(contactElement.distBow(:,1))
% minDistBulbElem = min(contactElement.distBulb(:,1))
% minDistBowNode = min(contactNode.distBow(:,1))
% minDistBulbNode = min(contactNode.distBulb(:,1))
% 
% if ship.bulbous == 0
%     minDist = min([minDistBowElem , minDistBowNode]);
% else
%     minDist = min([minDistBowElem , minDistBulbElem , minDistBowNode , minDistBulbNode]);
% end

vecTestMin = 0;
if listBowElem ~= 0
    vecTestMin = [vecTestMin ; contactElement.distBow(:,1)];
end
if listBulbElem ~= 0
    vecTestMin = [vecTestMin ; contactElement.distBulb(:,1)];
end
if listBowNode ~= 0
    vecTestMin = [vecTestMin ; contactNode.distBow(:,1)];
end
if listBulbNode ~= 0
    vecTestMin = [vecTestMin ; contactNode.distBulb(:,1)];
end
vecTestMin = vecTestMin(2:end);

minDist = min(vecTestMin);

contactElement.distBow(:,1) = contactElement.distBow(:,1) - minDist;
contactElement.distBulb(:,1) = contactElement.distBulb(:,1) - minDist;
contactNode.distBow(:,1) = contactNode.distBow(:,1) - minDist;
contactNode.distBulb(:,1) = contactNode.distBulb(:,1) - minDist;

M1 = contactElement.distBow;
M2 = contactElement.distBulb;
M3 = contactNode.distBow;
M4 = contactNode.distBulb;

%% Initialisation checks

contactNode.checkBow = zeros(size(contactNode.listBow,1),2);
contactNode.checkBulb = zeros(size(contactNode.listBulb,1),2);
contactElement.checkBow = zeros(size(contactElement.listBow,1),2);
contactElement.checkBulb = zeros(size(contactElement.listBulb,1),2);