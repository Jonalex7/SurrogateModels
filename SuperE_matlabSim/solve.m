function[properties,propertiesInit,contactElement,contactNode,dofs,dispTot,matrices,plasticity,cylinder,cylinderTot,ship,impact,punchProp,punchPropLeft,punchPropRight,basisJacket,Output] = solve(data,ship,properties,dofs,solveOption)

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% First calculations
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

%% Detection of potential impacted elements

[contactElement,contactNode] = automaticContact(data,properties,ship);

for i = 1:size(contactElement.distBow,1)
    if contactElement.distBow(i,1) == 0
        idFirstContactElement = contactElement.listBow(i);
        nbListContactElement = i;
    end
end

nbPossibleContacts = size(contactElement.listBow,1) + size(contactElement.listBulb,1) + size(contactNode.listBow,1) + size(contactNode.listBulb,1);

propertiesNoDiv = properties;
% nbPossibleContacts = 1;
%% Properties of local crushing

%%%% disp('Initialisation')

crushingModel = 1; %1 -> Loic - 2 -> Tim

if properties.defLoc == 1
    
    if crushingModel == 1 % Compute distance as done by Loïc
        
        [cylinder,ship,param] = dataImpact(properties,data,ship);
        ship.phibprime = atan((11.8938+4.8758)/(15.5461-5.6544));
        
        % Limit crushing
        d = limit(cylinder,ship);
        
        % Impact point
        impact = impactPoint(cylinder,ship,d,param);
        
    else
        
        idElem = idFirstContactElement;
        nbListElem = nbListContactElement;
        [cylinder,ship,param,impact] = impactPointTimInit(data,properties,ship,idElem,nbListElem,contactElement);
        ship.phibprime = atan((11.8938+4.8758)/(15.5461-5.6544));
        
    end
    
    % Parameters
    impact = parameters(cylinder,ship,param,impact);
end

%% Punching

[punchProp,punchPropLeft,punchPropRight] = punchingNode(properties,data);


%% Basis of jacket

basisJacket = basisJacketInit(properties,data);


%% Impact position

[data,properties,dofs,force,punchPropLeft,punchPropRight] = impactProperties(data,properties,dofs,punchPropLeft,punchPropRight,solveOption);
vShip = data.shipSpeed*[cos(data.shipTrajectory*pi/180) sin(data.shipTrajectory*pi/180)];
vTotShipMoins1 = norm(vShip);

propertiesInit = properties;

%%


% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% Initialisation of variables
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------


%% Initialisation

% Global movement

dofs.thirdNode = 0;

plasticity = zeros(3,dofs.nbElements);

forceElemStar  = zeros(12,dofs.nbElements);
dispTotStar    = zeros(dofs.nbDofs,1);
forceTotStar   = zeros(dofs.nbDofs,1);

deltaDispTot   = zeros(dofs.nbDofs,1);

dispNode = zeros(1,10^6);
forceNode = zeros(1,10^6);
timeVec = zeros(1,10^6);

dispGlobal = zeros(1,10^6);

Etot = zeros(1,10^6);
Ecin = zeros(1,10^6);
Ecin(1) = 0.5 * data.shipWeight * norm(vShip)^2 * 10^(-6);

RX1 = zeros(1,10^6);
RX2 = zeros(1,10^6);
RX3 = zeros(1,10^6);
RX4 = zeros(1,10^6);
RY1 = zeros(1,10^6);
RY2 = zeros(1,10^6);
RY3 = zeros(1,10^6);
RY4 = zeros(1,10^6);
RTot1 = zeros(1,10^6);
RTot2 = zeros(1,10^6);
RTot3 = zeros(1,10^6);
RTot4 = zeros(1,10^6);

% Local crushing : vertical part
resultVert.mode = 0;
resultVert.E90 = 0;
resultVert.d90 = 0;
resultVert.d90l = 0;
resultVert.d90r = 0;
resultVert.Cv = 0;
resultVert.Cvl = 0;
resultVert.Cvr = 0;
resultVert.gamma = 0;
resultVert.dispLoc = 0;
resultVertStar = resultVert;
forceLocalVert = zeros(1,10^6);

% Local crushing : horizontal part
resultHor.mode = 0;
resultHor.t = 0;
resultHor.dt = 0;
resultHor.Eb0 = 0;
resultHor.Ebt = 0;
resultHor.Em0 = 0;
resultHor.Emt = 0;
resultHor.Cht = 0;
resultHor.Ch = 0;
resultHor.Chl = 0;
resultHor.Chr = 0;
resultHor.R1 = 0;
resultHor.R2 = 0;
resultHor.dR1 = 0;
resultHor.dR2 = 0;
resultHor.E0 = 0;
resultHor.d0 = 0;
resultHor.d0l = 0;
resultHor.d0r = 0;
resultHor.dispLoc = 0;
resultHorStar = resultHor;
forceLocalHor = zeros(1,10^6);

if properties.defLoc == 1
    force.valueDisp = [0 0];
    resultHor = solveHor(properties,data,cylinder,ship,param,impact,force,resultHor,0);
end

% Local crushing : total
sequenceCrushing = zeros(nbPossibleContacts,3);
cylinderTot = zeros(nbPossibleContacts,12);
shipTot = zeros(nbPossibleContacts,12);
impactTot = zeros(nbPossibleContacts,33);
prevDispLocTot = zeros(nbPossibleContacts,1);
prevDispLocTotStar = prevDispLocTot;
resultVertTot = zeros(nbPossibleContacts,9);
resultVertTotStar = resultVertTot;
resultHorTot = zeros(nbPossibleContacts,19);
resultHorTotStar = resultHorTot;
forceLocalMax = 0;
forceLocalTot = zeros(1,10^6);
dispLocal = zeros(1,10^6);
nbContactPoints = 0;
prevContactPoints = 0;

% Punching
punchPropLeft.punchModeNodes = zeros(size(punchPropLeft.leg,1),4);%mode ; dP ; xi
punchPropRight.punchModeNodes = zeros(size(punchPropRight.leg,1),4);
forcePunchTot = zeros(1,10^6);
forcePunchTot(1:3) = ones(1,3)*25;
punchProp.maxDe = 1;
punchProp.redSect = ones(dofs.nbNodes,4);
punchProp.redSectStar = punchProp.redSect;
dispPunch = zeros(1,10^6);
EPunch = zeros(1,10^6);

% Basis of jacket
forceBasisTot = zeros(1,10^6);
dispBasis = zeros(1,10^6);

% Resultant
forceResult = zeros(1,10^6);
dispResult = zeros(1,10^6);
defoMode = zeros(1,10^6);

% ForceBrace
FBrace = zeros(1,10^6);
numberBrace = 38;

% Test values
MFoot1 = zeros(1,10^6);
MFoot2 = zeros(1,10^6);
NFoot = zeros(1,10^6);

thetaFoot = zeros(1,10^6);

vshipXAxis = zeros(1,10^6);
vshipYAxis = zeros(1,10^6);
vshipTotAxis = zeros(1,10^6);

% Second contact
secondContact = 0;
dispSecondContact = 0;

Eglob = 0;
Ecrush = 0;
Epunch = 0;
Ebase = 0;

% Points
goOn = 1;
time = solveOption.deltaT;
point = 1;
pointMaxLocDef = 1;
plotNumber = 1;

%%

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% Full computation
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

%% Solve

while goOn == 1
    
    Ecin(point+1) = 0.5 * data.shipWeight * norm(vShip)^2 * 10^(-6);
    
    % Displacement control
    force.valueDisp = vShip*solveOption.deltaT;
    vshipXAxis(point) = vShip(1);
    vshipYAxis(point) = vShip(2);
    vshipTotAxis(point) = norm(vShip);
    
    %% Check contact points
    
    % Element bow
    
    for i = 1:size(contactElement.listBow,1)
        
        if contactElement.checkBow(i,2) == 0
            distInit = contactElement.distBow(i,1);
            posVertImpact = contactElement.distBow(i,2);
            
            nodeIn = properties.nodein(contactElement.listBow(i));
            posVertIn = properties.nodes(nodeIn,3);
            nodeOut = properties.nodeout(contactElement.listBow(i));
            posVertOut = properties.nodes(nodeOut,3);
            rapportVert = abs((posVertImpact-posVertIn)/(posVertOut-posVertIn));
            
            meanDispX = (properties.nodes(nodeIn,1) - propertiesInit.nodes(nodeIn,1))*(1-rapportVert) + (properties.nodes(nodeOut,1) - propertiesInit.nodes(nodeOut,1))*rapportVert;
            meanDispY = (properties.nodes(nodeIn,2) - propertiesInit.nodes(nodeIn,2))*(1-rapportVert) + (properties.nodes(nodeOut,2) - propertiesInit.nodes(nodeOut,2))*rapportVert;
            distIteration = distInit + meanDispX*cos(data.shipTrajectory*pi/180) + meanDispY*sin(data.shipTrajectory*pi/180);
            contactElement.checkBow(i,1) = distIteration;
            
            if dispResult(point)+norm(force.valueDisp) > distIteration
                nbContactPoints = nbContactPoints+1;
                
                contactElement.checkBow(i,2) = nbContactPoints;
                
                sequenceCrushing(nbContactPoints,1) = 0;
                sequenceCrushing(nbContactPoints,2) = 1;
                sequenceCrushing(nbContactPoints,3) = i;
            end
            
        end
        
    end
    
    % Element bulb
    
    if contactElement.listBulb ~= 0
        
        for i = 1:size(contactElement.listBulb,1)
            
            if contactElement.checkBulb(i,2) == 0
                distInit = contactElement.distBulb(i,1);
                posVertImpact = contactElement.distBulb(i,2);
                
                nodeIn = properties.nodein(contactElement.listBulb(i));
                posVertIn = properties.nodes(nodeIn,3);
                nodeOut = properties.nodeout(contactElement.listBulb(i));
                posVertOut = properties.nodes(nodeOut,3);
                rapportVert = abs((posVertImpact-posVertIn)/(posVertOut-posVertIn));
                
                meanDispX = (properties.nodes(nodeIn,1) - propertiesInit.nodes(nodeIn,1))*(1-rapportVert) + (properties.nodes(nodeOut,1) - propertiesInit.nodes(nodeOut,1))*rapportVert;
                meanDispY = (properties.nodes(nodeIn,2) - propertiesInit.nodes(nodeIn,2))*(1-rapportVert) + (properties.nodes(nodeOut,2) - propertiesInit.nodes(nodeOut,2))*rapportVert;
                distIteration = distInit + meanDispX*cos(data.shipTrajectory*pi/180) + meanDispY*sin(data.shipTrajectory*pi/180);
                contactElement.checkBulb(i,1) = distIteration;
                
                if dispResult(point)+norm(force.valueDisp) > distIteration
                    nbContactPoints = nbContactPoints+1;
                    
                    contactElement.checkBulb(i,2) = nbContactPoints;
                    
                    sequenceCrushing(nbContactPoints,1) = 1;
                    sequenceCrushing(nbContactPoints,2) = 1;
                    sequenceCrushing(nbContactPoints,3) = i;
                end
                
            end
            
        end
        
    end
    
    % Node bow
    
    if contactNode.listBow ~= 0
        
        for i = 1:size(contactNode.listBow,1)
            
            if contactNode.checkBow(i,2) == 0
                nodeCheck = contactNode.listBow(i);
                distInit = contactNode.distBow(i);
                
                distIteration = distInit + (properties.nodes(nodeCheck,1) - propertiesInit.nodes(nodeCheck,1))*cos(data.shipTrajectory*pi/180) + (properties.nodes(nodeCheck,2) - propertiesInit.nodes(nodeCheck,2))*sin(data.shipTrajectory*pi/180);
                contactNode.checkBow(i,1) = distIteration;
                
                if dispResult(point)+norm(force.valueDisp) > distIteration
                    nbContactPoints = nbContactPoints+1;
                    
                    contactNode.checkBow(i,2) = nbContactPoints;
                    
                    sequenceCrushing(nbContactPoints,1) = 0;
                    sequenceCrushing(nbContactPoints,2) = 0;
                    sequenceCrushing(nbContactPoints,3) = i;
                end
            end
            
        end
        
    end
    
    % Node bulb
    
    if contactNode.listBulb ~= 0
        
        for i = 1:size(contactNode.listBulb,1)
            
            if contactNode.checkBulb(i,2) == 0
                nodeCheck = contactNode.listBulb(i);
                distInit = contactNode.distBulb(i);
                
                distIteration = distInit + (properties.nodes(nodeCheck,1) - propertiesInit.nodes(nodeCheck,1))*cos(data.shipTrajectory*pi/180) + (properties.nodes(nodeCheck,2) - propertiesInit.nodes(nodeCheck,2))*sin(data.shipTrajectory*pi/180);
                contactNode.checkBulb(i,1) = distIteration;
                
                if dispResult(point)+norm(force.valueDisp) > distIteration
                    nbContactPoints = nbContactPoints+1;
                    
                    contactNode.checkBulb(i,2) = nbContactPoints;
                    
                    sequenceCrushing(nbContactPoints,1) = 1;
                    sequenceCrushing(nbContactPoints,2) = 0;
                    sequenceCrushing(nbContactPoints,3) = i;
                end
            end
            
        end
        
    end
    
    impact.sequenceCrushing = sequenceCrushing;
%     sequenceCrushing
    
    
    %% Crushing properties
    
    if nbContactPoints > prevContactPoints
        for contactPoint = prevContactPoints+1:nbContactPoints
            if sequenceCrushing(contactPoint,2) == 0 || sequenceCrushing(contactPoint,2) == 1
                
                nbListElem = sequenceCrushing(contactPoint,3);
                if sequenceCrushing(contactPoint,1) == 0
                    idElem = contactElement.listBow(nbListElem);
                else
                    idElem = contactElement.listBulb(nbListElem);
                end
                
                [cylinder2,ship2,param2,impact2,cylinderVec,shipVec,impactVec] = impactPointTim(data,properties,propertiesNoDiv,ship,idElem,nbListElem,contactElement,contactNode,sequenceCrushing,contactPoint);
                cylinderTot(contactPoint,:) = cylinderVec;
                shipTot(contactPoint,:) = shipVec;
                impactTot(contactPoint,:) = impactVec;
                
            end
        end
    end
    
    
    %% Local crushing
    
    forceLocal2 = 0;
    
    for contactPoint = 1:nbContactPoints
        
        cylinderVec = cylinderTot(contactPoint,:);
        shipVec = shipTot(contactPoint,:);
        impactVec = impactTot(contactPoint,:);
        
        prevDispLocTot = prevDispLocTotStar;
        resultVertTot = resultVertTotStar;
        resultHorTot = resultHorTotStar;
        
        prevDispLoc = prevDispLocTot(contactPoint);
        resultVertVec = resultVertTot(contactPoint,:);
        resultHorVec = resultHorTot(contactPoint,:);
        
        % Vert
        [resultVert,resultVertVec] = solveVert2(properties,data,cylinderVec,shipVec,param,impactVec,force,resultVertVec,prevDispLoc);
        resultVertTot(contactPoint,:) = resultVertVec;
        forceLocalVert2 = abs(resultVert.Ptot90);
        
        % Hor
        [resultHor,resultHorVec] = solveHor2(properties,data,cylinderVec,shipVec,param,impactVec,force,resultHorVec,prevDispLoc);
        resultHorTot(contactPoint,:) = resultHorVec;
        forceLocalHor2 = abs(resultHor.Ptot0);
        
        % Tot
        PX = (abs(resultHor.PXtot0) * (1 - 2 * cylinderVec(3) / pi) + abs(resultVert.PXtot90) * 2 * cylinderVec(3) / pi) * param.k1; %cylinderVec(3) = cylinder.dzeta
        PY = (abs(resultHor.PYtot0) * (1 - 2 * cylinderVec(3) / pi) + abs(resultVert.PYtot90) * 2 * cylinderVec(3) / pi) * param.k2;
        
        forceLocal2 = forceLocal2 + (PX^2 + PY^2)^0.5;
        
        % Update
        prevDispLocTot(contactPoint) = resultVert.delta;
        
    end
    forceLocal2 = forceLocal2;
    forceLocalTot(point+1) = forceLocal2;%/1000;
    
    if properties.defLoc == 1
        prevDispLoc = dispLocal(point);
        
        resultVert = resultVertStar;
        resultHor = resultHorStar;
        
        if prevDispLoc > cylinder.maxLocalCrushingNode && data.impactedZone == 0 %0.95*data.De(1)
            
            forceLocalTot(point+1) = 4 * forceResult(pointMaxLocDef);
        
        else
        
            resultVert = solveVert(properties,data,cylinder,ship,param,impact,force,resultVert,prevDispLoc);
            forceLocalVert(point+1) = resultVert.Ptot90;
            forceVertabc = resultVert.Ptot90;
            
            resultHor = solveHor(properties,data,cylinder,ship,param,impact,force,resultHor,prevDispLoc);
            forceLocalHor(point+1) = resultHor.Ptot0;
            forceHorabc = resultHor.Ptot0;
            
            PX = (resultHor.PXtot0 * (1 - 2 * cylinder.dzeta / pi) + resultVert.PXtot90 * 2 * cylinder.dzeta / pi) * param.k1;
            PY = (resultHor.PYtot0 * (1 - 2 * cylinder.dzeta / pi) + resultVert.PYtot90 * 2 * cylinder.dzeta / pi) * param.k2;
            
            forceLocalTot(point+1) = (PX^2 + PY^2)^0.5;
            pointMaxLocDef = point+1;
            forceLocalMax = (PX^2 + PY^2)^0.5;
            
        end
    else
        forceLocalTot(point+1) = 0;
    end
    
    
    %% Global movement
    
%     data.fy = [ones(4,1)*255*10^6 ; 255*10^9];
    
    forceElem  = forceElemStar;
    dispTot    = dispTotStar;
    forceTot   = forceTotStar;
    
    % Initialisation for convergence
    error = 1;
    step = 0;
    
    deltaForceTotimoins1 = zeros(dofs.nbDofs,1);
    
    while error > solveOption.epsilon && step < 7
    
        % Stiffness matrices
        deltaForceVec = zeros(dofs.nbDofs,1);
        [matrices,deltaForceVec] = stiffness(solveOption,properties,forceElem,deltaForceVec,plasticity,dofs,punchProp);
        
        % Displacement control
        Kmod = matrices.K;
        for j = 1:length(force.dofsForce)
            Kmod(:,force.dofsForce(j)) = 0;
            Kmod(force.dofsForce(j),:) = 0;
            Kmod(force.dofsForce(j),force.dofsForce(j)) = 1;
            
            for i = 1:length(deltaForceVec)
                deltaForceVec(i) = deltaForceVec(i) - matrices.K(i,force.dofsForce(j))*force.valueDisp(j);
            end
        end
        
        for j = 1:length(force.dofsForce)
            deltaForceVec(force.dofsForce(j)) = force.valueDisp(j);
        end
        
        % Solve (variational)
        if rcond(Kmod(dofs.freeDofs,dofs.freeDofs)) < 10^(-25) || isnan(rcond(Kmod(dofs.freeDofs,dofs.freeDofs)))
            if rcond(Kmod(dofs.freeDofs,dofs.freeDofs)) < 10^(-25)
                disp(['rcond = ',num2str(rcond(Kmod(dofs.freeDofs,dofs.freeDofs)))])
            else
                disp('NaN')
            end
            time = 1000;
            break
        end
        
        deltaDispTot(dofs.freeDofs) = Kmod(dofs.freeDofs,dofs.freeDofs)\deltaForceVec(dofs.freeDofs);
        
        % Internal forces (variational)
        deltaForceElem = internalForces(properties,matrices,deltaDispTot,dofs);
        
        % Total forces (variational)
        deltaForceTot = matrices.K*deltaDispTot;
        
        % Total values
        dispTot   = dispTotStar   + deltaDispTot;
        forceElem = forceElemStar + deltaForceElem;
        forceTot  = forceTotStar  + deltaForceTot;
        
        % Plasticity
        plasticity = plasticityTest(properties,dofs,matrices,dispTot,forceElem,plasticity,punchProp);
        
        % Convergence
        sumError = 0;
        nbItem = 0;
        takeAccount = max(abs(deltaForceTot))/1000;
        for i = 1:dofs.nbDofs
            if abs(forceTot(i)) > takeAccount
                sumError = sumError + abs(deltaForceTot(i) - deltaForceTotimoins1(i))/abs(deltaForceTot(i));
                nbItem = nbItem + 1;
            end
        end
        error = sumError/nbItem;
        deltaForceTotimoins1 = deltaForceTot;
        
        % Update of properties
        for i = 1:size(properties.nodes,1)
            for j = 1:3
                properties.nodes(i,j) = propertiesInit.nodes(i,j) + dispTot(properties.nodes(i,j+3));
            end
        end
        [properties,dofs] = autoCalculation(data,properties,dofs,solveOption,0);
        
        % If needed, third node in element
        thirdHinge = 0;
        elemHinge = zeros(dofs.nbElements,1);
        for i = 1:dofs.nbElements-2*dofs.thirdNode
            if plasticity(1,i) == 1
                thirdHinge = thirdHinge + 1;
                elemHinge(thirdHinge) = i;
            end
        end
        if thirdHinge ~= 0
            elemHinge = elemHinge(1:thirdHinge);
            dofs.thirdNode = dofs.thirdNode + thirdHinge;
            break
        end
        
        % Incrementation
        step = step + 1;
    
    end
    
%     data.fy = [ones(4,1)*317*10^6 ; 317*10^9];
    
%     step
    
    %% Punching
    
    if point == 2
        
        [punchPropLeft,punchPropRight] = punchingNode2(propertiesInit,data,forceElem,punchProp,punchPropLeft,punchPropRight);
        
        punchPropLeft.behaviourNodesStar = punchPropLeft.behaviourNodes;
        punchPropRight.behaviourNodesStar = punchPropRight.behaviourNodes;
        
    elseif point >= 3
        
        punchPropLeft.behaviourNodes = punchPropLeft.behaviourNodesStar;
        punchPropRight.behaviourNodes = punchPropRight.behaviourNodesStar;
        
        [punchProp,punchPropLeft,punchPropRight,resultPunch] = punchingScenarios2(properties,data,force,forceElem,punchProp,punchPropLeft,punchPropRight);
        ddeltaShip = norm(force.valueDisp);
        dAddLocPunch = 0;
        if norm(resultPunch.FPunch) == 0
            forcePunchTot(point+1) = 2*forceResult(point);
        else
            if resultPunch.FLeft == 0
                resultPunch.FLeft = forceLocalMax;
                dAddLocPunch = ddeltaShip * min(cos(data.shipTrajectory*pi/180),sin(data.shipTrajectory*pi/180));
            end
            if resultPunch.FRight == 0
                resultPunch.FRight = forceLocalMax;
                dAddLocPunch = ddeltaShip * min(cos(data.shipTrajectory*pi/180),sin(data.shipTrajectory*pi/180));
            end
            FPunchX = resultPunch.FLeft*cos(punchPropLeft.angle) + resultPunch.FRight*cos(punchPropRight.angle);
            FPunchY = resultPunch.FLeft*sin(punchPropLeft.angle) + resultPunch.FRight*sin(punchPropRight.angle);
            forcePunchTot(point+1) = FPunchY*sin(data.shipTrajectory*pi/180) + FPunchX*cos(data.shipTrajectory*pi/180);
%             forcePunchTot(point+1) = sqrt((resultPunch.FLeft)^2 + (resultPunch.FRight)^2);
        end
    
    end
    
    forcePunchTot(point+1) = forcePunchTot(point+1)*1.5;
    
    %% Basis of jacket
    
    prevDispBasis = dispBasis(point);
    resultBasisJacket = solveBasisJacket(basisJacket,force,data,prevDispBasis);
    forceBasisTot(point+1) = resultBasisJacket.FAxisX*cos(data.shipTrajectory*pi/180) + resultBasisJacket.FAxisY*sin(data.shipTrajectory*pi/180);
    
    %% Post-treatment
    
    if time == 1000
        break
    end
    
    if thirdHinge ~= 0
        
        [properties,propertiesInit,dofs,forceElemStar,dispTotStar,forceTotStar,plasticity,punchProp,punchPropLeft,punchPropRight] = thirdNode(properties,propertiesInit,data,dofs,forceElemStar,dispTotStar,forceTotStar,plasticity,punchProp,punchPropLeft,punchPropRight,elemHinge,solveOption);
% % % %         M1 = draw(propertiesInit,dispTotStar,dofs,solveOption,plasticity,time);
% % % %         M(:,plotNumber) = M1;
        plotNumber = plotNumber + 1;
        time = time - solveOption.deltaT;
        
    else
        
        % Draw
        
        if abs(time - plotNumber*solveOption.deltaTFigures) < solveOption.deltaT/10
% % % %             M1 = draw(propertiesInit,dispTot,dofs,solveOption,plasticity,time);
% % % %             M(:,plotNumber) = M1;
            plotNumber = plotNumber+1;
        end
        
        % Save values
        dispNode(point+1) = sqrt((dispTot(force.dofsForce(1)))^2 + (dispTot(force.dofsForce(2)))^2);
        forceNode(point+1) = sqrt((forceTot(force.dofsForce(1)))^2 + (forceTot(force.dofsForce(2)))^2) * 10^(-6);
        timeVec(point+1) = time;
        
%         RX1(point+1) = -matrices.K(1,:)*dispTot * 10^(-6);
%         RX2(point+1) = -matrices.K(7,:)*dispTot * 10^(-6);
%         RX3(point+1) = -matrices.K(13,:)*dispTot * 10^(-6);
%         RX4(point+1) = -matrices.K(19,:)*dispTot * 10^(-6);
%         
%         RY1(point+1) = -matrices.K(2,:)*dispTot * 10^(-6);
%         RY2(point+1) = -matrices.K(8,:)*dispTot * 10^(-6);
%         RY3(point+1) = -matrices.K(14,:)*dispTot * 10^(-6);
%         RY4(point+1) = -matrices.K(20,:)*dispTot * 10^(-6);
        
        % Comparison of forces
        
        if forceLocalTot(point+1) <= forceNode(point+1) && forceLocalTot(point+1) <= forcePunchTot(point+1) && forceLocalTot(point+1) <= forceBasisTot(point+1) && properties.defLoc == 1
%             disp('Local')
            defoMode(point+1) = 1;
            forceResult(point+1) = forceLocalTot(point+1);
            dispResult(point+1) = dispResult(point) + norm(force.valueDisp);
            dispLocal(point+1) = dispLocal(point) + norm(force.valueDisp);
            dispGlobal(point+1) = dispGlobal(point);
            dispPunch(point+1) = dispPunch(point);
            dispBasis(point+1) = dispBasis(point);
            
            if forceResult(point+1) < forceResult(point)
                forceResult(point+1) = forceResult(point);
            end
            
            RX1(point+1) = RX1(point);
            RX2(point+1) = RX2(point);
            RX3(point+1) = RX3(point);
            RX4(point+1) = RX4(point);
            
            RY1(point+1) = RY1(point);
            RY2(point+1) = RY2(point);
            RY3(point+1) = RY3(point);
            RY4(point+1) = RY4(point);
            
            resultVertStar = resultVert;
            resultHorStar = resultHor;
            
            prevDispLocTotStar = prevDispLocTot;
            resultVertTotStar = resultVertTot;
            resultHorTotStar = resultHorTot;
            
            forceAxis = forceResult(point+1)*[cos(data.shipTrajectory*pi/180) sin(data.shipTrajectory*pi/180)]*10^6;
            
            EPunch(point+1) = EPunch(point);
            
            Ecrush = Ecrush + (forceResult(point)+forceResult(point+1))/2 * (dispResult(point+1)-dispResult(point));
            
        elseif forcePunchTot(point+1) < forceNode(point+1) && forcePunchTot(point+1) < forceLocalTot(point+1) && forcePunchTot(point+1) < forceBasisTot(point+1)
%             disp('Punching')
            defoMode(point+1) = 2;
            forceResult(point+1) = forcePunchTot(point+1);
            dispResult(point+1) = dispResult(point) + norm(force.valueDisp);
            dispLocal(point+1) = dispLocal(point) + dAddLocPunch;
            dispGlobal(point+1) = dispGlobal(point);
            dispPunch(point+1) = dispPunch(point) + norm(force.valueDisp);
            dispBasis(point+1) = dispBasis(point);
            
            if forceResult(point+1) < forceResult(point)
                forceResult(point+1) = forceResult(point);
            end
            
            RX1(point+1) = RX1(point);
            RX2(point+1) = RX2(point);
            RX3(point+1) = RX3(point);
            RX4(point+1) = RX4(point);
            
            RY1(point+1) = RY1(point);
            RY2(point+1) = RY2(point);
            RY3(point+1) = RY3(point);
            RY4(point+1) = RY4(point);
            
            punchPropLeft.behaviourNodesStar = punchPropLeft.behaviourNodes;
            punchPropRight.behaviourNodesStar = punchPropRight.behaviourNodes;
            
            punchPropLeft.punchModeNodesStar = punchPropLeft.punchModeNodes;
            punchPropRight.punchModeNodesStar = punchPropRight.punchModeNodes;
            
            punchProp.redSectStar = punchProp.redSect;
            
            forceAxis = forceResult(point+1)*[cos(data.shipTrajectory*pi/180) sin(data.shipTrajectory*pi/180)]*10^6;
            
            EPunch(point+1) = EPunch(point) + (forcePunchTot(point)+forcePunchTot(point+1))/2 * (dispResult(point+1)-dispResult(point));
            
            Epunch = Epunch + (forceResult(point)+forceResult(point+1))/2 * (dispResult(point+1)-dispResult(point));
            
        elseif forceBasisTot(point+1) < forceNode(point+1) && forceBasisTot(point+1) < forceLocalTot(point+1) && forceBasisTot(point+1) < forcePunchTot(point+1)
%             disp('Basis of jacket')
            defoMode(point+1) = 4;
            forceResult(point+1) = forceBasisTot(point+1);
            dispResult(point+1) = dispResult(point) + norm(force.valueDisp);
            dispLocal(point+1) = dispLocal(point);
            dispGlobal(point+1) = dispGlobal(point);
            dispPunch(point+1) = dispPunch(point);
            dispBasis(point+1) = dispBasis(point) + norm(force.valueDisp);
            
            if forceResult(point+1) < forceResult(point)
                forceResult(point+1) = forceResult(point);
            end
            
            RX1(point+1) = RX1(point);
            RX2(point+1) = RX2(point);
            RX3(point+1) = RX3(point);
            RX4(point+1) = RX4(point);
            
            RY1(point+1) = RY1(point);
            RY2(point+1) = RY2(point);
            RY3(point+1) = RY3(point);
            RY4(point+1) = RY4(point);
            
            forceAxis = forceResult(point+1)*[cos(data.shipTrajectory*pi/180) sin(data.shipTrajectory*pi/180)]*10^6;
            
            EPunch(point+1) = EPunch(point);
            
            Ebase = Ebase + (forceResult(point)+forceResult(point+1))/2 * (dispResult(point+1)-dispResult(point));
            
        else
%             disp('Global')
            defoMode(point+1) = 3;
            forceResult(point+1) = forceNode(point+1);
            dispResult(point+1) = dispResult(point) + norm(force.valueDisp);
            dispGlobal(point+1) = dispGlobal(point) + norm(force.valueDisp);
            dispLocal(point+1) = dispLocal(point);
            dispPunch(point+1) = dispPunch(point);
            dispBasis(point+1) = dispBasis(point);
            
            if forceResult(point+1) < forceResult(point)
                forceResult(point+1) = forceResult(point);
            end
            
            forceElemStar  = forceElem;
            dispTotStar    = dispTot;
            forceTotStar   = forceTot;
            
            RX1(point+1) = RX1(point) - matrices.K(1,:)*deltaDispTot * 10^(-6);
            RX2(point+1) = RX2(point) - matrices.K(7,:)*deltaDispTot * 10^(-6);
            RX3(point+1) = RX3(point) - matrices.K(13,:)*deltaDispTot * 10^(-6);
            RX4(point+1) = RX4(point) - matrices.K(19,:)*deltaDispTot * 10^(-6);
            
            RY1(point+1) = RY1(point) - matrices.K(2,:)*deltaDispTot * 10^(-6);
            RY2(point+1) = RY2(point) - matrices.K(8,:)*deltaDispTot * 10^(-6);
            RY3(point+1) = RY3(point) - matrices.K(14,:)*deltaDispTot * 10^(-6);
            RY4(point+1) = RY4(point) - matrices.K(20,:)*deltaDispTot * 10^(-6);
            
            forceAxis = forceResult(point+1)*[cos(data.shipTrajectory*pi/180) sin(data.shipTrajectory*pi/180)]*10^6;
            
            EPunch(point+1) = EPunch(point);
            
            Eglob = Eglob + (forceResult(point)+forceResult(point+1))/2 * (dispResult(point+1)-dispResult(point));
            
        end
        
        if secondContact == 0 && point >= 3
            dispSecondContact = dispResult(point+1);
            h2Contact = data.shipHeight - properties.nodes(25,3) - 1.5;
            distNodeInit = h2Contact * (tan(pi/2 - ship.phibprime) - tan(basisJacket.mu));
            distNodeSupp = dispTotStar(6*25-5) + dispBasis(point+1);
            distNodeSupp1 = dispTotStar(6*25-5);
            distNodeSupp2 = punchPropRight.behaviourNodesStar(3,3);
            distNodeSupp3 = dispBasis(point+1);
            diffDist = distNodeInit + distNodeSupp - dispResult(point+1);
            if dispResult(point+1) > distNodeInit + distNodeSupp
                secondContact = 1;
            end
        end
        
        for nbCylinder = 1:length(properties.L)-2*dofs.thirdNode
            epsilon = (properties.L(nbCylinder) - propertiesInit.L(nbCylinder)) / (propertiesInit.L(nbCylinder));
            if epsilon > 0.2
                disp(nbCylinder)
                disp(epsilon)
            end
        end
        
        MFoot1(point+1) = forceElemStar(5,5);
        MFoot2(point+1) = forceElemStar(6,5);
        NFoot(point+1) = forceElemStar(1,5);
        
        dispFoot = dispTotStar(6*10-5);
        thetaFoot(point+1) = atan(dispFoot/1.5);
    
        Etot(point+1) = Etot(point) + (forceResult(point)+forceResult(point+1))/2 * (dispResult(point+1)-dispResult(point));
        
        RTot1(point+1) = sqrt((RX1(point+1))^2+(RY1(point+1))^2);
        RTot2(point+1) = sqrt((RX2(point+1))^2+(RY2(point+1))^2);
        RTot3(point+1) = sqrt((RX3(point+1))^2+(RY3(point+1))^2);
        RTot4(point+1) = sqrt((RX4(point+1))^2+(RY4(point+1))^2);
        
        % Force in considered brace
        FBrace(point+1) = forceElemStar(1,numberBrace) * 10^(-6);
        
        % Update of speed
        if data.VInit == 1
            vShip = vShip - forceAxis./data.shipWeight*solveOption.deltaT;
            vTotShip = norm(vShip);
            
            if vTotShip > vTotShipMoins1
                goOn = 0;
                time = time - solveOption.deltaT;
            end
            vTotShipMoins1 = vTotShip;
            
        else
            if norm(dispTot(force.dofsForce)) >= data.dispmax
                goOn = 0;
            end
            if dispResult(point+1) >= data.dispmax
                goOn = 0;
            end
            disp(['Time = ',num2str(time)])
        end
        
        % Step correction
        
        if goOn == 1
            time = time + solveOption.deltaT;
            point = point + 1;
        end
        
        prevContactPoints = nbContactPoints;
    
    end

end

forceResultFinal = zeros(1,point);
for i = 1:point
    forceResultFinal(i) = forceResult(i);
    if i >= 2 && forceResultFinal(i) < forceResultFinal(i-1) && defoMode(i) == 2
        forceResultFinal(i) = forceResultFinal(i-1);
    end
end

%%

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% Graphs and results
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

%% Results

% disp(['Time before ship stops is equal to ',num2str(time),' seconds'])
% % % % disp(['2nd contact for crushing of ',num2str(dispSecondContact),' meters'])
% % % % M1 = draw(propertiesInit,dispTot,dofs,solveOption,plasticity,time); %%time
% % % % M(:,plotNumber) = M1;

%%%%% To plot the 3D Jacket at specific time
% % % % M1 = draw(propertiesInit,dispTot,dofs,solveOption,plasticity,0.31); %%time
% % % % M(:,plotNumber) = M1;

ab = Eglob;
cd = Ecrush;
ef = Epunch;
gh = Ebase;

% slopeEl = forceNode(2)/dispNode(2);
% dispEl = [0 max(dispNode(1:point-1))];
% forceEl = slopeEl*dispEl;
% 
% figure;
% hold on
% plot(dispNode(1:point),forceNode(1:point),'x-','markerSize',16)
% plot(dispEl,forceEl,'r--')
% xlabel('Displacement [m]','FontSize',solveOption.fontsizeLabel)
% ylabel('Force [MN]','FontSize',solveOption.fontsizeLabel)
% set(gca,'FontSize',solveOption.fontsizeAxis)
% hold off
% 
% figure;
% hold on
% plot(timeVec(1:point),forceNode(1:point))
% xlabel('Time [s]','FontSize',solveOption.fontsizeLabel)
% ylabel('Force [MN]','FontSize',solveOption.fontsizeLabel)
% set(gca,'FontSize',solveOption.fontsizeAxis)
% hold off
% 
% figure;
% hold on
% plot(timeVec(1:point),dispNode(1:point))
% xlabel('Time [s]','FontSize',solveOption.fontsizeLabel)
% ylabel('Displacement [m]','FontSize',solveOption.fontsizeLabel)
% set(gca,'FontSize',solveOption.fontsizeAxis)
% hold off

% % % % figure;
% % % % hold on
% % % % subplot(1,3,1)
% % % % hold on
% % % % plot(timeVec(1:point),RX1(1:point),'b')
% % % % plot(timeVec(1:point),RX2(1:point),'r')
% % % % plot(timeVec(1:point),RX3(1:point),'c')
% % % % plot(timeVec(1:point),RX4(1:point),'g')
% % % % plot(timeVec(1:point),RX1(1:point)+RX2(1:point)+RX3(1:point)+RX4(1:point),'k')
% % % % title('Axis X','FontSize',solveOption.fontsizeLabel)
% % % % xlabel('Time [s]','FontSize',solveOption.fontsizeLabel)
% % % % ylabel('Reaction forces [MN]','FontSize',solveOption.fontsizeLabel)
% % % % set(gca,'FontSize',solveOption.fontsizeAxis)
% % % % hold off
% % % % subplot(1,3,2)
% % % % hold on
% % % % plot(timeVec(1:point),RY1(1:point),'b')
% % % % plot(timeVec(1:point),RY2(1:point),'r')
% % % % plot(timeVec(1:point),RY3(1:point),'c')
% % % % plot(timeVec(1:point),RY4(1:point),'g')
% % % % plot(timeVec(1:point),RY1(1:point)+RY2(1:point)+RY3(1:point)+RY4(1:point),'k')
% % % % title('Axis Y','FontSize',solveOption.fontsizeLabel)
% % % % xlabel('Time [s]','FontSize',solveOption.fontsizeLabel)
% % % % ylabel('Reaction forces [MN]','FontSize',solveOption.fontsizeLabel)
% % % % set(gca,'FontSize',solveOption.fontsizeAxis)
% % % % hold off
% % % % subplot(1,3,3)
% % % % hold on
% % % % plot(timeVec(1:point),RTot1(1:point),'b')
% % % % plot(timeVec(1:point),RTot2(1:point),'r')
% % % % plot(timeVec(1:point),RTot3(1:point),'c')
% % % % plot(timeVec(1:point),RTot4(1:point),'g')
% % % % plot(timeVec(1:point),RTot1(1:point)+RTot2(1:point)+RTot3(1:point)+RTot4(1:point),'k')
% % % % title('Total','FontSize',solveOption.fontsizeLabel)
% % % % xlabel('Time [s]','FontSize',solveOption.fontsizeLabel)
% % % % ylabel('Reaction forces [MN]','FontSize',solveOption.fontsizeLabel)
% % % % set(gca,'FontSize',solveOption.fontsizeAxis)
% % % % legendBoundary = legend('RX1','RX2','RX3','RX4','Total');
% % % % set(legendBoundary,'FontSize',solveOption.fontsizeLegend)
% % % % hold off
% % % % 
% % % % figure;
% % % % hold on
% % % % plot(dispLocal(1:point),forceLocalVert(1:point),'b--')
% % % % plot(dispLocal(1:point),forceLocalHor(1:point),'r--')
% % % % plot(dispLocal(1:point),forceLocalTot(1:point),'k')
% % % % xlabel('Displacement [m]','FontSize',solveOption.fontsizeLabel)
% % % % ylabel('Local force [MN]','FontSize',solveOption.fontsizeLabel)
% % % % set(gca,'FontSize',solveOption.fontsizeAxis)
% % % % hold off
% % % % 
% % % % figure;
% % % % hold on
% % % % plot(timeVec(1:point),dispGlobal(1:point),'b')
% % % % plot(timeVec(1:point),dispLocal(1:point),'r')
% % % % plot(timeVec(1:point),dispPunch(1:point),'g')
% % % % plot(timeVec(1:point),dispBasis(1:point),'m')
% % % % plot(timeVec(1:point),dispResult(1:point),'k')
% % % % title('Global, local and total displacements','FontSize',solveOption.fontsizeTitle)
% % % % xlabel('Time [s]','FontSize',solveOption.fontsizeLabel)
% % % % ylabel('Displacement [m]','FontSize',solveOption.fontsizeLabel)
% % % % legendComp = legend('Global','Local','Punching','Basis','Total');
% % % % set(legendComp,'FontSize',solveOption.fontsizeLegend)
% % % % set(gca,'FontSize',solveOption.fontsizeAxis)
% % % % hold off
% % % % 
% % % % figure;
% % % % hold on
% % % % plot(timeVec(1:point),FBrace(1:point))
% % % % xlabel('Time [s]','FontSize',solveOption.fontsizeLabel)
% % % % ylabel('Force in brace [MN]','FontSize',solveOption.fontsizeLabel)
% % % % set(gca,'FontSize',solveOption.fontsizeAxis)
% % % % hold off
% % % % 
% % % % figure;         %%GRAPH OF FORCES CONTRIBUTION
% % % % hold on
% % % % plot(timeVec(1:point),forceNode(1:point),'b')
% % % % plot(timeVec(1:point),forceLocalTot(1:point),'r')
% % % % plot(timeVec(1:point),forcePunchTot(1:point),'g')
% % % % plot(timeVec(1:point),forceBasisTot(1:point),'m')
% % % % xlabel('Time [s]','FontSize',solveOption.fontsizeLabel)
% % % % ylabel('Forces [MN]','FontSize',solveOption.fontsizeLabel)
% % % % legendComp = legend('Global','Local','Punching');
% % % % set(legendComp,'FontSize',solveOption.fontsizeLegend)
% % % % set(gca,'FontSize',solveOption.fontsizeAxis)
% % % % hold off
% % % % 
% % % % figure;
% % % % hold on
% % % % plot(timeVec(1:point),Etot(1:point),'b')
% % % % plot(timeVec(1:point),Ecin(1:point),'r')
% % % % plot(timeVec(1:point),Etot(1:point)+Ecin(1:point),'k')
% % % % legendComp = legend('E_{int}','E_{cin}','E_{tot}');
% % % % set(legendComp,'FontSize',solveOption.fontsizeLegend)
% % % % xlabel('Time [s]','FontSize',solveOption.fontsizeLabel)
% % % % ylabel('Energy [MJ]','FontSize',solveOption.fontsizeLabel)
% % % % set(gca,'FontSize',solveOption.fontsizeAxis)
% % % % hold off

% % % % figure;
% % % % hold on
% % % % plot(dispResult(1:point),forceResult(1:point))
% % % % xlabel('Total displacement [m]','FontSize',solveOption.fontsizeLabel)
% % % % ylabel('Resistant force [MN]','FontSize',solveOption.fontsizeLabel)
% % % % set(gca,'FontSize',solveOption.fontsizeAxis)
% % % % hold off

% % % % figure;
% % % % hold on
% % % % plot(dispResult(1:point),Etot(1:point),'b')
% % % % plot(dispResult(1:point),Ecin(1:point),'r')
% % % % plot(dispResult(1:point),Etot(1:point)+Ecin(1:point),'k')
% % % % legendComp = legend('E_{int}','E_{cin}','E_{tot}');
% % % % set(legendComp,'FontSize',solveOption.fontsizeLegend)
% % % % xlabel('Total displacement [m]','FontSize',solveOption.fontsizeLabel)
% % % % ylabel('Energy [MJ]','FontSize',solveOption.fontsizeLabel)
% % % % set(gca,'FontSize',solveOption.fontsizeAxis)
% % % % hold off

% % % % figure;
% % % % hold on
% % % % plot(timeVec(1:point),forcePunchTot(1:point))
% % % % xlabel('Time [s]','FontSize',solveOption.fontsizeLabel)
% % % % ylabel('Resistant force punching [MN]','FontSize',solveOption.fontsizeLabel)
% % % % set(gca,'FontSize',solveOption.fontsizeAxis)
% % % % hold off

% % % % % figure;         %%GRAPH OF TOTAL RESISTANCE FORCE
% % % % % hold on
% % % % % plot(dispResult(1:point),forceResultFinal)
% % % % % xlabel('Total displacement [m]','FontSize',solveOption.fontsizeLabel)
% % % % % ylabel('Final Resistance Force [MN]','FontSize',solveOption.fontsizeLabel)
% % % % % set(gca,'FontSize',solveOption.fontsizeAxis)
% % % % % hold off

%JONATHAN (OUTPUT)---------------------------------------------------------

%saving maximum PENETRATION and FORCE
MaxForce = max(forceResultFinal);
index = find(forceResultFinal==MaxForce); %index may find several similar values CHECK
Penetration = dispResult(index(1));    %as expected similar index, only first one is considered
timeIndex = timeVec(index(1));      %as expected similar index, only first one is considered
Output = [timeIndex Penetration MaxForce];
%save('Output.mat','timeIndex','Penetration','MaxForce')

%JONATHAN (OUTPUT)---------------------------------------------------------

% figure;
% hold on
% plot(dispResult(1:point),MFoot1(1:point),'b')
% plot(dispResult(1:point),MFoot2(1:point),'r')
% plot(dispResult(1:point),NFoot(1:point),'g')
% hold off

% % % % figure;
% % % % hold on
% % % % plot(dispResult(1:point),thetaFoot(1:point),'b')
% % % % hold off
% % % % 
% % % % figure;
% % % % hold on
% % % % plot(dispResult(1:point),EPunch(1:point),'b')
% % % % xlabel('Total displacement [m]','FontSize',solveOption.fontsizeLabel)
% % % % ylabel('Punching energy [MJ]','FontSize',solveOption.fontsizeLabel)
% % % % hold off
% % % % 
% % % % figure;
% % % % hold on
% % % % plot(dispResult(1:point),vshipXAxis(1:point),'b')
% % % % plot(dispResult(1:point),vshipYAxis(1:point),'r')
% % % % plot(dispResult(1:point),vshipTotAxis(1:point),'k')
% % % % xlabel('Total displacement [m]','FontSize',solveOption.fontsizeLabel)
% % % % ylabel('v_{ship} [m/s]','FontSize',solveOption.fontsizeLabel)
% % % % set(gca,'FontSize',solveOption.fontsizeAxis)
% % % % hold off

% % % % timeSave = timeVec(1:point)';
% % % % ESave = Etot(1:point)';
% % % % forceSave = forceResultFinal';
% % % % dispSave = dispResult(1:point)';
% % % % save('ResultsJacket_A0_0001s.mat','timeSave','ESave','forceSave','dispSave')
% 
% dispSave = dispResult(1:point)';
% EPunchSave = EPunch(1:point)';
% save('PunchData45deg.mat','dispSave','EPunchSave')

if properties.defLoc == 0
    cylinder = 1;
    ship = 1;
    impact = 1;
end

% movie2avi(M,'structMovie.avi','fps',4)