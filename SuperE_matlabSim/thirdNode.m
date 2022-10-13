function[properties,propertiesInit,dofs,forceElem,dispTot,forceTot,plasticity,punchProp,punchPropLeft,punchPropRight] = thirdNode(properties,propertiesInit,data,dofs,forceElem,dispTot,forceTot,plasticity,punchProp,punchPropLeft,punchPropRight,elemHinge,solveOption)

for elem = 1:length(elemHinge)
    i = elemHinge(elem) - elem+1;
    
    % new node - initial structure
    
    L = propertiesInit.L(i);
    vecDir = propertiesInit.angle(1,:,i);
    
    xNew = propertiesInit.nodes(properties.nodein(i),1) + vecDir(1)*L/2;
    yNew = propertiesInit.nodes(properties.nodein(i),2) + vecDir(2)*L/2;
    zNew = propertiesInit.nodes(properties.nodein(i),3) + vecDir(3)*L/2;
    
    propertiesInit.nodes = [propertiesInit.nodes(:,:) ; xNew yNew zNew dofs.nbDofs+1 dofs.nbDofs+2 dofs.nbDofs+3 dofs.nbDofs+4 dofs.nbDofs+5 dofs.nbDofs+6];
    
    propertiesInit.type = [propertiesInit.type(1:i-1) ; propertiesInit.type(i+1:end) ; propertiesInit.type(i) ; propertiesInit.type(i)];
    propertiesInit.nodein = [propertiesInit.nodein(1:i-1) ; propertiesInit.nodein(i+1:end) ; propertiesInit.nodein(i) ; size(propertiesInit.nodes,1)];
    propertiesInit.nodeout = [propertiesInit.nodeout(1:i-1) ; propertiesInit.nodeout(i+1:end) ; size(propertiesInit.nodes,1) ; properties.nodeout(i)];
    
    [propertiesInit,dofs] = autoCalculation(data,propertiesInit,dofs,solveOption,0);
    if data.impactedZone == 1
        properties.impactedElemModel = properties.impactedElemModel - 1;
    end
    
    % new node - deformed structure
    
    T = zeros(12,12);
    T(1:3,1:3) = properties.rotation(:,:,i);
    T(4:6,4:6) = properties.rotation(:,:,i);
    T(7:9,7:9) = properties.rotation(:,:,i);
    T(10:12,10:12) = properties.rotation(:,:,i);
    
%     dispTot(properties.relatedDofs(i,:))
    dispLoc = T*dispTot(properties.relatedDofs(i,:));
    
    xDispLoc = (dispLoc(1) + dispLoc(7))/2;
    
    midSpan = properties.L(i)/2;
    h2 = 2*midSpan^3/(properties.L(i)^3) - 3*midSpan^2/(properties.L(i)^2) + 1;
    h3 = -(properties.L(i)*(midSpan^3/(properties.L(i)^3) - 2*midSpan^2/properties.L(i)^2 + midSpan/properties.L(i)));
    h5 = -2*midSpan^3/(properties.L(i)^3) + 3*midSpan^2/(properties.L(i)^2);
    h6 = -(properties.L(i)*(midSpan^3/(properties.L(i)^3) - midSpan^2/(properties.L(i)^2)));
    
    yDispLoc = (dispLoc(2) + dispLoc(8))/2 + h2*dispLoc(2) + h3*dispLoc(6) + h5*dispLoc(8) + h6*dispLoc(12);
    zDispLoc = (dispLoc(3) + dispLoc(9))/2 + h2*dispLoc(3) + h3*dispLoc(5) + h5*dispLoc(9) + h6*dispLoc(11);
    
    dispLocThirdNode = [xDispLoc ; yDispLoc ; zDispLoc];
    dispGlobThirdNode = properties.rotation(:,:,i).'*dispLocThirdNode;
    
    coordGlobThirdNode = [xNew+dispGlobThirdNode(1) ; yNew+dispGlobThirdNode(2) ; zNew+dispGlobThirdNode(3)];
    
    properties.nodes = [properties.nodes(:,:) ; coordGlobThirdNode(1) coordGlobThirdNode(2) coordGlobThirdNode(3) propertiesInit.nodes(end,4:9)];
    
    properties.type = propertiesInit.type;
    properties.nodein = propertiesInit.nodein;
    properties.nodeout = propertiesInit.nodeout;
    
    [properties,dofs] = autoCalculation(data,properties,dofs,solveOption,0);
    dofs.freeDofs = freeDofs(dofs);
    
    % Total displacement
    
    dispTot = [dispTot(:) ; dispGlobThirdNode(1) ; dispGlobThirdNode(2) ; dispGlobThirdNode(3) ; 0 ; 0 ; 0];
    
    % Efforts
    
%     forceElem(:,i)
    
    N = forceElem(1,i);
    Vy = forceElem(2,i);
    Vz = forceElem(3,i);
    T = forceElem(4,i);
    yDefl = h3*dispLoc(6) + h6*dispLoc(12);
    zDefl = h3*dispLoc(5) + h6*dispLoc(11);
    Myl = forceElem(6,i);
    Mzl = forceElem(5,i);
    Myr = forceElem(12,i);
    Mzr = forceElem(11,i);
    My = (Myl - Myr)/2 - abs(N*yDefl)*sign(yDefl);
    Mz1 = (Mzl - Mzr)/2;
    Mz2 = abs(N*zDefl)*sign(zDefl);
    Mz = (Mzl - Mzr)/2 - abs(N*zDefl)*sign(zDefl);
    
    forceElem = [forceElem(:,1:i-1) forceElem(:,i+1:end) zeros(12,2)];
    forceElem(:,end-1) = [N ; Vy ; Vz ; T ; Mzl ; Myl ; -N ; -Vy ; -Vz ; -T ; -Mz ; -My];
    forceElem(:,end) = [N ; Vy ; Vz ; T ; Mz ; My ; -N ; -Vy ; -Vz ; -T ; Mzr ; Myr];
    
    Npl = properties.Npl(end);
    Mpl = properties.Mpl(end);
    crit1 = (N/Npl)^2;
    crit2 = (My/Mpl)^2;
    crit3 = (Mz/Mpl)^2;
    crit = (N/Npl)^2 + (My/Mpl)^2 + (Mz/Mpl)^2;
    
%     forceElem(:,end-1)
%     forceElem(:,end)
    
    plasticityElem = plasticity(:,i);
    plasticity = [plasticity(:,1:i-1) plasticity(:,i+1:end) zeros(3,2)];
    plasticity(:,end-1) = [-1 ; plasticityElem(2) ; 1];
    plasticity(:,end) = [-1 ; 1 ; plasticityElem(3)];
    
    % Force
    
    forceTot = [forceTot(:) ; zeros(6,1)];
    
    % Nodes for puching
    
    for j = 1:size(punchPropLeft.leg,1)
        if punchPropLeft.leg(j,2) == i
            punchPropLeft.leg(j,2) = size(properties.nodein,1);
        elseif punchPropLeft.leg(j,2) > i
            punchPropLeft.leg(j,2) = punchPropLeft.leg(j,2) - 1;
        end
        
        if punchPropLeft.leg(j,4) == i
            punchPropLeft.leg(j,4) = size(properties.nodein,1);
        elseif punchPropLeft.leg(j,4) > i
            punchPropLeft.leg(j,4) = punchPropLeft.leg(j,4) - 1;
        end
    end
    
    for j = 1:size(punchPropRight.leg,1)
        if punchPropRight.leg(j,2) == i
            punchPropRight.leg(j,2) = size(properties.nodein,1);
        elseif punchPropRight.leg(j,2) > i
            punchPropRight.leg(j,2) = punchPropRight.leg(j,2) - 1;
        end
        
        if punchPropRight.leg(j,4) == i
            punchPropRight.leg(j,4) = size(properties.nodein,1);
        elseif punchPropRight.leg(j,4) > i
            punchPropRight.leg(j,4) = punchPropRight.leg(j,4) - 1;
        end
    end
    
    punchProp.redSect = [punchProp.redSect ;
        1 1 1 1];
    punchProp.redSectStar = [punchProp.redSectStar ;
        1 1 1 1];
        
end