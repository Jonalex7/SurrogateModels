function[data,properties,dofs,force,punchPropLeft,punchPropRight] = impactProperties(data,properties,dofs,punchPropLeft,punchPropRight,solveOption)

if data.impactedZone == 0
    force.dofsForce  = [data.impactedNode*6-5 data.impactedNode*6-4];
else
    zShip = data.shipHeight;
    
    x0 = properties.nodes(properties.nodein(data.impactedElement),1);
    y0 = properties.nodes(properties.nodein(data.impactedElement),2);
    z0 = properties.nodes(properties.nodein(data.impactedElement),3);
    
    vecDir = properties.angle(1,:,data.impactedElement);
    
    t = (zShip-z0)/vecDir(3);
    
    x = vecDir(1)*t + x0;
    y = vecDir(2)*t + y0;
    
    properties.nodes = [properties.nodes(:,:) ; x y zShip dofs.nbDofs+1 dofs.nbDofs+2 dofs.nbDofs+3 dofs.nbDofs+4 dofs.nbDofs+5 dofs.nbDofs+6];
    
    properties.type = [properties.type(1:data.impactedElement-1) ; properties.type(data.impactedElement+1:end) ; properties.type(data.impactedElement) ; properties.type(data.impactedElement)];
    properties.nodein = [properties.nodein(1:data.impactedElement-1) ; properties.nodein(data.impactedElement+1:end) ; properties.nodein(data.impactedElement) ; size(properties.nodes,1)];
    properties.nodeout = [properties.nodeout(1:data.impactedElement-1) ; properties.nodeout(data.impactedElement+1:end) ; size(properties.nodes,1) ; properties.nodeout(data.impactedElement)];
    
    [properties,dofs] = autoCalculation(data,properties,dofs,solveOption,0);
    dofs.freeDofs = freeDofs(dofs);
    
    force.dofsForce  = [dofs.nbDofs-5 dofs.nbDofs-4];     % Vector of numbers of Dofs where force is applied
    properties.impactedElemModel = [length(properties.type)-1 length(properties.type)];
    
    for i = 1:size(punchPropLeft.leg,1)
        if punchPropLeft.leg(i,2) > data.impactedElement
            punchPropLeft.leg(i,2) = punchPropLeft.leg(i,2) - 1;
        end
        if punchPropLeft.leg(i,4) > data.impactedElement
            punchPropLeft.leg(i,4) = punchPropLeft.leg(i,4) - 1;
        end
        
        if punchPropRight.leg(i,2) > data.impactedElement
            punchPropRight.leg(i,2) = punchPropRight.leg(i,2) - 1;
        end
        if punchPropRight.leg(i,4) > data.impactedElement
            punchPropRight.leg(i,4) = punchPropRight.leg(i,4) - 1;
        end
    end
end