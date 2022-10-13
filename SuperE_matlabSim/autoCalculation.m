function[properties,dofs] = autoCalculation(data,properties,dofs,solveOption,drawInit)

%% Data - auto-calculation

dofs.nbElements = length(properties.nodein); % Number of elements
dofs.nbNodes    = size(properties.nodes,1);  % Number of nodes
dofs.nbDofs     = dofs.nbNodes * 6;          % Number of Dofs

properties.origin      = zeros(dofs.nbElements,3); % Coordinates of origin of elements
properties.end         = zeros(dofs.nbElements,3); % Coordinates of end of element
properties.relatedDofs = zeros(dofs.nbElements,12); % Dofs related to the elemnt
for i = 1:dofs.nbElements
    properties.origin(i,:)      = properties.nodes(properties.nodein(i),1:3);
    properties.end(i,:)         = properties.nodes(properties.nodeout(i),1:3);
    properties.relatedDofs(i,:) = [properties.nodes(properties.nodein(i),4:9) properties.nodes(properties.nodeout(i),4:9)];
end

X = [1 0 0];
Y = [0 1 0];
Z = [0 0 1];

properties.L     = zeros(dofs.nbElements,1); % Length
properties.angle = zeros(3,3,dofs.nbElements); % vector of local axis
properties.rotation = zeros(3,3,dofs.nbElements);
for i = 1:dofs.nbElements
    Deltax = properties.end(i,1) - properties.origin(i,1);
    Deltay = properties.end(i,2) - properties.origin(i,2);
    Deltaz = properties.end(i,3) - properties.origin(i,3);
    properties.L(i) = sqrt(Deltax^2 + Deltay^2 + Deltaz^2);
    
    properties.angle(1,:,i) = [Deltax Deltay Deltaz];
    properties.angle(1,:,i) = properties.angle(1,:,i) / norm(properties.angle(1,:,i));
    if Deltaz ~= 0
        properties.angle(2,:,i) = [1 1 -(Deltax+Deltay)/(Deltaz)];
    elseif Deltay ~= 0
        properties.angle(2,:,i) = [1 -(Deltax+Deltaz)/(Deltay) 1];
    else
        properties.angle(2,:,i) = [-(Deltay+Deltaz)/(Deltax) 1 1];
    end
    properties.angle(2,:,i) = properties.angle(2,:,i) / norm(properties.angle(2,:,i));
    properties.angle(3,:,i) = cross(properties.angle(1,:,i),properties.angle(2,:,i));
    
    properties.rotation(:,:,i) = [X*properties.angle(1,:,i)' Y*properties.angle(1,:,i)' Z*properties.angle(1,:,i)' ;
        X*properties.angle(2,:,i)' Y*properties.angle(2,:,i)' Z*properties.angle(2,:,i)' ;
        X*properties.angle(3,:,i)' Y*properties.angle(3,:,i)' Z*properties.angle(3,:,i)'];
end

properties.De  = zeros(dofs.nbElements,1); % External diameter - [m]
properties.t   = zeros(dofs.nbElements,1); % Thickness - [m]
properties.fy  = zeros(dofs.nbElements,1);    % Yield limit - [N/m^2]
properties.E   = zeros(dofs.nbElements,1); % Young's modulus - [N/m^2]
properties.rho = zeros(dofs.nbElements,2); % Fixity factors
for i = 1:dofs.nbElements
    properties.De(i)    = data.De(properties.type(i));
    properties.t(i)     = data.t(properties.type(i));
    properties.fy(i)    = data.fy(properties.type(i));
    properties.E(i)     = data.E(properties.type(i));
    properties.rho(i,:) = data.rho(properties.type(i),:);
    
    if data.impactedZone == 1
        if i == properties.impactedElemModel(1) || i == properties.impactedElemModel(2)
            properties.fy(i) = properties.fy(i) * 10^3;
            properties.E(i)  = properties.E(i)  * 10^3;
        end
    end
end

properties.A   = zeros(dofs.nbElements,1); % Area
properties.I   = zeros(dofs.nbElements,1); % Inertia
properties.J   = zeros(dofs.nbElements,1); % Inertia of torsion
properties.Npl = zeros(dofs.nbElements,1); % Plastic axial resistance
properties.Mpl = zeros(dofs.nbElements,1); % Plastic bending moment
properties.Vpl = zeros(dofs.nbElements,1); % Plastic shear resistance
properties.Pcr = zeros(dofs.nbElements,1); % Critical axial load
for i = 1:dofs.nbElements
    Di = properties.De(i) - 2*properties.t(i);
    properties.A(i)   = pi/4*(properties.De(i)^2 - Di^2);
    properties.I(i)   = pi/64*(properties.De(i)^4 - Di^4);
    properties.J(i)   = pi/32*(properties.De(i)^4 - Di^4);
    properties.Npl(i) = properties.A(i) * properties.fy(i);
    properties.Mpl(i) = 4*properties.fy(i) * ((properties.De(i)/2)^3 - (Di/2)^3)/3;
    properties.Vpl(i) = 0.6*properties.A(i)*properties.fy(i)/sqrt(3);
    K = (1 + 0.145*((1-properties.rho(i,1)) + (1-properties.rho(i,2))) - 0.265*(1-properties.rho(i,1))*(1-properties.rho(i,2)))...
        / (2 - 0.364*((1-properties.rho(i,1)) + (1-properties.rho(i,2))) - 0.247*(1-properties.rho(i,1))*(1-properties.rho(i,2)));
    properties.Pcr(i) = pi^2 * properties.E(i) * properties.I(i) / ((K*properties.L(i))^2);
end

% properties.defoInit = properties.L.*properties.sensDefoInit;

%% Draw

if drawInit == 1
    figure;
    hold on
    for i = 1:dofs.nbElements
        plot3([properties.origin(i,1) properties.end(i,1)],[properties.origin(i,2) properties.end(i,2)],[properties.origin(i,3) properties.end(i,3)]);
    end
    view(solveOption.az,solveOption.el)
    axis equal
    axis([min(properties.nodes(:,1))-1 max(properties.nodes(:,1))+1 min(properties.nodes(:,2))-1 max(properties.nodes(:,2))+1 min(properties.nodes(:,3))-1 max(properties.nodes(:,3))+1])
    hold off
end