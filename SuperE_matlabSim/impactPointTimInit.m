function[cylinder,ship,param,impact] = impactPointTimInit(data,properties,ship,idElem,nbListElem,contactElement)

%% Parameters

param.k1 = 1;
param.k2 = 1;
param.Tol = 10^(-9);
param.ksi = 1;

%% Data of impact

deltax = properties.end(idElem,1) - properties.origin(idElem,1);
deltay = properties.end(idElem,2) - properties.origin(idElem,2);
if deltax == 0
    ship.alpha = 90 - data.shipTrajectory;
else
    ship.alpha = atan(deltay/deltax) * 180/pi - data.shipTrajectory; % degrees
end

%% Properties of element

% Data

cylinder.D = properties.De(idElem); % m
cylinder.t = properties.t(idElem); % m
deltax = properties.end(idElem,1) - properties.origin(idElem,1);
deltay = properties.end(idElem,2) - properties.origin(idElem,2);
deltaz = properties.end(idElem,3) - properties.origin(idElem,3);
if sqrt(deltax^2+deltay^2) == 0
    cylinder.dzeta = 90;
else
    cylinder.dzeta = atan(abs(deltaz)/sqrt(deltax^2+deltay^2)) * 180/pi; % degrees
end
cylinder.L = properties.L(idElem) - 2; % à réfléchir le "-2"
cylinder.fy = properties.fy(idElem); % N/m²
cylinder.rearImpact = 0; % 0 -> no ; 1 -> yes

% Mechanical properties

cylinder.R = (cylinder.D - cylinder.t)/2;
cylinder.dzeta = cylinder.dzeta*pi/180;

A = pi/4*(cylinder.D^2 - (cylinder.D-2*cylinder.t)^2);
cylinder.Np = A*cylinder.fy;
cylinder.Mp = 4*cylinder.fy * ((cylinder.D/2)^3 - ((cylinder.D-2*cylinder.t)/2)^3)/3;
cylinder.n0 = cylinder.fy*cylinder.t;
cylinder.m0 = cylinder.fy*cylinder.t^2/4;

ship.phib = ship.phib*pi/180;
ship.psib = ship.psib*pi/180;
ship.alpha = ship.alpha*pi/180;

%% Impact point

impact.xsI = contactElement.distBow(nbListElem,3);
impact.ysI = contactElement.distBow(nbListElem,4);
impact.zsI = contactElement.distBow(nbListElem,5);

reportHeight = abs((contactElement.distBow(nbListElem,2) - properties.origin(idElem,3) - 1*sin(cylinder.dzeta)) / (properties.end(idElem,3) - properties.origin(idElem,3) - 2*sin(cylinder.dzeta)));
impact.L1 = cylinder.L * reportHeight;
impact.L2 = cylinder.L - impact.L1;

impact.YI = impact.L1 * cos(cylinder.dzeta);
impact.ZI = impact.L1 * sin(cylinder.dzeta);

%% Maximum crushing if impact on node

alpha = data.shipTrajectory * pi/180;
Rl = data.De(1)/2;
Rb = data.De(4)/2;
if alpha == 0
    deltaMax = Rl + sqrt(Rl^2 - Rb^2);
else
    deltaMax = abs(sqrt(Rl^2 - Rb^2)/(tan(alpha)) - Rb + Rl/(sin(alpha))) * sin(alpha);
end
cylinder.maxLocalCrushingNode = deltaMax;