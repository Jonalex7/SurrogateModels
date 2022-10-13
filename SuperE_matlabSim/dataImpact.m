function[cylinder,ship,param] = dataImpact(properties,data,ship)

% Parameters

param.k1 = 1;
param.k2 = 1;
param.Tol = 10^(-9);
param.ksi = 1;

% Cylinder

if data.impactedZone == 0
    
    % Cylinder
    cylinder.D = properties.De(2);
    cylinder.t = properties.t(2);
    Deltax = properties.end(2,1) - properties.origin(2,1);
    Deltay = properties.end(2,2) - properties.origin(2,2);
    Deltaz = properties.end(2,3) - properties.origin(2,3);
    cylinder.dzeta = atan(abs(Deltaz)/sqrt(Deltax^2+Deltay^2)) * 180/pi; % degrees
    cylinder.L = 100 * cylinder.D; % m
    cylinder.fy = properties.fy(2); % MPa
    
    % Impact location
    ship.XP = 0;
    ship.YP = 50*cylinder.D * cos(cylinder.dzeta*pi/180) / tan(cylinder.dzeta*pi/180);
    ship.ZS = 50*cylinder.D * cos(cylinder.dzeta*pi/180);
    ship.alpha = abs(45 - data.shipTrajectory);
    
elseif data.impactedZone == 1
    
    % Cylinder
    cylinder.D = properties.De(data.impactedElement); % m
    cylinder.t = properties.t(data.impactedElement); % m
    Deltax = properties.end(data.impactedElement,1) - properties.origin(data.impactedElement,1);
    Deltay = properties.end(data.impactedElement,2) - properties.origin(data.impactedElement,2);
    Deltaz = properties.end(data.impactedElement,3) - properties.origin(data.impactedElement,3);
    if sqrt(Deltax^2+Deltay^2) == 0
        cylinder.dzeta = 90;
    else
        cylinder.dzeta = atan(abs(Deltaz)/sqrt(Deltax^2+Deltay^2)) * 180/pi; % degrees
    end
    cylinder.L = properties.L(data.impactedElement) - 2; % m
    cylinder.fy = properties.fy(data.impactedElement); % MPa
    
    % Impact location
    ship.XP = 0; % m
    ship.YP = abs((data.shipHeight - properties.origin(data.impactedElement,3) - 1)/tan(cylinder.dzeta*pi/180)); % m
    ship.ZS = abs(data.shipHeight - properties.origin(data.impactedElement,3) - 1); % m
    if Deltax == 0
        ship.alpha = 90 - data.shipTrajectory;
    else
        ship.alpha = atan(Deltay/Deltax) * 180/pi - data.shipTrajectory; % degrees
    end
    
end

cylinder.rearImpact = 0; % 0 -> no ; 1 -> yes

%% Change of units - autocalculation

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

%% Correction if dzeta < 0

if cylinder.dzeta < 0
    if abs(ship.alpha) <= pi/4
        cylinder.dzeta = -cylinder.dzeta;
        ship.alpha = -ship.alpha;
        ship.YP = cylinder.L * cos(cylinder.dzeta) - ship.YP;
        ship.ZS = ship.ZS + cylinder.L * sin(cylinder.dzeta);
    end
    if abs(ship.alpha) > pi/4 && abs(ship.alpha) <= pi/2
        cylinder.dzeta = -cylinder.dzeta;
        ship.alpha = -ship.alpha;
        ship.XP = -ship.XP;
        ship.ZS = ship.ZS + cylinder.L * sin(cylinder.dzeta);
    end
    param.k2 = -1;
end

%% Correction if rear impact

if cylinder.rearImpact == 1
    if abs(ship.alpha) <= pi/4
        ship.alpha = -ship.alpha;
    end
    if abs(ship.alpha) > pi/4 && abs(ship.alpha) <= pi/2
        ship.alpha = -ship.alpha;
        ship.XP = -ship.XP;
    end
    param.k1 = -1;
end

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