function[cylinder,ship,param,impact,cylinderVec,shipVec,impactVec] = impactPointTim(data,properties,propertiesNoDiv,ship,idElem,nbListElem,contactElement,contactNode,sequenceCrushing,contactPoint)

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

if sequenceCrushing(contactPoint,2) == 0
    
    cylinder.D = properties.De(2);
    cylinder.t = properties.t(2);
    Deltax = properties.end(2,1) - properties.origin(2,1);
    Deltay = properties.end(2,2) - properties.origin(2,2);
    Deltaz = properties.end(2,3) - properties.origin(2,3);
    cylinder.dzeta = atan(abs(Deltaz)/sqrt(Deltax^2+Deltay^2)) * 180/pi; % degrees
    cylinder.L = 100 * cylinder.D; % m
    cylinder.fy = properties.fy(2); % MPa
    cylinder.rearImpact = 0;
    
else
    
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
    cylinder.L = propertiesNoDiv.L(idElem) - 2; % à réfléchir le "-2"
    cylinder.fy = properties.fy(idElem); % N/m²
    cylinder.rearImpact = 0; % 0 -> no ; 1 -> yes
    
end

% Mechanical properties

cylinder.R = (cylinder.D - cylinder.t)/2;
cylinder.dzeta = cylinder.dzeta*pi/180;

A = pi/4*(cylinder.D^2 - (cylinder.D-2*cylinder.t)^2);
cylinder.Np = A*cylinder.fy;
cylinder.Mp = 4*cylinder.fy * ((cylinder.D/2)^3 - ((cylinder.D-2*cylinder.t)/2)^3)/3;
cylinder.n0 = cylinder.fy*cylinder.t;
cylinder.m0 = cylinder.fy*cylinder.t^2/4;

ship.phib = ship.phib;
ship.psib = ship.psib;
ship.alpha = ship.alpha*pi/180;

%% Impact point

if sequenceCrushing(contactPoint,2) == 0 % impact on node
    disp('node')
    if sequenceCrushing(contactPoint,1) == 0
        
        impact.xsI = contactElement.distBow(nbListElem,3);
        impact.ysI = contactElement.distBow(nbListElem,4);
        impact.zsI = contactElement.distBow(nbListElem,5);
        
        impact.L1 = cylinder.L / 2;
        impact.L2 = cylinder.L - impact.L1;
        
    else
        
        impact.xsI = contactElement.distBulb(nbListElem,3);
        impact.ysI = contactElement.distBulb(nbListElem,4);
        impact.zsI = contactElement.distBulb(nbListElem,5);
        
        impact.L1 = cylinder.L / 2;
        impact.L2 = cylinder.L - impact.L1;
        
    end
    
else % impact on element
    %%%%%disp('element')
    if sequenceCrushing(contactPoint,1) == 0
        
        impact.xsI = contactElement.distBow(nbListElem,3);
        impact.ysI = contactElement.distBow(nbListElem,4);
        impact.zsI = contactElement.distBow(nbListElem,5);
        
        abc = contactElement.distBow(nbListElem,2);
        bcd = propertiesNoDiv.origin(idElem,3);
        cde = propertiesNoDiv.end(idElem,3);
        
%         reportHeight = abs((contactElement.distBow(nbListElem,2) - properties.origin(idElem,3) - 1*sin(cylinder.dzeta)) / (properties.end(idElem,3) - properties.origin(idElem,3) - 2*sin(cylinder.dzeta)));
        reportHeight = abs((contactElement.distBow(nbListElem,2) - propertiesNoDiv.origin(idElem,3)) / (propertiesNoDiv.end(idElem,3) - propertiesNoDiv.origin(idElem,3)));
        impact.L1 = cylinder.L * reportHeight;
        impact.L2 = cylinder.L - impact.L1;
        
    else
        
        impact.xsI = contactElement.distBulb(nbListElem,3);
        impact.ysI = contactElement.distBulb(nbListElem,4);
        impact.zsI = contactElement.distBulb(nbListElem,5);
        
        abc2 = contactElement.distBulb(nbListElem,2);
        bcd2 = propertiesNoDiv.origin(idElem,3);
        cde2 = propertiesNoDiv.end(idElem,3);
        idElem;
        def2 = propertiesNoDiv.origin;
        
%         reportHeight = abs((contactElement.distBulb(nbListElem,2) - properties.origin(idElem,3) - 1*sin(cylinder.dzeta)) / (properties.end(idElem,3) - properties.origin(idElem,3) - 2*sin(cylinder.dzeta)));
        reportHeight = abs((contactElement.distBulb(nbListElem,2) - propertiesNoDiv.origin(idElem,3)) / (propertiesNoDiv.end(idElem,3) - propertiesNoDiv.origin(idElem,3)));
        impact.L1 = cylinder.L * reportHeight;
        impact.L2 = cylinder.L - impact.L1;
        
    end
    
end

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

%% Vector

% Cylinder
cylinderVec = zeros(1,12);
cylinderVec(1) = cylinder.D;
cylinderVec(2) = cylinder.t;
cylinderVec(3) = cylinder.dzeta;
cylinderVec(4) = cylinder.L;
cylinderVec(5) = cylinder.fy;
cylinderVec(6) = cylinder.rearImpact;
cylinderVec(7) = cylinder.R;
cylinderVec(8) = cylinder.Np;
cylinderVec(9) = cylinder.Mp;
cylinderVec(10) = cylinder.n0;
cylinderVec(11) = cylinder.m0;
cylinderVec(12) = deltaMax;

% Ship
if sequenceCrushing(contactPoint,1) == 1
    if sequenceCrushing(contactPoint,2) == 0
        zBulb = contactNode.distBulb(nbListElem,5);
    else
        zBulb = contactElement.distBulb(nbListElem,5);
    end
    
    xBulb = sqrt(1-(zBulb/(ship.hBulb/2))^2) * ship.pBulb;
    if zBulb < 0.05*ship.hBulb
        ship.phib = 90;
    else
        ship.phib = atan(xBulb/zBulb);
    end
    
    yBulb = sqrt(1-(zBulb/(ship.hBulb/2))^2) * ship.bBulb/2;
    if zBulb < 0.05*ship.hBulb
        ship.psib = 90;
    else
        ship.psib = atan(yBulb/zBulb);
    end
end
shipVec = zeros(1,12);
shipVec(1) = ship.bulbous;
shipVec(2) = ship.p;
shipVec(3) = ship.q;
shipVec(4) = ship.phib;
shipVec(5) = ship.psib;
shipVec(6) = ship.hb;
shipVec(7) = ship.hBulb;
shipVec(8) = ship.bBulb;
shipVec(9) = ship.pBulb;
shipVec(10) = ship.shiftCentreBulb;
shipVec(11) = ship.alpha;
shipVec(12) = ship.phibprime;

% Impact
impact = parameters(cylinder,ship,param,impact);
impactVec = zeros(1,33);
impactVec(1) = impact.xsI;
impactVec(2) = impact.ysI;
impactVec(3) = impact.zsI;
impactVec(4) = impact.L1;
impactVec(5) = impact.L2;
impactVec(6) = impact.YI;
impactVec(7) = impact.ZI;
impactVec(8) = impact.pv;
impactVec(9) = impact.qv;
impactVec(10) = impact.hbv;
impactVec(11) = impact.xvI;
impactVec(12) = impact.yvI;
impactVec(13) = impact.zvI;
impactVec(14) = impact.xvS;
impactVec(15) = impact.yvS;
impactVec(16) = impact.zvS;
impactVec(17) = impact.betav;
impactVec(18) = impact.df;
impactVec(19) = impact.nv;
impactVec(20) = impact.psiv0;
impactVec(21) = impact.psivf;
impactVec(22) = impact.eta1;
impactVec(23) = impact.eta2;
impactVec(24) = impact.xhI;
impactVec(25) = impact.yhI;
impactVec(26) = impact.zhI;
impactVec(27) = impact.alphah;
impactVec(28) = impact.xhS;
impactVec(29) = impact.yhS;
impactVec(30) = impact.zhS;
impactVec(31) = impact.betah;
impactVec(32) = impact.nh;
impactVec(33) = impact.psih0;