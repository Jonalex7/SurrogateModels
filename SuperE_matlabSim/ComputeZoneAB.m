function[M] = ComputeZoneAB(data,Rp,h,theta,delta)

%% Read data

dloc = delta;

alphaP = 0;
psi0 = 1/3*pi;%Jing-Ru
nhp = 1;
k = 1;%Loïc

fy = data.fy;
t = data.tLeg;
m0 = fy*t^2/4;
n0 = fy*t;

Rleg = data.RLeg;
Dbrace = data.DBrace;

Rstar = Rp*cos(theta);

%% Computation - Partie comprimée

psi = psi0 + (k * pi - psi0) * (dloc * cos(alphaP) / 2 / Rleg) ^ nhp;
dpsi = nhp * (k * pi - psi0) * dloc ^ (nhp - 1) * (cos(alphaP) / 2 / Rleg) ^ nhp;

% -----
% Rings
% -----

% Calculate the parameters of section
R1 = Rleg + dloc * cos(alphaP) * (psi - sin(psi)) / (pi * (1 - cos(psi)) - 2 * (psi - sin(psi)));
R2 = Rleg - dloc * cos(alphaP) * (pi - psi + sin(psi)) / (pi * (1 - cos(psi)) - 2 * (psi - sin(psi)));
dR1 = pi * dloc * (2 * (1 - cos(psi)) - psi * sin(psi)) / (pi * (1 - cos(psi)) - 2 * (psi - sin(psi))) * dpsi + psi - sin(psi);
dR1 = cos(alphaP) / (pi * (1 - cos(psi)) - 2 * (psi - sin(psi))) * dR1;
dR2 = pi * dloc * (pi - psi) * sin(psi) / (pi * (1 - cos(psi)) - 2 * (psi - sin(psi))) * dpsi - (pi - psi + sin(psi));
dR2 = cos(alphaP) / (pi * (1 - cos(psi)) - 2 * (psi - sin(psi))) * dR2;

% Calculate the speed
VB = (R1 - R2) * dpsi - (pi - psi) * dR1 - psi * dR2;
VF = R1 * dpsi - (pi - psi) * dR1;

% Calculate the resistance
Eb = 2 * m0 * (abs(VB) / R2 + abs(VF) * (1 / R2 - 1 / R1) + psi / R2 * abs(dR2) + (pi - psi) / R1 * abs(dR1));

% ----------
% Generators
% ----------

jmax = 100;
dEm = zeros(1,jmax+1);

for j = 0:jmax
    beta = j*pi/jmax;
    
    % Zone 1
    if Rleg * beta <= 1 / 3 * Dbrace
        dEmz = (dloc * cos(alphaP)) * cos(alphaP);
    elseif 1 / 3 * Dbrace < Rleg * beta && Rleg * beta <= (R1 - R2) * sin(psi)
        dEmz = (dloc * cos(alphaP) - Rleg * (1 - cos(beta))) * cos(alphaP);
        
        % Zone 2
    elseif Rleg * beta > (R1 - R2) * sin(psi) && Rleg * beta < (R1 - R2) * sin(psi) + R2 * psi
        
        % Calculation of phi
        phi = (Rleg * beta - (R1 - R2) * sin(psi)) / R2;
        dphi = -(R2 * sin(psi) * dR1 + (R1 - R2) * R2 * cos(psi) * dpsi + (Rleg * beta - R1 * sin(psi)) * dR2) / R2 ^ 2;
        
        % Calculation of dEm
        dEmz = (2 * R1 * (1 + cos(psi)) - R2 * (1 + cos(psi) - cos(phi) - cos(psi - phi)) - Rleg * (1 + cos(beta) + cos(psi) + cos(psi - beta))) * dR1;
        dEmz = dEmz + (2 * R2 * (1 - cos(psi - phi)) - R1 * (1 + cos(psi) - cos(phi) - cos(psi - phi)) - Rleg * (cos(phi - beta) + cos(phi) - cos(psi - beta) - cos(psi))) * dR2;
        dEmz = dEmz + (-R1 ^ 2 * sin(psi) + R2 ^ 2 * sin(psi - phi) - R1 * R2 * (sin(psi - phi) - sin(psi)) + (R1 - R2) * Rleg * (sin(psi) + sin(psi - beta))) * dpsi;
        dEmz = dEmz - R2 * (R2 * sin(psi - phi) - R1 * (sin(psi - phi) - sin(phi)) - Rleg * (sin(phi) + sin(phi - beta))) * dphi;
        
        % Zone 3
    elseif Rleg * beta >= (R1 - R2) * sin(psi) + R2 * psi
        
        % Calculation of phi
        phi = (Rleg * beta - (R1 - R2) * sin(psi) - R2 * psi) / R1;
        dphi = ((R2 * (psi - sin(psi)) - Rleg * beta) * dR1 + R1 * (sin(psi) - psi) * dR2 - R1 * ((R1 - R2) * cos(psi) + R2) * dpsi) / R1 ^ 2;
        
        % Calculation of dEm
        dEmz = 2 * R1 * dR1 * (1 + cos(phi + psi)) - R1 ^ 2 * sin(psi + phi) * (dpsi + dphi);
        
    end
    
    dEm(j+1) = dEmz;
end

dbeta = pi/jmax;
Emsum = 0;
for j = 1:jmax
    Emsum = Emsum + n0*Rleg*(dEm(j)+dEm(j+1))/2*dbeta;
end

Em = Emsum;

ksi = min((2 * Em / Eb) ^ 0.5,h);
FComp = (Em * (1 / ksi + 1 / Rstar) + Eb * (ksi + Rstar) / 2);

%% Computation - partie tendue

grandRayon = Rleg + delta/2;
petitRayon = sqrt(Rleg^2-Rleg*delta-delta^2/4);
exc = sqrt(grandRayon^2-petitRayon^2)/grandRayon;

betamax = pi;
nbDbeta = 100;
dbeta = betamax/nbDbeta;

% Anneaux

beta = 0;
xidotdsSegment = zeros(1,nbDbeta+1);
for j = 1:nbDbeta
    if exc == 0
        dxiddeltab = -(grandRayon^2*(Rleg+delta/2) + grandRayon*petitRayon^2)/(grandRayon^4) * (1 - exc^2*(cos(beta))^2)^(-3/2);
    else
        dxiddeltab = -(grandRayon^2*(Rleg+delta/2) + grandRayon*petitRayon^2)/(grandRayon^4) * (1 - exc^2*(cos(beta))^2)^(-3/2)...
            + 3/2*(petitRayon/grandRayon)^2 * (1 - exc^2*(cos(beta))^2)^(-5/2) * 2*exc*(cos(beta))^2 * ((grandRayon+Rleg+delta/2)/(2*sqrt(grandRayon^2-petitRayon^2))*grandRayon - sqrt(grandRayon^2-petitRayon^2)/2)/(grandRayon^2);
    end
    ds = sqrt((grandRayon*cos(beta))^2 + (petitRayon*sin(beta))^2);
    xidotdsSegment(j) = abs(dxiddeltab) * ds;
    beta = beta+dbeta;
end
Edotb = 0;
for j = 1:nbDbeta
    Edotb = Edotb + 2*m0 * (xidotdsSegment(j)+xidotdsSegment(j+1))*dbeta/2;
end

% Génératrices

beta = 0;
wdwddeltaSegment = zeros(1,nbDbeta+1);
for j = 1:nbDbeta
    if exc == 0
        delta = Rleg*10^(-2);
    end
    xA = Rleg*cos(beta);
    zA = Rleg*sin(beta);
    xB = delta/2 + grandRayon*cos(beta);
    zB = petitRayon*sin(beta);
    w = sqrt((xA-xB)^2 + (zA-zB)^2);
    dwddelta = 1/(2*w)*((xB-xA)*(1+cos(beta)) - (zB-zA)*(Rleg-delta/2)/petitRayon * sin(beta));
    wdwddeltaSegment(j) = w * dwddelta;
    beta = beta+dbeta;
end
Edotm = 0;
for j = 1:nbDbeta
    Edotm = Edotm + 2*Rleg*n0 * 2*(wdwddeltaSegment(j)+wdwddeltaSegment(j+1))*dbeta/2;
end

% abc = Edotm
% def = Edotb

ksi = min((2 * Edotm / Edotb) ^ 0.5,h);
FTrac = (Em * (1 / ksi + 1 / Rstar) + Eb * (ksi + Rstar) / 2);

%% Total

M = (FComp + FTrac) * Rp*cos(theta);