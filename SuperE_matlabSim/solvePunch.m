function[resultPunch,punchProp] = solvePunch(punchProp,i,delta,ddelta)

dP = punchProp.punchModeNodes(i,2);
CP = punchProp.punchModeNodes(i,3);

%% Data

% Leg
punching.Rleg = punchProp.nodeProp(i,1); % m
punching.Lleg = punchProp.nodeProp(i,2); % m
punching.tleg = punchProp.nodeProp(i,3); % m
punching.alpha = punchProp.nodeProp(i,5); % degrees
punching.ksi = punchProp.nodeProp(i,6); % -

% Leg
punching.Dbrace = punchProp.nodeProp(i,7); % m
punching.L1 = punchProp.nodeProp(i,8); % m
punching.L2 = punchProp.nodeProp(i,9); % m
punching.NB = punchProp.nodeProp(i,10);
punching.gap = punchProp.nodeProp(i,11); % m

punching.fy = punchProp.nodeProp(i,4); % MPa

punching.alphaP = 0;
% punching.psi0 = 3/4*pi;%Loïc
punching.psi0 = 1/3*pi;%Jing-Ru
punching.nhp = 1;
punching.k = 1;%Loïc
% punching.k = 0.5;%Jing-Ru
punching.a = 0.35*punching.Dbrace / cos(punching.alpha);
punching.Qb = 1;

% Properties
punching.D = 2*punching.Rleg + punching.tleg;
punching.Acs = pi/4*(punching.D^2 - (punching.D-2*punching.tleg)^2);
punching.Np = punching.Acs*punching.fy;
punching.Mp = 4*punching.fy * ((punching.D/2)^3 - ((punching.D-2*punching.tleg)/2)^3)/3;
punching.n0 = punching.fy*punching.tleg;
punching.m0 = punching.fy*punching.tleg^2/4;

Qb = punching.Qb;

Epunch = 0;

%% Local mode

% if mode == 0
    
dloc = delta;

psi = punching.psi0 + (punching.k * pi - punching.psi0) * (dloc * cos(punching.alphaP) / 2 / punching.Rleg) ^ punching.nhp;
dpsi = punching.nhp * (punching.k * pi - punching.psi0) * dloc ^ (punching.nhp - 1) * (cos(punching.alphaP) / 2 / punching.Rleg) ^ punching.nhp;

% -----
% Rings
% -----

% Calculate the parameters of section
R1 = punching.Rleg + dloc * cos(punching.alphaP) * (psi - sin(psi)) / (pi * (1 - cos(psi)) - 2 * (psi - sin(psi)));
R2 = punching.Rleg - dloc * cos(punching.alphaP) * (pi - psi + sin(psi)) / (pi * (1 - cos(psi)) - 2 * (psi - sin(psi)));
dR1 = pi * dloc * (2 * (1 - cos(psi)) - psi * sin(psi)) / (pi * (1 - cos(psi)) - 2 * (psi - sin(psi))) * dpsi + psi - sin(psi);
dR1 = cos(punching.alphaP) / (pi * (1 - cos(psi)) - 2 * (psi - sin(psi))) * dR1;
dR2 = pi * dloc * (pi - psi) * sin(psi) / (pi * (1 - cos(psi)) - 2 * (psi - sin(psi))) * dpsi - (pi - psi + sin(psi));
dR2 = cos(punching.alphaP) / (pi * (1 - cos(psi)) - 2 * (psi - sin(psi))) * dR2;

% Calculate the speed
VB = (R1 - R2) * dpsi - (pi - psi) * dR1 - psi * dR2;
VF = R1 * dpsi - (pi - psi) * dR1;

% Calculate the resistance
Eb = 2 * punching.m0 * (abs(VB) / R2 + abs(VF) * (1 / R2 - 1 / R1) + psi / R2 * abs(dR2) + (pi - psi) / R1 * abs(dR1));
Ebt = Eb;
dt = dloc;
if dloc == 0
    Eb0 = Eb;
end

% ----------
% Generators
% ----------

jmax = 100;
dEm = zeros(1,jmax+1);

for j = 0:jmax
    theta = j*pi/jmax;
    
    % Zone 1
    if punching.Rleg * theta <= 1 / 3 * punching.Dbrace
        dEmz = (dloc * cos(punching.alphaP)) * cos(punching.alphaP);
    elseif 1 / 3 * punching.Dbrace < punching.Rleg * theta && punching.Rleg * theta <= (R1 - R2) * sin(psi)
        dEmz = (dloc * cos(punching.alphaP) - punching.Rleg * (1 - cos(theta))) * cos(punching.alphaP);
        
        % Zone 2
    elseif punching.Rleg * theta > (R1 - R2) * sin(psi) && punching.Rleg * theta < (R1 - R2) * sin(psi) + R2 * psi
        
        % Calculation of phi
        phi = (punching.Rleg * theta - (R1 - R2) * sin(psi)) / R2;
        dphi = -(R2 * sin(psi) * dR1 + (R1 - R2) * R2 * cos(psi) * dpsi + (punching.Rleg * theta - R1 * sin(psi)) * dR2) / R2 ^ 2;
        
        % Calculation of dEm
        dEmz = (2 * R1 * (1 + cos(psi)) - R2 * (1 + cos(psi) - cos(phi) - cos(psi - phi)) - punching.Rleg * (1 + cos(theta) + cos(psi) + cos(psi - theta))) * dR1;
        dEmz = dEmz + (2 * R2 * (1 - cos(psi - phi)) - R1 * (1 + cos(psi) - cos(phi) - cos(psi - phi)) - punching.Rleg * (cos(phi - theta) + cos(phi) - cos(psi - theta) - cos(psi))) * dR2;
        dEmz = dEmz + (-R1 ^ 2 * sin(psi) + R2 ^ 2 * sin(psi - phi) - R1 * R2 * (sin(psi - phi) - sin(psi)) + (R1 - R2) * punching.Rleg * (sin(psi) + sin(psi - theta))) * dpsi;
        dEmz = dEmz - R2 * (R2 * sin(psi - phi) - R1 * (sin(psi - phi) - sin(phi)) - punching.Rleg * (sin(phi) + sin(phi - theta))) * dphi;
        
        % Zone 3
    elseif punching.Rleg * theta >= (R1 - R2) * sin(psi) + R2 * psi
        
        % Calculation of phi
        phi = (punching.Rleg * theta - (R1 - R2) * sin(psi) - R2 * psi) / R1;
        dphi = ((R2 * (psi - sin(psi)) - punching.Rleg * theta) * dR1 + R1 * (sin(psi) - psi) * dR2 - R1 * ((R1 - R2) * cos(psi) + R2) * dpsi) / R1 ^ 2;
        
        % Calculation of dEm
        dEmz = 2 * R1 * dR1 * (1 + cos(phi + psi)) - R1 ^ 2 * sin(psi + phi) * (dpsi + dphi);
        
    end
    
    dEm(j+1) = dEmz;
end

dtheta = pi/jmax;
Emsum = 0;
for j = 1:jmax
    Emsum = Emsum + punching.n0*punching.Rleg*(dEm(j)+dEm(j+1))/2*dtheta;
end

Em = Emsum;
if dloc == 0
    Em0 = Em;
end

% ----------------
% Local resistance
% ----------------

ksi1 = min((punching.NB * punching.a + punching.gap) + (2 * Em / Eb) ^ 0.5, punching.L1);
ksi2 = min((punching.NB * punching.a + punching.gap) + (2 * Em / Eb) ^ 0.5, punching.L2);

PX = (Em * (1 / (ksi1 - (punching.NB * punching.a + punching.gap)) + 1 / (ksi2 - (punching.NB * punching.a + punching.gap))) + Eb * ((punching.NB * punching.a + punching.gap) + (ksi1 + ksi2) / 2)) * cos(punching.alphaP);
PY = (Em * (1 / (ksi1 - (punching.NB * punching.a + punching.gap)) + 1 / (ksi2 - (punching.NB * punching.a + punching.gap))) + Eb * ((punching.NB * punching.a + punching.gap) + (ksi1 + ksi2) / 2)) * sin(punching.alphaP);
Epunch = Epunch + (Em * (1 / (ksi1 - (punching.NB * punching.a + punching.gap)) + 1 / (ksi2 - (punching.NB * punching.a + punching.gap))) + Eb * ((punching.NB * punching.a + punching.gap) + (ksi1 + ksi2) / 2)) * ddelta;

P = (PX ^ 2 + PY ^ 2) ^ 0.5;

% Calcul de la résistance globale initiale
beta = punching.Dbrace / 2 / punching.Rleg;
if beta > 0.6
    Qb = 0.3 / (beta * (1 - 0.833 * beta));
end

CP = Qb ^ 0.5 * beta * 0.5 * ((dloc * cos(punching.alphaP) / 2 / punching.Rleg) ^ 2 - 1) * (dloc * cos(punching.alphaP) / 2 / punching.Rleg - 2); % reduction factor
x = 0.5 * ((dloc * cos(punching.alphaP) / 2 / punching.Rleg) ^ 2 - 1) * (dloc * cos(punching.alphaP) / 2 / punching.Rleg - 2);

PG = (CP + punching.ksi) * punching.Mp * (punching.Lleg - 2*punching.NB*punching.a) / ((punching.L1 - punching.NB*punching.a) * (punching.L2 - punching.NB*punching.a));

if P > PG
    dP = dloc;
    mode = 1;
end
    
% end

% mode = 0;

%% Global mode

% if mode == 1
%     % Enfoncement global
%     dglob = delta - dP;
%     
%     % Calculation of resistance
%     N = (punching.Lleg - (punching.NB * 2 * punching.a + punching.gap)) / ((punching.L1 - (punching.NB * punching.a + punching.gap / 2)) * (punching.L2 - (punching.NB * punching.a + punching.gap / 2)) ^ 2 + (punching.L2 - (punching.NB * punching.a + punching.gap / 2)) * (punching.L1 - (punching.NB * punching.a + punching.gap / 2)) ^ 2);
%     N = (punching.L1 - (punching.NB * punching.a + punching.gap / 2)) * (punching.L2 - (punching.NB * punching.a + punching.gap / 2)) * N;
%     N = punching.Np ^ 2 / 2 / (CP + punching.ksi) / punching.Mp * dglob * cos(punching.alphaP) * (1 + N) / 2;
%     N = min(N, punching.Np);
%     P = (punching.Lleg - (punching.NB * 2 * punching.a + punching.gap)) / (punching.L1 - (punching.NB * punching.a + punching.gap / 2)) / (punching.L2 - (punching.NB * punching.a + punching.gap / 2));
%     P = P * (punching.Mp * (CP + punching.ksi) * (1 - N ^ 2 / punching.Np ^ 2) + N * dglob * cos(punching.alphaP));
%     PX = P * cos(punching.alphaP);
%     PY = P * sin(punching.alphaP);
%     Epunch = Epunch + P * cos(punching.alphaP) * ddelta;
%     
% end

%% Inertia of reduced section

Rmax = punching.Rleg+punching.tleg/2;
Rmin = punching.Rleg-punching.tleg/2;
Iinit = pi/64 * ((2*Rmax)^4 - (2*Rmin)^4);
Atot = pi*(Rmax^2 - Rmin^2);

% Static moment
R1down = R1-punching.tleg/2;
R1top = R1+punching.tleg/2;
theta1down = psi;
theta1top = 2*pi-psi;
SC1 = -2*sin(psi) * (R1top^3 - R1down^3)/3 + R1*(R1top^2 - R1down^2)/2*(theta1top-theta1down) - punching.Rleg * (R1top^2 - R1down^2)/2*(theta1top-theta1down);

R2down = R2-punching.tleg/2;
R2top = R2+punching.tleg/2;
theta2top = psi;
theta2down = 0;
SC2 = 2* (sin(psi)*(R2top^3 - R2down^3)/3 + (R1 - (R1-R2)*cos(pi-psi) - punching.Rleg) * (R2top^2-R2down^2)/2 * (theta2top-theta2down));

SStraight = 2*(R1-R2)*sin(psi)*punching.tleg*(punching.Rleg-dloc);

Stot = SC1 + SC2 + SStraight;

xCG = Stot/Atot;

% Inertia
% Curve C1
IC1part1 = 0.5 * (R1top^4-R1down^4)/4 * ((0.5*sin(2*theta1top)+theta1top) - (0.5*sin(2*theta1down)+theta1down));
IC1part2 = 2 * (R1top^3-R1down^3)/3 * (sin(theta1top)-sin(theta1down)) * (R1 - punching.Rleg - xCG);
IC1part3 = (R1top^2-R1down^2)/2 * (theta1top-theta1down) * (R1 - punching.Rleg - xCG)^2;

IC1 = IC1part1 + IC1part2 + IC1part3;

% Curve C2
IC2part1 = 0.5 * (R2top^4-R2down^4)/4 * (0.5*sin(2*theta2top) + theta2top);
IC2part2 = 2 * (R2top^3-R2down^3)/3 * sin(theta2top) * (R1 - (R1-R2)*cos(pi-psi) - punching.Rleg - xCG);
IC2part3 = (R2top^2-R2down^2)/2 * (theta2top-theta2down) * (R1 - (R1-R2)*cos(pi-psi) - punching.Rleg - xCG)^2;

IC2 = 2 * (IC2part1 + IC2part2 + IC2part3);

% Straight part
IStraightOwn = (R1-R2)*sin(psi)*punching.tleg^3/12;
IStraightTransp = (R1-R2)*sin(psi)*punching.tleg * (punching.Rleg-dloc-xCG)^2;

IStraight = IStraightOwn + IStraightTransp;

% Total
IRedSect = IC1 + IC2 + IStraight;

% Comparison
rapportI = IRedSect/Iinit;

%% Results

resultPunch.Ptotpunch = P * 10^(-6);
resultPunch.PXtotpunch = PX * 10^(-6);
resultPunch.PYtotpunch = PY * 10^(-6);
resultPunch.Etotpunch = Epunch * 10^(-6);

punchProp.punchModeNodes(i,2) = delta;
punchProp.punchModeNodes(i,3) = x;
punchProp.punchModeNodes(i,4) = rapportI;