function[resultVert,resultVertVec] = solveVert2(properties,data,cylinderVec,shipVec,param,impactVec,force,resultVertVec,prevDispLoc)

mode = resultVertVec(1); % 0-> local ; 1-> global
ddelta = norm(force.valueDisp);
% delta = resultVert.dispLoc + ddelta;
delta = prevDispLoc + ddelta;

E90 = resultVertVec(2);
d90 = resultVertVec(3);
d90l = resultVertVec(4);
d90r = resultVertVec(5);
Cv = resultVertVec(6);
Cvl = resultVertVec(7);
Cvr = resultVertVec(8);
gamma = resultVertVec(9);

% Cylinder
cylinder.D = cylinderVec(1);
cylinder.t = cylinderVec(2);
cylinder.dzeta = cylinderVec(3);
cylinder.L = cylinderVec(4);
cylinder.fy = cylinderVec(5);
cylinder.rearImpact = cylinderVec(6);
cylinder.R = cylinderVec(7);
cylinder.Np = cylinderVec(8);
cylinder.Mp = cylinderVec(9);
cylinder.n0 = cylinderVec(10);
cylinder.m0 = cylinderVec(11);
deltaMax = cylinderVec(12);

% Ship
ship.bulbous = shipVec(1);
ship.p = shipVec(2);
ship.q = shipVec(3);
ship.phib = shipVec(4);
ship.psib = shipVec(5);
ship.hb = shipVec(6);
ship.hBulb = shipVec(7);
ship.bBulb = shipVec(8);
ship.pBulb = shipVec(9);
ship.shiftCentreBulb = shipVec(10);
ship.alpha = shipVec(11);
ship.phibprime = shipVec(12);

% Impact
impact.xsI = impactVec(1);
impact.ysI = impactVec(2);
impact.zsI = impactVec(3);
impact.L1 = impactVec(4);
impact.L2 = impactVec(5);
impact.YI = impactVec(6);
impact.ZI = impactVec(7);
impact.pv = impactVec(8);
impact.qv = impactVec(9);
impact.hbv = impactVec(10);
impact.xvI = impactVec(11);
impact.yvI = impactVec(12);
impact.zvI = impactVec(13);
impact.xvS = impactVec(14);
impact.yvS = impactVec(15);
impact.zvS = impactVec(16);
impact.betav = impactVec(17);
impact.df = impactVec(18);
impact.nv = impactVec(19);
impact.psiv0 = impactVec(20);
impact.psivf = impactVec(21);
impact.eta1 = impactVec(22);
impact.eta2 = impactVec(23);
impact.xhI = impactVec(24);
impact.yhI = impactVec(25);
impact.zhI = impactVec(26);
impact.alphah = impactVec(27);
impact.xhS = impactVec(28);
impact.yhS = impactVec(29);
impact.zhS = impactVec(30);
impact.betah = impactVec(31);
impact.nh = impactVec(32);
impact.psih0 = impactVec(33);

%% Local mode

if mode == 0
    
    dloc = delta;
    
    psi = impact.psiv0 + (impact.psivf - impact.psiv0) * (dloc / impact.df) ^ impact.nv;
    dpsi = impact.nv * (impact.psivf - impact.psiv0) * dloc ^ (impact.nv - 1) / impact.df ^ impact.nv;
    
    % -----
    % Rings
    % -----
    
    % Impact excentré
    if impact.xvS ~= 0
        
        % Calcul des paramètres de section
        xA = (impact.qv ^ 2 / impact.pv ^ 2 * impact.xvS - (dloc - impact.zvS) * cot(impact.betav) - sign(impact.xvS) * impact.qv / impact.pv * (impact.qv ^ 2 + impact.pv ^ 2 * cot(impact.betav) ^ 2 - (impact.xvS * cot(impact.betav) + dloc - impact.zvS) ^ 2) ^ 0.5) / (impact.qv ^ 2 / impact.pv ^ 2 + cot(impact.betav) ^ 2);
        gamma = atan(impact.qv / impact.pv * (impact.xvS - xA) / (impact.pv ^ 2 - (xA - impact.xvS) ^ 2) ^ 0.5);
        R1 = cylinder.R + (psi - sin(psi)) * (cylinder.R - xA * (cot(impact.betav) * cos(gamma) + sin(gamma))) / (pi * (1 - cos(psi)) - 2 * (psi - sin(psi)));
        R2 = (cylinder.R + xA * (cot(impact.betav) * cos(gamma) + sin(gamma)) - R1 * (1 + cos(psi))) / (1 - cos(psi));
        dxA = (impact.qv / impact.pv * sign(impact.xvS) * (impact.xvS * cot(impact.betav) + dloc - impact.zvS) / (impact.qv ^ 2 + impact.pv ^ 2 * cot(impact.betav) ^ 2 - (impact.xvS * cot(impact.betav) + dloc - impact.zvS) ^ 2) ^ 0.5 - cot(impact.betav)) / (impact.qv ^ 2 / impact.pv ^ 2 + cot(impact.betav) ^ 2);
        dgamma = -impact.pv ^ 3 * impact.qv / (impact.pv ^ 4 + (impact.qv ^ 2 - impact.pv ^ 2) * (xA - impact.xvS) ^ 2) / (impact.pv ^ 2 - (xA - impact.xvS) ^ 2) ^ 0.5 * dxA;
        
        % Calcul de dR1
        dR1_xA = -(cot(impact.betav) * cos(gamma) + sin(gamma)) * (psi - sin(psi)) / (pi * (1 - cos(psi)) - 2 * (psi - sin(psi)));
        dR1_psi = (cylinder.R - xA * (cot(impact.betav) * cos(gamma) + sin(gamma))) * (2 * (1 - cos(psi)) - psi * sin(psi)) / (pi * (1 - cos(psi)) - 2 * (psi - sin(psi))) ^ 2 * pi;
        dR1_gamma = -xA * (cos(gamma) - cot(impact.betav) * sin(gamma)) * (psi - sin(psi)) / (pi * (1 - cos(psi)) - 2 * (psi - sin(psi)));
        dR1 = dR1_psi * dpsi + dR1_xA * dxA + dR1_gamma * dgamma;
        
        % Caclul de dR2
        dR2_xA = (cot(impact.betav) * cos(gamma) + sin(gamma)) / (1 - cos(psi));
        dR2_psi = (2 * R1 - cylinder.R - xA * (cot(impact.betav) * cos(gamma) + sin(gamma))) * sin(psi) / (1 - cos(psi)) ^ 2;
        dR2_gamma = xA * (cos(gamma) - cot(impact.betav) * sin(gamma)) / (1 - cos(psi));
        dR2 = dR2_psi * dpsi + dR2_xA * dxA + dR2_gamma * dgamma - (1 + cos(psi)) / (1 - cos(psi)) * dR1;
        
        % Calcul des vitesses
        VB = (R1 - R2) * dpsi - (pi - psi) * dR1 - psi * dR2;
        VF = R1 * dpsi - (pi - psi) * dR1;
        
        % Calcul de la résistance
        Eb = min(2 * cylinder.m0 * (abs(VB) / R2 + abs(VF) * (1 / R2 - 1 / R1) + psi / R2 * abs(dR2) + (pi - psi) / R1 * abs(dR1)), 8 * cylinder.m0 / cylinder.R);
        
        
        % Impact centré
    else
        
        % Calcul des paramètres de section
        xA = 0;
        gamma = 0;
        R1 = cylinder.R + dloc * (psi - sin(psi)) / (pi * (1 - cos(psi)) - 2 * (psi - sin(psi)));
        R2 = (2 * cylinder.R - dloc - R1 * (1 + cos(psi))) / (1 - cos(psi));
        dxA = 0;
        dgamma = 0;
        
        % Calcul de dR1
        dR1_xA = (psi - sin(psi)) / (pi * (1 - cos(psi)) - 2 * (psi - sin(psi)));
        dR1_psi = dloc * pi * (2 * (1 - cos(psi)) - psi * sin(psi)) / ((pi * (1 - cos(psi)) - 2 * (psi - sin(psi))) ^ 2);
        dR1 = dR1_psi * dpsi + dR1_xA;
        
        % Calcul de dR2
        dR2_xA = -1 / (1 - cos(psi));
        dR2_psi = sin(psi) * (2 * R1 - 2 * cylinder.R + dloc) / ((1 - cos(psi)) ^ 2);
        dR2 = dR2_psi * dpsi + dR2_xA - (1 + cos(psi)) / (1 - cos(psi)) * dR1;
        
        % Calcul des vitesses
        VB = (R1 - R2) * dpsi - (pi - psi) * dR1 - psi * dR2;
        VF = R1 * dpsi - (pi - psi) * dR1;
        
        % Calcul de la résistance
        Eb = min(2 * cylinder.m0 * (abs(VB) / R2 + abs(VF) * (1 / R2 - 1 / R1) + psi / R2 * abs(dR2) + (pi - psi) / R1 * abs(dR1)), 8 * cylinder.m0 / cylinder.R);
        
    end
    
    % ----------
    % Generators
    % ----------
    
    jmax = 100;
    dEm = zeros(1,jmax+1);
    
    if impact.xvS ~= 0
        xD = xA * sin(gamma) / sin(impact.betav) * cos(gamma - impact.betav);
        dxD = sin(gamma) * cos(gamma - impact.betav) / sin(impact.betav) * dxA + xA * cos(2 * gamma - impact.betav) / sin(impact.betav) * dgamma;
    else
        xD = 0;
        dxD = 0;
    end
    
    for j = 0:jmax
        theta = j*pi/jmax;
        
        % Zone 1
        if cylinder.R * theta <= (R1 - R2) * sin(psi)
            
            % Impact excentré
            if impact.xvS ~= 0
                dEm1 = 1 / sin(impact.betav) * (xA * sin(gamma - impact.betav) * dgamma - cos(gamma - impact.betav) * dxA) * (cylinder.R * cos(theta) - xD / sin(gamma));
                dEm1 = dEm1 + cylinder.R * (cylinder.R * theta * (1 - cos(theta)) - (cylinder.R - xD / sin(gamma)) * sin(theta)) * dgamma;
                dEm2 = 1 / sin(impact.betav) * (xA * sin(gamma - impact.betav) * dgamma - cos(gamma - impact.betav) * dxA) * (cylinder.R * cos(theta) - xD / sin(gamma));
                dEm2 = dEm2 - cylinder.R * (cylinder.R * theta * (1 - cos(theta)) - (cylinder.R - xD / sin(gamma)) * sin(theta)) * dgamma;
                
                % Impact centré
            else
                dEm1 = cylinder.R * cos(theta) - (cylinder.R - dloc);
                dEm2 = cylinder.R * cos(theta) - (cylinder.R - dloc);
            end
            
            % Zone 2
        elseif cylinder.R * theta > (R1 - R2) * sin(psi) && cylinder.R * theta < (R1 - R2) * sin(psi) + R2 * psi
            
            % Calcul de phi
            phi = (cylinder.R * theta - (R1 - R2) * sin(psi)) / R2;
            dphi1 = -(R2 * sin(psi) * dR1 + (R1 - R2) * R2 * cos(psi) * dpsi + (cylinder.R * theta - R1 * sin(psi)) * dR2 - cylinder.R * R2 * dgamma) / R2 ^ 2;
            dphi2 = -(R2 * sin(psi) * dR1 + (R1 - R2) * R2 * cos(psi) * dpsi + (cylinder.R * theta - R1 * sin(psi)) * dR2 + cylinder.R * R2 * dgamma) / R2 ^ 2;
            
            % Calcul de dEm
            dEm1 = (2 * R1 * (1 + cos(psi)) - R2 * (1 + cos(psi) - cos(phi) - cos(psi - phi)) - cylinder.R * (1 + cos(theta) + cos(psi) + cos(psi - theta))) * dR1;
            dEm1 = dEm1 + (2 * R2 * (1 - cos(psi - phi)) - R1 * (1 + cos(psi) - cos(phi) - cos(psi - phi)) - cylinder.R * (cos(phi - theta) + cos(phi) - cos(psi - theta) - cos(psi))) * dR2;
            dEm1 = dEm1 + (-R1 ^ 2 * sin(psi) + R2 ^ 2 * sin(psi - phi) - R1 * R2 * (sin(psi - phi) - sin(psi)) + (R1 - R2) * cylinder.R * (sin(psi) + sin(psi - theta))) * dpsi;
            dEm1 = dEm1 - R2 * (R2 * sin(psi - phi) - R1 * (sin(psi - phi) - sin(phi)) - cylinder.R * (sin(phi) + sin(phi - theta))) * dphi1;
            dEm1 = dEm1 - cylinder.R * ((cylinder.R - R1) * sin(theta) + (R1 - R2) * sin(psi - theta) + R2 * sin(phi - theta)) * dgamma;
            dEm2 = (2 * R1 * (1 + cos(psi)) - R2 * (1 + cos(psi) - cos(phi) - cos(psi - phi)) - cylinder.R * (1 + cos(theta) + cos(psi) + cos(psi - theta))) * dR1;
            dEm2 = dEm2 + (2 * R2 * (1 - cos(psi - phi)) - R1 * (1 + cos(psi) - cos(phi) - cos(psi - phi)) - cylinder.R * (cos(phi - theta) + cos(phi) - cos(psi - theta) - cos(psi))) * dR2;
            dEm2 = dEm2 + (-R1 ^ 2 * sin(psi) + R2 ^ 2 * sin(psi - phi) - R1 * R2 * (sin(psi - phi) - sin(psi)) + (R1 - R2) * cylinder.R * (sin(psi) + sin(psi - theta))) * dpsi;
            dEm2 = dEm2 - R2 * (R2 * sin(psi - phi) - R1 * (sin(psi - phi) - sin(phi)) - cylinder.R * (sin(phi) + sin(phi - theta))) * dphi2;
            dEm2 = dEm2 + cylinder.R * ((cylinder.R - R1) * sin(theta) + (R1 - R2) * sin(psi - theta) + R2 * sin(phi - theta)) * dgamma;
            
            % Zone 3
        else
            
            % Calcul de phi
            phi = (cylinder.R * theta - (R1 - R2) * sin(psi) - R2 * psi) / R1;
            dphi1 = ((R2 * (psi - sin(psi)) - cylinder.R * theta) * dR1 + R1 * (sin(psi) - psi) * dR2 - R1 * ((R1 - R2) * cos(psi) + R2) * dpsi - cylinder.R * R1 * dgamma) / R1 ^ 2;
            dphi2 = ((R2 * (psi - sin(psi)) - cylinder.R * theta) * dR1 + R1 * (sin(psi) - psi) * dR2 - R1 * ((R1 - R2) * cos(psi) + R2) * dpsi + cylinder.R * R1 * dgamma) / R1 ^ 2;
            
            % Calcul de dEm
            dEm1 = 2 * R1 * dR1 * (1 + cos(phi + psi)) - R1 ^ 2 * sin(psi + phi) * (dpsi + dphi1);
            dEm1 = dEm1 - cylinder.R ^ 2 * sin(theta) * dgamma;
            dEm2 = 2 * R1 * dR1 * (1 + cos(phi + psi)) - R1 ^ 2 * sin(psi + phi) * (dpsi + dphi2);
            dEm2 = dEm2 + cylinder.R ^ 2 * sin(theta) * dgamma;
            
        end
        
        dEm(j+1) = abs(dEm1) + abs(dEm2);
    end
    
    dtheta = pi/jmax;
    Emsum = 0;
    for j = 1:jmax
        Emsum = Emsum + cylinder.n0*cylinder.R*(dEm(j)+dEm(j+1))/2*dtheta;
    end
    
    Em = Emsum;
    
    % ----------------
    % Local resistance
    % ----------------
    
    if impact.xvS ~= 0
        da = -cylinder.R / impact.xvI * dxA;
    else
        da = 1;
    end
    
    if data.impactedZone == 1 %Elements
        if properties.type(data.impactedElement) == 1 %Leg
            
            ksi1 = (2 * Em / Eb) ^ 0.5;
            ksi2 = (2 * Em / Eb) ^ 0.5;
            
            theta = pi/2 - gamma;
            if delta == 0
                d90l = 0;
                d90r = 0;
            else
                if impact.xvS ~= 0
                    dW0 = 1 / sin(impact.betav) * (xA * sin(gamma - impact.betav) * dgamma - cos(gamma - impact.betav) * dxA) * (cylinder.R * cos(gamma - impact.betav) - xD / sin(gamma));
                    dW0 = dW0 + cylinder.R * (cylinder.R * theta * (1 - cos(theta)) - (cylinder.R - xD / sin(gamma)) * sin(theta)) * dgamma;
                    dW0 = dW0 / ((cylinder.R * cos(gamma - impact.betav) - xD / sin(gamma)) ^ 2 + cylinder.R ^ 2 * (gamma - impact.betav - sin(gamma - impact.betav)) ^ 2) ^ 0.5;
                else % Impact centré
                    dW0 = 1;
                end
            end
            
            % Calcul des déplacements
            for j = 0:1
                y = j * (impact.L1 + impact.L2);
                if y > impact.L1 - ksi1 && y <= impact.L1
                    d90l = d90l + dW0 * (1 - (impact.L1 - y) / ksi1) * ddelta;
                elseif y > impact.L1 && y < impact.L1 + ksi2
                    d90r = d90r + dW0 * (1 + (impact.L1 - y) / ksi2) * ddelta;
                end
            end
            
            CvLeft = 0.5 * ((d90l/(2*cylinder.R))^2 - 1) * (d90l/(2*cylinder.R) - 2);
            CvRight = 0.5 * ((d90r/(2*cylinder.R))^2 - 1) * (d90r/(2*cylinder.R) - 2);
            
        else
            
            ksi1 = min((2 * Em / Eb) ^ 0.5, impact.L1);
            ksi2 = min((2 * Em / Eb) ^ 0.5, impact.L2);
            
            CvLeft = 1;
            CvRight = 1;
            
        end
        
    else %Nodes
        
        ksi1 = (2 * Em / Eb) ^ 0.5;
        ksi2 = (2 * Em / Eb) ^ 0.5;
        
    end
    
    if impact.xvS ~= 0
        PX90 = (Em * (1 / ksi1 + 1 / ksi2) + Eb * (ksi1 + ksi2) / 2) / da / cos(gamma - impact.betav) * cos(gamma + ship.alpha);
        PY90 = (Em * (1 / ksi1 + 1 / ksi2) + Eb * (ksi1 + ksi2) / 2) / da / cos(gamma - impact.betav) * sin(gamma + ship.alpha);
        E90 = E90 + (Em * (1 / ksi1 + 1 / ksi2) + Eb * (ksi1 + ksi2) / 2) * ddelta;
    else
        PX90 = Em * (1 / ksi1 + 1 / ksi2) + Eb * (ksi1 + ksi2) / 2;
        PY90 = 0;
        E90 = E90 + (Em * (1 / ksi1 + 1 / ksi2) + Eb * (ksi1 + ksi2) / 2) * ddelta;
    end
    
    P90 = (PX90 ^ 2 + PY90 ^ 2) ^ 0.5;
    
    if data.impactedZone == 1
        % Calcul de la résistance globale initiale
        if impact.xvS ~= 0
            aloc = abs(xA - impact.xvI) * (1 + cot(impact.betav) ^ 2) ^ 0.5;
        else
            aloc = dloc;
        end
        Cv = 0.5 * ((cylinder.R - (cylinder.R - aloc) * cos(gamma - impact.betav)) ^ 2 / 4 / cylinder.R ^ 2 - 1) * ((cylinder.R - (cylinder.R - aloc) * cos(gamma - impact.betav)) / 2 / cylinder.R - 2);
        
        PG90 = cylinder.Mp / (impact.L1 * impact.L2) * (impact.L1*(Cv+CvRight) + impact.L2*(Cv+CvLeft));
        
        % Transition verticale
        if P90 > PG90
            mode = 1;
            d90 = dloc;
        end
    end
    
end

%% Global mode

if mode == 1
    % Enfoncement global
    dglob = delta - d90;
    
    % Calcul de la résistance globale
    N = cylinder.Np ^ 2 / 2 / cylinder.Mp / (Cv + param.ksi) * dglob * cos(gamma);
    N = min(N, cylinder.Np);
    P90 = ((Cv + param.ksi) * cylinder.Mp * (1 - N ^ 2 / cylinder.Np ^ 2) + N * dglob * cos(gamma)) * (impact.L1 + impact.L2) / impact.L1 / impact.L2;
    PX90 = P90 * cos(ship.alpha + gamma);
    PY90 = P90 * sin(ship.alpha + gamma);
    E90 = E90 + (PX90 * cos(ship.alpha) + PY90 * sin(ship.alpha)) * ddelta;
end

%% Total

Ptot90 = P90*10^(-6);
PXtot90 = PX90*10^(-6);
PYtot90 = PY90*10^(-6);
Etot90 = E90*10^(-6);

%% Results

resultVert.delta = delta;
resultVert.Ptot90 = Ptot90;
resultVert.PXtot90 = PXtot90;
resultVert.PYtot90 = PYtot90;
resultVert.Etot90 = Etot90;

resultVertVec = zeros(1,9);
resultVertVec(1) = mode;
resultVertVec(2) = E90;
resultVertVec(3) = d90;
resultVertVec(4) = d90l;
resultVertVec(5) = d90r;
resultVertVec(6) = Cv;
resultVertVec(7) = Cvl;
resultVertVec(8) = Cvr;
resultVertVec(9) = gamma;