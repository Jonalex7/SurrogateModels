function[resultHor] = solveHor(properties,data,cylinder,ship,param,impact,force,resultHor,prevDispLoc)

mode = resultHor.mode; % 0-> local ; 1-> global
ddelta = norm(force.valueDisp);
% delta = resultHor.dispLoc + ddelta;
delta = prevDispLoc + ddelta;

t = resultHor.t;
dt = resultHor.dt;
Eb0 = resultHor.Eb0;
Ebt = resultHor.Ebt;
Em0 = resultHor.Em0;
Emt = resultHor.Emt;
Cht = resultHor.Cht;
Ch = resultHor.Ch;
Chl = resultHor.Chl;
Chr = resultHor.Chr;

R1 = resultHor.R1;
R2 = resultHor.R2;
dR1 = resultHor.dR1;
dR2 = resultHor.dR2;

E0 = resultHor.E0;
d0 = resultHor.d0;
d0l = resultHor.d0l;
d0r = resultHor.d0r;

%% Local mode

if mode == 0
    
    dloc = delta;
    
    psi = impact.psih0 + (pi - impact.psih0) * (dloc * cos(impact.alphah) * sin(impact.betah) / 2 / cylinder.R) ^ impact.nh;
    dpsi = impact.nh * (pi - impact.psih0) * dloc ^ (impact.nh - 1) * (cos(impact.alphah) * sin(impact.betah) / 2 / cylinder.R) ^ impact.nh;
    
    % -----
    % Rings
    % -----
    
    if t == 1
        Eb = Ebt + (Ebt - Eb0) /dt * (dloc - dt);
    elseif t == 0
        % Calcul des paramètres de section
        R1 = cylinder.R + dloc * cos(impact.alphah) * sin(impact.betah) * (psi - sin(psi)) / (pi * (1 - cos(psi)) - 2 * (psi - sin(psi)));
        R2 = cylinder.R - dloc * cos(impact.alphah) * sin(impact.betah) * (pi - psi + sin(psi)) / (pi * (1 - cos(psi)) - 2 * (psi - sin(psi)));
        dR1 = pi * dloc * (2 * (1 - cos(psi)) - psi * sin(psi)) / (pi * (1 - cos(psi)) - 2 * (psi - sin(psi))) * dpsi + psi - sin(psi);
        dR1 = cos(impact.alphah) * sin(impact.betah) / (pi * (1 - cos(psi)) - 2 * (psi - sin(psi))) * dR1;
        dR2 = pi * dloc * (pi - psi) * sin(psi) / (pi * (1 - cos(psi)) - 2 * (psi - sin(psi))) * dpsi - (pi - psi + sin(psi));
        dR2 = cos(impact.alphah) * sin(impact.betah) / (pi * (1 - cos(psi)) - 2 * (psi - sin(psi))) * dR2;
        
        % Calcul des vitesses
        VB = (R1 - R2) * dpsi - (pi - psi) * dR1 - psi * dR2;
        VF = R1 * dpsi - (pi - psi) * dR1;
        
        % Calcul de la résistance
        Eb = 2 * cylinder.m0 * (abs(VB) / R2 + abs(VF) * (1 / R2 - 1 / R1) + psi / R2 * abs(dR2) + (pi - psi) / R1 * abs(dR1));
        Ebt = Eb;
        dt = dloc;
        if dloc == 0
            Eb0 = Eb;
        end
    end
    
    % ----------
    % Generators
    % ----------
    
    if t == 1
        Em = Emt + (Emt - Em0) / dt * (dloc - dt);
    elseif t == 0
        
        jmax = 100;
        dEmVec = zeros(1,jmax+1);
        
        for j = 0:jmax
            theta = j*pi/jmax;
            
            % Zone 1
            if cylinder.R * theta <= (R1 - R2) * sin(psi)
                dEm = (dloc * cos(impact.alphah) * sin(impact.betah) - cylinder.R * (1 - cos(theta))) * cos(impact.alphah) * sin(impact.betah);
                
                % Zone 2
            elseif cylinder.R * theta > (R1 - R2) * sin(psi) && cylinder.R * theta < (R1 - R2) * sin(psi) + R2 * psi
                % Calcul de phi
                phi = (cylinder.R * theta - (R1 - R2) * sin(psi)) / R2;
                dphi = -(R2 * sin(psi) * dR1 + (R1 - R2) * R2 * cos(psi) * dpsi + (cylinder.R * theta - R1 * sin(psi)) * dR2) / R2 ^ 2;
                
                % Calcul de dEm
                dEm = (2 * R1 * (1 + cos(psi)) - R2 * (1 + cos(psi) - cos(phi) - cos(psi - phi)) - cylinder.R * (1 + cos(theta) + cos(psi) + cos(psi - theta))) * dR1;
                dEm = dEm + (2 * R2 * (1 - cos(psi - phi)) - R1 * (1 + cos(psi) - cos(phi) - cos(psi - phi)) - cylinder.R * (cos(phi - theta) + cos(phi) - cos(psi - theta) - cos(psi))) * dR2;
                dEm = dEm + (-R1 ^ 2 * sin(psi) + R2 ^ 2 * sin(psi - phi) - R1 * R2 * (sin(psi - phi) - sin(psi)) + (R1 - R2) * cylinder.R * (sin(psi) + sin(psi - theta))) * dpsi;
                dEm = dEm - R2 * (R2 * sin(psi - phi) - R1 * (sin(psi - phi) - sin(phi)) - cylinder.R * (sin(phi) + sin(phi - theta))) * dphi;
                
                % Zone 3
            elseif cylinder.R * theta >= (R1 - R2) * sin(psi) + R2 * psi
                % Calcul de phi
                phi = (cylinder.R * theta - (R1 - R2) * sin(psi) - R2 * psi) / R1;
                dphi = ((R2 * (psi - sin(psi)) - cylinder.R * theta) * dR1 + R1 * (sin(psi) - psi) * dR2 - R1 * ((R1 - R2) * cos(psi) + R2) * dpsi) / R1 ^ 2;
                
                % Calcul de dEm
                dEm = 2 * R1 * dR1 * (1 + cos(phi + psi)) - R1 ^ 2 * sin(psi + phi) * (dpsi + dphi);
                
            end
            
            dEmVec(j+1) = abs(dEm);
        end
        
        dtheta = pi/jmax;
        Emsum = 0;
        for j = 1:jmax
            Emsum = Emsum + 2*cylinder.n0*cylinder.R*(dEmVec(j)+dEmVec(j+1))/2*dtheta;
        end
        
        Em = Emsum;
        
        Emt = Em;
        if dloc == 0
            Em0 = Em;
        end
        
    end
    
    % ----------------
    % Local resistance
    % ----------------
    
    if data.impactedZone == 1
        if properties.type(data.impactedElement) == 1 %Leg
            
            ksi1 = (2 * Em / Eb) ^ 0.5;
            ksi2 = (2 * Em / Eb) ^ 0.5;
            
            if delta == 0
                d0l = 0;
                d0r = 0;
            else
                dW0 = cos(impact.alphah);
            end
            
            % Calcul des déplacements
            for j = 0:1
                y = j * (impact.L1 + impact.L2);
                if y > impact.yhI - dloc * sin(impact.alphah) - ksi1 && y <= impact.yhI - dloc * sin(impact.alphah)
                    d0l = d0l + dW0 * (1 - (impact.yhI - dloc * sin(impact.alphah) - y) / ksi1) * ddelta;
                elseif y > impact.yhI - dloc * sin(impact.alphah) && y < impact.yhI - dloc * sin(impact.alphah) + ksi2
                    d0r = d0r + dW0 * (1 + (impact.yhI - dloc * sin(impact.alphah) - y) / ksi2) * ddelta;
                end
            end
            
            ChLeft = 0.5 * ((d0l/(2*cylinder.R))^2 - 1) * (d0l/(2*cylinder.R) - 2);
            ChRight = 0.5 * ((d0r/(2*cylinder.R))^2 - 1) * (d0r/(2*cylinder.R) - 2);
            
        else %Brace and horizontal elements
            
            ksi1 = min((2 * Em / Eb) ^ 0.5, impact.yhI - dloc * sin(impact.alphah));
            ksi2 = min((2 * Em / Eb) ^ 0.5, cylinder.L - impact.yhI + dloc * sin(impact.alphah));
            
            ChLeft = 1;
            ChRight = 1;
            
        end
    
    else
        
        ksi1 = (2 * Em / Eb) ^ 0.5;
        ksi2 = (2 * Em / Eb) ^ 0.5;
        
    end
    
    PX0 = (Em * (1 / ksi1 + 1 / ksi2) + Eb * (ksi1 + ksi2) / 2) / cos(impact.alphah) * cos(ship.alpha - impact.eta1);
    PY0 = (Em * (1 / ksi1 + 1 / ksi2) + Eb * (ksi1 + ksi2) / 2) / cos(impact.alphah) * sin(ship.alpha - impact.eta1);
    E0 = E0 + (Em * (1 / ksi1 + 1 / ksi2) + Eb * (ksi1 + ksi2) / 2) * ddelta;
    
    P0 = (PX0 ^ 2 + PY0 ^ 2) ^ 0.5;
    
    if data.impactedZone == 1
        % Calcul de la résistance globale initiale
        if t == 1
            Ch = 1 + (Cht - 1) * dloc / dt;
        elseif t == 0
            Ch = 0.5 * ((dloc * cos(impact.alphah) * sin(impact.betah) / 2 / cylinder.R) ^ 2 - 1) * (dloc * cos(impact.alphah) * sin(impact.betah) / 2 / cylinder.R - 2);
            Cht = Ch;
        end
        
        l1 = impact.L1 - dloc * sin(impact.alphah);
        l2 = impact.L2 + dloc * sin(impact.alphah);
        PG0 = (cylinder.Mp)/(l1*l2) * (l1*(Ch+ChRight) + l2*(ChLeft+Ch));
        
        % Test transition débordement
        if (cylinder.R - dloc * cos(impact.alphah) * sin(impact.betah)) * cos(impact.betah) - (R1 - R2) * sin(psi) < impact.zhS - ship.hb || (cylinder.R - dloc * cos(impact.alphah) * sin(impact.betah)) * cos(impact.betah) + (R1 - R2) * sin(psi) > impact.zhS
            if dt ~= 0
                t = 1;
            else
                t = 0;
                dt = 0;
            end
        end
        
        % Transition horizontale
        if P0 > PG0
            mode = 1;
            d0 = dloc;
        end
    end
    
end

%% Global mode

if mode == 1
    
    % Enfoncement global
    dglob = delta - d0;
    
    % Calcul de la résistance
    N = (impact.L1 + impact.L2) / ((impact.L1 - d0 * sin(impact.alphah)) * (impact.L2 + delta * sin(impact.alphah)) ^ 2 + (impact.L2 + d0 * sin(impact.alphah)) * (impact.L1 - delta * sin(impact.alphah)) ^ 2);
    N = (impact.L1 - delta * sin(impact.alphah)) * (impact.L2 + delta * sin(impact.alphah)) * N;
    N = cylinder.Np ^ 2 / 2 / (Ch + param.ksi) / cylinder.Mp * dglob * cos(impact.alphah) * sin(impact.betah) * (1 + N) / 2;
    N = min(N, cylinder.Np);
    P0 = (impact.L1 + impact.L2) / (impact.L1 - delta * sin(impact.alphah)) / (impact.L2 + delta * sin(impact.alphah)) * sin(impact.betah);
    P0 = P0 * (cylinder.Mp * (Ch + param.ksi) * (1 - N ^ 2 / cylinder.Np ^ 2) + N * dglob * cos(impact.alphah) * sin(impact.betah));
    PX0 = P0 * cos(ship.alpha - impact.eta1);
    PY0 = P0 * sin(ship.alpha - impact.eta1);
    E0 = E0 + P0 * cos(impact.alphah) * ddelta;
    
end

%% Total

Ptot0 = P0*10^(-6);
PXtot0 = PX0*10^(-6);
PYtot0 = PY0*10^(-6);
Etot0 = E0*10^(-6);

%% Results

resultHor.delta = delta;
resultHor.Ptot0 = Ptot0;
resultHor.PXtot0 = PXtot0;
resultHor.PYtot0 = PYtot0;
resultHor.Etot0 = Etot0;

resultHor.dispLoc = delta;
resultHor.mode = mode;
resultHor.E0 = E0;
resultHor.d0 = d0;
resultHor.d0l = d0l;
resultHor.d0r = d0r;

resultHor.t = t;
resultHor.dt = dt;
resultHor.Eb0 = Eb0;
resultHor.Ebt = Ebt;
resultHor.Em0 = Em0;
resultHor.Emt = Emt;
resultHor.Cht = Cht;
resultHor.Ch = Ch;
resultHor.Chl = Chl;
resultHor.Chr = Chr;

resultHor.R1 = R1;
resultHor.R2 = R2;
resultHor.dR1 = dR1;
resultHor.dR2 = dR2;