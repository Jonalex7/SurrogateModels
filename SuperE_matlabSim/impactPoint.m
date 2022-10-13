function[impact] = impactPoint(cylinder,ship,d,param)

dmin = d(1);
dmax = d(2);

% Initialisation

M = 1;
d = dmax;
XI = -100;
YI = -100;
ZI = -100;
dis_min = 10^20;

% Recherche du point de contact

while dis_min > param.Tol && M <= 500

    d = (dmax + dmin) / 2;

    if abs(ship.alpha) > pi / 4 && abs(ship.alpha) <= pi / 2
        XS = ship.XP + d * cos(ship.alpha);
        YS = d * sin(ship.alpha);
    elseif abs(ship.alpha) <= pi / 4
        XS = d * cos(ship.alpha);
        YS = ship.YP + d * sin(ship.alpha);
    end
    

    % Cas 1: dzeta <= 45 deg
    if cylinder.dzeta <= pi / 4
        [XI,YI,ZI,cont,dis_min] = contact1(XS,YS,XI,YI,ZI,dis_min);

    % Cas 2: dzeta > 45 deg
    else
        [XI,YI,ZI,cont,dis_min] = contact2(XS,YS,XI,YI,ZI,dis_min);
    end

    % Incrémentation
    if cont == 0
        dmax = d;
    elseif cont == 1
        dmin = d;
    end
    M = M + 1;

end

%Position initiale

if abs(ship.alpha) > pi / 4 && abs(ship.alpha) <= pi / 2
    XS = ship.XP + d * cos(ship.alpha);
    YS = d * sin(ship.alpha);
elseif abs(ship.alpha) <= pi / 4
    XS = d * cos(ship.alpha);
    YS = ship.YP + d * sin(ship.alpha);
end

% Coordonnées de I dans les axes liés au navire

impact.xsI = (XI - XS) * cos(ship.alpha) + (YI - YS) * sin(ship.alpha);
impact.ysI = -(XI - XS) * sin(ship.alpha) + (YI - YS) * cos(ship.alpha);
impact.zsI = ZI - ship.ZS;

impact.XI = XI;
impact.YI = YI;
impact.ZI = ZI;

impact.L1 = YI * cos(cylinder.dzeta) + ZI * sin(cylinder.dzeta);
impact.L2 = cylinder.L - YI * cos(cylinder.dzeta) - ZI * sin(cylinder.dzeta);

% Contrôle

if XI >= XS + (ship.p + (ZI - ship.ZS) / tan(ship.psib)) * abs(sin(ship.alpha))
    disp('Erreur : arrêt du calcul')
end


%% function contact 1

    function[XI,YI,ZI,cont,dis_min] = contact1(XS,YS,XI,YI,ZI,dis_min)
        % Valeurs initiales
        
        cont = 0;
        R1 = ship.q + (min(ship.ZS, cylinder.L * sin(cylinder.dzeta) + cylinder.R * cos(cylinder.dzeta)) - ship.ZS) / tan(ship.phib);
        R2 = ship.p + (min(ship.ZS, cylinder.L * sin(cylinder.dzeta) + cylinder.R * cos(cylinder.dzeta)) - ship.ZS) / tan(ship.psib);
        
        % Limites étrave
        if ship.alpha >= 0
            Y1 = YS - ((ship.q - ship.hb / tan(ship.phib)) ^ 2 * sin(ship.alpha) ^ 2 + (ship.p - ship.hb / tan(ship.psib)) ^ 2 * cos(ship.alpha) ^ 2) ^ 0.5;
            Y2 = YS - (R1 ^ 2 * sin(ship.alpha) ^ 2 + R2 ^ 2 * cos(ship.alpha) ^ 2) ^ 0.5;
            Y3 = YS + (ship.p - ship.hb / tan(ship.psib)) * cos(ship.alpha);
            Y4 = YS + R2 * cos(ship.alpha);
            Yinf_e = Y2;
            Ysup_e = Y4;
        else
            Y1 = YS + ((ship.q - ship.hb / tan(ship.phib)) ^ 2 * sin(ship.alpha) ^ 2 + (ship.p - ship.hb / tan(ship.psib)) ^ 2 * cos(ship.alpha) ^ 2) ^ 0.5;
            Y2 = YS + (R1 ^ 2 * sin(ship.alpha) ^ 2 + R2 ^ 2 * cos(ship.alpha) ^ 2) ^ 0.5;
            Y3 = YS - (ship.p - ship.hb / tan(ship.psib)) * cos(ship.alpha);
            Y4 = YS - R2 * cos(ship.alpha);
            Yinf_e = Y4;
            Ysup_e = Y2;
        end
        
        % Limites cylindre
        Yinf_c = -cylinder.R * sin(cylinder.dzeta);
        Ysup_c = cylinder.L * cos(cylinder.dzeta) + cylinder.R * sin(cylinder.dzeta);
        
        % Limites Y
        Yinf = max(Yinf_e, Yinf_c) + param.Tol;
        Ysup = min(Ysup_e, Ysup_c) - param.Tol;
        
        % Boucle sur Y
        for i = 0:100
            
            Y = Yinf + (Ysup - Yinf) * i / 100;
            
            % Calcul des limites sur Z
            
            % Limites étrave
            if ship.alpha >= 0
                if Y2 <= Y && Y < Y1
                    Zinf_e = (Y - YS) ^ 2 * (cot(ship.phib) ^ 2 * sin(ship.alpha) ^ 2 + cot(ship.psib) ^ 2 * cos(ship.alpha) ^ 2) - (ship.p * cot(ship.phib) - ship.q * cot(ship.psib)) ^ 2 * sin(ship.alpha) ^ 2 * cos(ship.alpha) ^ 2;
                    Zinf_e = Zinf_e ^ 0.5 - (ship.q * cot(ship.phib) * sin(ship.alpha) ^ 2 + ship.p * cot(ship.psib) * cos(ship.alpha) ^ 2);
                    Zinf_e = ship.ZS + Zinf_e / (cot(ship.phib) ^ 2 * sin(ship.alpha) ^ 2 + cot(ship.psib) ^ 2 * cos(ship.alpha) ^ 2);
                end
                if Y1 <= Y && Y <= Y3
                    Zinf_e = ship.ZS - ship.hb;
                elseif Y3 < Y && Y <= Y4
                    Zinf_e = ship.ZS + (ship.p * cos(ship.alpha) + YS - Y) / cos(ship.alpha) / cot(ship.psib);
                end
                
            else
                
                if Y4 <= Y && Y <= Y3
                    Zinf_e = ship.ZS + (ship.p * cos(ship.alpha) + Y - YS) / cos(ship.alpha) / cot(ship.psib);
                elseif Y3 <= Y && Y <= Y2
                    Zinf_e = ship.ZS - ship.hb;
                elseif Y2 < Y && Y <= Y1
                    Zinf_e = (Y - YS) ^ 2 * (cot(ship.phib) ^ 2 * sin(ship.alpha) ^ 2 + cot(ship.psib) ^ 2 * cos(ship.alpha) ^ 2) - (ship.p * cot(ship.phib) - ship.q * cot(ship.psib)) ^ 2 * sin(ship.alpha) ^ 2 * cos(ship.alpha) ^ 2;
                    Zinf_e = Zinf_e ^ 0.5 - (ship.q * cot(ship.phib) * sin(ship.alpha) ^ 2 + ship.p * cot(ship.psib) * cos(ship.alpha) ^ 2);
                    Zinf_e = ship.ZS + Zinf_e / (cot(ship.phib) ^ 2 * sin(ship.alpha) ^ 2 + cot(ship.psib) ^ 2 * cos(ship.alpha) ^ 2);
                end
                
            end
            
            Zsup_e = ship.ZS;
            
            % Limites cylindre
            if cylinder.R * sin(cylinder.dzeta) <= Y && Y <= cylinder.L * cos(cylinder.dzeta) - cylinder.R * sin(cylinder.dzeta)
                Zsup_c = (Y * sin(cylinder.dzeta) + cylinder.R) / cos(cylinder.dzeta);
                Zinf_c = (Y * sin(cylinder.dzeta) - cylinder.R) / cos(cylinder.dzeta);
            end
            
            if cylinder.dzeta ~= 0 && -cylinder.R * sin(cylinder.dzeta) <= Y && Y < cylinder.R * sin(cylinder.dzeta)
                Zsup_c = (Y * sin(cylinder.dzeta) + cylinder.R) / cos(cylinder.dzeta);
                Zinf_c = -Y / tan(cylinder.dzeta);
            end
            
            if cylinder.dzeta ~= 0 && cylinder.L * cos(cylinder.dzeta) - cylinder.R * sin(cylinder.dzeta) < Y && Y <= cylinder.L * cos(cylinder.dzeta) + cylinder.R * sin(cylinder.dzeta)
                Zsup_c = (cylinder.L - Y * cos(cylinder.dzeta)) / sin(cylinder.dzeta);
                Zinf_c = (Y * sin(cylinder.dzeta) - cylinder.R) / cos(cylinder.dzeta);
            end
            
            Zinf = max(Zinf_c, Zinf_e) + param.Tol;
            Zsup = min(Zsup_c, Zsup_e) - param.Tol;
            
            % Boucle sur Z
            
            for j = 0:100
                
                Z = Zinf + (Zsup - Zinf) * j / 100;
                dis = 10 ^ 20;
                RX = ship.q + (Z - ship.ZS) / tan(ship.phib);
                RY = ship.p + (Z - ship.ZS) / tan(ship.psib);
                
                % Contact impossible
                if Zinf >= Zsup
                    dis = 10 ^ 20;
                    
                    % Contact possible
                elseif Zinf < Zsup
                    
                    % Calcul de Xc+
                    Xc1 = cylinder.R * (1 - (Z - Y * tan(cylinder.dzeta)) ^ 2 / (cylinder.R / cos(cylinder.dzeta)) ^ 2) ^ 0.5;
                    
                    % Calcul distances
                    Xe2 = RX ^ 2 * sin(ship.alpha) ^ 2 + RY ^ 2 * cos(ship.alpha) ^ 2 - (Y - YS) ^ 2;
                    Xe1 = (Y - YS) * (RX ^ 2 - RY ^ 2) * sin(ship.alpha) * cos(ship.alpha);
                    Xe1 = Xe1 - RX * RY * (RX ^ 2 * sin(ship.alpha) ^ 2 + RY ^ 2 * cos(ship.alpha) ^ 2 - (Y - YS) ^ 2) ^ 0.5;
                    Xe1 = Xe1 / (RX ^ 2 * sin(ship.alpha) ^ 2 + RY ^ 2 * cos(ship.alpha) ^ 2) + XS;
                    dis = abs(Xe1 - Xc1);
                    
                    % Test intrusion
                    if Xe1 <= Xc1
                        cont = 1;
                    end
                end
                
                % Test contact
                if dis < dis_min
                    XI = Xc1;
                    YI = Y;
                    ZI = Z;
                    dis_min = dis;
                end
                
                if cont == 1
                    break
                end
                
            end
            
            if cont == 1
                break
            end
            
        end
        
    end

%% function contact 2

    function[XI,YI,ZI,cont,dis_min] = contact2(XS,YS,XI,YI,ZI,dis_min)
        % Valeurs initiales
        
        cont = 0;
        
        % Limites cylindre
        Zinf_c = -cylinder.R * cos(cylinder.dzeta);
        Zsup_c = cylinder.L * sin(cylinder.dzeta) + cylinder.R * cos(cylinder.dzeta);
        
        % Limites étrave
        Zinf_e = ship.ZS - ship.hb;
        Zsup_e = ship.ZS;
        
        % Limites Z
        Zinf = max(Zinf_c, Zinf_e) + param.Tol;
        Zsup = min(Zsup_c, Zsup_e) - param.Tol;
        
        % Boucle sur Z
        for i = 0:100
            
            Z = Zinf + (Zsup - Zinf) * i / 100;
            RX = ship.q + (Z - ship.ZS) / tan(ship.phib);
            RY = ship.p + (Z - ship.ZS) / tan(ship.psib);
            
            % Calcul des limites sur Y
            
            % Limites étrave
            if ship.alpha >= 0
                Yinf_e = YS - (RX ^ 2 * sin(ship.alpha) ^ 2 + RY ^ 2 * cos(ship.alpha) ^ 2) ^ 0.5;
                Ysup_e = YS + RY * cos(ship.alpha);
            else
                Yinf_e = YS - RY * cos(ship.alpha);
                Ysup_e = YS + (RX ^ 2 * sin(ship.alpha) ^ 2 + RY ^ 2 * cos(ship.alpha) ^ 2) ^ 0.5;
            end
            
            % Limites cylindre
            if cylinder.dzeta ~= pi / 2 && -cylinder.R * cos(cylinder.dzeta) <= Z && Z < cylinder.R * cos(cylinder.dzeta)
                Yinf_c = -Z * tan(cylinder.dzeta);
                Ysup_c = (Z * cos(cylinder.dzeta) + cylinder.R) / sin(cylinder.dzeta);
            elseif cylinder.R * cos(cylinder.dzeta) <= Z && Z <= cylinder.L * sin(cylinder.dzeta) - cylinder.R * cos(cylinder.dzeta)
                Yinf_c = (Z * cos(cylinder.dzeta) - cylinder.R) / sin(cylinder.dzeta);
                Ysup_c = (Z * cos(cylinder.dzeta) + cylinder.R) / sin(cylinder.dzeta);
            elseif cylinder.dzeta ~= pi / 2 && cylinder.L * sin(cylinder.dzeta) - cylinder.R * cos(cylinder.dzeta) < Z && Z <= cylinder.L * sin(cylinder.dzeta) + cylinder.R * cos(cylinder.dzeta)
                Yinf_c = (Z * cos(cylinder.dzeta) - cylinder.R) / sin(cylinder.dzeta);
                Ysup_c = (cylinder.L - Z * sin(cylinder.dzeta)) / cos(cylinder.dzeta);
            end
            
            % Limites Y
            Yinf = max(Yinf_c, Yinf_e) + param.Tol;
            Ysup = min(Ysup_c, Ysup_e) - param.Tol;
            
            % Boucle sur Y
            
            for j = 0:100
                
                Y = Yinf + (Ysup - Yinf) * j / 100;
                
                % Contact impossible
                if Yinf >= Ysup
                    dis = 10 ^ 20;
                    
                    % Contact possible
                else
                    
                    % Calcul Xc- et Xc+
                    Xc2 = 1 - (Y - Z * cot(cylinder.dzeta)) ^ 2 / (cylinder.R / sin(cylinder.dzeta)) ^ 2;
                    Xc1 = -cylinder.R * (1 - (Y - Z * cot(cylinder.dzeta)) ^ 2 / (cylinder.R / sin(cylinder.dzeta)) ^ 2) ^ 0.5;
                    Xc2 = cylinder.R * (1 - (Y - Z * cot(cylinder.dzeta)) ^ 2 / (cylinder.R / sin(cylinder.dzeta)) ^ 2) ^ 0.5;
                    
                    % Calcul distances
                    if ship.alpha >= 0
                        
                        if Y < YS - RY * cos(ship.alpha)
                            Xe1 = (Y - YS) * (RX ^ 2 - RY ^ 2) * sin(ship.alpha) * cos(ship.alpha);
                            Xe1 = Xe1 - RX * RY * (RX ^ 2 * sin(ship.alpha) ^ 2 + RY ^ 2 * cos(ship.alpha) ^ 2 - (Y - YS) ^ 2) ^ 0.5;
                            Xe1 = Xe1 / (RX ^ 2 * sin(ship.alpha) ^ 2 + RY ^ 2 * cos(ship.alpha) ^ 2) + XS;
                            Xe2 = (Y - YS) * (RX ^ 2 - RY ^ 2) * sin(ship.alpha) * cos(ship.alpha);
                            Xe2 = Xe2 + RX * RY * (RX ^ 2 * sin(ship.alpha) ^ 2 + RY ^ 2 * cos(ship.alpha) ^ 2 - (Y - YS) ^ 2) ^ 0.5;
                            Xe2 = Xe2 / (RX ^ 2 * sin(ship.alpha) ^ 2 + RY ^ 2 * cos(ship.alpha) ^ 2) + XS;
                        else
                            Xe1 = (Y - YS) * (RX ^ 2 - RY ^ 2) * sin(ship.alpha) * cos(ship.alpha);
                            Xe1 = Xe1 - RX * RY * (RX ^ 2 * sin(ship.alpha) ^ 2 + RY ^ 2 * cos(ship.alpha) ^ 2 - (Y - YS) ^ 2) ^ 0.5;
                            Xe1 = Xe1 / (RX ^ 2 * sin(ship.alpha) ^ 2 + RY ^ 2 * cos(ship.alpha) ^ 2) + XS;
                            if ship.alpha ~= pi / 2
                                Xe2 = XS + (Y - YS) * tan(ship.alpha) + 2 * RY * sin(ship.alpha);
                            else
                                Xe2 = XS + RY;
                            end
                        end
                        
                    else
                        
                        if YS + RY * cos(ship.alpha) < Y
                            Xe1 = (Y - YS) * (RX ^ 2 - RY ^ 2) * sin(ship.alpha) * cos(ship.alpha);
                            Xe1 = Xe1 - RX * RY * (RX ^ 2 * sin(ship.alpha) ^ 2 + RY ^ 2 * cos(ship.alpha) ^ 2 - (Y - YS) ^ 2) ^ 0.5;
                            Xe1 = Xe1 / (RX ^ 2 * sin(ship.alpha) ^ 2 + RY ^ 2 * cos(ship.alpha) ^ 2) + XS;
                            Xe2 = (Y - YS) * (RX ^ 2 - RY ^ 2) * sin(ship.alpha) * cos(ship.alpha);
                            Xe2 = Xe2 + RX * RY * (RX ^ 2 * sin(ship.alpha) ^ 2 + RY ^ 2 * cos(ship.alpha) ^ 2 - (Y - YS) ^ 2) ^ 0.5;
                            Xe2 = Xe2 / (RX ^ 2 * sin(ship.alpha) ^ 2 + RY ^ 2 * cos(ship.alpha) ^ 2) + XS;
                        else
                            Xe1 = (Y - YS) * (RX ^ 2 - RY ^ 2) * sin(ship.alpha) * cos(ship.alpha);
                            Xe1 = Xe1 - RX * RY * (RX ^ 2 * sin(ship.alpha) ^ 2 + RY ^ 2 * cos(ship.alpha) ^ 2 - (Y - YS) ^ 2) ^ 0.5;
                            Xe1 = Xe1 / (RX ^ 2 * sin(ship.alpha) ^ 2 + RY ^ 2 * cos(ship.alpha) ^ 2) + XS;
                            if ship.alpha ~= -pi / 2
                                Xe2 = XS + (Y - YS) * tan(ship.ship.alpha) - 2 * RY * sin(ship.ship.alpha);
                            elseif ship.ship.alpha == -pi / 2
                                Xe2 = XS + RY;
                            end
                        end
                        
                    end
                    
                    dis = min([abs(Xe1 - Xc1), abs(Xe2 - Xc2), abs(Xe1 - Xc2), abs(Xe2 - Xc1)]);
                    
                    % Test intrusion
                    if Xe1 <= Xc1 && Xc1 <= Xe2
                        cont = 1;
                    elseif Xe1 <= Xc2 && Xc2 <= Xe2
                        cont = 1;
                    end
                    
                end
                
                % Test contact
                if dis < dis_min
                    if dis == abs(Xe1 - Xc1) || dis == abs(Xe2 - Xc1)
                        XI = Xc1;
                    elseif dis == abs(Xe1 - Xc2) || dis == abs(Xe2 - Xc2)
                        XI = Xc2;
                    end
                    YI = Y;
                    ZI = Z;
                    dis_min = dis;
                end
                
                if cont == 1
                    break
                end
                
            end
            
            if cont == 1
                break
            end
            
        end
        
    end

end