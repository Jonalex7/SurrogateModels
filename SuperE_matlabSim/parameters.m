function[impact] = parameters(cylinder,ship,param,impact)

%% PARAMETRES POUR UN CYLINDRE VERTICAL

% Repositionnement

% Propriétés de l'impacteur fictif
impact.pv = ship.p + impact.zsI / tan(ship.psib);
impact.qv = ship.q + impact.zsI / tan(ship.phib);
impact.hbv = ship.hb + impact.zsI;

% Coordonnées dans les axes du cylindre
impact.xvI = -impact.qv * cylinder.R * impact.ysI / (ship.p ^ 2 * (ship.p ^ 2 - impact.ysI ^ 2) + ship.q ^ 2 * impact.ysI ^ 2) ^ 0.5;
impact.yvI = impact.YI * cos(cylinder.dzeta) + impact.ZI * sin(cylinder.dzeta);
impact.zvI = impact.pv * cylinder.R * (impact.pv ^ 2 - impact.ysI ^ 2) ^ 0.5 / (impact.pv ^ 2 * (impact.pv ^ 2 - impact.ysI ^ 2) + impact.qv ^ 2 * impact.ysI ^ 2) ^ 0.5;

% Repositionnement du navire
impact.xvS = -impact.ysI * (1 + impact.qv * cylinder.R / (impact.pv ^ 2 * (impact.pv ^ 2 - impact.ysI ^ 2) + impact.qv ^ 2 * impact.ysI ^ 2) ^ 0.5);
if abs(impact.xvS) < param.Tol
    impact.xvS = 0;
end
impact.yvS = impact.yvI;
impact.zvS = (impact.pv ^ 2 - impact.ysI ^ 2) ^ 0.5 * (impact.qv / impact.pv + impact.pv * cylinder.R / (impact.pv ^ 2 * (impact.pv ^ 2 - impact.ysI ^ 2) + impact.qv ^ 2 * impact.ysI ^ 2) ^ 0.5);

% Calcul de beta

if impact.xvS == 0
    impact.betav = 0;

else
    
    % Impact dégénéré à droite
    if impact.zvI == 0 && abs(impact.xvI - cylinder.R) < param.Tol
        impact.betav = pi / 2;
    % Impact dégénéré à gauche
    elseif impact.zvI == 0 && bs(impact.xvI + cylinder.R) < param.Tol
        impact.betav = -pi / 2;
    % Cas général
    else
        impact.betav = atan(impact.xvI / impact.zvI);
    end
    
end

% Calcul de l'enfoncement maximum df

impact.df = maxCrush(cylinder,impact,impact.pv,impact.qv);

% Paramètres de résistance

impact.nv = 1;
impact.psiv0 = 3 * pi / 4;

if impact.xvS > 0 && impact.xvS - impact.pv <= -cylinder.R
    impact.psivf = pi;
elseif impact.xvS > 0 && impact.xvS - impact.pv > -cylinder.R
    impact.psivf = impact.psiv0 - (pi - impact.psiv0) * (impact.xvS - ship.p - cylinder.R) / 2 / cylinder.R;
elseif impact.xvS < 0 && impact.xvS + impact.pv >= cylinder.R
    impact.psivf = pi;
elseif impact.xvS < 0 && impact.xvS + impact.pv < cylinder.R
    impact.psivf = impact.psiv0 + (pi - impact.psiv0) * (impact.xvS + ship.p + cylinder.R) / 2 / cylinder.R;
elseif impact.xvS == 0
    impact.psivf = pi;
end


%% PARAMETRES POUR UN CYLINDRE HORIZONTAL

% Repositionnement

% Calcul des angles
RX = ship.q + impact.zsI * cot(ship.phib);
RY = ship.p + impact.zsI * cot(ship.psib);
gy = RX / RY * impact.ysI / (RY ^ 2 - impact.ysI ^ 2) ^ 0.5;
gz = -(RY ^ 2 - impact.ysI ^ 2) ^ 0.5 / RY * cot(ship.phib) - RX / RY ^ 2 * impact.ysI ^ 2 / (RY ^ 2 - impact.ysI ^ 2) ^ 0.5 * cot(ship.psib);
impact.eta1 = atan(gy);
impact.eta2 = atan(gz / (1 + gy ^ 2) ^ 0.5);

% Coordonnées dans les axes du cylindre
impact.xhI = cylinder.R * cos(impact.eta2);
impact.yhI = impact.YI * cos(cylinder.dzeta) + impact.ZI * sin(cylinder.dzeta);
impact.zhI = -cylinder.R * sin(impact.eta2);

% Repositionnement si eta1 < 0
impact.alphah = abs(impact.eta1);
if impact.eta1 < 0
    impact.yhI = cylinder.L - impact.yhI;
end

% Repositionnement du navire
impact.xhS = impact.xhI - impact.xsI * cos(impact.alphah) + impact.ysI * sin(impact.alphah);
impact.yhS = impact.yhI - impact.xsI * sin(impact.alphah) - impact.ysI * cos(impact.alphah);
impact.zhS = impact.zhI - impact.zsI;

% Calcul de beta

impact.betah = atan(impact.xhI / impact.zhI);

impact.nh = 1;
impact.psih0 = 3 * pi / 4;



%% Calculation of df

    function[df] = maxCrush(cylinder,impact,pv,qv)
        % Calcul de d1
        
        if impact.xvS == 0
            d1 = 2 * cylinder.R;
        
        else
            
            % Impact excentré à droite
            if impact.xvS > 0 && impact.xvS - pv < -cylinder.R
                xv1 = -cylinder.R;
                xv2 = 0;
                Err = 1;
            end
            
            % Situation impossible à droite
            if impact.xvS > 0 && impact.xvS - pv > -cylinder.R
                xvJ = 0;
                zvJ = 0;
                d1 = 10 ^ 20;
                Err = 0;
            end
            
            % Situation tangente en -cylinder.R
            if impact.xvS > 0 && abs(impact.xvS - pv + cylinder.R) < param.Tol
                d1 = impact.zvS;
                xvJ = -cylinder.R;
                zvJ = 0;
                Err = 0;
            end
            
            % Impact dégénéré à droite
            if impact.xvS > 0 && abs(impact.xvS - pv - cylinder.R) < param.Tol
                d1 = 0;
                xvJ = cylinder.R;
                zvJ = 0;
                Err = 0;
            end
            
            % Impact excentré à gauche
            if impact.xvS < 0 && impact.xvS + pv > cylinder.R
                xv1 = 0;
                xv2 = cylinder.R;
                Err = 1;
            end
            
            % Situation impossible à gauche
            if impact.xvS < 0 && impact.xvS + pv < cylinder.R
                xvJ = 0;
                zvJ = 0;
                d1 = 10 ^ 20;
                Err = 0;
            end
            
            % Situation tangente en R
            if abs(impact.xvS + pv - cylinder.R) < param.Tol
                d1 = impact.zvS;
                xvJ = cylinder.R;
                zvJ = 0;
                Err = 0;
            end
            
            % Impact dégénéré à gauche
            if impact.xvS < 0 && abs(impact.xvS + pv + cylinder.R) < param.Tol
                d1 = 0;
                xvJ = -cylinder.R;
                zvJ = 0;
                Err = 0;
            end
            
            while Err == 1 && abs(xv1 - xv2) > param.Tol
                
                xv0 = (xv1 + xv2) / 2;
                f0 = xv0 * (pv ^ 2 - (xv0 - impact.xvS) ^ 2) ^ 0.5 + qv / pv * (impact.xvS - xv0) * (cylinder.R ^ 2 - xv0 ^ 2) ^ 0.5;
                f1 = xv1 * (pv ^ 2 - (xv1 - impact.xvS) ^ 2) ^ 0.5 + qv / pv * (impact.xvS - xv1) * (cylinder.R ^ 2 - xv1 ^ 2) ^ 0.5;
                f2 = xv2 * (pv ^ 2 - (xv2 - impact.xvS) ^ 2) ^ 0.5 + qv / pv * (impact.xvS - xv2) * (cylinder.R ^ 2 - xv2 ^ 2) ^ 0.5;
                
                if f0 > 0
                    if f1 > 0 && f2 < 0
                        xv1 = xv0;
                    elseif f1 < 0 && f2 > 0
                        xv2 = xv0;
                    elseif f1 >= 0 && f2 >= 0
                        Err = 0;
                    elseif f1 <= 0 && f2 <= 0
                        Err = 0;
                    end
                else
                    if f1 > 0 && f2 < 0
                        xv2 = xv0;
                    elseif f1 < 0 && f2 > 0
                        xv1 = xv0;
                    elseif f1 >= 0 && f2 >= 0
                        Err = 0;
                    elseif f1 <= 0 && f2 <= 0
                        Err = 0;
                    end
                end
                
            end
            
            if Err == 1
                xvJ = xv0;
                zvJ = -(cylinder.R ^ 2 - xvJ ^ 2) ^ 0.5;
                d1 = impact.zvS - zvJ - qv / pv * (pv ^ 2 - (impact.xvS - xvJ) ^ 2) ^ 0.5;
            end
            
        end
        
        % Calcul de d2
        
        if impact.xvS ~= 0
            d2 = (qv ^ 2 + pv ^ 2 * cot(impact.betav) ^ 2) ^ 0.5 - (1 - qv ^ 2 / pv ^ 2) * (impact.xvS - impact.xvI) * cot(impact.betav);
        else
            d2 = 10 ^ 20;
        end
        
        % Calcul de d3
        
        if impact.xvS == 0
            d3 = 10 ^ 20;
        elseif impact.xvS > 0
            d3 = impact.zvS - (impact.xvS - pv) * cot(impact.betav);
        elseif impact.xvS < 0
            d3 = impact.zvS - (impact.xvS + pv) * cot(impact.betav);
        end
        
        % Calcul de df
        
        df = min([d1 d2 d3]);
        
    end
end