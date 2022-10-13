function[d] = limit(cylinder,ship)

R1 = ship.q + (min(ship.ZS, cylinder.L * sin(cylinder.dzeta) + cylinder.R * cos(cylinder.dzeta)) - ship.ZS) / tan(ship.phib);
R2 = ship.p + (min(ship.ZS, cylinder.L * sin(cylinder.dzeta) + cylinder.R * cos(cylinder.dzeta)) - ship.ZS) / tan(ship.psib);

% -pi/2 <= alpha < -pi /4

if ship.alpha >= -pi / 2 && ship.alpha < -pi / 4
    if ship.alpha ~= -pi / 2
        dmin1 = -(cylinder.R + ship.XP - R2 * sin(ship.alpha)) / cos(ship.alpha);
        dmin2 = (cylinder.L * cos(cylinder.dzeta) + cylinder.R * sin(cylinder.dzeta) + R2 * cos(ship.alpha)) / sin(ship.alpha);
        dmax1 = (cylinder.R + (R1 ^ 2 * cos(ship.alpha) ^ 2 + R2 ^ 2 * sin(ship.alpha) ^ 2) ^ 0.5 - ship.XP) / cos(ship.alpha);
        dmax2 = -(cylinder.R * sin(cylinder.dzeta) + (R1 ^ 2 * sin(ship.alpha) ^ 2 + R2 ^ 2 * cos(ship.alpha) ^ 2) ^ 0.5) / sin(ship.alpha);
    else
        dmin = -(cylinder.L * cos(cylinder.dzeta) + cylinder.R * sin(cylinder.dzeta));
        dmax = cylinder.R * sin(cylinder.dzeta) + R1;
    end
end

% -pi/4 <= alpha < 0

if ship.alpha >= -pi / 4 && ship.alpha < 0
    dmin1 = (R2 * sin(ship.alpha) - cylinder.R) / cos(ship.alpha);
    dmin2 = (cylinder.L * cos(cylinder.dzeta) + cylinder.R * sin(cylinder.dzeta) + R2 * cos(ship.alpha) - ship.YP) / sin(ship.alpha);
    dmax1 = (cylinder.R + (R1 ^ 2 * cos(ship.alpha) ^ 2 + R2 ^ 2 * sin(ship.alpha) ^ 2) ^ 0.5) / cos(ship.alpha);
    dmax2 = -(cylinder.R * sin(cylinder.dzeta) + ship.YP + (R1 ^ 2 * sin(ship.alpha) ^ 2 + R2 ^ 2 * cos(ship.alpha) ^ 2) ^ 0.5) / sin(ship.alpha);
end

% 0 <= alpha <= pi/4

if ship.alpha >= 0 && ship.alpha <= pi / 4
    if ship.alpha ~= 0
        dmin1 = -(R2 * sin(ship.alpha) - cylinder.R) / cos(ship.alpha);
        dmin2 = -(R2 * cos(ship.alpha) + cylinder.R * sin(cylinder.dzeta) + ship.YP) / sin(ship.alpha);
        dmax1 = (cylinder.R + (R1 ^ 2 * cos(ship.alpha) ^ 2 + R2 ^ 2 * sin(ship.alpha) ^ 2) ^ 0.5) / cos(ship.alpha);
        dmax2 = (cylinder.L * cos(cylinder.dzeta) + cylinder.R * sin(cylinder.dzeta) - ship.YP + (R1 ^ 2 * sin(ship.alpha) ^ 2 + R2 ^ 2 * cos(ship.alpha) ^ 2) ^ 0.5) / sin(ship.alpha);
    else
        dmin = -cylinder.R;
        dmax = cylinder.R + R1;
    end
end

% pi/4 < alpha <= pi /2

if ship.alpha > pi / 4 && ship.alpha <= pi / 2
    if ship.alpha ~= pi / 2
        dmin1 = -(cylinder.R + R2 * sin(ship.alpha) + ship.XP) / cos(ship.alpha);
        dmin2 = -(R2 * cos(ship.alpha) + cylinder.R * sin(cylinder.dzeta)) / sin(ship.alpha);
        dmax1 = (cylinder.R + (R1 ^ 2 * cos(ship.alpha) ^ 2 + R2 ^ 2 * sin(ship.alpha) ^ 2) ^ 0.5 - ship.XP) / cos(ship.alpha);
        dmax2 = (cylinder.L * cos(cylinder.dzeta) + cylinder.R * sin(cylinder.dzeta) + (R1 ^ 2 * sin(ship.alpha) ^ 2 + R2 ^ 2 * cos(ship.alpha) ^ 2) ^ 0.5) / sin(ship.alpha);
    else
        dmin = -cylinder.R * sin(cylinder.dzeta);
        dmax = cylinder.L * cos(cylinder.dzeta) + cylinder.R * sin(cylinder.dzeta) + R1;
    end
end

% Choix de dmin et dmax

if ship.alpha ~= -pi / 2 && ship.alpha ~= 0 && ship.alpha ~= pi / 2
    if abs(dmax1 - dmin1) <= abs(dmax2 - dmin2)
        dmax = dmax1;
        dmin = dmin1;
    else
        dmax = dmax2;
        dmin = dmin2;
    end
end

d(1) = dmin;
d(2) = dmax;