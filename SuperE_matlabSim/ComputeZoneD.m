function[M] = ComputeZoneD(data,theta)

Mpl = 4*data.fy*((data.DeLeg/2)^3 - (data.DiLeg/2)^3)/3;

%% Read data

fy = data.fy;
D = data.DLeg;
R = data.RLeg;
t = data.tLeg;

m0 = fy*t^2/4;
n0 = fy*t;

H = sqrt(D*t)/1.5;

thetamax = atan(2*H/R)/2;
nbDtheta = 500;
dtheta = thetamax/nbDtheta;


%% Computation

if theta == 0
    
    M = Mpl;
    
else
    
    betamax = pi/2;
    nbDbeta = 100;
    dbeta = betamax/nbDbeta;
    
    % Circle 1
    beta = 0;
    dalphaDthetaH1Segment = zeros(1,nbDbeta+1);
    for j = 1:nbDbeta
        dalphaDthetaH1Segment(j) = R*cos(beta)/(sqrt(1-(1-R*cos(beta)/(2*H)*tan(theta))^2)) / (2*H * (cos(theta))^2) + cos(beta);
        beta = beta+dbeta;
    end
    EdotH1 = 0;
    for j = 1:nbDbeta
        EdotH1 = EdotH1 + 2*m0*R * (dalphaDthetaH1Segment(j)+dalphaDthetaH1Segment(j+1))*dbeta/2;
    end
    
    if EdotH1 > Mpl
        EdotH1 = Mpl;
    end
    
    % Circle 3
    beta = 0;
    dalphaDthetaH3Segment = zeros(1,nbDbeta+1);
    for j = 1:nbDbeta
        dalphaDthetaH3Segment(j) = R*cos(beta)/(sqrt(1-(1-R*cos(beta)/(2*H)*tan(theta))^2)) / (2*H * (cos(theta))^2);
        beta = beta+dbeta;
    end
    EdotH3 = 0;
    for j = 1:nbDbeta
        EdotH3 = EdotH3 + 2*m0*R * (dalphaDthetaH3Segment(j)+dalphaDthetaH3Segment(j+1))*dbeta/2;
    end
    
    if EdotH3 > Mpl
        EdotH3 = Mpl;
    end
    
    % Circle 2
    
    beta = 0;
    dalphaDthetaH2Segment = zeros(1,nbDbeta+1);
    for j = 1:nbDbeta
        alpha = acos(1 - R*cos(beta)/(2*H)*tan(theta));
        dalphaDbeta = -R*sin(beta)/(sqrt(1-(1-R*cos(beta)/(2*H)*tan(theta))^2)) / (2*H) * tan(theta);
        ds = sqrt((H*cos(alpha)*dalphaDbeta)^2 + (R+H*sin(alpha))^2);
        dalphaDthetaH2Segment(j) = dalphaDthetaH3Segment(j)*ds;
        beta = beta+dbeta;
    end
    EdotH2 = 0;
    for j = 1:nbDbeta
        EdotH2 = EdotH2 + 4*m0 * (dalphaDthetaH2Segment(j)+dalphaDthetaH2Segment(j+1))*dbeta/2;
    end
    
    if EdotH2 > Mpl
        EdotH2 = Mpl;
    end
    
    % Circles 1 - behind
    beta = 0;
    dalphaDthetaH1bSegment = zeros(1,nbDbeta+1);
    for j = 1:nbDbeta
        dalphaDthetaH1bSegment(j) = cos(beta);
        beta = beta+dbeta;
    end
    EdotH1b = 0;
    for j = 1:nbDbeta
        EdotH1b = EdotH1b + 2*m0*R * (dalphaDthetaH1bSegment(j)+dalphaDthetaH1bSegment(j+1))*dbeta/2;
    end
    
    if EdotH1b > Mpl
        EdotH1b = Mpl;
    end
    
    % Membran effect
    
    nbDH = 100;
    dz = H/nbDH;
    
    EmBeforeDerivation = zeros(1,2);
    
    thetaLoop = theta - dtheta;
    for j = 1:2
        
        z = 0;
        EmSegment = zeros(1,nbDH+1);
        for k = 1:nbDH+1
            
            beta = 0;
            LSegment = zeros(1,nbDbeta+1);
            for l = 1:nbDbeta
                alpha = acos(1 - R*cos(beta)/(2*H)*tan(thetaLoop));
                dalphaDbeta = -R*sin(beta)/(sqrt(1-(1-R*cos(beta)/(2*H)*tan(thetaLoop))^2)) / (2*H) * tan(thetaLoop);
                LSegment(l) = sqrt((z*cos(alpha)*dalphaDbeta)^2 + (R+z*sin(alpha))^2);
                
                beta = beta+dbeta;
            end
            LSegment(end) = R;
            LSum = 0;
            for l = 1:nbDbeta
                LSum = LSum + (LSegment(l)+LSegment(l+1))*dbeta/2;
            end
            
            EmSegment(k) = (2*LSum)^2;
            z = z+dz;
        end
        EmIntHSum = 0;
        for k = 1:nbDH
            EmIntHSum = EmIntHSum + (EmSegment(k)+EmSegment(k+1))*dz/2;
        end
        
        EmBeforeDerivation(j) = EmIntHSum;
        thetaLoop = thetaLoop + 2*dtheta;
    end
    EdotM = 4*n0/(pi*R) * (EmBeforeDerivation(2)-EmBeforeDerivation(1))/dtheta;
    
    if EdotM > 4*Mpl
        EdotM = 4*Mpl;
    end
    
    % Total
    
    M = EdotH1 + EdotH1b + EdotH2 + EdotH3 + EdotM;
    if M > Mpl
        M = Mpl;
    end
    
end