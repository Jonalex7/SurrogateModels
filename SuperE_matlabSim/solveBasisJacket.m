function[resultBasisJacket] = solveBasisJacket(basisJacket,force,data,prevDispBasis)

ddeltaprimeddelta = 1;%hstar/H;

%% Axis X

delta = (prevDispBasis+norm(force.valueDisp))*cos(data.shipTrajectory*pi/180);
deltaprime = delta*ddeltaprimeddelta;

thetaLeft = asin(cos(basisJacket.mu)*deltaprime/(sqrt((basisJacket.h)^2+deltaprime^2+2*basisJacket.h*deltaprime*sin(basisJacket.mu))));
thetaRight = asin(cos(basisJacket.mu)*deltaprime/(sqrt((basisJacket.h)^2+deltaprime^2-2*basisJacket.h*deltaprime*sin(basisJacket.mu))));

dthetaLeftddeltaprime = cos(basisJacket.mu)/sqrt((basisJacket.h)^2+deltaprime^2+2*basisJacket.h*deltaprime*sin(basisJacket.mu)-(cos(basisJacket.mu)*deltaprime)^2)...
    * (1 - deltaprime*(deltaprime+basisJacket.h*sin(basisJacket.mu))/((basisJacket.h)^2+deltaprime^2+2*basisJacket.h*deltaprime*sin(basisJacket.mu)));
dthetaRightddeltaprime = cos(basisJacket.mu)/sqrt((basisJacket.h)^2+deltaprime^2-2*basisJacket.h*deltaprime*sin(basisJacket.mu)-(cos(basisJacket.mu)*deltaprime)^2)...
    * (1 - deltaprime*(deltaprime-basisJacket.h*sin(basisJacket.mu))/((basisJacket.h)^2+deltaprime^2-2*basisJacket.h*deltaprime*sin(basisJacket.mu)));

% -------
% Traction left foot
% -------

hprime = sqrt((basisJacket.h)^2 + deltaprime^2 + 2*basisJacket.h*deltaprime*sin(basisJacket.mu));
Npl = basisJacket.ALeg * basisJacket.fy;
FTrac = Npl * (deltaprime+basisJacket.h*sin(basisJacket.mu))/hprime * ddeltaprimeddelta * 10^(-6);

% -------
% Rotation left foot
% -------

n0 = basisJacket.fy * basisJacket.tLeg;
MRotLeftLeg = 4*n0 * (basisJacket.RLeg)^2 * cos(thetaLeft) / (basisJacket.h * cos(basisJacket.mu));
FRotLeftLeg = MRotLeftLeg * dthetaLeftddeltaprime * ddeltaprimeddelta * 10^(-6);

% -------
% Zone A
% -------

Rp = basisJacket.RBrace/cos(basisJacket.mu);
deltaA = Rp * sin(thetaLeft);
MLocZoneA = ComputeZoneAB(basisJacket,Rp,basisJacket.h,thetaLeft,deltaA);
FZoneA = MLocZoneA * dthetaLeftddeltaprime * ddeltaprimeddelta * 10^(-6);

% -------
% Zone B
% -------

Rp = basisJacket.RBrace/cos(basisJacket.mu);
deltaB = Rp * sin(thetaRight);
MLocZoneB = ComputeZoneAB(basisJacket,Rp,basisJacket.h,thetaRight,deltaB);
FZoneB = MLocZoneB * dthetaRightddeltaprime * ddeltaprimeddelta * 10^(-6);

% -------
% Zone C
% -------

Mpl = 4*basisJacket.fy*((basisJacket.DeLeg/2)^3 - (basisJacket.DiLeg/2)^3)/3;
FZoneC = Mpl * dthetaLeftddeltaprime * ddeltaprimeddelta * 10^(-6);

% -------
% Zone D
% -------

MLocZoneD = ComputeZoneD(basisJacket,thetaRight);
FZoneD = MLocZoneD * dthetaRightddeltaprime * ddeltaprimeddelta * 10^(-6);

% -------
% Total
% -------

FBasisJacket = FZoneA + FZoneB + FZoneC + FZoneD + FTrac + FRotLeftLeg;
resultBasisJacket.FAxisX = FBasisJacket;

%% Axis Y

delta = (prevDispBasis+norm(force.valueDisp))*sin(data.shipTrajectory*pi/180);
deltaprime = delta*ddeltaprimeddelta;

thetaLeft = asin(cos(basisJacket.mu)*deltaprime/(sqrt((basisJacket.h)^2+deltaprime^2+2*basisJacket.h*deltaprime*sin(basisJacket.mu))));
thetaRight = asin(cos(basisJacket.mu)*deltaprime/(sqrt((basisJacket.h)^2+deltaprime^2-2*basisJacket.h*deltaprime*sin(basisJacket.mu))));

dthetaLeftddeltaprime = cos(basisJacket.mu)/sqrt((basisJacket.h)^2+deltaprime^2+2*basisJacket.h*deltaprime*sin(basisJacket.mu)-(cos(basisJacket.mu)*deltaprime)^2)...
    * (1 - deltaprime*(deltaprime+basisJacket.h*sin(basisJacket.mu))/((basisJacket.h)^2+deltaprime^2+2*basisJacket.h*deltaprime*sin(basisJacket.mu)));
dthetaRightddeltaprime = cos(basisJacket.mu)/sqrt((basisJacket.h)^2+deltaprime^2-2*basisJacket.h*deltaprime*sin(basisJacket.mu)-(cos(basisJacket.mu)*deltaprime)^2)...
    * (1 - deltaprime*(deltaprime-basisJacket.h*sin(basisJacket.mu))/((basisJacket.h)^2+deltaprime^2-2*basisJacket.h*deltaprime*sin(basisJacket.mu)));

% -------
% Traction left foot
% -------

hprime = sqrt((basisJacket.h)^2 + deltaprime^2 + 2*basisJacket.h*deltaprime*sin(basisJacket.mu));
Npl = basisJacket.ALeg * basisJacket.fy;
FTrac = Npl * (deltaprime+basisJacket.h*sin(basisJacket.mu))/hprime * ddeltaprimeddelta * 10^(-6);

% -------
% Rotation left foot
% -------

n0 = basisJacket.fy * basisJacket.tLeg;
MRotLeftLeg = 4*n0 * (basisJacket.RLeg)^2 * cos(thetaLeft) / (basisJacket.h * cos(basisJacket.mu));
FRotLeftLeg = MRotLeftLeg * dthetaLeftddeltaprime * ddeltaprimeddelta * 10^(-6);

% -------
% Zone A
% -------

Rp = basisJacket.RBrace/cos(basisJacket.mu);
deltaA = Rp * sin(thetaLeft);
MLocZoneA = ComputeZoneAB(basisJacket,Rp,basisJacket.h,thetaLeft,deltaA);
FZoneA = MLocZoneA * dthetaLeftddeltaprime * ddeltaprimeddelta * 10^(-6);

% -------
% Zone B
% -------

Rp = basisJacket.RBrace/cos(basisJacket.mu);
deltaB = Rp * sin(thetaRight);
MLocZoneB = ComputeZoneAB(basisJacket,Rp,basisJacket.h,thetaRight,deltaB);
FZoneB = MLocZoneB * dthetaRightddeltaprime * ddeltaprimeddelta * 10^(-6);

% -------
% Zone C
% -------

Mpl = 4*basisJacket.fy*((basisJacket.DeLeg/2)^3 - (basisJacket.DiLeg/2)^3)/3;
FZoneC = Mpl * dthetaLeftddeltaprime * ddeltaprimeddelta * 10^(-6);

% -------
% Zone D
% -------

MLocZoneD = ComputeZoneD(basisJacket,thetaRight);
FZoneD = MLocZoneD * dthetaRightddeltaprime * ddeltaprimeddelta * 10^(-6);

% -------
% Total
% -------

FBasisJacket = FZoneA + FZoneB + FZoneC + FZoneD + FTrac + FRotLeftLeg;
resultBasisJacket.FAxisY = FBasisJacket;