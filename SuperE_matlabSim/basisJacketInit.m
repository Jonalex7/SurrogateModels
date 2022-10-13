function[basisJacket] = basisJacketInit(properties,data)

basisJacket.DeLeg = data.De(1);
basisJacket.tLeg = data.t(1);
basisJacket.DLeg = basisJacket.DeLeg - basisJacket.tLeg;
basisJacket.DiLeg = basisJacket.DeLeg - 2*basisJacket.tLeg;
basisJacket.RLeg = basisJacket.DLeg/2;
basisJacket.ALeg = pi * ((basisJacket.DeLeg/2)^2 - (basisJacket.DiLeg/2)^2);

basisJacket.DeBrace = data.De(4);
basisJacket.tBrace = data.t(4);
basisJacket.DBrace = basisJacket.DeBrace - basisJacket.tBrace;
basisJacket.DiBrace = basisJacket.DeBrace - 2*basisJacket.tBrace;
basisJacket.RBrace = basisJacket.DBrace/2;

basisJacket.fy = data.fy(1);

basisJacket.xTop = data.b(2)/2;
basisJacket.xBottom = data.b(1)/2;
basisJacket.totalHeight = data.h(end);
basisJacket.mu = atan((basisJacket.xBottom - basisJacket.xTop)/basisJacket.totalHeight);

basisJacket.hstar = 2.5;
basisJacket.h = basisJacket.hstar/(cos(basisJacket.mu));
if data.impactedZone == 0 % node
    basisJacket.H = properties.nodes(data.impactedNode,3);
else % element
    basisJacket.H = data.shipHeight;
end

basisJacket.psi0 = 3/4*pi;% Loïc
% basisJacket.psi0 = 1/3*pi;% Jing-Ru
basisJacket.k = 1;
basisJacket.nhp = 1;
basisJacket.alphaP = 0;