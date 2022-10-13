punchProp.nodeProp = [0.625 34.428 0.05 255000000 0.594 1 0.65 20.04 14.38 1 0 1.5 4.3378000000];
i = 1;

% delta = 0.378;
% delta = 0.3838;
% delta = 0.35;
delta = 0;

% ddelta = 0.0057;
ddelta = 0;

punchProp.punchModeNodes = [0 0 0];

punchProp.CP = 0;
punchProp.dP = 0;

resultPunch = solvePunch(punchProp,i,delta,ddelta);

resultPunch.Ptotpunch