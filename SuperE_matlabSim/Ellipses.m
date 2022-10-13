function[contactEllipses,solveX,solveY] = Ellipses(dataCyl)

nbPoints = 100; %Points to define the resolution of the elipses

%% Ellipse 1

a1 = dataCyl.a1;
b1 = dataCyl.b1;

shiftx1 = dataCyl.shiftx1;
shifty1 = dataCyl.shifty1;
alphaDeg1 = dataCyl.alphaDeg1;

x1Vec = -a1:2*a1/nbPoints:a1;
y1Vec = zeros(1,length(x1Vec));

for i = 1:nbPoints+1
    
    x = x1Vec(i);
    y1Vec(i) = sqrt((1-(x/a1)^2)*b1^2);
    
end

x1Vec = [x1Vec x1Vec];
y1Vec = [y1Vec -y1Vec];

alpha1 = alphaDeg1*pi/180;

X1Vec = zeros(1,length(x1Vec));
Y1Vec = zeros(1,length(y1Vec));

for i = 1:length(X1Vec)
    
    X1Vec(i) = x1Vec(i)*cos(alpha1) - y1Vec(i)*sin(alpha1);
    Y1Vec(i) = x1Vec(i)*sin(alpha1) + y1Vec(i)*cos(alpha1);
    
end

X1Vec = X1Vec+shiftx1;
Y1Vec = Y1Vec+shifty1;

%% Ellipse 2

a2 = dataCyl.a2;
b2 = dataCyl.b2;

shiftx2 = dataCyl.shiftx2;
shifty2 = dataCyl.shifty2;
alphaDeg2 = dataCyl.alphaDeg2;

x2Vec = -a2:2*a2/nbPoints:a2;
y2Vec = zeros(1,length(x2Vec));

for i = 1:nbPoints+1
    
    x = x2Vec(i);
    y2Vec(i) = sqrt((1-(x/a2)^2)*b2^2);
    
end

x2Vec = [x2Vec x2Vec];
y2Vec = [y2Vec -y2Vec];

alpha2 = alphaDeg2*pi/180;

X2Vec = zeros(1,length(x2Vec));
Y2Vec = zeros(1,length(y2Vec));

for i = 1:length(X2Vec)
    
    X2Vec(i) = x2Vec(i)*cos(alpha2) - y2Vec(i)*sin(alpha2);
    Y2Vec(i) = x2Vec(i)*sin(alpha2) + y2Vec(i)*cos(alpha2);
    
end

X2Vec = X2Vec+shiftx2;
Y2Vec = Y2Vec+shifty2;

%% Affichage

% figure;
% hold on
% plot(x1Vec,y1Vec,'b')
% plot(X1Vec,Y1Vec,'r')
% plot(X2Vec,Y2Vec,'m')
% axis equal
% grid on
% hold off

%% Résolution

% solveX = shiftx1 - a1*cos(alpha2);
% solveY = shifty1 - a1*sin(alpha2);
% 
% for i = 1:20
%     
%     f1 = (((solveX-shiftx1)*cos(alpha1)+(solveY-shifty1)*sin(alpha1))/a1)^2 + ((-(solveX-shiftx1)*sin(alpha1)+(solveY-shifty1)*cos(alpha1))/b1)^2 - 1;
%     f2 = (((solveX-shiftx2)*cos(alpha2)+(solveY-shifty2)*sin(alpha2))/a2)^2 + ((-(solveX-shiftx2)*sin(alpha2)+(solveY-shifty2)*cos(alpha2))/b2)^2 - 1;
%     
%     df1dX = 2*(((solveX-shiftx1)*cos(alpha1)+(solveY-shifty1)*sin(alpha1))/a1) * cos(alpha1)/a1 - 2*((-(solveX-shiftx1)*sin(alpha1)+(solveY-shifty1)*cos(alpha1))/b1) * sin(alpha1)/b1;
%     df2dX = 2*(((solveX-shiftx2)*cos(alpha2)+(solveY-shifty2)*sin(alpha2))/a2) * cos(alpha2)/a2 - 2*((-(solveX-shiftx2)*sin(alpha2)+(solveY-shifty2)*cos(alpha2))/b2) * sin(alpha2)/b2;
%     df1dY = 2*(((solveX-shiftx1)*cos(alpha1)+(solveY-shifty1)*sin(alpha1))/a1) * sin(alpha1)/a1 + 2*((-(solveX-shiftx1)*sin(alpha1)+(solveY-shifty1)*cos(alpha1))/b1) * cos(alpha1)/b1;
%     df2dY = 2*(((solveX-shiftx2)*cos(alpha2)+(solveY-shifty2)*sin(alpha2))/a2) * sin(alpha2)/a2 + 2*((-(solveX-shiftx2)*sin(alpha2)+(solveY-shifty2)*cos(alpha2))/b2) * cos(alpha2)/b2;
%     
%     matrice = [df1dX df1dY;
%         df2dX df2dY];
%     vecteur = [f1 ; f2];
%     
%     resultat = matrice\vecteur;
%     
%     solveX = solveX - resultat(1);
%     solveY = solveY - resultat(2);
% 
% end

solveX = shiftx2 - a2*cos(alpha2);
solveY = shifty2 - a2*sin(alpha2);

for i = 1:20
    
    f1 = (((solveX-shiftx1)*cos(alpha1)+(solveY-shifty1)*sin(alpha1))/a1)^2 + ((-(solveX-shiftx1)*sin(alpha1)+(solveY-shifty1)*cos(alpha1))/b1)^2 - 1;
    f2 = (((solveX-shiftx2)*cos(alpha2)+(solveY-shifty2)*sin(alpha2))/a2)^2 + ((-(solveX-shiftx2)*sin(alpha2)+(solveY-shifty2)*cos(alpha2))/b2)^2 - 1;
    
    df1dX = 2*(((solveX-shiftx1)*cos(alpha1)+(solveY-shifty1)*sin(alpha1))/a1) * cos(alpha1)/a1 - 2*((-(solveX-shiftx1)*sin(alpha1)+(solveY-shifty1)*cos(alpha1))/b1) * sin(alpha1)/b1;
    df2dX = 2*(((solveX-shiftx2)*cos(alpha2)+(solveY-shifty2)*sin(alpha2))/a2) * cos(alpha2)/a2 - 2*((-(solveX-shiftx2)*sin(alpha2)+(solveY-shifty2)*cos(alpha2))/b2) * sin(alpha2)/b2;
    df1dY = 2*(((solveX-shiftx1)*cos(alpha1)+(solveY-shifty1)*sin(alpha1))/a1) * sin(alpha1)/a1 + 2*((-(solveX-shiftx1)*sin(alpha1)+(solveY-shifty1)*cos(alpha1))/b1) * cos(alpha1)/b1;
    df2dY = 2*(((solveX-shiftx2)*cos(alpha2)+(solveY-shifty2)*sin(alpha2))/a2) * sin(alpha2)/a2 + 2*((-(solveX-shiftx2)*sin(alpha2)+(solveY-shifty2)*cos(alpha2))/b2) * cos(alpha2)/b2;
    
    matrice = [df1dX df1dY;
        df2dX df2dY];
    vecteur = [f1 ; f2];
    
    check = rcond(matrice);
    
    if check < 10^(-12)
        
        solveX = shiftx2 - a2*cos(alpha2) + 0.1;
        solveY = shifty2 - a2*sin(alpha2);
        
%         disp('je change')
        
    else
        
        resultat = matrice\vecteur;
        
        solveX = solveX - resultat(1);
        solveY = solveY - resultat(2);
        
    end

end

%% Résultats

if abs(f1) < 10^(-6) && abs(f2) < 10^(-6)
    contactEllipses = 1;
else
    contactEllipses = 0;
end