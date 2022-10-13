function[deltaFSurf] = normalSurf4(properties,forceFunct,rhoV,i,xiMply,xiMplz)

% stop = 0;

Mpl = properties.Mpl(i);
Npl = properties.Npl(i);

Np = abs(forceFunct(1));
Myp = abs(forceFunct(5));
Mzp = abs(forceFunct(6));
rho = rhoV;

if imag(Np+Myp+Mzp) ~= 0
    disp('imaginaireAvant')
    disp(['MyAvant = ',num2str(xiMply*Myp)])
    disp(['MzAvant = ',num2str(xiMplz*Mzp)])
    disp(['NAvant = ',num2str(Np)])
end

N = Np;
My = Myp;
Mz = Mzp;

% Mn = Mpl*(1 - (N/Npl)^2);
% interactionCrit = (My/((1-rho)*Mn))^2 + (Mz/((1-rho)*Mn))^2;

t = 1;
for boucle = 1:5
    F1 = (My/((1-rho)*xiMply*Mpl*(1 - (N/Npl)^2)))^2 + (Mz/((1-rho)*xiMplz*Mpl*(1 - (N/Npl)^2)))^2 - 1;
    dF1dN = 4*N/(Npl^2*(1-(N/Npl)^2)^3) * ((My/xiMply)^2+(Mz/xiMplz)^2)/(((1-rho)*Mpl)^2);
    dF1dMy = 2*My/(((1-rho)*xiMply*Mpl*(1-(N/Npl)^2))^2);
    dF1dMz = 2*Mz/(((1-rho)*xiMplz*Mpl*(1-(N/Npl)^2))^2);
    dF1dt = 0;
    
    F2 = dF1dN*t + Np - N;
    dF2dN = ((My/xiMply)^2+(Mz/xiMplz)^2)/((Npl*(1-rho)*Mpl)^2) * (4/((1-(N/Npl)^2)^3) + 24*N^2/(Npl^2*(1-(N/Npl)^2)^4))*t - 1;
    dF2dMy = 4*N/(Npl^2*(1-(N/Npl)^2)^3) * 2*My/(((1-rho)*xiMply*Mpl)^2) * t;
    dF2dMz = 4*N/(Npl^2*(1-(N/Npl)^2)^3) * 2*Mz/(((1-rho)*xiMplz*Mpl)^2) * t;
    dF2dt = dF1dN;
    
    F3 = dF1dMy*t + Myp - My;
    dF3dN = 2*My/(((1-rho)*xiMply*Mpl)^2) * 4*N/(Npl^2 * (1-(N/Npl)^2)^3) * t;
    dF3dMy = 2/(((1-rho)*xiMply*Mpl*(1-(N/Npl)^2))^2) * t - 1;
    dF3dMz = 0;
    dF3dt = dF1dMy;
    
    F4 = dF1dMz*t + Mzp - Mz;
    dF4dN = 2*Mz/(((1-rho)*xiMplz*Mpl)^2) * 4*N/(Npl^2 * (1-(N/Npl)^2)^3) * t;
    dF4dMy = 0;
    dF4dMz = 2/(((1-rho)*xiMplz*Mpl*(1-(N/Npl)^2))^2) * t - 1;
    dF4dt = dF1dMz;
    
    jacobienne = [dF1dN dF1dMy dF1dMz dF1dt ;
        dF2dN dF2dMy dF2dMz dF2dt ;
        dF3dN dF3dMy dF3dMz dF3dt ;
        dF4dN dF4dMy dF4dMz dF4dt];
    
    vec = [F1 ; F2 ; F3 ; F4];
    
    deltaSol = jacobienne\vec;
    
    N = N - deltaSol(1);
    My = My - deltaSol(2);
    Mz = Mz - deltaSol(3);
    t = t - deltaSol(4);
    
end

% N
% My
% Mz
% t
% 
% F1 = (My/((1-rho)*Mpl*(1 - (N/Npl)^2)))^2 + (Mz/((1-rho)*Mpl*(1 - (N/Npl)^2)))^2

if imag(deltaSol) ~= 0
    disp('imaginaire')
    disp(['Mpl = ',num2str(Mpl)])
    disp(['Npl = ',num2str(Npl)])
    disp(['MyPoint = ',num2str(My)])
    disp(['MzPoint = ',num2str(Mz)])
    disp(['NPoint = ',num2str(N)])
    disp(jacobienne)
    disp(vec)
    disp(['i = ',num2str(i)])
    abcd = efgh;
%     disp(['contrMx = ',num2str(contrMx)])
%     disp(['contrMy = ',num2str(contrMy)])
%     disp(['contrN = ',num2str(contrN)])
%     disp(['rhoV = ',num2str(rhoV)])
%     stop = 1;
end

deltaFSurf = zeros(6,1);
if forceFunct(1) <= 0
    deltaFSurf(1) = Np - N;
else
    deltaFSurf(1) = -(Np - N);
end

if forceFunct(5) <= 0
    deltaFSurf(5) = Myp - My;
else
    deltaFSurf(5) = -(Myp - My);
end

if forceFunct(6) <= 0
    deltaFSurf(6) = Mzp - Mz;
else
    deltaFSurf(6) = -(Mzp - Mz);
end