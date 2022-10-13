function[matrices,forceVec] = stiffness(solveOption,properties,forceElem,forceVec,plasticity,dofs,punchProp)

K     = zeros(dofs.nbDofs,dofs.nbDofs);
Kelem = zeros(12,12,dofs.nbElements);
T     = zeros(12,12,dofs.nbElements);

for i = 1:dofs.nbElements
    %% Correction for fictive elements
    
    if properties.type(i) == 5;
        forceElem(:,i) = zeros(12,1);
    end
    
    %% Properties
    
    E = properties.E(i);
    G = E/2.6;
    A = properties.A(i);
    Av = 0.6*A;
    I = properties.I(i);
    J = properties.J(i);
    L = properties.L(i);
%     e0d = properties.defoInit(i);
    rho1 = properties.rho(i,1);
    rho2 = properties.rho(i,2);
    
    F = forceElem(:,i);
    
%     w = 8*F(1)*e0d/L;
    
    %% Effect of punching on inertia
    
    modInertia = 0;
    
    % Reduction of stiffness based on plastic moment reduction
    aa = punchProp.redSectStar(properties.nodein(i),1); % axis X
    bb = punchProp.redSectStar(properties.nodeout(i),1);
    cc = punchProp.redSectStar(properties.nodein(i),2); % axis Y
    dd = punchProp.redSectStar(properties.nodeout(i),2);
    
    % Reduction of stiffness based on inertia reduction
    a = punchProp.redSectStar(properties.nodein(i),3); % axis X
    b = punchProp.redSectStar(properties.nodeout(i),3);
    c = punchProp.redSectStar(properties.nodein(i),4); % axis Y
    d = punchProp.redSectStar(properties.nodeout(i),4);
    
    if a ~= 1 || b ~= 1
        modInertia = 1;
    elseif c ~= 1 || d ~= 1
        modInertia = 1;
    end
    
    if modInertia == 0
        Ix = I;
        Iy = I;
        
        xiMplyLeft = 1;
        xiMplzLeft = 1;
        xiMplyRight = 1;
        xiMplzRight = 1;
        
    elseif modInertia == 1
        
        vecDir1 = properties.angle(1,:,i);
        dir = -vecDir1(1)/vecDir1(3);
        axisXElement = [1 0 dir];
        axisXElement = axisXElement/norm(axisXElement);
        vecDir2 = properties.angle(2,:,i);
        alpha = acos(axisXElement * vecDir2');
        IprimeX = I*(a+b)/2;
        IprimeY = I*(c+d)/2;
        Ix = (IprimeX * abs(cos(alpha)) + IprimeY * abs(sin(alpha))) / (abs(cos(alpha)) + abs(sin(alpha)));
        Iy = (IprimeX * abs(sin(alpha)) + IprimeY * abs(cos(alpha))) / (abs(cos(alpha)) + abs(sin(alpha)));
        
        xiMplyLeft = (aa * abs(cos(alpha)) + cc * abs(sin(alpha))) / (abs(cos(alpha)) + abs(sin(alpha)));
        xiMplzLeft = (aa * abs(sin(alpha)) + cc * abs(cos(alpha))) / (abs(cos(alpha)) + abs(sin(alpha)));
        xiMplyRight = (bb * abs(cos(alpha)) + dd * abs(sin(alpha))) / (abs(cos(alpha)) + abs(sin(alpha)));
        xiMplzRight = (bb * abs(sin(alpha)) + dd * abs(cos(alpha))) / (abs(cos(alpha)) + abs(sin(alpha)));
    end
    
    %% Elementary stiffness matrices (elastic - large displacement)
    
    u = sqrt(abs(F(1)/(E*I/(L^2))));
    if imag(u) ~= 0
        disp('probleme sur u')
    end
    v = u;
    
    if F(1) > 0 %compression
        
        % Axis y
        phi = 12*E*Iy/(G*Av*L^2);
        k55 = (1+phi/4)/(1+phi)*(3*rho1*(1 - rho2)*u^2 + 9*rho1*rho2*(1 - u/tan(u)))...
            / ((1 - rho1)*(1 - rho2)*u^2 + 3*(rho1 + rho2 - 2*rho1*rho2)*(1 - u/tan(u)) + 9*rho1*rho2*(tan(u/2)/(u/2) - 1)) * E*Iy/L;
        k115 = (1-phi/2)/(1+phi)*(9*rho1*rho2*(u/sin(u) - 1))...
            / ((1 - rho1)*(1 - rho2)*u^2 + 3*(rho1 + rho2 - 2 *rho1*rho2)*(1 - u/tan(u)) + 9*rho1*rho2*(tan(u/2)/(u/2) - 1)) * E*Iy/L;
        k1111 = (1+phi/4)/(1+phi)*(3*rho2*(1 - rho1)*u^2 + 9*rho1*rho2*(1 - u/tan(u)))...
            / ((1 - rho1)*(1 - rho2)*u^2 + 3*(rho1 + rho2 - 2*rho1*rho2)*(1 - u/tan(u)) + 9*rho1*rho2*(tan(u/2)/(u/2) - 1)) * E*Iy/L;
        k35 = (k55+k115)/L;
        k59 = -k35;
        k113 = (k1111+k115)/L;
        k119 = -(k1111+k115)/L;
        k33 = (k55+k1111+2*k115-F(1)*L)/(L^2);
        k39 = -(k55+k1111+2*k115-F(1)*L)/(L^2);
        k99 = k33;
        
        % Axis x
        phi = 12*E*Ix/(G*Av*L^2);
        k66 = (1+phi/4)/(1+phi)*(3*rho1*(1 - rho2)*u^2 + 9*rho1*rho2*(1 - u/tan(u)))...
            / ((1 - rho1)*(1 - rho2)*u^2 + 3*(rho1 + rho2 - 2*rho1*rho2)*(1 - u/tan(u)) + 9*rho1*rho2*(tan(u/2)/(u/2) - 1)) * E*Ix/L;
        k126 = (1-phi/2)/(1+phi)*(9*rho1*rho2*(u/sin(u) - 1))...
            / ((1 - rho1)*(1 - rho2)*u^2 + 3*(rho1 + rho2 - 2 *rho1*rho2)*(1 - u/tan(u)) + 9*rho1*rho2*(tan(u/2)/(u/2) - 1)) * E*Ix/L;
        k1212 = (1+phi/4)/(1+phi)*(3*rho2*(1 - rho1)*u^2 + 9*rho1*rho2*(1 - u/tan(u)))...
            / ((1 - rho1)*(1 - rho2)*u^2 + 3*(rho1 + rho2 - 2*rho1*rho2)*(1 - u/tan(u)) + 9*rho1*rho2*(tan(u/2)/(u/2) - 1)) * E*Ix/L;
        k26 = (k66+k126)/L;
        k68 = -k26;
        k122 = (k1212+k126)/L;
        k128 = -(k1212+k126)/L;
        k22 = (k66+k1212+2*k126-F(1)*L)/(L^2);
        k28 = -(k66+k1212+2*k126-F(1)*L)/(L^2);
        k88 = k22;
        
        % Compression - torsion
        Hy = v*(F(5)^2 + F(11)^2)*(1/tan(v) + v/(sin(v)^2)) - 2*(F(5) + F(11))^2 + 2*v*F(5)*F(11)*(1 + v/tan(v))/sin(v);
        Hz = v*(F(6)^2 + F(12)^2)*(1/tan(u) + u/(sin(u)^2)) - 2*(F(6) + F(12))^2 + 2*u*F(6)*F(12)*(1 + u/tan(u))/sin(u);
        s1 = 1 / (1 + (Hy+Hz)*E*A/(4*F(1)^3*L^2));
        k11 = s1*E*A/L;
        k44 = G*J/L;
        
    elseif F(1) == 0
        
        % Axis y
        phi = 12*E*Iy/(G*Av*L^2);
        k55 = (1+phi/4)/(1+phi)*12*rho1 / (4 - rho1*rho2) * E*Iy/L;
        k115 = (1-phi/2)/(1+phi)*6*rho1*rho2/(4 - rho1*rho2) * E*Iy/L;
        k1111 = (1+phi/4)/(1+phi)*12*rho2 / (4 - rho1*rho2) * E*Iy/L;
        k35 = (k55+k115)/L;
        k59 = -k35;
        k113 = (k1111+k115)/L;
        k119 = -(k1111+k115)/L;
        k33 = (k55+k1111+2*k115-F(1)*L)/(L^2);
        k39 = -(k55+k1111+2*k115-F(1)*L)/(L^2);
        k99 = k33;
        
        % Axis x
        phi = 12*E*Ix/(G*Av*L^2);
        k66 = (1+phi/4)/(1+phi)*12*rho1 / (4 - rho1*rho2) * E*Ix/L;
        k126 = (1-phi/2)/(1+phi)*6*rho1*rho2/(4 - rho1*rho2) * E*Ix/L;
        k1212 = (1+phi/4)/(1+phi)*12*rho2 / (4 - rho1*rho2) * E*Ix/L;
        k26 = (k66+k126)/L;
        k68 = -k26;
        k122 = (k1212+k126)/L;
        k128 = -(k1212+k126)/L;
        k22 = (k66+k1212+2*k126-F(1)*L)/(L^2);
        k28 = -(k66+k1212+2*k126-F(1)*L)/(L^2);
        k88 = k22;
        
        % Compression - torsion
        k11 = E*A/L;
        k44 = G*J/L;
        
    else %tension
        
        % Axis y
        phi = 12*E*Iy/(G*Av*L^2);
        k55 = (1+phi/4)/(1+phi)*(-3*rho1*(1 - rho2)*u^2 + 9*rho1*rho2*(1 - u/tanh(u)))...
            / (-(1 - rho1)*(1 - rho2)*u^2 + 3*(rho1 + rho2 - 2*rho1*rho2)*(1 - u/tanh(u)) + 9*rho1*rho2*(tanh(u/2)/(u/2) - 1)) * E*Iy/L;
        k115 = (1-phi/2)/(1+phi)*(9*rho1*rho2*(u/sinh(u) - 1))...
            / (-(1 - rho1)*(1 - rho2)*u^2 + 3*(rho1 + rho2 - 2 *rho1*rho2)*(1 - u/tanh(u)) + 9*rho1*rho2*(tanh(u/2)/(u/2) - 1)) * E*Iy/L;
        k1111 = (1+phi/4)/(1+phi)*(-3*rho2*(1 - rho1)*u^2 + 9*rho1*rho2*(1 - u/tanh(u)))...
            / (-(1 - rho1)*(1 - rho2)*u^2 + 3*(rho1 + rho2 - 2*rho1*rho2)*(1 - u/tanh(u)) + 9*rho1*rho2*(tanh(u/2)/(u/2) - 1)) * E*Iy/L;
        k35 = (k55+k115)/L;
        k59 = -k35;
        k113 = (k1111+k115)/L;
        k119 = -(k1111+k115)/L;
        k33 = (k55+k1111+2*k115-F(1)*L)/(L^2);
        k39 = -(k55+k1111+2*k115-F(1)*L)/(L^2);
        k99 = k33;
        
        % Axis x
        phi = 12*E*Ix/(G*Av*L^2);
        k66 = (1+phi/4)/(1+phi)*(-3*rho1*(1 - rho2)*u^2 + 9*rho1*rho2*(1 - u/tanh(u)))...
            / (-(1 - rho1)*(1 - rho2)*u^2 + 3*(rho1 + rho2 - 2*rho1*rho2)*(1 - u/tanh(u)) + 9*rho1*rho2*(tanh(u/2)/(u/2) - 1)) * E*Ix/L;
        k126 = (1-phi/2)/(1+phi)*(9*rho1*rho2*(u/sinh(u) - 1))...
            / (-(1 - rho1)*(1 - rho2)*u^2 + 3*(rho1 + rho2 - 2 *rho1*rho2)*(1 - u/tanh(u)) + 9*rho1*rho2*(tanh(u/2)/(u/2) - 1)) * E*Ix/L;
        k1212 = (1+phi/4)/(1+phi)*(-3*rho2*(1 - rho1)*u^2 + 9*rho1*rho2*(1 - u/tanh(u)))...
            / (-(1 - rho1)*(1 - rho2)*u^2 + 3*(rho1 + rho2 - 2*rho1*rho2)*(1 - u/tanh(u)) + 9*rho1*rho2*(tanh(u/2)/(u/2) - 1)) * E*Ix/L;
        k26 = (k66+k126)/L;
        k68 = -k26;
        k122 = (k1212+k126)/L;
        k128 = -(k1212+k126)/L;
        k22 = (k66+k1212+2*k126-F(1)*L)/(L^2);
        k28 = -(k66+k1212+2*k126-F(1)*L)/(L^2);
        k88 = k22;
        
        % Compression - torsion
        Hy = v*(F(5)^2 + F(11)^2)*(1/tanh(v) + v/(sinh(v)^2)) - 2*(F(5) + F(11))^2 + 2*v*F(5)*F(11)*(1 + v/tanh(v))/sinh(v);
        Hz = v*(F(6)^2 + F(12)^2)*(1/tanh(u) + u/(sinh(u)^2)) - 2*(F(6) + F(12))^2 + 2*u*F(6)*F(12)*(1 + u/tanh(u))/sinh(u);
        s1 = 1 / (1 + (Hy+Hz)*E*A/(4*F(1)^3*L^2));
        k11 = s1*E*A/L;
        k44 = G*J/L;
        
    end
    
    Kelem(:,:,i) = [k11   0     0     0     0     0   -k11    0     0     0     0      0;
                     0   k22    0     0     0    k26    0    k28    0     0     0    k122;
                     0    0    k33    0    k35    0     0     0    k39    0   k113     0;
                     0    0     0    k44    0     0     0     0     0   -k44    0      0;
                     0    0    k35    0    k55    0     0     0    k59    0   k115     0;
                     0   k26    0     0     0    k66    0    k68    0     0     0    k126;
                   -k11   0     0     0     0     0    k11    0     0     0     0      0;
                     0   k28    0     0     0    k68    0    k88    0     0     0    k128;
                     0    0    k39    0    k59    0     0     0    k99    0   k119     0;
                     0    0     0   -k44    0     0     0     0     0    k44    0      0;
                     0    0   k113    0   k115    0     0     0   k119    0   k1111    0;
                     0  k122    0     0     0   k126    0   k128    0     0     0    k1212];
    
    %% Rotation matrices
    
    T(1:3,1:3,i) = properties.rotation(:,:,i);
    T(4:6,4:6,i) = properties.rotation(:,:,i);
    T(7:9,7:9,i) = properties.rotation(:,:,i);
    T(10:12,10:12,i) = properties.rotation(:,:,i);
    
    %% Effect of plastic hinges
    
    if solveOption.behaviour == 1 && i ~= 5
        
        % Initialisation
        
        k11 = Kelem(1:6,1:6,i);
        k12 = Kelem(1:6,7:12,i);
        k21 = Kelem(7:12,1:6,i);
        k22 = Kelem(7:12,7:12,i);
        
        k11h = zeros(6,6);
        k22h = zeros(6,6);
        C11 = zeros(6,6);
        C22 = zeros(6,6);
        
        % Gradient of plastic surface
        
        Npl = properties.Npl(i);
        Mpl = properties.Mpl(i);
        Vpl = properties.Vpl(i);
        
        N1 = forceElem(1,i);
        My1 = forceElem(5,i);
        Mz1 = forceElem(6,i);
        
        rhoV1 = 0;
        g10 = 0;
        g11 = 0;
        if sqrt(forceElem(2,i)^2+forceElem(3,i)^2) > 0.5*Vpl
            rhoV1 = (2*sqrt((forceElem(2,i)^2)+(forceElem(3,i)^2))/Vpl - 1)^2;
            if imag(rhoV1) ~= 0
                disp('1')
                ab = forceElem(2,i)
                cd = forceElem(3,i)
                ef = Vpl
                Fdisplay = F
            end
%             g10 = -forceElem(6,i)/(1.04*properties.Mpl(i))*4*(2*forceElem(2,i)/properties.Vpl(i) - 1)/((1-rhoV1)^2);
%             g11 = -forceElem(5,i)/(1.04*properties.Mpl(i))*4*(2*forceElem(3,i)/properties.Vpl(i) - 1)/((1-rhoV1)^2);
        end
        dF1dN1 = 4*N1/(Npl^2*(1-(N1/Npl)^2)^3) * ((My1/xiMplyLeft)^2+(Mz1/xiMplzLeft)^2)/(((1-rhoV1)*Mpl)^2);
        dF1dMy1 = 2*My1/(((1-rhoV1)*xiMplyLeft*Mpl*(1-(N1/Npl)^2))^2);
        dF1dMz1 = 2*Mz1/(((1-rhoV1)*xiMplzLeft*Mpl*(1-(N1/Npl)^2))^2);
        g1 = [dF1dN1 ; g10 ; g11 ; 0 ; dF1dMy1 ; dF1dMz1];
        
        N2 = forceElem(7,i);
        My2 = forceElem(11,i);
        Mz2 = forceElem(12,i);
        
        rhoV2 = 0;
        g20 = 0;
        g21 = 0;
        if sqrt(forceElem(8,i)^2+forceElem(9,i)^2) > 0.5*Vpl
            rhoV2 = (2*sqrt((forceElem(8,i)^2)+(forceElem(9,i)^2))/Vpl - 1)^2;
%             g20 = -forceElem(12,i)/(1.04*properties.Mpl(i))*4*(2*forceElem(8,i)/properties.Vpl(i) - 1)/((1-rhoV2)^2);
%             g21 = -forceElem(11,i)/(1.04*properties.Mpl(i))*4*(2*forceElem(9,i)/properties.Vpl(i) - 1)/((1-rhoV2)^2);
        end
        dF1dN2 = 4*N2/(Npl^2*(1-(N2/Npl)^2)^3) * ((My2/xiMplyRight)^2+(Mz2/xiMplzRight)^2)/(((1-rhoV2)*Mpl)^2);
        dF1dMy2 = 2*My2/(((1-rhoV2)*xiMplyRight*Mpl*(1-(N2/Npl)^2))^2);
        dF1dMz2 = 2*Mz2/(((1-rhoV2)*xiMplzRight*Mpl*(1-(N2/Npl)^2))^2);
        g2 = [dF1dN2 ; g20 ; g21 ; 0 ; dF1dMy2 ; dF1dMz2];
        
        % Stiffness
        
        if plasticity(2,i) == 1 && plasticity(3,i) == 1
            det1 = (g1.'*(k11+k11h+C11)*g1)*(g2.'*(k22+k22h+C22)*g2) - (g1.'*k12*g2)*(g2.'*k21*g1);
            h11T = 1/det1*((g2.'*(k22+k22h+C22)*g2)*g1.'*k11 - (g1.'*k12*g2)*g2.'*k21);
            h12T = 1/det1*((g2.'*(k22+k22h+C22)*g2)*g1.'*k12 - (g1.'*k12*g2)*g2.'*k22);
            h21T = 1/det1*((g1.'*(k11+k11h+C11)*g1)*g2.'*k21 - (g2.'*k21*g1)*g1.'*k11);
            h22T = 1/det1*((g1.'*(k11+k11h+C11)*g1)*g2.'*k22 - (g2.'*k21*g1)*g1.'*k12);
            
        elseif plasticity(2,i) == 1
            det2 = g1.'*(k11+k11h+C11)*g1;
            h11T = 1/det2*g1.'*k11;
            h12T = 1/det2*g1.'*k12;
            h21T = zeros(1,6);
            h22T = zeros(1,6);
            
        elseif plasticity(3,i) == 1
            det3 = g2.'*(k22+k22h+C22)*g2;
            h11T = zeros(1,6);
            h12T = zeros(1,6);
            h21T = 1/det3*g2.'*k21;
            h22T = 1/det3*g2.'*k22;
            
        else
            h11T = zeros(1,6);
            h12T = zeros(1,6);
            h21T = zeros(1,6);
            h22T = zeros(1,6);
            
        end
        
        k11ep = k11 - k11*g1*h11T - k12*g2*h21T;
        k12ep = k12 - k11*g1*h12T - k12*g2*h22T;
        k21ep = k21 - k21*g1*h11T - k22*g2*h21T;
        k22ep = k22 - k21*g1*h12T - k22*g2*h22T;
        
        Kelem(:,:,i) = [k11ep k12ep ;
            k21ep k22ep];
        
        % Plastic correction
        
        if plasticity(2,i) == 1 && plasticity(3,i) == 1
            forceFunct1 = forceElem(1:6,i);
            deltaFSurf1 = normalSurf4(properties,forceFunct1,rhoV1,i,xiMplyLeft,xiMplzLeft);
            forceFunct2 = forceElem(7:12,i);
            deltaFSurf2 = normalSurf4(properties,forceFunct2,rhoV2,i,xiMplyRight,xiMplzRight);
            
            t11 = 1/det1*(g2.'*(k22ep+k22h+C22)*g2)*g1.';
            t12 = 1/det1*(-g1.'*k12ep*g2)*g2.';
            t21 = 1/det1*(-g2.'*k21ep*g1)*g1.';
            t22 = 1/det1*(g1.'*(k11ep+k11h+C11)*g1)*g2.';
            
%             t11 = 1/det1*(g2.'*(k22+k22h+C22)*g2)*g1.';
%             t12 = 1/det1*(-g2.'*k21*g1)*g1.';
%             t21 = 1/det1*(-g1.'*k12*g2)*g2.';
%             t22 = 1/det1*(g1.'*(k11+k11h+C11)*g1)*g2.';
            
        elseif plasticity(2,i) == 1
            forceFunct1 = forceElem(1:6,i);
            deltaFSurf1 = normalSurf4(properties,forceFunct1,rhoV1,i,xiMplyLeft,xiMplzLeft);
            deltaFSurf2 = zeros(6,1);
            
            t11 = 1/det2*g1.';
            t12 = zeros(1,6);
            t21 = zeros(1,6);
            t22 = zeros(1,6);
            
        elseif plasticity(3,i) == 1
            deltaFSurf1 = zeros(6,1);
            forceFunct2 = forceElem(7:12,i);
            deltaFSurf2 = normalSurf4(properties,forceFunct2,rhoV2,i,xiMplyRight,xiMplzRight);
            
            t11 = zeros(1,6);
            t12 = zeros(1,6);
            t21 = zeros(1,6);
            t22 = 1/det3*g2.';
            
        else
            deltaFSurf1 = zeros(6,1);
            deltaFSurf2 = zeros(6,1);
            
            t11 = zeros(1,6);
            t12 = zeros(1,6);
            t21 = zeros(1,6);
            t22 = zeros(1,6);
            
        end
        
        t11ep = -k11ep*g1*t11 - k12ep*g2*t21;
        t12ep = -k11ep*g1*t12 - k12ep*g2*t22;
        t21ep = -k21ep*g1*t11 - k22ep*g2*t21;
        t22ep = -k21ep*g1*t12 - k22ep*g2*t22;
        
%         t11ep = -k11*g1*t11 - k12*g2*t21;
%         t12ep = -k11*g1*t12 - k12*g2*t22;
%         t21ep = -k21*g1*t11 - k22*g2*t21;
%         t22ep = -k21*g1*t12 - k22*g2*t22;

    end
    
    %% Global matrix
    
    relatedDofs = properties.relatedDofs(i,:);
    
    K(relatedDofs,relatedDofs) = K(relatedDofs,relatedDofs) + T(:,:,i).'*Kelem(:,:,i)*T(:,:,i);
    
    tep = [t11ep t12ep;
        t21ep t22ep];
    deltaFSurf = [deltaFSurf1 ; deltaFSurf2];
    
    forceVec(relatedDofs) = forceVec(relatedDofs) - T(:,:,i).'*(tep*deltaFSurf);

end

matrices.K = K;
matrices.Kelem = Kelem;
matrices.T = T;