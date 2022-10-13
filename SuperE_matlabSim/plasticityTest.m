function[plasticity] = plasticityTest(properties,dofs,matrices,dispTot,forceElem,plasticity,punchProp)

plasticitym1 = plasticity;

plasticity = zeros(3,dofs.nbElements);

for i = 1:dofs.nbElements
    %% Properties
    
    Npl = properties.Npl(i);
    Mpl = properties.Mpl(i);
    
    %% Plastic hinge at midspan
    midSpan = properties.L(i)/2;
%     h2 = 2*midSpan^3/(properties.L(i)^3) - 3*midSpan^2/(properties.L(i)^2) + 1;
    h3 = properties.L(i)*(midSpan^3/(properties.L(i)^3) - 2*midSpan^2/properties.L(i)^2 + midSpan/properties.L(i));
%     h5 = -2*midSpan^3/(properties.L(i)^3) + 3*midSpan^2/(properties.L(i)^2);
    h6 = properties.L(i)*(midSpan^3/(properties.L(i)^3) - midSpan^2/(properties.L(i)^2));
    
    dispLocal = matrices.T(:,:,i)*dispTot(properties.relatedDofs(i,:));
    
    vMidSpany = h3*dispLocal(5) + h6*dispLocal(11);
    vMidSpanz = h3*dispLocal(6) + h6*dispLocal(12);
    
    NEd = abs(forceElem(1,i));
    MEdy = abs((forceElem(5,i) - forceElem(11,i))/2 + forceElem(1,i)*vMidSpany);
    MEdz = abs((forceElem(6,i) - forceElem(12,i))/2 + forceElem(1,i)*vMidSpanz);
    
%     NEd = forceElem(1,i);
%     MEd = (forceElem(3,i) - forceElem(6,i))/2 + forceElem(1,i)*vMidSpan;
    
    rho = 0;
    if sqrt(forceElem(2,i)^2+forceElem(3,i)^2) > 0.5*properties.Vpl(i)
        rho = (2*sqrt(forceElem(2,i)^2+forceElem(3,i)^2)/properties.Vpl(i) - 1)^2;
    end
    interactionCritEl = (MEdy/((1-rho)*Mpl*(1 - (NEd/Npl)^2)))^2 + (MEdz/((1-rho)*Mpl*(1 - (NEd/Npl)^2)))^2;
    
    if interactionCritEl < 1
        if forceElem(1,i) > 0 % Compression
            plasticity(1,i) = -1;
        else                  % Tension
            plasticity(1,i) = -2;
        end
    else
        if forceElem(1,i) > 0 % Compression
            plasticity(1,i) = 1;
        else                  % Tension
            plasticity(1,i) = 2;
        end
    end
    
    %% Plastic hinge at extremity
    
    % Reduction of stiffness based on plastic moment reduction
    a = punchProp.redSectStar(properties.nodein(i),1); % axis X
    b = punchProp.redSectStar(properties.nodeout(i),1);
    c = punchProp.redSectStar(properties.nodein(i),2); % axis Y
    d = punchProp.redSectStar(properties.nodeout(i),2);
    
    vecDir1 = properties.angle(1,:,i);
    
    % Bending moment on the left part
    
    if vecDir1(3) == 0
        Mply = Mpl;
        Mplz = Mpl;
    else
        dir = -vecDir1(1)/vecDir1(3);
        axisXElement = [1 0 dir];
        axisXElement = axisXElement/norm(axisXElement);
        vecDir2 = properties.angle(2,:,i);
        alpha = acos(axisXElement * vecDir2');
        xiMplyLeft = (a * abs(cos(alpha)) + c * abs(sin(alpha))) / (abs(cos(alpha)) + abs(sin(alpha)));
        xiMplzLeft = (a * abs(sin(alpha)) + c * abs(cos(alpha))) / (abs(cos(alpha)) + abs(sin(alpha)));
        Mply = xiMplyLeft * Mpl;
        Mplz = xiMplzLeft * Mpl;
    end
    
    rhoV = 0;
    if sqrt(forceElem(2,i)^2+forceElem(3,i)^2) > 0.5*properties.Vpl(i)
        rhoV = (2*sqrt(forceElem(2,i)^2+forceElem(3,i)^2)/properties.Vpl(i) - 1)^2;
    end
    MEdy = abs(forceElem(5,i));
    MEdz = abs(forceElem(6,i));
    interactionCritLeft = (MEdy/((1-rhoV)*Mply*(1 - (NEd/Npl)^2)))^2 + (MEdz/((1-rhoV)*Mplz*(1 - (NEd/Npl)^2)))^2;
%     if i == 5% || i == 4
%         disp('Left')
% %         disp(['rho = ',num2str(rhoV)])
% %         disp(['V tot = ',num2str(sqrt(forceElem(2,i)^2+forceElem(3,i)^2))])
% %         disp(['Rapport V = ',num2str(sqrt(forceElem(2,i)^2+forceElem(3,i)^2) / (0.5*properties.Vpl(i)))])
%         disp(interactionCritLeft)
% %         disp(['contrN = ',num2str((NEd/properties.Npl(i))^2)])
% %         disp(['contrMy = ',num2str((MEdy/((1-rhoV)*properties.Mpl(i)))^2)])
% %         disp(['contrMz = ',num2str((MEdz/((1-rhoV)*properties.Mpl(i)))^2)])
%     end
    
    if interactionCritLeft > 1 
        plasticity(2,i) = 1;
%         disp('oui')
    end
    
    % Bending moment on the right part
    
    if vecDir1(3) == 0
        Mply = Mpl;
        Mplz = Mpl;
    else
        dir = -vecDir1(1)/vecDir1(3);
        axisXElement = [1 0 dir];
        axisXElement = axisXElement/norm(axisXElement);
        vecDir2 = properties.angle(2,:,i);
        alpha = acos(axisXElement * vecDir2');
        xiMplyRight = (b * abs(cos(alpha)) + d * abs(sin(alpha))) / (abs(cos(alpha)) + abs(sin(alpha)));
        xiMplzRight = (b * abs(sin(alpha)) + d * abs(cos(alpha))) / (abs(cos(alpha)) + abs(sin(alpha)));
        Mply = xiMplyRight * Mpl;
        Mplz = xiMplzRight * Mpl;
    end
    
    rhoV = 0;
    if sqrt(forceElem(8,i)^2+forceElem(9,i)^2) > 0.5*properties.Vpl(i)
        rhoV = (2*sqrt(forceElem(8,i)^2+forceElem(9,i)^2)/properties.Vpl(i) - 1)^2;
    end
    MEdy = abs(forceElem(11,i));
    MEdz = abs(forceElem(12,i));
    interactionCritRight = (MEdy/((1-rhoV)*Mply*(1 - (NEd/Npl)^2)))^2 + (MEdz/((1-rhoV)*Mplz*(1 - (NEd/Npl)^2)))^2;
%     if i == 5% || i == 4
% %         disp('Right')
% %         disp(['rho = ',num2str(rhoV)])
% %         disp(['Rapport V = ',num2str(sqrt(forceElem(8,i)^2+forceElem(9,i)^2) / (0.5*properties.Vpl(i)))])
%         disp(interactionCritRight)
% %         disp(['contrN = ',num2str((NEd/properties.Npl(i))^2)])
% %         disp(['contrMy = ',num2str((MEdy/((1-rhoV)*properties.Mpl(i)))^2)])
% %         disp(['contrMz = ',num2str((MEdz/((1-rhoV)*properties.Mpl(i)))^2)])
% %         crit = interactionCritRight;
%     end
    
    if interactionCritRight > 1 
        plasticity(3,i) = 1;
%         disp('oui')
    end    
    
end

%% If plasticity reached, stays in plasticity

% for i = 1:dofs.nbElements
%     for j = 1:3
%         if j == 1
%             if plasticitym1(j,i) > 0 && interactionCritEl > 0.95
%                 plasticity(j,i) = plasticitym1(j,i);
%             end
%         elseif j == 2
%             if plasticitym1(j,i) == 1 && interactionCritLeft > 0.95
%                 plasticity(j,i) = 1;
%             end
%         elseif j == 3
%             if plasticitym1(j,i) == 1 && interactionCritRight > 0.95
%                 plasticity(j,i) = 1;
%             end
%         end
%     end
% end