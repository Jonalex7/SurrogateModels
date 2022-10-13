function[forceElem] = internalForces(properties,matrices,dispTot,dofs)

forceElem = zeros(12,dofs.nbElements);

for i = 1:dofs.nbElements
%     if dofs.thirdNode ~= 0
%         i
%         matrices.Kelem(:,:,i)
%         matrices.T(:,:,i)
%         properties.relatedDofs(i,:)
%         dispTot(properties.relatedDofs(i,:))
%     end
    forceElem(:,i) = matrices.Kelem(:,:,i)*matrices.T(:,:,i)*dispTot(properties.relatedDofs(i,:));
end