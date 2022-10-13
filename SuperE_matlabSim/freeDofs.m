function[freeDofs] = freeDofs(dofs)

freeDofs = zeros(1,dofs.nbDofs-length(dofs.blokedDofs));

incBD = 1;
incFD = 1;
for i = 1:dofs.nbDofs
    if incBD <= length(dofs.blokedDofs) && dofs.blokedDofs(incBD) == i;
        incBD = incBD + 1;
    else
        freeDofs(incFD) = i;
        incFD = incFD + 1;
    end
end