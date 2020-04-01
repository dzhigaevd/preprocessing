function [out] = centerofmass(X,Y,Z,V)
%COM Summary of this function goes here
%   X, Y, Z are the coordinates (have to be the same size as V and each other)
%   V are the values

% X
d = 1;
for i=1:size(X,1)
    for j=1:size(X,2)
        for k=1:size(X,3)
            dCOMx(d) = X(i,j,k)*V(i,j,k);
            d = d +1;
        end
    end
end
COMx = sum(dCOMx)/sum(V,'all');

% Y
d = 1;
for i=1:size(Y,1)
    for j=1:size(Y,2)
        for k=1:size(Y,3)
            dCOMy(d) = Y(i,j,k)*V(i,j,k);
            d = d +1;
        end
    end
end
COMy = sum(dCOMy)/sum(V,'all');

% Z
d = 1;
for i=1:size(Z,1)
    for j=1:size(Z,2)
        for k=1:size(Z,3)
            dCOMz(d) = Z(i,j,k)*V(i,j,k);
            d = d +1;
        end
    end
end
COMz = sum(dCOMz)/sum(V,'all');

out = [COMx, COMy, COMz];

end

