
%  First deriv in z direction for boundary layer (2nd order numerics)
%  F is on UVP nodes and dFdz is on w nodes, 2:Nz

function [dfdz] = ddz_uvp(f,dz)
    dfdz = zeros([size(f,1), size(f,2), size(f,3)-1]);
    for k=1:size(f,3)-1
        dfdz(:,:,k)=(f(:,:,k+1)-f(:,:,k))./dz;
    end 
end
