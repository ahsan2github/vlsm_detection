function [divlamb] = calc_Div_lamb_vector(omegax,omegay,omegaz,u,v,w,Nx,Nz,z_i,l_r,dz)
%NOTE: output is on w-nodes


    % get u/v at w nodes
    for k=size(u,3):-1:2
        u(:,:,k)=(u(:,:,k)+u(:,:,k-1))/2;
        v(:,:,k)=(v(:,:,k)+v(:,:,k-1))/2;    
    end
    u(:,:,1) = zeros([size(u,1) size(u,2)]);
    v(:,:,1) = zeros([size(v,1) size(v,2)]);
    % get d/dz(omega_y) and d/dz(omega_x) at uvp-nodes
    ddz_omega_y = ddz_w(omegay,dz);
    ddz_omega_x = ddz_w(omegax,dz);
    % interpolate ddz_omega_x, ddz_omega_y to w-nodes
    for i =2:Nz
       ddz_omega_y(:,:,i) = (ddz_omega_y(:,:,i) + ddz_omega_y(:,:,i-1)) ./ 2.0;
       ddz_omega_x(:,:,i) = (ddz_omega_x(:,:,i) + ddz_omega_x(:,:,i-1)) ./ 2.0;
    end
    
    % get ddx_omega_z, ddy_omega_z at w nodes
    [ddx_omega_z,ddy_omega_z]=ddx(omegaz,z_i,l_r,Nx);
    % get ddx_omega_y and ddy_omega_x at w nodes
    [ddx_omega_y, ~] = ddx(omegay,z_i,l_r,Nx); 
    [~, ddy_omega_x] = ddx(omegax,z_i,l_r,Nx);
    
    divlamb= u .* (ddy_omega_z - ddz_omega_y) + v .* (ddx_omega_z - ddz_omega_x) + ...
             w .* (ddx_omega_y - ddy_omega_x) - omegax.*omegax - omegay.*omegay - omegaz.*omegaz;
    
end

    
