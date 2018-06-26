function [omega1, omega2, omega3]=calc_vorticity_vector(uu,vv,ww,z_i,l_r,Nx,Ny,Nz,dz)
%NOTE: output is on w-nodes

    % calculate gradients ddx, ddy of  u,v at u-v-p nodes
    [~,dudy]=ddx(uu,z_i,l_r,Nx);
    [dvdx,~]=ddx(vv,z_i,l_r,Ny);
    % interpolate ddy of u,v at w nodes (2:Nz)
    dudy = 0.50 .* (dudy(:,:,1:Nz-1) + dudy(:,:,2:Nz));
    dvdx = 0.50 .* (dvdx(:,:,1:Nz-1) + dvdx(:,:,2:Nz));
    % calculate ddx of w at w-nodes (2:Nz)
    [dwdx,dwdy]=ddx(ww,z_i,l_r,Nx); dwdx = dwdx(:,:,2:Nz);dwdy = dwdy(:,:,2:Nz);
    % calculate ddz of u, v at w nodes (2:Nz)
    dudz = ddz_uvp(uu, dz); dvdz = ddz_uvp(vv, dz);


    omega1 = zeros([Nx, Ny, Nz]); omega1(:,:,2:Nz) = dwdy - dvdz;
    omega2 = zeros([Nx, Ny, Nz]); omega2(:,:,2:Nz) = dudz - dwdx;
    omega3 = zeros([Nx, Ny, Nz]); omega3(:,:,2:Nz) = dvdx - dudy;
    omega1(:,:,1) = -vv(:,:,1)/dz/2.0; omega2(:,:,1) = uu(:,:,1)/dz/2.0;

    
%     omega3=zeros(size(dudy,1),size(dudy,2),size(dudy,3));
%     for k=size(dudy,3):-1:2
%           omega3(:,:,k)=(dvdx(:,:,k)+dvdx(:,:,k-1))/2-(dudy(:,:,k)+dudy(:,:,k-1))/2;
%     end
%     omega{3}=omega3;


return
end
