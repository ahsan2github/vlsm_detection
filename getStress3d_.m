function [horizontal_shear_stress, uw_total,uv_total, vw_total, u_star_mean]= ...
    getStress3d_(parentDir, Nx, Ny, Nz, dx, dy, dz, u_star, Ugal, Vgal,zo, ...
    l_r, z_i,frameVec)
uw_total = zeros([Nx, Ny, Nz]);
uv_total = zeros([Nx, Ny, Nz]);
vw_total = zeros([Nx, Ny, Nz]);

tau12 = zeros([Nx, Ny, Nz-1]); 
tau13 = zeros([Nx, Ny, Nz-1]);
tau23 = zeros([Nx, Ny, Nz-1]);     
u_mean_nz = zeros([Nz, 1]);
v_mean_nz = zeros([Nz, 1]);
w_mean_nz = zeros([Nz, 1]);
tau_reynolds_s = 0.0;
u_star_mean = 0.0;
horizontal_shear_stress = zeros([Nx, Ny, Nz]);
delta = (dx*dy*dz)^(1.0/3.0);
framecnt = 0;

for frameNo = frameVec
    frameStr = sprintf('%4.4i',frameNo);
    fn = [parentDir '/output/u_frame/u_frame'  frameStr   '.bin'];
    fh = fopen(fn, 'r');
    uu = fread(fh,'double');
    uu = reshape(uu, [Nx,Ny,Nz]).*u_star + Ugal;     
    fclose(fh);    

    fn = [parentDir '/output/v_frame/v_frame'  frameStr   '.bin'];
    fh = fopen(fn, 'r');
    vv = fread(fh,'double');
    vv = reshape(vv, [Nx,Ny,Nz]).*u_star + Vgal; 
    fclose(fh);     

    fn = [parentDir '/output/w_frame/w_frame'  frameStr   '.bin'];
    fh = fopen(fn, 'r');
    ww = fread(fh,'double');
    ww = reshape(ww, [Nx,Ny,Nz]).*u_star; 
    fclose(fh); 

    fn = [parentDir '/output/Cs2_frame/Cs2_frame'  frameStr   '.bin'];
    fh = fopen(fn, 'r');
    Cs2 = fread(fh,'double');
    Cs2 = reshape(Cs2, [Nx,Ny,Nz]); 
    fclose(fh); 
    % calculate mean at u,v,p nodes    
    for jj = 1:Nz                
        u_mean_nz(jj) = mean(mean(uu(:,:,jj)));
        v_mean_nz(jj) = mean(mean(vv(:,:,jj)));
        w_mean_nz(jj) = mean(mean(ww(:,:,jj)));
    end
    uw_nxnynz = zeros([Nx, Ny, Nz]);
    uv_nxnynz = zeros([Nx, Ny, Nz]);
    vw_nxnynz = zeros([Nx, Ny, Nz]);
    % calculate reynolds stress at uvp nodes, dimension [Nx, Ny, Nz]      
    for i = 1:Nz
        uw_nxnynz(:,:,i) = uw_nxnynz(:,:,i) + (uu(:,:,i)-u_mean_nz(i)).*(ww(:,:,i)-w_mean_nz(i));
        uv_nxnynz(:,:,i) = uv_nxnynz(:,:,i) + (uu(:,:,i)-u_mean_nz(i)).*(vv(:,:,i)-v_mean_nz(i));
        vw_nxnynz(:,:,i) = vw_nxnynz(:,:,i) + (vv(:,:,i)-v_mean_nz(i)).*(ww(:,:,i)-w_mean_nz(i));            
    end
    % calculate the surface shear stress from actual velocity field,
    % similarity theory
    u_r = sqrt((uu(:,:,1)).^2+vv(:,:,1).^2);
    tau_xzs = -(u_r(:,:)./log(dz/2/zo)).^2 .*(uu(:,:,1)./...
        u_r(:,:));
    tau_yzs= -(u_r(:,:)./log(dz/2/zo)).^2 .*(vv(:,:,1)./...
        u_r(:,:));

    tau_reynolds_s = tau_reynolds_s + (tau_xzs.^2 + tau_yzs.^2).^0.5;
    u_star_mean = u_star_mean + mean(mean(((tau_xzs.^2 + tau_yzs.^2).^0.25)));
    % interpolate uw, uv, vw at w-nodes, 2:Nz w nodes
    uw = uw_nxnynz(:,:,1:Nz-1) + uw_nxnynz(:,:,2:Nz);
    uv = uv_nxnynz(:,:,1:Nz-1) + uv_nxnynz(:,:,2:Nz);
    vw = vw_nxnynz(:,:,1:Nz-1) + vw_nxnynz(:,:,2:Nz);
    clearvars uw_nxnynz uv_nxnynz vw_nxnynz;
    % calculate gradients ddx, ddy of  u,v at u-v-p nodes
    [dudx,dudy]=ddx(uu,z_i,l_r,Nx);
    [dvdx,dvdy]=ddx(vv,z_i,l_r,Nx);
    % interpolate ddx, ddy of u,v at w nodes (2:Nz)
    dudx = 0.50 .* (dudx(:,:,1:Nz-1) + dudx(:,:,2:Nz));
    dudy = 0.50 .* (dudy(:,:,1:Nz-1) + dudy(:,:,2:Nz));
    dvdx = 0.50 .* (dvdx(:,:,1:Nz-1) + dvdx(:,:,2:Nz));
    dvdy = 0.50 .* (dvdy(:,:,1:Nz-1) + dvdy(:,:,2:Nz));

    % calculate ddx of w at w-nodes (2:Nz)
    [dwdx,dwdy]=ddx(ww,z_i,l_r,Nx); dwdx = dwdx(:,:,2:Nz);dwdy = dwdy(:,:,2:Nz);
    % calculate ddz of u, v at w nodes (2:Nz)
    dudz = ddz_uvp(uu, dz); dvdz = ddz_uvp(vv, dz);
    % calculate ddz of w at uvp nodes, will have a vertical dimension Nz-1
    dwdz = ddz_w(ww, dz); 
    % interpolate dwdz at w-nodes, this will leave dwdz with a vertical
    % dimension of Nz-2, 2:Nz-1
    dwdz = dwdz(:,:,1:Nz-2)+dwdz(:,:,2:Nz-1);                
    % calculate strain rates at w nodes, will have vertical dimensiom
    % of Nz-2, 2:Nz-1 w-nodes
    S11 = 0.5.*(dudx(:,:,1:end-1) + dudx(:,:,1:end-1)); 
    S12 = 0.5.*(dudy(:,:,1:end-1) + dvdx(:,:,1:end-1)); 
    S13 = 0.5.*(dudz(:,:,1:end-1) + dwdx(:,:,1:end-1)); 
    S22 = 0.5.*(dvdy(:,:,1:end-1) + dvdy(:,:,1:end-1)); 
    S23 = 0.5.*(dvdz(:,:,1:end-1) + dwdy(:,:,1:end-1)); 
    S33 = 0.5.*(dwdz(:,:,1:end) + dwdz(:,:,1:end));
    abs_S = sqrt(2.0*(S11.*S11 + S12.*S12 + S13.*S13 + S12.*S12 + S22.*S22 + ...
            S23.*S23 + S13.*S13 + S23.*S23 + S33.*S33));
	
	clearvars dudx dudy dudz dvdx dvdy dvdz dwdx dwdy dwdz;
    % at this point all quantities are available at nodes 2:Nz-1   
    % if mom_nodes == 0 in LESinputs, the first CS2 value is calculated
    % at dz/2 location, then subsequent Cs values are calculated at
    % regular w-nodes starting from 2nd w-node
    tau12(:,:,1) = 0.0; tau13(:,:,1) = tau_xzs; 
    tau23(:,:,1) = tau_yzs;
    for i = 2:Nz-1
        tau12(:,:,i) = ( -2.0 .* (delta .* sqrt(Cs2(:,:,i)))^2 .* abs_S(:,:,i-1) .* S12(:,:,i-1) );
        tau13(:,:,i) = ( -2.0 .* (delta .* sqrt(Cs2(:,:,i)))^2 .* abs_S(:,:,i-1) .* S13(:,:,i-1) );
        tau23(:,:,i) = ( -2.0 .* (delta .* sqrt(Cs2(:,:,i)))^2 .* abs_S(:,:,i-1) .* S23(:,:,i-1) ) ;
    end        
    clearvars S11 S12 S13 S22 S23 S33;
    % add turb stress to resolved stress, stresses are at 2:Nz-1 w-nodes
    uw_total(:,:,1) = uw_total(:,:,1) + tau13(:,:,1);
    uv_total(:,:,1) = uv_total(:,:,1) + 0.0;
    vw_total(:,:,1) = vw_total(:,:,1) + tau23(:,:,1);
    uw_total(:,:,2:Nz-1) = uw_total(:,:,2:Nz-1) + uw(:,:,1:end-1) + tau13(:,:,2:end); 
    uv_total(:,:,2:Nz-1) = uv_total(:,:,2:Nz-1) + uv(:,:,1:end-1) + tau12(:,:,2:end); 
    vw_total(:,:,2:Nz-1) = vw_total(:,:,2:Nz-1) + vw(:,:,1:end-1) + tau23(:,:,2:end);
    horizontal_shear_stress(:,:,1) = horizontal_shear_stress(:,:,1)+...
        (tau13(:,:,1).^2 + tau23(:,:,1).^2).^(0.5);
    horizontal_shear_stress(:,:,2:end-1) = horizontal_shear_stress(:,:,2:end-1) + ...
        ((uw(:,:,1:end-1) + tau13(:,:,2:end)).^2 + ...
        (vw(:,:,1:end-1) + tau23(:,:,2:end)).^2).^(0.5);
framecnt = framecnt + 1;    
end
u_star_mean = u_star_mean/framecnt;
uw_total = uw_total./framecnt;
uv_total = uv_total./framecnt;
vw_total = vw_total./framecnt;
horizontal_shear_stress = horizontal_shear_stress/framecnt;
end

