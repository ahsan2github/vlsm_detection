%%
clear; clc; close all;
tic;

frameVec=[9] %,10,11,12,13,14,15,16,17,18,19];
bxw = 15; byw = 3.0; bzh = 1;
width = 3.0; height = width / 1.1; formatEng = '-depsc';
SS ='';
parentDir = '/raid/home/mohammad/NeutralChannel_2048x2048x64_free';
%parentDir = '/glade/scratch/khanm/NeutralChannel_2048x2048x64_free';
%parentDir = '/scratch/ibrix/chpc_gen/u0699351/NeutralChannel_2048x2048x64_free';
outputDir = [parentDir '/output/'];
readinputs(parentDir);
u_star_in = u_star; clear u_star;
BLH = 1500.0;
heightVec = [2.0] %, 4.0, 7.0,10, 13, 17, 20, ...
    %23, 26, 29, 32, 36, 38, 42, 45]; 
avg_stress_slcth  = zeros([Nx, Ny, length(heightVec)]);


% get shear stress for frames and save them as binary
%for frameNo = frameVec 
%    stresFieldName = sprintf('chnl_shearStress%4.4i.bin',frameNo);
%    fid = fopen(stresFieldName, 'w');    
%    [shear_stress,~,~,~, u_star]= ...
%    getStress3d_(parentDir, Nx, Ny, Nz, dx, dy, dz, u_star_in, Ugal, Vgal,zo, l_r, z_i,frameNo);
%    fwrite(fid,shear_stress,'double');
%    fclose(fid);
%end

% get Divergence of the lamb vector for frames and save them as binary
for frameNo = frameVec 
    frameStr = sprintf('%4.4i',frameNo);
    fn = [parentDir '/output/u_frame/u_frame'  frameStr   '.bin'];
    fh = fopen(fn, 'r');
    uu = fread(fh,'double');
    uu = reshape(uu, [Nx,Ny,Nz]).*u_star_in + Ugal;
    fclose(fh);
    frameStr = sprintf('%4.4i',frameNo);
    fn = [parentDir '/output/v_frame/v_frame'  frameStr   '.bin'];
    fh = fopen(fn, 'r');
    vv = fread(fh,'double');
    vv = reshape(vv, [Nx,Ny,Nz]).*u_star_in + Vgal;  
    fclose(fh);      
    frameStr = sprintf('%4.4i',frameNo);
    fn = [parentDir '/output/w_frame/w_frame'  frameStr   '.bin'];
    fh = fopen(fn, 'r');
    ww = fread(fh,'double');
    ww = reshape(ww, [Nx,Ny,Nz]).*u_star_in;
    fclose(fh); 
    
    [omegax, omegay, omegaz]=calc_vorticity_vector(uu,vv,ww,z_i,l_r,Nx,Ny,Nz,dz);
    omegaFName_1 = sprintf('chnl_omega_x_%4.4i.bin',frameNo);
    fid = fopen(omegaFName_1, 'w');
    fwrite(fid,omegax,'double');
    fclose(fid);
    omegaFName_2 = sprintf('chnl_omega_y_%4.4i.bin',frameNo);       
    fid = fopen(omegaFName_2, 'w');
    fwrite(fid,omegay,'double');
    fclose(fid);
    omegaFName_3 = sprintf('chnl_omega_z_%4.4i.bin',frameNo);       
    fid = fopen(omegaFName_3, 'w');
    fwrite(fid,omegaz,'double');
    fclose(fid);
    
    divlamb = calc_Div_lamb_vector(omegax,omegay,omegaz,uu,vv,ww,Nx,z_i,l_r,dz);
    divLambFName = sprintf('chnl_DivLambVec%4.4i.bin',frameNo);
    fid = fopen(divLambFName, 'w');
    fwrite(fid,divlamb,'double');
    fclose(fid);
end
msg = ['Omega and Divergence of the lamb  vector calculation done!'];

disp(msg);
% create a box just to get box dimension
frameStr = sprintf('%4.4i',frameVec(1));
fn = [parentDir '/output/u_frame/u_frame'  frameStr   '.bin'];
fh = fopen(fn, 'r');
uu = fread(fh,'double');
uu = reshape(uu, [Nx,Ny,Nz]).*u_star_in + Ugal; 
fclose(fh);

tic
[tbox, vpcnt, pxw, pyw] = getAvgBoxBotm2Blh(Nx, Ny, Nz, z_i, BLH, l_z,...
	bxw, byw,round(Nx/2), round(Ny/2), uu);

toc

[bnx, bny, bnz] = size(tbox);
uBox = zeros([bnx, bny, bnz]);
vBox = zeros([bnx, bny, bnz]);
wBox = zeros([bnx, bny, bnz]);
omega1Box = zeros([bnx, bny, bnz]);
omega2Box = zeros([bnx, bny, bnz]);
omega3Box = zeros([bnx, bny, bnz]);
stressBox = zeros([bnx, bny, bnz]);
pBox = zeros([bnx, bny, bnz]);
divLambBox = zeros([bnx, bny, bnz]);
strucFound = 0;
total_frames_analyzed = 0;
ccc = 1;
for frameNo = frameVec 
    frameStr = sprintf('%4.4i',frameNo);
    fn = [parentDir '/output/u_frame/u_frame'  frameStr   '.bin'];
    fh = fopen(fn, 'r');
    uu = fread(fh,'double');
    uu = reshape(uu, [Nx,Ny,Nz]).*u_star_in + Ugal;
    uu=removeLevelMean(uu);
    fclose(fh);  
    
    frameStr = sprintf('%4.4i',frameNo);
    fn = [parentDir '/output/v_frame/v_frame'  frameStr   '.bin'];
    fh = fopen(fn, 'r');
    vv = fread(fh,'double');
    vv = reshape(vv, [Nx,Ny,Nz]).*u_star_in + Vgal;
    vv =removeLevelMean(vv);    
    fclose(fh);  
    
    frameStr = sprintf('%4.4i',frameNo);
    fn = [parentDir '/output/w_frame/w_frame'  frameStr   '.bin'];
    fh = fopen(fn, 'r');
    ww = fread(fh,'double');
    ww = reshape(ww, [Nx,Ny,Nz]).*u_star_in;
    ww =removeLevelMean(ww);    
    fclose(fh); 
         
    omegaFName = sprintf('chnl_omega_x_%4.4i.bin',frameNo);
    ff = fopen(omegaFName, 'r');
    omega1 = fread(ff,'double');
    fclose(ff); 
    omega1 = reshape(omega1, [Nx Ny Nz]);
    
    omegaFName = sprintf('chnl_omega_y_%4.4i.bin',frameNo);
    ff = fopen(omegaFName, 'r');
    omega2 = fread(ff,'double');
    fclose(ff); 
    omega2 = reshape(omega2, [Nx Ny Nz]);
    
    omegaFName = sprintf('chnl_omega_z_%4.4i.bin',frameNo);
    ff = fopen(omegaFName, 'r');
    omega3 = fread(ff,'double');
    fclose(ff); 
    omega3 = reshape(omega3, [Nx Ny Nz]);
    
    frameStr = sprintf('%4.4i',frameNo);
    fn = [parentDir '/output/p_frame/p_frame'  frameStr   '.bin'];
    fh = fopen(fn, 'r');
    pp= fread(fh,'double');
    pp = reshape(pp, [Nx,Ny,Nz]);   
    fclose(fh); 

    stresFieldName = sprintf('chnl_shearStress%4.4i.bin',frameNo);         
    ff = fopen(stresFieldName, 'r');
    data_shear_stress = fread(ff,'double');
    fclose(ff);
    data_shear_stress = reshape(data_shear_stress, [Nx, Ny, Nz]);
    shear_stress_field = data_shear_stress(:,:,2);   
    
    divLambFName = sprintf('chnl_DivLambVec%4.4i.bin',frameNo);
    ff = fopen(divLambFName, 'r');
    div_lamb_vec = fread(ff,'double');
    fclose(ff);    
    div_lamb_vec = reshape( div_lamb_vec, [Nx, Ny, Nz]);
    div_lamb_vec_at_lvl_2 = div_lamb_vec(:,:,2);
    
    [N, binc] = hist(div_lamb_vec_at_lvl_2(:)./max(div_lamb_vec_at_lvl_2(:)),20);
    figure(); plot(binc, N./sum(N(:)));

    % calculation for high stress structures
    lowDivLambBin = div_lamb_vec_at_lvl_2(:) <0; % gives a binary array 
    lowDivLambBin  = reshape(lowDivLambBin, size(squeeze(div_lamb_vec_at_lvl_2)));
    
    figure(); imagesc(div_lamb_vec_at_lvl_2); hold on; colormap parula; colorbar;    

    avg_stress_slcth  = avg_stress_slcth  + data_shear_stress(:,:,heightVec);

    props = regionprops(lowDivLambBin, 'Centroid',...
        'MajorAxisLength','MinorAxisLength', 'Orientation');

    mAxisLengths = [props(1:end).MajorAxisLength]; 
    % select major axis length ranges 
    thresMajorAxisLength_low = floor(Nx * (BLH*1)/(2*pi*z_i)); %
    thresMajorAxisLength_high = floor(Nx * (BLH*24)/(2*pi*z_i)); %
    c1 = mAxisLengths >= thresMajorAxisLength_low;
    c2 = mAxisLengths <= thresMajorAxisLength_high;   
    crit = find(c1 .* c2);
    props = props(crit);
    count = length(props);
    stressMax = zeros(count,1);
    vpcnt = zeros(count,1);
    for indx = 1:count
        [box_stress, vpcnt, pxw, pyw] = getAvgBoxBotm2Blh(Nx, Ny, Nz, z_i, BLH, l_z,...
            bxw, byw, props(indx).Centroid(2), props(indx).Centroid(1), data_shear_stress);
        stressMax(indx) = max(max(box_stress(:,:,1)));
        vpcnt(indx) = vpcnt;
    end
    divLambMin = zeros(count,1);
    for indx = 1:count
        [box_divLamb, vpcnt, pxw, pyw] = getAvgBoxBotm2Blh(Nx, Ny, Nz, z_i, BLH, l_z,...
            bxw, byw, props(indx).Centroid(2), props(indx).Centroid(1), div_lamb_vec);
        divLambMin(indx) = (min(min(box_divLamb(:,:,2))));
    end
    tracker = 1;
    while(~isempty(props))
        if(tracker > length(props)), break; end;
        cxCondition = zeros([length(props) 1]);
        cyCondition = zeros([length(props) 1]);        
        for counter = 1:length(props)
            cxCondition(counter) = abs(props(counter).Centroid(2) - props(tracker).Centroid(2)) < pxw*2;
            cyCondition(counter) = abs(props(counter).Centroid(1) - props(tracker).Centroid(1)) < pyw*2;
        end
        % find which elements are connectected to elem(tracker)
        bothCondition = cxCondition .* cyCondition;
        conditionIsTrue = find(bothCondition); 
        divLamb_temp = divLambMin(conditionIsTrue);  
        [~, iindx] = min(divLamb_temp);
        tmp_elem = props(conditionIsTrue(iindx));       
        % remove the overlapped boxes
        for test2 = length(conditionIsTrue):-1:1
           props(conditionIsTrue(test2)) = [];       
        end
        props(tracker) = tmp_elem;
        tracker = tracker + 1;
    end 
    strucFound  = strucFound + tracker-1;
    for counter = length(props):-1:1
        if(isempty(props(counter).Centroid)), props(counter)=[]; end
    end
%   plot to check sanity
%     for indx = 1:length(props)
%         plot(round(props(indx).Centroid(1)), round(props(indx).Centroid(2)),'or','MarkerFaceColor','r');        
%     end    
    
%     for indx = 1:length(props)
%         if(round(props(indx).Centroid(2)-pxw) > 0 &&  round(props(indx).Centroid(2)-pxw) < Nx && ...
%                 round(props(indx).Centroid(1)-pyw) > 0 && round(props(indx).Centroid(1)-pyw)  < Ny)
%             rectangle('Position',[round(props(indx).Centroid(1)-pyw), round(props(indx).Centroid(2)-pxw),...
%                 round(pyw*2), round(pxw*2)], 'EdgeColor','r','LineWidth',2);
%         end
%     end
%      get the quantities required from no overlapped boxes     
    
    for counter = 1:length(props)
    
        [box_stress_, ~, ~, ~] = getAvgBoxBotm2Blh(Nx, Ny, Nz, z_i, BLH, l_z,...
         bxw, byw,props(counter).Centroid(2), props(counter).Centroid(1), data.shear_stress);
         stressBox =  stressBox + box_stress_;
         
        [box_u, vpcnt, pxw, pyw] = getAvgBoxBotm2Blh(Nx, Ny, Nz, z_i, BLH, l_z,...
         bxw, byw,props(counter).Centroid(2), props(counter).Centroid(1), uu);
         uBox = uBox + box_u;
         
        [box_v, ~, ~, ~] = getAvgBoxBotm2Blh(Nx, Ny, Nz, z_i, BLH, l_z,...
         bxw, byw,props(counter).Centroid(2), props(counter).Centroid(1), vv);
         vBox =  vBox + box_v;
         
        [box_w, ~, ~, ~] = getAvgBoxBotm2Blh(Nx, Ny, Nz, z_i, BLH, l_z,...
         bxw, byw,props(counter).Centroid(2), props(counter).Centroid(1), ww);
         wBox =  wBox + box_w;
         
        [box_omega1, ~, ~, ~] = getAvgBoxBotm2Blh(Nx, Ny, Nz, z_i, BLH, l_z,...
         bxw, byw,props(counter).Centroid(2), props(counter).Centroid(1), omega1);
         omega1Box =  omega1Box + box_omega1;
         
        [box_omega2, ~, ~, ~] = getAvgBoxBotm2Blh(Nx, Ny, Nz, z_i, BLH, l_z,...
         bxw, byw,props(counter).Centroid(2), props(counter).Centroid(1), omega1);
         omega2Box =  omega2Box + box_omega2; 

         [box_omega3, ~, ~, ~] = getAvgBoxBotm2Blh(Nx, Ny, Nz, z_i, BLH, l_z,...
         bxw, byw,props(counter).Centroid(2), props(counter).Centroid(1), omega3);
         omega3Box =  omega3Box + box_omega3;  
         
         [box_pressure, ~, ~, ~] = getAvgBoxBotm2Blh(Nx, Ny, Nz, z_i, BLH, l_z,...
         bxw, byw,props(counter).Centroid(2), props(counter).Centroid(1), pp);
         pBox =  pBox + box_pressure;      
         
         [box_, ~, ~, ~] = getAvgBoxBotm2Blh(Nx, Ny, Nz, z_i, BLH, l_z,...
         bxw, byw,props(counter).Centroid(2), props(counter).Centroid(1), div_lamb_vec);
         divLambBox =  divLambBox + box_;                  
           
    end
%     hold off;
%     ylabel('Stream-wise'); xlabel('Span-wise');
    total_frames_analyzed = total_frames_analyzed + 1;
    ccc = ccc + 1;
end
avg_stress_slcth  = avg_stress_slcth ./total_frames_analyzed;
uBox = uBox./strucFound;
vBox = vBox./strucFound;
wBox = wBox./strucFound;
omega1Box = omega1Box./strucFound;
omega2Box = omega2Box./strucFound;
omega3Box = omega3Box./strucFound;
stressBox = stressBox./strucFound;
pBox = pBox./strucFound;
divLambBox = divLambBox./strucFound;
msg = ['Done, calculating'];
disp(msg);

save('chnl_divLambVecTrigger_box.mat','uBox','vBox','wBox','stressBox','pBox','omega1Box',...
    'omega2Box','omega3Box','divLambBox','strucFound','vpcnt','avg_stress_slcth','heightVec','frameVec','-v7.3');

toc

