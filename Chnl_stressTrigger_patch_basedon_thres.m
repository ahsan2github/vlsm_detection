%%
clear; clc; close all;
tic;

frameVec=[10,11,12,13] %,14,15,16,17,18,19];
bw_horz = 1.5; bw_vert = 10; bzh = 1;
width = 3.0; height = width / 1.1; formatEng = '-depsc';
SS ='';
parentDir = '/raid/home/mohammad/NeutralChannel_2048x2048x64_free';
%parentDir = '/glade/scratch/khanm/NeutralChannel_2048x2048x64_free';
%parentDir = '/scratch/ibrix/chpc_gen/u0699351/NeutralChannel_2048x2048x64_free';
outputDir = [parentDir '/output/'];
readinputs(parentDir);
u_star_in = u_star; clear u_star;
BLH = 1500.0;
heightVec = [2.0, 4.0, 7.0,10, 13, 17, 20, ...
    23, 26, 29, 32, 36, 38, 42, 45]; 
avg_stress_slcth  = zeros([Nx, Ny, length(heightVec)]);
%%%=================================================================================================
%% get stress fields for frames and save them 
%for frameNo = frameVec 
%   stresFieldName = sprintf('shearStress%4.4i.mat',frameNo);
%   [shear_stress,~,~,~, u_star]= ...
%   getStress3d_(parentDir, Nx, Ny, Nz, dx, dy, dz, u_star_in, Ugal, Vgal,zo, l_r, z_i,frameNo);
%   save(stresFieldName','shear_stress','u_star','-v7.3');
%end
%%%=================================================================================================
%disp(['Stress Calculation Done !']);

% create a box just to get box dimension
frameStr = sprintf('%4.4i',frameVec(1));
fn = [parentDir '/output/u_frame/u_frame'  frameStr   '.bin'];
fh = fopen(fn, 'r');
uu = fread(fh,'double');
uu = reshape(uu, [Nx,Ny,Nz]).*u_star_in + Ugal; 
fclose(fh);

tic
[box, vpcnt, pxw, pyw, ~, ~, ~, ~] = getAvgBoxBotm2Blh(Nx, Ny, Nz, z_i, BLH, l_z,...
	bw_horz, bw_vert,round(Nx/2), round(Ny/2), uu);

toc

[bnx, bny, bnz] = size(box);
uBox = zeros([bnx, bny, bnz]);
vBox = zeros([bnx, bny, bnz]);
wBox = zeros([bnx, bny, bnz]);
omega1Box = zeros([bnx, bny, bnz]);
omega2Box = zeros([bnx, bny, bnz]);
omega3Box = zeros([bnx, bny, bnz]);
stressBox = zeros([bnx, bny, bnz]);
pBox = zeros([bnx, bny, bnz]);
lambda2Box = zeros([bnx, bny, bnz]);
divLambVecBox = zeros([bnx, bny, bnz]);
meanPackingFactor = 0;
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
    
    frameStr = sprintf('%4.4i',frameNo);
	fn = [parentDir '/output/lambda2_frame/lambda2_frame'  frameStr   '.bin'];
	fh = fopen(fn, 'r');
	lambda2 = fread(fh,'double');
	lambda2 = reshape(lambda2, [Nx,Ny,Nz]);
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
        
    divLambFName = sprintf('chnl_DivLambVec%4.4i.bin',frameNo);
    ff = fopen(divLambFName, 'r');
    div_lamb_vec = fread(ff,'double');
    fclose(ff);    
    div_lamb_vec = reshape( div_lamb_vec, [Nx, Ny, Nz]);    

    stresFieldName = sprintf('shearStress%4.4i.mat',frameNo);
	data = load(stresFieldName);
	u_star = data.u_star;
	shear_stress_field = data.shear_stress(:,:,2);    

	surf_pressure = pp(:,:,1);	
%   [N, binc] = hist(surf_pressure(:)./max(surf_pressure(:)),20);
% 	figure(); plot(binc, N./sum(N(:)));
	% calculation for high stress structures
	highStressBin = shear_stress_field(:) > mean(shear_stress_field(:)); %...
        % + 1.0 * std(shear_stress_field(:)); % gives a binary array 
	highStressBin = reshape(highStressBin, size(squeeze(shear_stress_field)));
    
	figure(); imagesc(highStressBin); hold on; colormap gray; hold on;   

    avg_stress_slcth  = avg_stress_slcth  + data.shear_stress(:,:,heightVec);

	props = regionprops(highStressBin, 'Centroid',...
		'MajorAxisLength','MinorAxisLength', 'Orientation');

    mAxisLengths = [props(1:end).MajorAxisLength]; 
    % select major axis length ranges 
    thresMajorAxisLength_low = floor(Nx * (BLH*7)/(2*pi*z_i)); %
    thresMajorAxisLength_high = floor(Nx * (BLH*13)/(2*pi*z_i)); %
    c1 = mAxisLengths >= thresMajorAxisLength_low;
    c2 = mAxisLengths <= thresMajorAxisLength_high;   
    crit = find(c1 .* c2);
    props = props(crit);
    count = length(props);
    stressMax = zeros(count,1);
    vpcnt = zeros(count,1);
    for indx = 1:count
        [box_stress, vpcnt, pxw, pyw, ~,~,~,~] = getAvgBoxBotm2Blh(Nx, Ny, Nz, z_i, BLH, l_z,...
            bw_horz, bw_vert, props(indx).Centroid(1), props(indx).Centroid(2), data.shear_stress);
        stressMax(indx) = mean(mean(box_stress(:,:,2)));
        vpcnt(indx) = vpcnt;
    end
    pressMax = zeros(count,1);
    for indx = 1:count
        [box_press, vpcnt, pxw, pyw,~,~,~,~] = getAvgBoxBotm2Blh(Nx, Ny, Nz, z_i, BLH, l_z,...
            bw_horz, bw_vert, props(indx).Centroid(1), props(indx).Centroid(2), pp);
        pressMax(indx) = abs(mean(mean(box_press(:,:,2))));
    end
    tracker = 1;
    while(~isempty(props))
        if(tracker > length(props)), break; end;
        cxCondition = zeros([length(props) 1]);
        cyCondition = zeros([length(props) 1]);        
        for counter = 1:length(props)
            cxCondition(counter) = abs(props(counter).Centroid(1) - props(tracker).Centroid(1)) < pxw*2;
            cyCondition(counter) = abs(props(counter).Centroid(2) - props(tracker).Centroid(2)) < pyw*2;
        end
        % find which elements are connectected to elem(tracker)
        bothCondition = cxCondition .* cyCondition;
        conditionIsTrue = find(bothCondition); 
        stressmax_temp = stressMax(conditionIsTrue);  
        [~, iindx] = max(stressmax_temp);
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
%% ============== plot to check sanity  ===================================
    phi = linspace(0,2*pi,50);
    cosphi = cos(phi);
    sinphi = sin(phi);
    totalShearBox = zeros([length(props) 1]);    
    totalShearEllipse = zeros([length(props) 1]);
    areaFracBox = zeros([length(props) 1]);
    areaFracEllipse = zeros([length(props) 1]);  
    packingFactor = zeros([length(props) 1]);
    for indx = 1:length(props)
        xbar = props(indx).Centroid(1);
        ybar = props(indx).Centroid(2);

        a = props(indx).MajorAxisLength/2;
        b = props(indx).MinorAxisLength/2;
        theta = pi*props(indx).Orientation/180;    
        R = [ cos(theta)   sin(theta)
             -sin(theta)   cos(theta)];

        xy = [a*cosphi; b*sinphi];
        xy = R*xy;

        x = xy(1,:) + xbar;
        y = xy(2,:) + ybar;
        plot(x,y,'g','LineWidth',2);                    
        [avg, detectedElemNo, totalElem, vector_x, vector_y] = avgInEllipse(props(indx).Centroid(1),...
        props(indx).Centroid(2), a, b,...
        props(indx).Orientation, highStressBin, highStressBin);
        plot(vector_y, vector_x,'.b');
        plot(xbar, ybar,'+g'); 
        
        [avgStress, nDataPixels, nTotalPixels , ~, ~] = ...
                avgInEllipse(props(indx).Centroid(1),  props(indx).Centroid(2), props(indx).MajorAxisLength/2, ...
                props(indx).MinorAxisLength/2, props(indx).Orientation,...
                highStressBin,  shear_stress_field);
        totalShearEllipse(indx) =  nDataPixels * avgStress;
        packingFactor(indx) = nDataPixels/nTotalPixels;
        areaFracEllipse(indx) = nDataPixels/(Nx*Ny);        

    end
    for indx = 1:length(props)
        plot(round(props(indx).Centroid(1)), round(props(indx).Centroid(2)),'or','MarkerFaceColor','r');   
        txt = num2str(indx);
        text(round(props(indx).Centroid(1)),round(props(indx).Centroid(2)),txt,'HorizontalAlignment','right',...
            'Color','yellow');     
    end    
     
    
    for indx = 1:length(props)
        if(round(props(indx).Centroid(1)-pxw) > 0 &&  round(props(indx).Centroid(1)-pxw) < Ny && ...
            round(props(indx).Centroid(2)-pyw) > 0 && round(props(indx).Centroid(2)-pyw)  < Nx)
            
            [box_stress, vpcnt, pxw, pyw,xst,xed,yst,yed] = getAvgBoxBotm2Blh(Nx, Ny, Nz, z_i, BLH, l_z,...
                bw_horz, bw_vert, props(indx).Centroid(1), props(indx).Centroid(2), data.shear_stress);
            totalShearBox(indx) = sum(sum(box_stress(:,:,2)));
            areaFracBox(indx) = numel(box_stress(:,:,2))/(Nx*Ny) * packingFactor(indx);            
            
%             rectangle('Position',[round(props(indx).Centroid(1)-pxw), round(props(indx).Centroid(2)-pyw),...
%                 round(pxw*2), round(pyw*2)], 'EdgeColor','r','LineWidth',2);   
            rectangle('Position',[xst, yst,...
                xed-xst, yed-yst], 'EdgeColor','r','LineWidth',2);
        end
    end   


    
%% =======================================================================    
     
    for ii=1:strucFound
        [avgStress, nDataPixels, nTotalPixels , ~, ~] = ...
        avgInEllipse(props(indx).Centroid(1),  props(indx).Centroid(2), props(indx).MajorAxisLength/2, ...
        props(indx).MinorAxisLength/2, props(indx).Orientation,...
        highStressBin,  shear_stress_field);
        meanPackingFactor = meanPackingFactor +(nDataPixels/nTotalPixels);
        
    end

% %      get the quantities required from non overlapped boxes     
    
    for counter = 1:length(props)
        
        [box_stress_, ~, ~, ~,~,~,~,~] = getAvgBoxBotm2Blh(Nx, Ny, Nz, z_i, BLH, l_z,...
         bw_horz, bw_vert,props(counter).Centroid(1), props(counter).Centroid(2), data.shear_stress);
         stressBox =  stressBox + box_stress_;
         
        [box_u, vpcnt, pxw, pyw,~,~,~,~] = getAvgBoxBotm2Blh(Nx, Ny, Nz, z_i, BLH, l_z,...
         bw_horz, bw_vert,props(counter).Centroid(1), props(counter).Centroid(2), uu);
         uBox = uBox + box_u;
         
        [box_v, ~, ~, ~,~,~,~,~] = getAvgBoxBotm2Blh(Nx, Ny, Nz, z_i, BLH, l_z,...
         bw_horz, bw_vert,props(counter).Centroid(1), props(counter).Centroid(2), vv);
         vBox =  vBox + box_v;
         
        [box_w, ~, ~, ~,~,~,~,~] = getAvgBoxBotm2Blh(Nx, Ny, Nz, z_i, BLH, l_z,...
         bw_horz, bw_vert,props(counter).Centroid(1), props(counter).Centroid(2), ww);
         wBox =  wBox + box_w;
         
        [box_omega1, ~, ~, ~,~,~,~,~] = getAvgBoxBotm2Blh(Nx, Ny, Nz, z_i, BLH, l_z,...
         bw_horz, bw_vert,props(counter).Centroid(1), props(counter).Centroid(2), omega1);
         omega1Box =  omega1Box + box_omega1;
         
        [box_omega2, ~, ~, ~,~,~,~,~] = getAvgBoxBotm2Blh(Nx, Ny, Nz, z_i, BLH, l_z,...
         bw_horz, bw_vert,props(counter).Centroid(1), props(counter).Centroid(2), omega1);
         omega2Box =  omega2Box + box_omega2; 
      
         [box_omega3, ~, ~, ~,~,~,~,~] = getAvgBoxBotm2Blh(Nx, Ny, Nz, z_i, BLH, l_z,...
         bw_horz, bw_vert,props(counter).Centroid(1), props(counter).Centroid(2), omega3);
         omega3Box =  omega3Box + box_omega3;  
         
         [box_pressure, ~, ~, ~,~,~,~,~] = getAvgBoxBotm2Blh(Nx, Ny, Nz, z_i, BLH, l_z,...
         bw_horz, bw_vert,props(counter).Centroid(1), props(counter).Centroid(2), pp);
         pBox =  pBox + box_pressure;             
               
         [box_, ~, ~, ~,~,~,~,~] = getAvgBoxBotm2Blh(Nx, Ny, Nz, z_i, BLH, l_z,...
         bw_horz, bw_vert,props(counter).Centroid(1), props(counter).Centroid(2), lambda2);
         lambda2Box =  lambda2Box + box_;   
         
         [box_, ~, ~, ~,~,~,~,~] = getAvgBoxBotm2Blh(Nx, Ny, Nz, z_i, BLH, l_z,...
         bw_horz, bw_vert,props(counter).Centroid(1), props(counter).Centroid(2), div_lamb_vec);
         divLambVecBox =  divLambVecBox + box_;            
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
lambda2Box = lambda2Box./strucFound;
divLambVecBox = divLambVecBox./strucFound;
meanPackingFactor = meanPackingFactor/strucFound;
msg = ['Done, calculating'];
disp(msg);

save('chnl_stressTrigger_gratherThanMean.mat','uBox','vBox','wBox','stressBox','pBox','omega1Box','divLambVecBox',...
    'omega2Box','omega3Box','lambda2Box','strucFound','vpcnt','avg_stress_slcth','heightVec','frameVec',...
    'meanPackingFactor','-v7.3');

toc

