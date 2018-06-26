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
heightVec = [2.0, 4.0, 7.0,10, 13, 17, 20, ...
    23, 26, 29, 32, 36, 38, 42, 45]; 
avg_stress_slcth  = zeros([Nx, Ny, length(heightVec)]);
%%%==========================================================================
%% get stress fields for frames and save them 
%for frameNo = frameVec 
%   stresFieldName = sprintf('shearStress%4.4i.mat',frameNo);
%   [shear_stress,~,~,~, u_star]= ...
%   getStress3d_(parentDir, Nx, Ny, Nz, dx, dy, dz, u_star_in, Ugal, Vgal,zo, l_r, z_i,frameNo);
%   save(stresFieldName','shear_stress','u_star','-v7.3');
%end
%%%==========================================================================
%disp(['Stress Calculation Done !']);

% create a box just to get box dimension
frameStr = sprintf('%4.4i',frameVec(1));
fn = [parentDir '/output/u_frame/u_frame'  frameStr   '.bin'];
fh = fopen(fn, 'r');
uu = fread(fh,'double');
uu = reshape(uu, [Nx,Ny,Nz]).*u_star_in + Ugal; 
fclose(fh);

tic
[box, vpcnt, pxw, pyw] = getAvgBoxBotm2Blh(Nx, Ny, Nz, z_i, BLH, l_z,...
	bxw, byw,round(Nx/2), round(Ny/2), uu);

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
    fname = sprintf('chnl_Omega_%4.4i.mat',frameNo);
    [omega1, omega2, omega3]=calc_vorticity_vector(uu,vv,ww,z_i,l_r,Nx,Ny,Nz,dz);
    save(fname, 'omega1', 'omega2','omega3','-v7.3');
    
    frameStr = sprintf('%4.4i',frameNo);
	fn = [parentDir '/output/p_frame/p_frame'  frameStr   '.bin'];
	fh = fopen(fn, 'r');
	pp= fread(fh,'double');
	pp = reshape(pp, [Nx,Ny,Nz]);   
    fclose(fh); 

    stresFieldName = sprintf('shearStress%4.4i.mat',frameNo);
	data = load(stresFieldName);
	u_star = data.u_star;
	shear_stress_field = data.shear_stress(:,:,1);    

	surf_pressure = pp(:,:,1);	
%   [N, binc] = hist(surf_pressure(:)./max(surf_pressure(:)),20);
% 	figure(); plot(binc, N./sum(N(:)));
	% calculation for high stress structures
	highStressBin = shear_stress_field(:) > mean(shear_stress_field(:)); %+ 3.0 * std(shear_stress_field(:)); % gives a binary array 
	highStressBin = reshape(highStressBin, size(squeeze(shear_stress_field)));
    
% 	figure(); imagesc(highPressBin); hold on; colormap parula;    

    avg_stress_slcth  = avg_stress_slcth  + data.shear_stress(:,:,heightVec);

	props = regionprops(highStressBin, 'Centroid',...
		'MajorAxisLength','MinorAxisLength', 'Orientation');
    count = length(props);
    stressMax = zeros(count,1);
    vpcnt = zeros(count,1);
    for indx = 1:count
        [box_stress, vpcnt, pxw, pyw] = getAvgBoxBotm2Blh(Nx, Ny, Nz, z_i, BLH, l_z,...
            bxw, byw, props(indx).Centroid(2), props(indx).Centroid(1), data.shear_stress);
        stressMax(indx) = max(max(box_stress(:,:,1)));
        vpcnt(indx) = vpcnt;
    end
    pressMax = zeros(count,1);
    for indx = 1:count
        [box_press, vpcnt, pxw, pyw] = getAvgBoxBotm2Blh(Nx, Ny, Nz, z_i, BLH, l_z,...
            bxw, byw, props(indx).Centroid(2), props(indx).Centroid(1), pp);
        pressMax(indx) = abs(max(max(box_press(:,:,1))));
    end
    tracker = 1;
    while(~isempty(props))
        if(tracker > length(props)), break; end;
        cxCondition = zeros([length(props) 1]);
        cyCondition = zeros([length(props) 1]);
        tmp_elem = props(tracker);
        for counter = 1:length(props)
            cxCondition(counter) = abs(props(counter).Centroid(2) - props(tracker).Centroid(2)) < pxw*2;
            cyCondition(counter) = abs(props(counter).Centroid(1) - props(tracker).Centroid(1)) < pyw*2;
        end
        % find which elements are connectected to elem(tracker)
        bothCondition = cxCondition .* cyCondition;
        conditionIsTrue = find(bothCondition);         
        for mm = conditionIsTrue                    
            if(abs(stressMax(mm)) > abs(stressMax(tracker))) 
                tmp_elem = props(mm);                            
            end        
        end        
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
msg = ['Done, calculating'];
disp(msg);

save('chnl_stressTrigger_gratherThanMean.mat','uBox','vBox','wBox','stressBox','pBox','omega1Box',...
    'omega2Box','omega3Box','strucFound','vpcnt','avg_stress_slcth','heightVec','frameVec','-v7.3');

toc

