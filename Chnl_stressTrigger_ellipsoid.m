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

vertStressFracProfile = zeros([length(heightVec) 1]);
vertAreaFracProfile = zeros([length(heightVec) 1]);
meanPackingFactor = 0;
strucFound = 0;
interrogationLvl = 2;

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
	shear_stress_field = data.shear_stress(:,:,:);    

	surf_pressure = pp(:,:,1);	
%   [N, binc] = hist(surf_pressure(:)./max(surf_pressure(:)),20);
% 	figure(); plot(binc, N./sum(N(:)));
	% calculation for high stress structures
	highStressBin = shear_stress_field(:) > mean(shear_stress_field(:)); %...
        % + 1.0 * std(shear_stress_field(:)); % gives a binary array 
	highStressBin = reshape(highStressBin, size(shear_stress_field));
    
	figure(); imagesc(highStressBin(:,:,interrogationLvl)); hold on; colormap gray; hold on;   

    avg_stress_slcth  = avg_stress_slcth  + data.shear_stress(:,:,heightVec);

	props = regionprops(highStressBin(:,:,interrogationLvl), 'Centroid',...
		'MajorAxisLength','MinorAxisLength', 'Orientation','BoundingBox');

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
        [avgStress, ~, ~, ~, ~,~] = getAvgEllipsoid(props(indx).Centroid(1),...
            props(indx).Centroid(2), props(indx).MajorAxisLength/2, props(indx).MinorAxisLength/2,...
            props(indx).Orientation, highStressBin, shear_stress_field, heightVec, Nx, Ny, Nz);
        stressMax(indx) = avgStress(interrogationLvl);    
    end
% % delete overlapped ellipses
    tracker = 1;
    while(~isempty(props))
        if(tracker > length(props)), break; end;
        cxCondition = zeros([length(props) 1]);
        cyCondition = zeros([length(props) 1]);
        for counter = 1:length(props)
            horz_HalfWidth_tracker = props(tracker).BoundingBox(3)/2;
            vert_HalfHeight_tracker = props(tracker).BoundingBox(4)/2;
            horz_HalfWidth_counter = props(counter).BoundingBox(3)/2;
            vert_HalfHeight_counter = props(counter).BoundingBox(4)/2;  
            horz_center_coord_tracker = props(tracker).BoundingBox(1)+ horz_HalfWidth_tracker;   
            vert_center_coord_tracker = props(tracker).BoundingBox(2)+ vert_HalfHeight_tracker;
            horz_center_coord_counter = props(counter).BoundingBox(1)+ horz_HalfWidth_counter;   
            vert_center_coord_counter = props(counter).BoundingBox(2)+ vert_HalfHeight_counter;                  
            cxCondition(counter) = abs(horz_center_coord_counter - horz_center_coord_tracker) < ...
                (horz_HalfWidth_tracker + horz_HalfWidth_counter);
            cyCondition(counter) = abs(vert_center_coord_counter - vert_center_coord_tracker ) < ...
                (vert_HalfHeight_tracker + vert_HalfHeight_counter);
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
           
        [~, detectedElemNo, totalElem, vector_lin, ~,~] = getAvgEllipsoid(props(indx).Centroid(1),...
            props(indx).Centroid(2), props(indx).MajorAxisLength/2, props(indx).MinorAxisLength/2,...
            props(indx).Orientation, highStressBin, shear_stress_field, heightVec, Nx, Ny, Nz);
        for kk = 1:length(heightVec)
            if(kk ==1) 
                veclin = vector_lin(1:sum(detectedElemNo(1:kk)));
            else
                veclin = vector_lin(sum(detectedElemNo(1:kk-1))+1:sum(detectedElemNo(1:kk)));
            end
            [vert_c horz_c page_c] = ind2sub(size(shear_stress_field),veclin);
                        
            if(heightVec(kk) == interrogationLvl), plot(horz_c, vert_c,'.b'); end
            if(heightVec(kk) == interrogationLvl) 
                page_c 
            end
        end
    end
    for indx = 1:length(props)
        plot(round(props(indx).Centroid(1)), round(props(indx).Centroid(2)),'or','MarkerFaceColor','r');   
        txt = num2str(indx);
        text(round(props(indx).Centroid(1)),round(props(indx).Centroid(2)),txt,'HorizontalAlignment','right',...
            'Color','yellow');     
    end    
     
   
%% =======================================================================    
    for kk = 1:length(props) 
        [avgStress, detectedElemNo, totalElem, vector_lin, dataVec,avgNoFilt] = getAvgEllipsoid(props(indx).Centroid(1),...
            props(indx).Centroid(2), props(indx).MajorAxisLength/2, props(indx).MinorAxisLength/2,...
            props(indx).Orientation, highStressBin, shear_stress_field, heightVec, Nx, Ny, Nz);
        vertStressFracProfile = vertStressFracProfile  +  (avgStress./avgNoFilt);
        vertAreaFracProfile = vertAreaFracProfile + (detectedElemNo./totalElem);
    end
    toc
end
vertStressFracProfile = vertStressFracProfile./length(frameVec);
vertAreaFracProfile = vertAreaFracProfile./length(frameVec);
msg = ['Done, calculating'];
disp(msg);
figure();
plot(vertAreaFracProfile, vertStressFracProfile);



