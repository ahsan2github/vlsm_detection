% % this version finds structures within band (0.75X-1.25X)delta
% % where X ranges from 2-24

% % change heightVector
% % change frameVec
% % change output mat file name 
  
clear; clc;
tic;
heightVec = [];  
envelope = 2:2:24;
frameVec=[];
areaFracVec_low = zeros(length(heightVec), length(envelope));
fluxFracVec_low = zeros(length(heightVec),length(envelope));
noRect_low = zeros(length(heightVec),length(envelope));
uu_frac_low = zeros(length(heightVec), length(envelope));
minorAxislength_low_avg= zeros(length(heightVec), length(envelope));
majorAxislength_low_avg = zeros(length(heightVec), length(envelope));
orientation_low_avg= zeros(length(heightVec), length(envelope));

%parentDir = '/raid/home/mohammad/NeutralChannel_2048x2048x64_free';
%parentDir = '/glade/scratch/khanm/NeutralChannel_2048x2048x64_free';
parentDir = '/scratch/ibrix/chpc_gen/u0699351/NeutralChannel_2048x2048x64_free';
outputDir = [parentDir '/output/'];
readinputs(parentDir);
BLH = 1500.0;
zo = 0.1;
u_star_in = u_star; clear u_star;
% get stress fields for frames and save them 
for frameNo = frameVec 
    stresFieldName = sprintf('chnl_shearStress%4.4i.bin',frameNo);
    fid = fopen(stresFieldName, 'w');    
    [shear_stress,~,~,~, u_star]= ...
    getStress3d_(parentDir, Nx, Ny, Nz, dx, dy, dz, u_star_in, Ugal, Vgal,zo, l_r, z_i,frameNo);
    fwrite(fid,shear_stress,'double');
    fclose(fid);
end
%==========================================================================
disp ('Stress Calculation Done');
% calculation for low momentum

for outerloop = 1:length(heightVec)   
    for slide = 1:length(envelope)
        total_frames_analyzed = 0;
        struc2dfound = 0;
        box_area_frac = 0;
        box_uu_sum = 0; 
        box_uu_frac = 0.0;   
        box_stress_frac = 0.0;
        minorAxisLength_avg = 0;
        majorAxisLength_avg = 0;
        orientation_avg = 0;
        clearvars elem;
        for frameNo = frameVec 
            frameStr = sprintf('%4.4i',frameNo);
            fn = [parentDir '/output/u_frame/u_frame'  frameStr   '.bin'];
            fh = fopen(fn, 'r');
            uu = fread(fh,'double');
            uu = reshape(uu, [Nx,Ny,Nz]).*u_star_in + Ugal;     
            fclose(fh);        
            stresFieldName = sprintf('chnl_shearStress%4.4i.bin',frameNo);
            ff = fopen(stresFieldName, 'r');
            data_shear_stress = fread(ff,'double');
            fclose(ff);
            data_shear_stress = reshape(data_shear_stress, [Nx, Ny, Nz]);
            shear_stress_field = data_shear_stress(:,:, heightVec(outerloop));
            turb_u_lvl = uu(:,:,heightVec(outerloop)) - mean(mean(uu(:,:,heightVec(outerloop))));
            %  Fourier Filtering       
            u_filtered = cuttoff_2D_multiF(turb_u_lvl, 2 * pi * z_i / heightVec(outerloop) / dz);
            u_filtered  = squeeze(u_filtered );
            % calculation for low momentum structures
            lowStreakBin_u = u_filtered(:) < (mean(u_filtered(:))-3.0*std(u_filtered(:))); % gives a binary array 
            lowStreakBin_u = reshape(lowStreakBin_u, size(u_filtered));        
            % identify the image regions with data
            props = regionprops(lowStreakBin_u, 'Centroid',...
                'MajorAxisLength','MinorAxisLength', 'Orientation');
            [maxLength, mxi] = max([props(:).MajorAxisLength]);
            mxcx = floor(props(mxi).Centroid(1)); mxcy = floor(props(mxi).Centroid(2));
            mAxisLengths = [props(1:end).MajorAxisLength]; 
            % select major axis length ranges 
            thresMajorAxisLength_low = floor(Nx * (BLH*envelope(slide)*0.75)/(2*pi*z_i)); %
            thresMajorAxisLength_high = floor(Nx * (BLH*envelope(slide)*1.25)/(2*pi*z_i)); %
            c1 = mAxisLengths >= thresMajorAxisLength_low;
            c2 = mAxisLengths <= thresMajorAxisLength_high;   
            crit = find(c1 .* c2);
            count = numel(crit);            
            if(count >= 1) 
                elmCount = 0;
                for indx = crit
                    cx = props(indx).Centroid(1); cy = props(indx).Centroid(2);
                    [avgStress, nDataPixels, nTotalPixels , ~, ~] = ...
                            avgInEllipse(cx, cy, props(indx).MajorAxisLength/2, ...
                            props(indx).MinorAxisLength/2, props(indx).Orientation,...
                            lowStreakBin_u,  shear_stress_field);
                    trueAreaFrac =  nDataPixels / nTotalPixels;   
                    [avgU, ~, ~, vector_x, vector_y] = avgInEllipse(cx, cy,...
                            props(indx).MajorAxisLength/2, props(indx).MinorAxisLength/2, props(indx).Orientation,...
                            lowStreakBin_u, uu(:,:,heightVec(outerloop)));                   
                    box_area_frac = box_area_frac + (nDataPixels/Nx/Ny);
                    box_stress_frac = box_stress_frac + ((avgStress* nDataPixels)/(mean(shear_stress_field(:))*Nx*Ny));
                    % get major / minor axis length as multiple of BLH
                    minorAxisLength_avg = minorAxisLength_avg + (props(indx).MinorAxisLength...
                        * 2* pi*z_i/size(lowStreakBin_u,2)/BLH);
                    majorAxisLength_avg = majorAxisLength_avg + (props(indx).MajorAxisLength...
                        * 2* pi*z_i/size(lowStreakBin_u,2)/BLH);
                    if(props(indx).Orientation < 0)
                        ang_ort = abs(props(indx).Orientation) + 90;
                    else
                        ang_ort = abs(props(indx).Orientation);
                    end    
                    orientation_avg = orientation_avg + ang_ort;
                    for tt = 1:length(vector_x)
                        box_uu_sum = box_uu_sum + turb_u_lvl(vector_x(tt), vector_y(tt))^2;
                    end
                    elmCount = elmCount + 1;
                end
                struc2dfound = struc2dfound + elmCount;
            else
                box_area_frac = box_area_frac + 0.0;
                box_stress_frac = box_stress_frac  + 0.0; 
                box_uu_sum = box_uu_sum + 0;
                minorAxisLength_avg = minorAxisLength_avg + 0.0;
                majorAxisLength_avg = majorAxisLength_avg + 0.0;
                orientation_avg = orientation_avg + 0.0;                
                msg =['Under given criteria no area was found, ' 'frame: ' ...
                        num2str(frameNo)  '  srucLength : '  num2str(envelope(slide)) ...
			', height:' num2str(heightVec(outerloop))];
	        disp(msg)
            end
            total_frames_analyzed = total_frames_analyzed + 1; 
            % convert box_uu_sum to fractional energy
            box_uu_frac = box_uu_frac + (box_uu_sum /sum(turb_u_lvl(:).^2));  
        end

        areaFracVec_low(outerloop, slide) = box_area_frac/total_frames_analyzed;
        fluxFracVec_low(outerloop, slide) = box_stress_frac/total_frames_analyzed; 
        uu_frac_low(outerloop, slide) = box_uu_frac / total_frames_analyzed;
        noRect_low(outerloop, slide)= struc2dfound/total_frames_analyzed;
        minorAxislength_low_avg(outerloop, slide) = minorAxisLength_avg/struc2dfound;
        majorAxislength_low_avg(outerloop, slide) = majorAxisLength_avg / struc2dfound;
        orientation_low_avg(outerloop, slide) = orientation_avg/struc2dfound;       

    end
    msg = ['Done, calculating heigh :' num2str(heightVec(outerloop))];
    disp(msg)
    toc
end

save('chnl_areaFrac_stressFrac_hxx_xx_band_3sig.mat','areaFracVec_low','fluxFracVec_low','noRect_low',...
        'majorAxislength_low_avg','minorAxislength_low_avg','orientation_low_avg',...
        'heightVec', 'uu_frac_low','envelope','frameVec','-v7.3');
toc
