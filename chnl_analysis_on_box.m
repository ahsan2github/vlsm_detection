  
clear; clc;
tic;
envelope = 2:1:8;
strucFound = 0;
width = 3.0; height = width / 1.1; formatEng = '-depsc';
SS ='';
parentDir = '/raid/home/mohammad/NeutralChannel_2048x2048x64_free';
%parentDir = '/glade/scratch/khanm/NeutralChannel_2048x2048x64_free';
%parentDir = '/scratch/ibrix/chpc_gen/u0699351/NeutralChannel_2048x2048x64_free';
outputDir = [parentDir '/output/'];
readinputs(parentDir);
BLH = 1500.0;
zo = 0.1;
delta  = (dx * dy * dz)^(1/3.0);
data = load('chnl_stressTrigger_gratherThanMean.mat');
heightVec = data.heightVec;
areaFracVec_low = zeros(length(heightVec), length(envelope));
fluxFracVec_low = zeros(length(heightVec),length(envelope));
uu_frac_low = zeros(length(heightVec), length(envelope));
minorAxislength_low_avg= zeros(length(heightVec), length(envelope));
majorAxislength_low_avg = zeros(length(heightVec), length(envelope));
noRect_low = zeros(length(heightVec),length(envelope));
nn = length(heightVec);

for outerloop = 1:nn    
    for slide = 1:length(envelope)
        box_area_frac = 0;
        box_uu_sum = 0; 
        box_uu_frac = 0.0;   
        box_stress_frac = 0.0;
        minorAxisLength_avg = 0;
        majorAxisLength_avg = 0;
        clearvars elem;
        uu = data.uBox;
        u_lvl = uu(:,:,heightVec(outerloop));
        lowStreakBin_u = u_lvl(:) < 0.0; % gives a binary array 
        lowStreakBin_u = reshape(lowStreakBin_u, size(u_lvl));        
        % identify the image regions with data
        props = regionprops(lowStreakBin_u, 'Centroid',...
            'MajorAxisLength','MinorAxisLength', 'Orientation');
        [maxLength, mxi] = max([props(:).MajorAxisLength]);
        mxcx = floor(props(mxi).Centroid(1)); mxcy = floor(props(mxi).Centroid(2));
        mAxisLengths = [props(1:end).MajorAxisLength]; 
        % select major axis length ranges 
        thresMajorAxisLength_low = floor(Nx * (BLH*envelope(slide)*0.5/(2*pi*z_i))); %
        thresMajorAxisLength_high = floor(Nx * (BLH*envelope(slide)*1.5/(2*pi*z_i))); %
        c1 = mAxisLengths >= thresMajorAxisLength_low;
        c2 = mAxisLengths <= thresMajorAxisLength_high;   
        crit = find(c1 .* c2);
        shear_stress_field = data.stressBox(:,:, heightVec(outerloop));
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
                box_area_frac = box_area_frac + (nDataPixels/numel(u_lvl));
                box_stress_frac = box_stress_frac + ((avgStress* nDataPixels)/(mean(shear_stress_field(:))*numel(u_lvl)));
                % get minor axis length as multiple of BLH
                minorAxisLength_avg = minorAxisLength_avg + (props(indx).MinorAxisLength...
                    * 2* pi*z_i/size(lowStreakBin_u,2)/BLH);
                majorAxisLength_avg = majorAxisLength_avg + (props(indx).MajorAxisLength...
                    * 2* pi*z_i/size(lowStreakBin_u,2)/BLH);
                for tt = 1:length(vector_x)
                    box_uu_sum = box_uu_sum + u_lvl(vector_x(tt), vector_y(tt))^2;
                end
                elmCount = elmCount + 1;
            end
            strucFound = strucFound + elmCount;
            box_uu_frac = box_uu_frac + box_uu_sum/sum(u_lvl(:).^2);
        else
            box_area_frac = box_area_frac + 0.0;
            box_stress_frac = box_stress_frac  + 0.0; 
            box_uu_sum = box_uu_sum + 0;
            minorAxisLength_avg = minorAxisLength_avg + 0.0;
            majorAxisLength_avg = majorAxisLength_avg + 0.0;
            box_uu_frac = box_uu_frac + 0.0;          
            msg =['Under given criteria no area was found, Structure Length, :'  ...
                    num2str(envelope(slide)), ', height:' num2str(heightVec(outerloop))];
            disp(msg);
        end
        
        minorAxisLength_avg = minorAxisLength_avg/strucFound;
        majorAxisLength_avg = majorAxisLength_avg /strucFound;
        areaFracVec_low(outerloop, slide) = box_area_frac;
        fluxFracVec_low(outerloop, slide) = box_stress_frac;
        uu_frac_low(outerloop, slide) = box_uu_frac;
        noRect_low(outerloop, slide)= strucFound;
        minorAxislength_low_avg(outerloop, slide) = minorAxisLength_avg;
        majorAxislength_low_avg(outerloop, slide) = majorAxisLength_avg;
        msg = ['Done, calculating heigh :' num2str(heightVec(outerloop))];
	    %disp(msg);
    end
end

%% overall low momentum stats
uu = data.uBox;
stressFrac_box_low = zeros([length(heightVec) 1]);
areaFrac_box_low = zeros([length(heightVec) 1]);
for outerloop = 1:length(heightVec)    
    u_lvl = uu(:,:,heightVec(outerloop));
    lowStreakBin_u = u_lvl(:) < 0.0; % gives a binary array 
    crit = find(lowStreakBin_u); % find the locations of TRUE
    shear_stress_field = data.stressBox(:,:, heightVec(outerloop));
    stressFrac_box_low(outerloop) = (mean(shear_stress_field(crit))*numel(crit))/...
        (mean(shear_stress_field(:))*numel(u_lvl));
    areaFrac_box_low(outerloop) = numel(crit)/numel(u_lvl);
end
save('chnl_analysis_on_box_stressTrig.mat','areaFracVec_low','fluxFracVec_low','noRect_low',...
    'majorAxislength_low_avg','minorAxislength_low_avg','heightVec', 'uu_frac_low',...
    'envelope','strucFound','stressFrac_box_low','areaFrac_box_low','-v7.3');
toc
