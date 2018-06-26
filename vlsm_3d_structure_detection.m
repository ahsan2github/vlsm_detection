%%
clear; clc; close all;
tic;

frameVec= [39,40,41,42,43]; 
width = 3.0; height = width / 1.1; formatSpec = '-depsc';
SS ='';
parentDir = '/raid/home/mohammad/Ekman_Ugeo2_2048x2048x64_Lz750';
outputDir = [parentDir '/output/'];
readinputs(parentDir);
BLH = 550.69;
u_star_in = u_star; clear u_star;
heightVec = 1:Nz; 

%%%=================================================================================================
%% get stress fields for frames and save them 
%for frameNo = frameVec 
%   stresFieldName = sprintf('shearStress%4.4i.mat',frameNo);
%   [shear_stress,~,~,~, u_star]= ...
%   getStress3d_(parentDir, Nx, Ny, Nz, dx, dy, dz, u_star_in, Ugal, Vgal,zo, l_r, z_i,frameNo);
%   save(stresFieldName','ek02shear_stress','u_star','-v7.3');
%end
%%%=================================================================================================
%disp(['Stress Calculation Done !']);

strucFound = 0;
interrogationLvl = round(BLH*0.5/dz);
% ======================================================================
fixedMajorAxisLength = floor(Nx * (BLH*15.0)/(2*pi*z_i));
fixedMinorAxisLength = floor(Nx * (BLH*5.0)/(2*pi*z_i));
% ======================================================================
a = round(fixedMajorAxisLength/2);
b = round(fixedMinorAxisLength/2);
dataBox_u_frame = zeros([length(-b:b) length(-a:a) length(heightVec)]);

for frameNo = frameVec 
	frameStr = sprintf('%4.4i',frameNo);
	fn = [parentDir '/output/u_frame/u_frame'  frameStr   '.bin'];
	fh = fopen(fn, 'r');
	uu = fread(fh,'double');
	uu = reshape(uu, [Nx,Ny,Nz]).*u_star_in + Ugal;
	uu=removeLevelMean(uu);
    fclose(fh);  

    stresFieldName = sprintf('ek02_shearStress%4.4i.mat',frameNo);
	data = load(stresFieldName);
	u_star = data.u_star;
	shear_stress_field = data.shear_stress(:,:,:);    
    % get 2D binary field
    criteriaLvl =  uu(:,:,interrogationLvl);
    criteriaBinLvl = criteriaLvl(:) <  mean(criteriaLvl(:)) - 0.75*std(criteriaLvl(:));
    criteriaBinLvl = reshape(criteriaBinLvl, squeeze(size(criteriaLvl)));
    % get 3D binary field
    criteriaBin = uu(:) < mean(criteriaLvl(:)) - 0.75*std(criteriaLvl(:));
    criteriaBin = reshape(criteriaBin, size(shear_stress_field));
    
	figure(); imagesc(criteriaBinLvl); hold on; colormap gray; hold on;   

	props = regionprops(criteriaBinLvl, 'Centroid',...
		'MajorAxisLength','MinorAxisLength', 'Orientation','BoundingBox');

    mAxisLengths = [props(1:end).MajorAxisLength]; 
    % ======================================================================
    % select major axis length ranges 
    thresMajorAxisLength_low = floor(Nx * (BLH*1)/(2*pi*z_i)); %
    thresMajorAxisLength_high = floor(Nx * (BLH*15)/(2*pi*z_i)); %
    % ======================================================================
    c1 = mAxisLengths >= thresMajorAxisLength_low;
    c2 = mAxisLengths <= thresMajorAxisLength_high;   
    crit = find(c1 .* c2);
    props = props(crit);
    if(isempty(props)), error('no structure found under given criteria'); end
    count = length(props);
    % delete overlapped ellipses
    % ======================================================================
    [props, horzcc, vertcc, horzDim, vertDim, horzBLC, vertBLC] = ...
        removeOverlappedBoxes(props,a,b,criteriaBinLvl, ....
            shear_stress_field(:,:,interrogationLvl))   ;
    disp(['Non overlapped Structure isolation done, elapsed time : ' num2str(toc)]);
    
    
%% ============== plot to check sanity  ===================================
    phi = linspace(0,2*pi,360);
    cosphi = cos(phi);
    sinphi = sin(phi);
    for indx = 1:length(props)
        xbar = props(indx).Centroid(1);
        ybar = props(indx).Centroid(2);

        a = fixedMajorAxisLength/2;
        b = fixedMinorAxisLength/2;
        theta = pi*props(indx).Orientation/180;    
        R = [ cos(theta)   sin(theta)
             -sin(theta)   cos(theta)];

        xy = [a*cosphi; b*sinphi];
        xy = R*xy;

        x = xy(1,:) + xbar;
        y = xy(2,:) + ybar;
        plot(x,y,'g','LineWidth',2);
        
        [avg, detectedElemNo, totalElem, vector_h, vector_v] = avgInEllipse1(props(indx).Centroid(1),...
           props(indx).Centroid(2), a, b,...
           props(indx).Orientation, criteriaBinLvl, criteriaBinLvl);
       
        [horz_org , vert_org,horz_dim, vert_dim] = getBBox(props(indx).Centroid(1),...
                props(indx).Centroid(2), a, b, props(indx).Orientation, Nx, Ny);        
        plot(vector_h, vector_v,'.b');
        plot(xbar, ybar,'+g');  
        rectangle('Position', [horz_org , vert_org,horz_dim, vert_dim], 'EdgeColor','r','LineWidth',2);
    end
    for indx = 1:length(props)
        plot(round(props(indx).Centroid(1)), round(props(indx).Centroid(2)),'or','MarkerFaceColor','r');   
        txt = num2str(indx);
        text(round(props(indx).Centroid(1)),round(props(indx).Centroid(2)),txt,'HorizontalAlignment','right',...
            'Color','yellow');     
    end    
     
   
%% =======================================================================    
    a = round(fixedMajorAxisLength/2);
    b = round(fixedMinorAxisLength/2);
    dataBox_u = zeros([length(-b:b) length(-a:a) length(heightVec)]);
    
    toc
    disp('data extraction starts ....');
    for kk=1:length(props)
        [dataBox_,detectedElemNo, totalElem, avg_nofilter] = getAvgBox(props(kk).Centroid(1),...
            props(kk).Centroid(2), a, b, props(kk).Orientation,...
            criteriaBin, uu, heightVec, Nx, Ny, Nz, false);
        packingfactor = numel(find(dataBox_u(:)))/numel(dataBox_u);
        dataBox_u = dataBox_u + dataBox_;
    end
    dataBox_u = dataBox_u./length(props);
    toc
    dataBox_u_frame = dataBox_u_frame + dataBox_u;
end
dataBox_u_frame = dataBox_u_frame./length(frameVec);

%%
pxw = size(dataBox_u,1)/2;
pyw = size(dataBox_u,2)/2;
X = linspace(-pxw,pxw, size(dataBox_u,1)); X = X .* dx ./ BLH;
Y = linspace(-pyw,pyw, size(dataBox_u,2)); Y = Y .* dy ./ BLH;
Z = linspace(1,size(dataBox_u,3),size(dataBox_u,3))-0.5; Z = Z .* dz ./ BLH;
[XX, YY, ZZ] = meshgrid(X, Y, Z);
smooth_box = smooth3(permute(dataBox_u, [2 1 3]),'gaussian',[9 9 9], 0.5);
figh = figure(); 
set(figh, 'Visible', 'on');
hold all; 
set(gcf, 'PaperUnits', 'inches', 'PaperPosition', [0 0 width height]);
thresVal = mean(smooth_box(:))-1.5*std(smooth_box(:));
p = patch(isosurface(X,Y,Z,smooth_box,thresVal));
set(p,'FaceColor','red');
set(p, 'EdgeColor','None');
% thresVal = mean(smooth_box(:))+0.50*std(smooth_box(:));
% p = patch(isosurface(X,Y,Z,smooth_box,thr esVal));
% set(p,'FaceColor','yellow');
% set(p, 'EdgeColor','None');
% thresVal = mean(smooth_box(:))-3.0*std(smooth_box(:));
% p = patch(isosurface(X,Y,Z,smooth_box,thresVal));
% set(p,'FaceColor','yellow');
% set(p, 'EdgeColor','None');
% axis tight; camlight headlight;
% lighting Gouraud; view(-160,52);
xlabel('$L_{minor}/\delta$','FontSize',14,'interpreter','latex');
ylabel('$ L_{major}/\delta$','FontSize',14,'interpreter','latex');
zlabel('$z/\delta$','FontSize',14,'interpreter','latex');
% set(gca, 'FontName','Helvetica');
axis tight; set(gca,'BoxStyle','full');
camlight('left'); view(-74,26); 
zlim([0 1]);

SPR = sprintf('%s%s','','vlsm_ek02.eps');  
%set(gcf, 'Color', 'w');
% print(formatSpec ,SPR);

%save('ek02_strucBox_u_15del.mat','dataBox_u_frame','dataBox_u', 'X','Y','Z','frameVec',...
%     'fixedMajorAxisLength','fixedMinorAxisLength','strucFound');

