
function [dataBox,detectedElemNo, totalElem, avg_nofilter] = getAvgBox(cx, cy, a, b,...
    orientation, bw, dataArray, heightVec, Nx, Ny, Nz, maskOptn)
    % Author : Mohammad A. Khan
    % ==================================  inputs =====================================
    % this function is to be used in conjunction with regionprops
    % cx, cy are the x and y coordinates of the centroid
    % a, b are the 1/2 lengths of the major and minor axes of the 
    % ellipse returned by regionprops, orientation is the orientation
    % returned by regionprops, which is the degree of angle between
    % ellipse's major axis and x-axis, bw is the binary array where
    % elements were detected, it can also be though of the mask array
    % against the dataArray, this dataArray is the actual scalar
    % array from which data is extracted
    % maskOptn is optional, if set to 0, data will NOT be masked from bw
    % =================================  outputs ======================================
    % avg is the average over box volume calculated from dataArray
    % totalElem is the total number of pixels in box
    % jvector_lin, jvector_data are java Vectors of linear coordinates
    % of all the interior points of the ellipsoid and jvector_data saves data
    % corresponding to mask points
    % detectedElemNo is an array that stores the number of detected
    % mask points at each horizontal level
    if nargin > 12
        error('getAvgBox:TooManyInputs', ...
            'requires at most 1 optional inputs');
    end
    switch nargin
        case 11
            maskOptn = true;
    end
    jvector_lin = java.util.ArrayList; 
    jvector_data = java.util.ArrayList; 
    % column direction aligned with major axis  
    % row direction aligned with minor axis
      
    errormsg = ['Array dimensions are not consistent'];
    if(size(bw,3) < length(heightVec)), error(errormsg); end
    if(size(dataArray,3) < length(heightVec)), error(errormsg); end
    if(prod(size(bw) == size(dataArray))== 0), error(errormsg); end
    phi = [0,pi/2, pi, 2*pi*3/4.0];
    cosphi = cos(phi);
    sinphi = sin(phi);
    % convrt to radian
    theta = pi* orientation/180;   
    R = [ cos(theta)   sin(theta)
         -sin(theta)   cos(theta)];
    % calculate the four end points of major and minor axes 
    xy = [a*cosphi; b*sinphi];
    % generate points on major axis
    xx = ceil(min(xy(1,:))):floor(max(xy(1,:)));
    minorex = floor(min(xy(2,:))):round(max(xy(2,:)));
    dataBox = zeros([length(minorex) length(xx) length(heightVec)]);
    avg_nofilter = zeros([length(heightVec) 1]);
    sumTotal = zeros([length(heightVec) 1]);
    detectedElemNo = zeros([length(heightVec) 1]);
    % for each point on major axis calculate all possible
    % y locations that fall within ellipse
    for jj = 1:length(heightVec)  
        lvlElemNo = 0;
        for i = xx
            yy = floor(min(xy(2,:))):round(max(xy(2,:)));
            xcopy = ones(1,length(yy)).*i;
            % rotate the ellipse that was generated around (0,0) point 
            % to match the oreintation given and translate the ellipse  to the
            % centroid location
            rtxy = R*[xcopy;yy];
            xcor = (rtxy(1,:)+cx);
            ycor = (rtxy(2,:)+cy);
            % establish bound check
            xout1 = find(xcor < 1);
            xout2 = find(xcor > size(bw,2));
            xcor([xout1, xout2]) = [];
            ycor([xout1, xout2]) = [];
            yout1 = find(ycor < 1);
            yout2 = find(ycor > size(bw,1));
            xcor([yout1, yout2]) = [];
            ycor([yout1, yout2]) = [];
            xcor = round(xcor);
            ycor = round(ycor);
            % delete all the points that fall on the background
            % of the binary mask image at different heights
            if ~isempty(xcor)
                if(maskOptn)
                    for ii = length(xcor):-1:1
                        if(~bw(ycor(ii),xcor(ii),heightVec(jj)))
                            xcor(ii) = [];
                            ycor(ii) = [];
                        end                
                    end
                end
            end        
            if ~isempty(xcor)  
                for ii = 1:length(ycor)
                    % xcor has been calculated so far in column direction
                    % ycor has been calcuated in row direction
                    % remeber graph x direction is actually column direction
                    % and graph y direction is row direction 
                    sumTotal(jj) = sumTotal(jj) + dataArray(ycor(ii),xcor(ii),heightVec(jj));
                    linIndx = sub2ind([Nx, Ny, Nz],ycor(ii), xcor(ii), heightVec(jj)); 
                    jvector_lin.add(linIndx);
                    jvector_data.add(dataArray(linIndx)); 
                    lvlElemNo = lvlElemNo + 1;
                end                
            end            
        end
        detectedElemNo(jj) = lvlElemNo; 
%         disp(['lvlElemNo :' num2str(lvlElemNo)]);  
%         disp(['jvector_lin.size :' num2str(jvector_lin.size)]);
%         disp(['jvector_data.size :' num2str(jvector_data.size)]);
%         disp(['sum(detectedElemNo(:)) :' num2str(sum(detectedElemNo(:)))]);
        avg_nofilter(jj) = mean(mean(dataArray(:,:,heightVec(jj))));
    end
    vector_lin = zeros([jvector_lin.size 1]);
    dataVec = zeros([jvector_lin.size 1]);
    for tt = 0:jvector_lin.size()-1
        vector_lin(tt+1) =  jvector_lin.get(tt); % linear  index
        dataVec(tt+1) = jvector_data.get(tt);
    end
    clear jvector_lin  jvector_data;
    % undo the previous rotation
    theta = -pi* orientation/180;   
    R = [ cos(theta)   sin(theta)
         -sin(theta)   cos(theta)]; 
    vector_hc = zeros([length(vector_lin) 1]);
    for kk = 1:length(vector_lin)
        [vc, hc, pc] = ind2sub([Nx, Ny, Nz], vector_lin(kk));
        vc = vc - cy; hc = hc - cx;
        res = R * [hc; vc];
        hc = res(1); vc = res(2);        
        hc = round(hc + length(xx)/2 + 1);
        vc = floor(vc + length(minorex)/2 + 1);       
        pc = find(heightVec == pc);
        dataBox(vc, hc, pc) = dataVec(kk);
        vector_hc(kk) = hc;
        
    end
    avg = sumTotal./detectedElemNo;
    totalElem = round(a*2*b*2*length(heightVec));
    disp(['=======================================================================']);
end
