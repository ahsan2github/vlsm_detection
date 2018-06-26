
function [avg, detectedElemNo, totalElem, vector_x, vector_y] = avgInEllipse(cx, cy, a, b,...
    orientation, bw, dataArray)
    % Author : Mohammad A. Khan
    % =======================  inputs ==========================
    % this function is to be used in conjunction with regionprops
    % cx, cy are the x and y coordinated of the centroid
    % a, b are the 1/2 lengths of the major and minor axes of the 
    % ellipse returned by regionprops, orientation is the orientation
    % returned by regionprops, which is the degree of angle between
    % ellipse's major axis and x-axis, bw is the binary array where
    % elements were detected, it can also be though of the mask array
    % against the dataArray, this dataArray is the actual scalar
    % array on which the averaging over ellipsoidal area is expected
    % ====================  outputs =============================
    % avg is the average over ellipse area calculated from dataArray
    % totalElem is the total number of pixels in ellipse
    % jvector_x, jvector_y are java Vectors of x,y coordinates
    % of all the interior points of the ellipse
    
    jvector_x = java.util.Vector;
    jvector_y = java.util.Vector;    
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
    totalElem = 0;
    sumTotal = 0;
    % for each point on major axis calculate all possible
    % y locations that fall within ellipse
    for i = xx
        proot = [-sqrt((1-i*i/a/a)*b*b), sqrt((1-i*i/a/a)*b*b)];
        yy = real((min(proot))):real((max(proot)));
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
        % of the binary mask image
        if ~isempty(xcor)
            for ii = length(xcor):-1:1
                if(~bw(ycor(ii),xcor(ii)))
                    xcor(ii) = [];
                    ycor(ii) = [];
                end                
            end
        end        
        if ~isempty(xcor)            
            totalElem = totalElem + numel(xcor);
            for ii = 1:length(xcor)
                % xcor has been calculated so far in column direction
                % ycor has been calcuated in row direction
                % remeber graph x direction is actually column direction
                % and graph y direction is row direction 
                sumTotal = sumTotal + dataArray(ycor(ii),xcor(ii));
                jvector_x.addElement(xcor(ii));
                jvector_y.addElement(ycor(ii));
            end

        end
    end
    vector_x = zeros([jvector_x.size 1]);
    vector_y = zeros([jvector_x.size 1]);
    for tt = 0:jvector_x.size()-1
        vector_x(tt+1) =  jvector_y.elementAt(tt);
        vector_y(tt+1) =  jvector_x.elementAt(tt);
    end
    avg = sumTotal/totalElem;
    detectedElemNo = jvector_x.size;
    totalElem = round(pi*a*b);
end
