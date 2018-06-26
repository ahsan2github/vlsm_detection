%%% ----------------- spectral filter and (optional) x and y direction 
%%% derivative calculation.  Note: f should be the ratio Lx/d 
%%% (where d is the filter scale) f should correspond to a vector of 
%%% integers the output will have the same size as x,f.  The first returned 
%%% variable is an array of filtered velocities, the 2nd and 3rd are the
%%% x and y derivatives, respectively.  The 4th and 5th, if requested, are
%%% the unfiltered x and y derivatives, respectively.  If only 1 output is 
%%% requested, only the filtered values are returned (and not the filtered
%%% derivs), if only 3 are requested, then only the filtered values and the
%%% and the filtered derivs are returned. Note setting f=size(x)/2 only
%%% will return unfiltered derivs
function [varargout] = cuttoff_2D_multiF(x,f)

nn=size(x);
x_h = fft2(x);
f=round(f);

u=zeros([size(x),length(f)]);
if(nargout > 1)
    ki_ii=zeros(nn);
    ki_jj=zeros(nn);
end

if(nargout > 1)
    
    dudxf=zeros(size(u));
    dudyf=zeros(size(u));
    
    for j=1:nn(2)
        jj=j-1;
        if(jj > nn(2)/2); jj=jj-nn(2); end
        
        for i=1:nn(1)
            ii=i-1;
            if(ii > nn(1)/2); ii=ii-nn(1); end
            
            %multiply by ki (wavenumber times imaginary number)
            if(jj == nn(2)/2 || ii == nn(1)/2)
                ki_ii(i,j) = 0;
                ki_jj(i,j) = 0;
            else
                ki_ii(i,j) = ii*sqrt(-1);
                ki_jj(i,j) = jj*sqrt(-1);
            end
        end
    end
    
end

for i=1:length(f)
    x_hat = x_h;
    x_hat(f(i)+1:nn(1)-(f(i)-1),:) = 0; %set k_x above cutt off zero
    x_hat(:,f(i)+1:nn(1)-(f(i)-1)) = 0; %set k_y above cutt off zero

    u(:,:,i) = ifft2(x_hat,'symmetric');
    if(nargout > 1)
        dudxf(:,:,i) = ifft2(x_hat.*ki_ii,'symmetric');
        dudyf(:,:,i) = ifft2(x_hat.*ki_jj,'symmetric');
    end
end

varargout(1)={u};
if(nargout > 1)
    varargout(2)={dudxf};
    varargout(3)={dudyf};
    if(nargout > 3)
        varargout(4) = {ifft2(x_h.*ki_ii,'symmetric')};
        varargout(5) = {ifft2(x_h.*ki_jj,'symmetric')};
    end
end
