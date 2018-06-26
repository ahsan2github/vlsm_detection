function [out]=removeLevelMean(inp)
	for i = 1:size(inp,3)
		out(:,:,i)=inp(:,:,i)-mean(mean(inp(:,:,i)));
	end
end
