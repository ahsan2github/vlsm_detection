function [varargout] = removeOverlappedBoxes(props,a,b,binArrayMask, dataMat)       
% ---------------- delete overlapped boxes ---------------------------------------
    strucFound = 0;    
    horzcc = zeros([length(props) 1]);
    vertcc = zeros([length(props) 1]);
    horzDim = zeros([length(props) 1]);
    vertDim = zeros([length(props) 1]);
    horzBLC = zeros([length(props) 1]);
    vertBLC = zeros([length(props) 1]);
    stressMax = zeros([length(props) 1]);
    for ss = 1:length(props)
        [horzBLC(ss),vertBLC(ss),horzDim(ss),vertDim(ss),horzcc(ss),vertcc(ss)] = ...
            getBBox(props(ss).Centroid(1), props(ss).Centroid(2),...
            a, b, props(ss).Orientation,size(binArrayMask,1), size(binArrayMask,2)); 
        [avgStress, ~, ~, ~, ~] = avgInEllipse1(props(ss).Centroid(1),...
             props(ss).Centroid(2), a, b,...
             props(ss).Orientation, binArrayMask, dataMat);
        stressMax(ss) = avgStress;       
    end      
    [res, si] = sort(stressMax,'descend');  
    stressMax = res;
    props = props(si);
    horzBLC = horzBLC(si);
    vertBLC = vertBLC(si);
    horzDim = horzDim(si);
    vertDim = vertDim(si);
    horzcc = horzcc(si);
    vertcc = vertcc(si);
    tracker = 1;
    while(~isempty(props))
        if(tracker > length(props)), break; end;
        CollisionCondition = zeros([length(props) 1]);
        tmp_elem = props(tracker);
        tmp_horzDim = horzDim(tracker);
        tmp_vertDim = vertDim(tracker);
        tmp_horzcc = horzcc(tracker);
        tmp_vertcc = vertcc(tracker);
        tmp_horzBLC = horzBLC(tracker);
        tmp_vertBLC = vertBLC(tracker);   
        tmp_stressMax = stressMax(tracker);    
        for counter = 1:length(props)
            CollisionCondition(counter) = boxCollide(horzcc(tracker), vertcc(tracker),...
                   horzDim(tracker), vertDim(tracker),horzcc(counter), vertcc(counter),...
                   horzDim(counter), vertDim(counter));
        end
        % find which elements are connectected to elem(tracker)
        conditionIsTrue = find(CollisionCondition);
        if(~isempty(conditionIsTrue))
            % remove the overlapped boxes
            conditionIsTrue = conditionIsTrue(2:end);    
            for test2 = length(conditionIsTrue):-1:1
               props(conditionIsTrue(test2)) = [];  
               stressMax(conditionIsTrue(test2)) = [];
               horzDim(conditionIsTrue(test2)) = [];
               vertDim(conditionIsTrue(test2)) = [];
               horzcc(conditionIsTrue(test2)) = [];
               vertcc(conditionIsTrue(test2)) = [];
               horzBLC(conditionIsTrue(test2)) = [];
               vertBLC(conditionIsTrue(test2)) = [];
            end            
        end
        props(tracker) = tmp_elem;
        horzDim(tracker) = tmp_horzDim;
        vertDim(tracker) = tmp_vertDim;
        horzcc(tracker) = tmp_horzcc;
        vertcc(tracker) = tmp_vertcc;
        horzBLC(tracker) = tmp_horzBLC;
        vertBLC(tracker) = tmp_vertBLC;
        stressMax(tracker) = tmp_stressMax;
        tracker = tracker + 1;
    end 
    strucFound  = strucFound + tracker-1;
    for counter = length(props):-1:1
        if(isempty(props(counter).Centroid)), props(counter)=[]; end
        if(isempty(props(counter).MajorAxisLength)), props(counter)=[]; end
    end
    varargout{1} = props;
    varargout{2} = horzcc;
    varargout{3} = vertcc;
    varargout{4} = horzDim;
    varargout{5} = vertDim;
    varargout{6} = horzBLC;
    varargout{7} = vertBLC;
    disp(['Overlapped Structure isolation done, found structures: ' num2str(strucFound )]);
end
