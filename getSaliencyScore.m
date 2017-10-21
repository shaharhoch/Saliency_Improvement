function [ saliency_score ] = getSaliencyScore( saliency, object_mask )
saliency_score = double(0);
for i=1:numel(object_mask)
    if(object_mask(i) == 1)
        saliency_score = saliency_score + double(saliency(i));
    end 
end
saliency_score = saliency_score/sum(object_mask(:));
end

