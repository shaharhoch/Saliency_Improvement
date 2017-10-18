function [ out ] = s_canvasSize( orig, sz, pn, options)
%parse input
if ~exist('options','var')
    options.loc = 'center'; %expand in all directions
end
if ~exist('pn','var')
    pn = orig(1,1,:);
end

p=double(sz(1)); q=double(sz(2)); m=size(orig,1); n=size(orig,2);
out = orig;
if (p~=m || q~=n)
    pval=-300;
    pq=double(p)*double(q);

    if (strcmp(options.loc, 'center'))
    orig_pad = padarray(double(orig), double([floor((p-m)/2) floor((q-n)/2)]), pval,'post');
    orig_pad = padarray(orig_pad, double([ceil((p-m)/2) ceil((q-n)/2)]), pval,'pre');
    elseif(strcmp(options.loc, 'first'))
    orig_pad = padarray(double(orig), double([(p-m) (q-n)]), pval,'post');
    end
    i = find(orig_pad(:,:,1)==-300);
    orig_pad(i) = double(pn(1));
    if (size(orig,3)== 3)
        orig_pad(i+pq) = double(pn(2));
        orig_pad(i+2*pq) = double(pn(3));
    end
    out=cast(orig_pad,class(orig));
end
end

