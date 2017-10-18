function [ ex_rgb,ex_lab] = s_imextrapolate( imrgb, imlab, options )
%S_IMEXTRAPOLATE perform extrapolation according to Criminisi inpainting
%algorithm

dbglvl =0;

%define source and target regions
mask = s_im2bin(imlab, 'lab'); % segment the image to background and fragment
source = s_scale_mask(mask,-options.contDistFromBoundary); %shrink fragment mask
target = logical(s_scale_mask(mask,options.extrapPixelsSize) - source); %expend mask for target region
%------------------------------

psz = options.psz;
% error check
if ~ismatrix(target) || ~islogical(target) ||...
    ~ismatrix(source) || ~islogical(source)
    error('Invalid mask'); 
end
if mod(psz,2)==0 
    error('Patch size psz must be odd.'); 
end

ex_rgb=imrgb;
ex_lab=imlab;

sz = size(source);
if ~isequal(size(imrgb),size(imlab))||...
    ~isequal(size(imrgb),[sz,3]) ||...
    ~isequal(sz,size(target))
    error('wrong input sizes');
end


[G,Gx,Gy] = h_grad(imlab);
G(target)=0; Gx(target)=0; Gy(target)=0;

G=G./max(G(:));
realIntensities = imgradient(Gx,Gy);
nzInten_i= find(realIntensities ~=0); 
Gx(nzInten_i)=Gx(nzInten_i)./realIntensities(nzInten_i);
Gy(nzInten_i)=Gy(nzInten_i)./realIntensities(nzInten_i);

[Gx,Gy]=deal(-Gy,Gx); % Rotate gradient 90 degrees

[imlab_rots,imrgb_rots, src_rots,G_rots,Gx_rots,Gy_rots] = h_rotate(options.rotations,imlab,imrgb,source,G,Gx,Gy);
srcp_rots=zeros([2,numel(source),size(src_rots,3)],'uint32');
src_pch_rots_starts_rc=zeros([2,numel(source),size(src_rots,3)],'uint32');
for ri=1:size(src_rots,3)
    srcp_rots(:,:,ri)=h_src2patches(src_rots(:,:,ri), psz);
    valid_srcp_i = srcp_rots(1,:,ri)~=0;
    [src_pch_rots_starts_rc(1,valid_srcp_i,ri),src_pch_rots_starts_rc(2,valid_srcp_i,ri)]=ind2sub(sz,srcp_rots(1,valid_srcp_i,ri));
end

%init
C=double(source);
prevdR=[];
Hp=[];
N=[];
% Seed 'rand' for reproducible results (good for testing)
rand('state',0);
iter=1;
% Loop until entire fill region has been covered
while any(target(:))
    % Find contour & normalized gradients of fill region
    dinvsrc = double(~source);
    dR = find(conv2(dinvsrc,[1,1,1;1,-8,1;1,1,1],'same')>0);
    if (~strcmp(options.Dmode,'g'))
        [Nx,Ny] = imgradientxy(dinvsrc);
        N = [Nx(dR(:)) Ny(dR(:))];
        N = normr(N);
        N(~isfinite(N))=0; % handle NaN and Inf
    end
    % Compute confidences along the fill front
    [C,Hp] = h_computeC(dR,prevdR,C,Hp,sz,psz,target);
    % Compute patch priorities = confidence term * data term
    D = h_computeD(G,Gx,Gy,dR,N,options.Dmode)+eps;
    priorities = C(dR).* D;
    
    % Find patch with maximum priority, Hp
    [~,mx_prio_i] = max(priorities);
    tar_p = dR(mx_prio_i);
    tar_pch_ind = Hp(:,:,mx_prio_i);
    toFill = target(tar_pch_ind);
    if sum(toFill(:))==0
        error('chose a mask with nothing to fill');
    end
    % Find exemplar that minimizes error, Hq
    [src_pch_ind,chosen_rot_i] = h_bestexemplar(imlab_rots,tar_pch_ind,toFill, source, srcp_rots, G_rots);
      
    % Copy image data from src patch to target patch
    chosen_rot_offset = (chosen_rot_i-1)*numel(imlab);
    tar_pch_ind_toFill = tar_pch_ind(toFill);
    src_pch_ind_toFill = src_pch_ind(toFill);
    tar_3di=s_ind2to3(tar_pch_ind_toFill,sz);
%     tar_3di=tar_3di(:);
    src_3di=s_ind2to3(src_pch_ind_toFill,sz) + chosen_rot_offset;
%     src_3di=src_3di(:);
    
    ex_lab(tar_3di) = imlab_rots(src_3di);
    prev_ex_rgb = ex_rgb;
    ex_rgb(tar_3di) = imrgb_rots(src_3di);
    imlab_rots(:,:,:,1) = ex_lab;
    
    target(tar_pch_ind_toFill) = false;
    source(tar_pch_ind_toFill) = true;
    
    % Propagate confidence & isophote values
    C(tar_pch_ind_toFill) = C(tar_p);
    if (options.isCalcGradEveryIter)
        [G,Gx,Gy] = h_grad(ex_lab, options.Gmode);
        G(target)=0; Gx(target)=0; Gy(target)=0;

        G=G./max(G(:));
        [Gx,Gy]=deal(-Gy,Gx); % Rotate gradient 90 degrees
    else
        chosen_rot_offset = (chosen_rot_i-1)*numel(G);
        src_i=src_pch_ind_toFill + chosen_rot_offset;
    
        Gx(tar_pch_ind_toFill) = Gx_rots(src_i);
        Gy(tar_pch_ind_toFill) = Gy_rots(src_i);
        G(tar_pch_ind_toFill) = G_rots(src_i);
    end
    
    d_showex(ex_rgb,src_pch_ind,tar_pch_ind,sz,dbglvl);%---------->Debug
    
    prevdR=dR;
    iter=iter+1;
end

end


%------------------------------
%       Help methods
%------------------------------
function [varargout] = h_rotate(angles, varargin)
    ROT_METHOD = 'nearest';
    ROT_SIZE = 'crop';    

    n_ang = length(angles);
    n_varin = length(varargin);
    if nargout~=n_varin
        error('wrong input and output');
    end
    if n_ang == 0
        varargout = varargin;
    else
        for k = 1:n_varin
           sz=size(varargin{k});
           nd = ndims(varargin{k});
           varargout{k}=zeros([sz,n_ang]);
           for ang_i = 1:n_ang
               tmp = imrotate(varargin{k},angles(ang_i),ROT_METHOD,ROT_SIZE);
               if nd==2
                   varargout{k}(:,:,ang_i)=tmp;
               elseif nd==3
                   varargout{k}(:,:,:,ang_i)=tmp;
               end
           end
        end
    end
end

function [g,gx,gy] = h_grad(imlab)
    [gx,gy] = imgradientxy(imlab(:,:,1));
    g = imgradient(gx, gy);
end

function [srcp] = h_src2patches(src, psz)
DEBUG_VERB = 0;

nel = numel(src);
srcp = GetValidPatches(size(src,1), size(src,2), psz, psz, nel, logical(src), uint32(1:nel));
srcp = uint32(srcp);
if (DEBUG_VERB >=1)
    dsrc = double(src);
    st = srcp(1,:);
    dsrc(st(st ~= 0)) = 0.5;
    figure(500); imshow(dsrc); title('src region with valid patches starts');
end
end

function [Hp,rows,cols] = h_getpatch(imsz,psz,pi)
    w=(psz-1)/2;
    [r,c]=ind2sub(imsz,pi);
    rows = (max(r-w,1):min(r+w,imsz(1)))';
    cols = (max(c-w,1):min(c+w,imsz(2)));
    mrows=repmat(rows,1,size(cols,2));
    mcols=repmat(cols,size(rows,1),1);
    Hp=sub2ind(imsz,mrows,mcols);
end

function [C,new_patches] = h_computeC(dR,prevdR,C,Hp,sz,psz,target)
    [~,dR_loc_in_prevdR]=ismember(dR,prevdR);
    new_patches=zeros(psz,psz,length(dR));
    for ki=1:length(dR)
        k=dR(ki);
        if dR_loc_in_prevdR(ki)~=0 % already calc in prev iteration
            kpatch_ind = Hp(:,:,dR_loc_in_prevdR(ki));
        else
            [kpatch_ind,~,~] = h_getpatch(sz,psz,k);
        end
        kpch_tar_mask = target(kpatch_ind);
        if any(kpch_tar_mask(:))
            src_ind = kpatch_ind(~kpch_tar_mask); 
            C(k) = sum(C(src_ind))/numel(kpatch_ind);
        else
            C(k)=0;
        end
        new_patches(:,:,ki)=kpatch_ind;
    end
end

function [D] = h_computeD(G,Gx,Gy,dR,N,mode)
    if strcmp(mode,'g')
        D=G(dR) + log(1+G(dR));
    end
    if strcmp(mode,'c')
        D=abs(Gx(dR).*N(:,1)+Gy(dR).*N(:,2));
    end
    if strcmp(mode,'all')
        D1=G(dR) + log(1+G(dR));
        D2 = abs(Gx(dR).*N(:,1)+Gy(dR).*N(:,2));    
        D=D1.*D2;
    end
    if any(D<0)
        error('Data term should not be negative');
    end
end

function [Hq, chosen_rot_i] = h_bestexemplar(imlab_rots,tar_pch_ind,toFill,source,srcp_rots, G_rots)
    smoothFactor=0.5;    
    psz=size(tar_pch_ind,1); mm=size(imlab_rots,1); nn=size(imlab_rots,2);
    
    n_imgs = size(imlab_rots,4);
    [errors, intensities] = bestexemplarhelper(mm,nn,psz,psz,n_imgs,...
                double(imlab_rots),uint32(tar_pch_ind), source,srcp_rots, G_rots);
    nerrors=sqrt(errors);
%     nerrors=nerrors./max(nerrors(:));
%     nintensities=intensities./numel(find(toFill));
    %errors : (mm*nn) X n_imgs, contains all the errors for every
    %patch in every image
    mixerror = (1 - smoothFactor)*nerrors + smoothFactor * intensities;
%     mixerror = nerrors;

    [bestErrorsVal,bestErrorsIdx]=min(mixerror); % output the val and idx of the minimal error for each image
    [~,chosen_rot_i]=min(bestErrorsVal); % find the image with the minimal error
    q = (srcp_rots(1, bestErrorsIdx(chosen_rot_i), chosen_rot_i )...
        + srcp_rots(2, bestErrorsIdx(chosen_rot_i ), chosen_rot_i )) / 2 ;
    [Hq,~, ~] = h_getpatch([mm nn],psz,q);
end

%------------------------------
%       Debug methods
%------------------------------
function d_showex(ex_rgb,src_pch_ind,tar_pch_ind,sz,dbglvl)
    if (dbglvl>=1)
        src_3di = s_ind2to3(src_pch_ind(:),sz);
        tar_3di = s_ind2to3(tar_pch_ind(:),sz);
        ex_rgb_wpatch=ex_rgb;
        ex_rgb_wpatch(src_3di) = 80;
        ex_rgb_wpatch(tar_3di) = 170;
        figure(306),imshow(ex_rgb_wpatch);
    end
end