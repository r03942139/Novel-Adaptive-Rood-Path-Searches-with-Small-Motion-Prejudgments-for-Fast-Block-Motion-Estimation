function [ssim_index,ssim_map]=GBIQA(img1,img2)
%     img1=rgb2ycbcr(img1,'forward');
%     img2=rgb2ycbcr(img2,'forward');
    
    Bsize=8;
    C4=1e-5;
    Kp=200;
    
    [M,N,c]=size(img1);
    g=zeros(size(img1));
    
    %% Gradient similarity
    for k=1:c
        gx=cal_g(img1(:,:,k));
        gy=cal_g(img2(:,:,k));
        R=abs(gx-gy)./max(gx,gy);
        K=Kp./max(gx,gy);
        
        g(:,:,k)=(2*(1-R)+K)./(1+(1-R).^2+K);
    end
    g=max(0,min(1,g));
    
    %% Luminance similarity
    L=255; % luminance levels
    
    e=1-((img1-img2)/L).^2;
    
    %% Overall weighting
    p=0.1; % positive weighting parameter
    W=p.*g;
    
    ssim_map=(1-W).*g+W.*e;
    ssim_index=mean(ssim_map(:));
end

function g=cal_g(img)
    M1=[0 0 0 0 0; 1 3 8 3 1; 0 0 0 0 0; -1 -3 -8 -3 -1; 0 0 0 0 0];
    M2=[0 0 1 0 0; 0 8 3 0 0; 1 3 0 -3 -1; 0 0 -3 -8 0; 0 0 -1 0 0];
    M3=fliplr(M2);
    M4=M1.';
    window=ones(size(M1))/numel(M1);
    
    g=zeros(size(img));
    g(:,:,1)=filter2(window,filter2(M1,img,'same'),'same');
    g(:,:,2)=filter2(window,filter2(M2,img,'same'),'same');
    g(:,:,3)=filter2(window,filter2(M3,img,'same'),'same');
    g(:,:,4)=filter2(window,filter2(M4,img,'same'),'same');
    
    g=max(abs(g),[],3);
    
end