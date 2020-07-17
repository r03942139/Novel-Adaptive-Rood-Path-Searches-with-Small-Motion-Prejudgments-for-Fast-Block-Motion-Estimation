function [SSIM_out, GBIQA_out] = testSSIMGBIQA(Seq,Seq_r)
    [M,N,C,K] = size(Seq);
    
    mssim = zeros(3,K);
    for j = 2:K
        [mssim(1,j), ssim_map] = find_ssim(Seq_r(:,:,1,j),Seq(:,:,1,j)); % Y
        [mssim(2,j), ssim_map] = find_ssim(Seq_r(:,:,2,j),Seq(:,:,2,j));% Cb
        [mssim(3,j), ssim_map] = find_ssim(Seq_r(:,:,3,j),Seq(:,:,3,j));% Cr
    end
    SSIM_out = (mean(mssim(1,:)) + mean(mssim(2,:)) + mean(mssim(3,:)) )/3;
    
    mGBIQA = zeros(3,K);
    for j = 2:K
        [mGBIQA(1,j), ssim_map] = GBIQA(Seq_r(:,:,1,j),Seq(:,:,1,j)); % Y
        [mGBIQA(2,j), ssim_map] = GBIQA(Seq_r(:,:,2,j),Seq(:,:,2,j));% Cb
        [mGBIQA(3,j), ssim_map] = GBIQA(Seq_r(:,:,3,j),Seq(:,:,3,j));% Cr
    end
    GBIQA_out = (mean(mGBIQA(1,:)) + mean(mGBIQA(2,:))+ mean(mGBIQA(3,:)) )/3;


end