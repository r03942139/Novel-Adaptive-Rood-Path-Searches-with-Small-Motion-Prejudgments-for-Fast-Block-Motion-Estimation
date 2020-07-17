% Computes the Mean Absolute Difference (MAD) for the given two blocks
% Input
%       currentBlk : The block for which we are finding the MAD
%       refBlk : the block w.r.t. which the MAD is being computed
%       n : the side of the two square blocks
%
% Output
%       cost : The MAD for the two blocks
%
% Written by Aroh Barjatya

%% It has been revised to the Sum absolute difference (SAD) in 2016.10.01 by Hung-Yi Chen

function cost = costFuncMAD(currentBlk,refBlk, n, edge_wei)
% SAD
% SAD        = sum(sum(abs(Block - Block_ref).*edge_wei));%/(My*Nx);
if nargin == 3
    err = sum(sum(abs(currentBlk- refBlk))); 
elseif nargin == 4
    err = sum(sum(abs(currentBlk- refBlk).*edge_wei));
end
% err = sum(sum(abs(currentBlk- refBlk)));
cost = err;
% cost = sum(sum(xor(currentBlk,refBlk)))/ (n*n);
% cost = err ;
