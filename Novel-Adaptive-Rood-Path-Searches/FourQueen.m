% 2004 A Hierarchical N-Queen Decimation Lattice 
% Input: N*N data
function Out4Queen = FourQueen(Input)

[s1,s2] = size(Input); % 16x16 expected
Out4Queen = zeros(8); % 8x8 pattern

for i = 1:4
    if i == 1
        I = Input(1:8,1:8); % Up-left
    elseif i ==2
        I = Input(1:8,9:16); % Up-right
    elseif i ==3
        I = Input(9:16,1:8); % Down-left
    else % i == 4
        I = Input(9:16,9:16); % Down-left
    end
    
Out4Queen(2*i-1:2*i,1:8) = [ I(3,1),I(1,2),I(4,3),I(2,4),I(3,5),I(1,6),I(4,7),I(2,8);
                             I(7,1),I(5,2),I(8,3),I(6,4),I(7,5),I(5,6),I(8,7),I(6,8)];
end





    



end