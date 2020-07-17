% Computes motion vectors using Adaptive Rood Pattern Search method
%
% Based on the paper by Hadi Amirpour1 ¡P Amir Mousavinia1
% SIViP (2016) 10:1393¡V1400
%
% Input
%   imgP : The image for which we want to find motion vectors
%   imgI : The reference image
%   mbSize : Size of the macroblock
%   p : Search parameter  (read literature to find what this means)
%
% Ouput
%   motionVect : the motion vectors for each integral macroblock in imgP
%   ARPScomputations: The average number of points searched for a macroblock
%
% Written by Aroh Barjatya


function [motionVect, ARPScomputations] = motionEstDSPMV(imgP, imgI, mbSize, p, MVpre)

[row,col] = size(imgI);

vectors = zeros(2,row*col/mbSize^2);
costs = ones(1,9) * 65537;

% The index points for Small Diamond search pattern
SDSP(1,:) = [ 0 -1];
SDSP(2,:) = [-1  0];
SDSP(3,:) = [ 0  0];
SDSP(4,:) = [ 1  0];
SDSP(5,:) = [ 0  1];

SDSP(6,:) = [ 1  1];
SDSP(7,:) = [ 1  -1];
SDSP(8,:) = [ -1  -1];
SDSP(9,:) = [ -1  1];

% We will be storing the positions of points where the checking has been
% already done in a matrix that is initialised to zero. As one point is
% checked, we set the corresponding element in the matrix to one. 

checkMatrix = zeros(2*p+1,2*p+1);

computations = 0;

% we start off from the top left of the image
% we will walk in steps of mbSize
% mbCount will keep track of how many blocks we have evaluated

mbCount = 1;
UP = col/mbSize;
for i = 1 : mbSize : row-mbSize+1
    for j = 1 : mbSize : col-mbSize+1
        
        % the Adapive Rood Pattern search starts
        % we are scanning in raster order
        
        x = j;
        y = i;
        
        ini = 1;
        LDSP = [];
            % have center
            LDSP(ini,:) = [0,0];
            ini = ini+1;
            if j > 1 % have left (1)
                LDSP(ini,:) = [vectors(2,mbCount-1), vectors(1,mbCount-1)];
                ini = ini+1;
            end
            if i > 1 % have up (2)
                LDSP(ini,:) = [vectors(2,mbCount- UP), vectors(1,mbCount- UP)];
                ini = ini+1;
                if j > 1 % have upper-left (3)
                    LDSP(ini,:) = [vectors(2,mbCount- UP-1), vectors(1,mbCount- UP-1)];
                    ini = ini+1;
                end
                if j < col % have upper-right (5) 
                    LDSP(ini,:) = [vectors(2,mbCount- UP+1), vectors(1,mbCount- UP+1)];
                    ini = ini+1;
                end  
            end
            % Temporal-Previous
            LDSP(ini,:) = [MVpre(2,mbCount), MVpre(1,mbCount)];
            maxIndex = ini;
            
        
        % do the DPS-MV
        for k = 1:maxIndex
            refBlkVer = y + LDSP(k,2);   % row/Vert co-ordinate for ref block
            refBlkHor = x + LDSP(k,1);   % col/Horizontal co-ordinate
            if ( refBlkVer < 1 || refBlkVer+mbSize-1 > row ...
                 || refBlkHor < 1 || refBlkHor+mbSize-1 > col)
                continue; % outside image boundary
            end
            if checkMatrix(LDSP(k,2) + p+1, LDSP(k,1) + p+1) == 1;
                continue
            end
            costs(k) = costFuncMAD(imgP(i:i+mbSize-1,j:j+mbSize-1), ...
                  imgI(refBlkVer:refBlkVer+mbSize-1, refBlkHor:refBlkHor+mbSize-1), mbSize);
            computations =  computations + 1;
            checkMatrix(LDSP(k,2) + p+1, LDSP(k,1) + p+1) = 1;
            
        end
        
        [cost, point] = min(costs);
        
%% The 'doneFlag' is set to '1' when the minimum is at the center of the diamond           
        x = x + LDSP(point, 1);
        y = y + LDSP(point, 2);
        costs = ones(1,9) * 65537;
        costs(3) = cost;
        
        if cost < 512
            vectors(1,mbCount) = y - i;    % row co-ordinate for the vector
            vectors(2,mbCount) = x - j;    % col co-ordinate for the vector            
            mbCount = mbCount +1;
            costs = ones(1,9) * 65537;
            checkMatrix = zeros(2*p+1,2*p+1);
            continue
        else % Next further step
        
            for k = 1:size(SDSP,1)
                refBlkVer = y + SDSP(k,2);   % row/Vert co-ordinate for ref block
                refBlkHor = x + SDSP(k,1);   % col/Horizontal co-ordinate
                if ( refBlkVer < 1 || refBlkVer+mbSize-1 > row ...
                      || refBlkHor < 1 || refBlkHor+mbSize-1 > col)
                      continue;
                end

                if (k == 3)
                    continue
                elseif (refBlkHor < j-p || refBlkHor > j+p || refBlkVer < i-p ...
                            || refBlkVer > i+p)
                        continue;
                elseif (checkMatrix(y-i+SDSP(k,2)+p+1 , x-j+SDSP(k,1)+p+1) == 1)
                    continue
                end
            
                costs(k) = costFuncMAD(imgP(i:i+mbSize-1,j:j+mbSize-1), ...
                             imgI(refBlkVer:refBlkVer+mbSize-1, ...
                                 refBlkHor:refBlkHor+mbSize-1), mbSize);
                checkMatrix(y-i+SDSP(k,2)+p+1, x-j+SDSP(k,1)+p+1) = 1;
                computations =  computations + 1;
                
  
            end
            
            [cost, point] = min(costs);
           
            x = x + SDSP(point, 1);
            y = y + SDSP(point, 2);
            costs = ones(1,9) * 65537;
%             end  % while loop ends here

        vectors(1,mbCount) = y - i;    % row co-ordinate for the vector
        vectors(2,mbCount) = x - j;    % col co-ordinate for the vector            
        mbCount = mbCount + 1;
        costs = ones(1,9) * 65537;
        checkMatrix = zeros(2*p+1,2*p+1);
        end
    end
end
    
motionVect = vectors;
ARPScomputations = computations/(mbCount-1) ; 