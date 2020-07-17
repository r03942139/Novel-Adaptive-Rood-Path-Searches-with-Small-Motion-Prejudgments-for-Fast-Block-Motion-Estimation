% Computes motion vectors using Proposed method
% By Hung-Yi Chen 
% 
% Input
%   imgP : The image for which we want to find motion vectors
%   imgI : The reference image
%   mbSize : Size of the macroblock
%   p : Search parameter  (read literature to find what this means)
%
% Ouput
%   motionVect : the motion vectors for each integral macroblock in imgP
% 
%   ARPScomputations: The average number of points searched for a macroblock


function [motionVect, Chencomputations] = MEChenARPS2(imgP, imgI, mbSize, p, MVpre)

[row,col] = size(imgI);

vectors = zeros(2,row*col/mbSize^2);

MVx = zeros(row/mbSize,col/mbSize); MVy = MVx;
% costs = ones(1,7) * 65537; % ARPS
costs = ones(1,11) * 65537; %A RPS
% costs = ones(1,7) * 65537; % Chen's ARPS

% The index points for Small Diamond search pattern
SDSP(1,:) = [ 0 -1];
SDSP(2,:) = [-1  0];
SDSP(3,:) = [ 0  0];
SDSP(4,:) = [ 1  0];
SDSP(5,:) = [ 0  1];

% The index points for Chen-based Diamond search pattern
ChenDSP(1,:) = [ 0 -1];
ChenDSP(2,:) = [-1  0];
ChenDSP(3,:) = [ 0  0];
ChenDSP(4,:) = [ 1  0];
ChenDSP(5,:) = [ 0  1];
            
ChenDSP(6,:) = [ 1  1];
ChenDSP(7,:) = [ -1  1];
ChenDSP(8,:) = [ 1  -1];
ChenDSP(9,:) = [ -1 -1];

ChenDSP(10,:) = [-2  0];
ChenDSP(11,:) = [ 2  0];
% ChenDSP(10,:) = [ 0 -2];
% ChenDSP(11,:) = [ 0  2];

% We will be storing the positions of points where the checking has been
% already done in a matrix that is initialised to zero. As one point is
% checked, we set the corresponding element in the matrix to one. 

checkMatrix = zeros(2*p+1,2*p+1); % [2*Search_Range+1,2*Search_Range+1]
computations = 0;

% we start off from the top left of the image
% we will walk in steps of mbSize
% mbCount will keep track of how many blocks we have evaluated
mbCount = 1;
for i = 1 : mbSize : row-mbSize+1
    for j = 1 : mbSize : col-mbSize+1
        
        % the Adapive Rood Pattern search starts
        % we are scanning in raster order
        x = j;
        y = i;
        
        costs(3) = costFuncMAD(imgP(i:i+mbSize-1,j:j+mbSize-1), ...
                                    imgI(i:i+mbSize-1,j:j+mbSize-1),mbSize);
        
%% =================== Stage 1. ===================    
        if costs(3) <= 512 % Zero Motion Prejudgment
            vectors(1,mbCount) = 0;    % row co-ordinate for the vector
            vectors(2,mbCount) = 0;    % col co-ordinate for the vector   
            MVy(round(i/mbSize)+1,round(j/mbSize)+1) = 0; % 
            MVx(round(i/mbSize)+1,round(j/mbSize)+1) = 0;
            mbCount = mbCount + 1;
            computations =  computations + 1;
%             costs = ones(1,7) * 65537;
            costs = ones(1,11) * 65537;
            checkMatrix = zeros(2*p+1,2*p+1);
            continue;
%% =================== Stage 2. ===================                
        elseif costs(3) <= 1024        
            if i>1 && i< row-mbSize+1 && j>1 && j< col-mbSize+1
                maxIndex = 9;
            else
                maxIndex = 5;
            end          
            for k = 1:maxIndex
                if k == 3 % center point already calculated
                    continue; % costs(3) [0,0] is already calculated
                end
                
                refBlkVer = y + ChenDSP(k,2);   % row/Vert co-ordinate for ref block
                refBlkHor = x + ChenDSP(k,1);   % col/Horizontal co-ordinate
                
                if ( refBlkVer < 1 || refBlkVer+mbSize-1 > row ...
                    || refBlkHor < 1 || refBlkHor+mbSize-1 > col)
                    continue; % outside image boundary
                end
                
                if k <= 5;
                    % General SAD
                    PPP = imgP(i:i+mbSize-1,j:j+mbSize-1); % Horizontal
                    III = imgI(refBlkVer:refBlkVer+mbSize-1, refBlkHor:refBlkHor+mbSize-1);
                    costs(k) = costFuncMAD(PPP, III, mbSize);
                    computations =  computations +1;

                else % k >5 [1,1] [-1,1] [1,-1] [-1 -1]
                    % ========= Two-Quincunx SAD =========
                    PPP = TwoQuincunx(imgP(i:i+mbSize-1,j:j+mbSize-1));
                    III = TwoQuincunx(imgI(refBlkVer:refBlkVer+mbSize-1, refBlkHor:refBlkHor+mbSize-1));
                    
                    costs(k) = costFuncMAD(PPP, III, mbSize)*2;
                    computations =  computations + 0.5;
                end
%                 costs(k) = costFuncMAD(imgP(i:i+mbSize-1,j:j+mbSize-1), ...
%                   imgI(refBlkVer:refBlkVer+mbSize-1, refBlkHor:refBlkHor+mbSize-1), mbSize);
                
                checkMatrix(ChenDSP(k,2) + p+1, ChenDSP(k,1) + p+1) = 1;
            end
            [cost, point] = min(costs);
            x = x + ChenDSP(point, 1);
            y = y + ChenDSP(point, 2);
            
            % Gradient-decent
            refBlkVer = y + ChenDSP(point,2) + ChenDSP(point,2);   % row/Vert co-ordinate for ref block
            refBlkHor = x + ChenDSP(point,1) + ChenDSP(point,1);   % col/Horizontal co-ordinate   
            if ( refBlkVer >= 1 && refBlkVer+mbSize-1 <= row ...
                && refBlkHor >= 1 && refBlkHor+mbSize-1 <= col)
                % outside image boundary
                PPP = TwoQuincunx(imgP(i:i+mbSize-1,j:j+mbSize-1));
                III = TwoQuincunx(imgI(refBlkVer:refBlkVer+mbSize-1, refBlkHor:refBlkHor+mbSize-1));
                cost_grad = costFuncMAD(PPP, III, mbSize)*2;
                computations =  computations + 0.5;
                if cost_grad < cost
                    x = x + ChenDSP(point, 1) + ChenDSP(point, 1);
                    y = y + ChenDSP(point, 2) + ChenDSP(point, 2);
                end
            end
            
            MVy(round(i/mbSize)+1,round(j/mbSize)+1) = y - i; 
            MVx(round(i/mbSize)+1,round(j/mbSize)+1) = x - j; 
            vectors(1,mbCount) = y - i;    % row co-ordinate for the vector
            vectors(2,mbCount) = x - j;    % col co-ordinate for the vector            
            mbCount = mbCount + 1;
%             costs = ones(1,7) * 65537;
            costs = ones(1,11) * 65537;
            checkMatrix = zeros(2*p+1,2*p+1);
            continue;
%% =================== Stage 3. ===================                    
        elseif costs(3) <= 1536     
            if i>1 && i< row-mbSize+1 && j>1 && j< col-mbSize+1
                maxIndex = 11;
            else
                maxIndex = 5;
            end          
            for k = 1:maxIndex
                if k == 3 % center point already calculated
                    continue; % costs(3) [0,0] is already calculated
                end
                
                refBlkVer = y + ChenDSP(k,2);   % row/Vert co-ordinate for ref block
                refBlkHor = x + ChenDSP(k,1);   % col/Horizontal co-ordinate
                
                if ( refBlkVer < 1 || refBlkVer+mbSize-1 > row ...
                    || refBlkHor < 1 || refBlkHor+mbSize-1 > col)
                    continue; % outside image boundary
                end
                
                if k <= 5;
                    % General SAD
                    PPP = imgP(i:i+mbSize-1,j:j+mbSize-1);
                    III = imgI(refBlkVer:refBlkVer+mbSize-1, refBlkHor:refBlkHor+mbSize-1);
                    costs(k) = costFuncMAD(PPP, III, mbSize);
                    computations =  computations +1;

                elseif k<=9   % [1,1] [-1,1] [1,-1] [-1 -1]
                    % ========= Two-Quincunx SAD =========
                    PPP = TwoQuincunx(imgP(i:i+mbSize-1,j:j+mbSize-1));
                    III = TwoQuincunx(imgI(refBlkVer:refBlkVer+mbSize-1, refBlkHor:refBlkHor+mbSize-1));

                    costs(k) = costFuncMAD(PPP, III, mbSize)*2;
                    computations =  computations + 0.5;
                else % k>9
                    % ========= Four-Queen SAD =========
                    PPP = FourQueen(imgP(i:i+mbSize-1,j:j+mbSize-1));
                    III = FourQueen(imgI(refBlkVer:refBlkVer+mbSize-1, refBlkHor:refBlkHor+mbSize-1));
                    
                    costs(k) = costFuncMAD(PPP, III, mbSize)*4;
                    computations =  computations + 0.25;
                end
%                 costs(k) = costFuncMAD(imgP(i:i+mbSize-1,j:j+mbSize-1), ...
%                   imgI(refBlkVer:refBlkVer+mbSize-1, refBlkHor:refBlkHor+mbSize-1), mbSize);
                
                checkMatrix(ChenDSP(k,2) + p+1, ChenDSP(k,1) + p+1) = 1;
            end
            [cost, point] = min(costs);
            x = x + ChenDSP(point, 1);
            y = y + ChenDSP(point, 2);
            
            MVy(round(i/mbSize)+1,round(j/mbSize)+1) = y - i; 
            MVx(round(i/mbSize)+1,round(j/mbSize)+1) = x - j; 
            vectors(1,mbCount) = y - i;    % row co-ordinate for the vector
            vectors(2,mbCount) = x - j;    % col co-ordinate for the vector            
            mbCount = mbCount + 1;
%             costs = ones(1,7) * 65537;
            costs = ones(1,11) * 65537;
            checkMatrix = zeros(2*p+1,2*p+1);
            continue;
 
        end          
%% =================== Stage 4. ===================                        
%% (1) If we are in the left most column then we have to make sure that
  % we just do the LDSP with stepSize = 2                     
        checkMatrix(p+1,p+1) = 1;
        computations =  computations + 1; 
        UP = col/mbSize;
        if j-1 < 1 
            if i-1 >= 1
                stepSize1 = max(abs(vectors(1,mbCount- UP)), abs(vectors(2,mbCount- UP)) ); % UP
            else
                stepSize1 = 2;
            end
            % Original CHEN-ARPS
            stepSize2 = max(abs(MVpre(1,mbCount)), abs(MVpre(2,mbCount))); % temporal previous
%             stepSize = round((stepSize1+stepSize2)/2);
            stepSize = max(stepSize1,stepSize2);
            if stepSize2 == 0
                maxIndex = 5;
            else
                maxIndex = 6;
                LDSP(6,:) = [ MVpre(2, mbCount)  MVpre(1, mbCount)]; % Previous Motion Vector
            end
        else 
                % Original CHEN-ARPS
            stepSize1 = max(abs(vectors(1,mbCount-1)), abs(vectors(2,mbCount-1)) );  % raster left
            stepSize2 = max(abs(MVpre(1,mbCount)), abs(MVpre(2,mbCount))); % temporal previous
%             stepSize = round((stepSize1+stepSize2)/2);
            stepSize = max(stepSize1,stepSize2);

            if vectors(1,mbCount-1) == MVpre(1,mbCount) && vectors(2,mbCount-1) == MVpre(2,mbCount)  % left = previous
% now we have to make sure that if the point due to motion vector is one of 
% the LDSP points then we don't calculate it again                
                if ( (abs(vectors(1,mbCount-1)) == stepSize && vectors(2,mbCount-1) == 0) ...
                    || (abs(vectors(2,mbCount-1)) == stepSize && vectors(1,mbCount-1) == 0) ) ...
%             if ( (abs(MVy(round(i/mbSize)+1,round(j/mbSize)) ) == stepSize && MVx(round(i/mbSize)+1,round(j/mbSize)) == 0) ...
%                  || (abs( MVx(round(i/mbSize)+1,round(j/mbSize)) ) == stepSize && MVy(round(i/mbSize)+1,round(j/mbSize)) == 0) ) ...     
                    maxIndex = 5; % we just have to check at the rood pattern 5 points             
                else
                    maxIndex = 6; % we have to check 6 points
                    LDSP(6,:) = [ vectors(2, mbCount-1)  vectors(1, mbCount-1)]; % Left-side Motion Vector
                end 
            else % vectors(:,mbCount-1) ~= MVpre(:,mbCount)  % left != previous
                if ( (abs(vectors(1,mbCount-1)) == stepSize && vectors(2,mbCount-1) == 0) ...
                    || (abs(vectors(2,mbCount-1)) == stepSize && vectors(1,mbCount-1) == 0) ) ...
                    maxIndex = 6; % we have to check 6 points
                    LDSP(6,:) = [ MVpre(2, mbCount)  MVpre(1, mbCount)]; % Previous Motion Vector
                    
                elseif ( (abs(MVpre(1,mbCount-1)) == stepSize && MVpre(2,mbCount-1) == 0) ...
                    || (abs(MVpre(2,mbCount-1)) == stepSize && MVpre(1,mbCount-1) == 0) ) ...
                    maxIndex = 6; % we have to check 6 points
                    LDSP(6,:) = [ vectors(2, mbCount-1)  vectors(1, mbCount-1) ]; % Left-side Motion Vector
                else
                    maxIndex = 7; % we have to check 7 points
                    LDSP(6,:) = [ vectors(2, mbCount-1)  vectors(1, mbCount-1) ]; % Left-side Motion Vector
                    LDSP(7,:) = [ MVpre(2, mbCount)  MVpre(1, mbCount)]; % Previous Motion Vector
                end
            end
        end
        % The index points for first and only Large Diamond search pattern
        LDSP(1,:) = [ 0 -stepSize];
        LDSP(2,:) = [-stepSize 0]; 
        LDSP(3,:) = [ 0  0];
        LDSP(4,:) = [stepSize  0];
        LDSP(5,:) = [ 0  stepSize];
        
        % do the LDSP
        for k = 1:maxIndex
            refBlkVer = y + LDSP(k,2);   % row/Vert co-ordinate for ref block
            refBlkHor = x + LDSP(k,1);   % col/Horizontal co-ordinate
            if ( refBlkVer < 1 || refBlkVer+mbSize-1 > row ...
                 || refBlkHor < 1 || refBlkHor+mbSize-1 > col)
             
                continue; % outside image boundary
            end

            if checkMatrix(LDSP(k,2) + p+1, LDSP(k,1) + p+1) == 1
                continue;
            end
            
            if k == 3 || stepSize == 0
                continue; % center point already calculated
            end
            costs(k) = costFuncMAD(imgP(i:i+mbSize-1,j:j+mbSize-1), ...
                  imgI(refBlkVer:refBlkVer+mbSize-1, refBlkHor:refBlkHor+mbSize-1), mbSize);
            computations =  computations + 1;
            checkMatrix(LDSP(k,2) + p+1, LDSP(k,1) + p+1) = 1; 
        end
        
            [cost, point] = min(costs);
            x = x + LDSP(point, 1);
            y = y + LDSP(point, 2);


% The 'doneFlag' is set to '1' when the minimum is at the center of the diamond           
        
%         costs = ones(1,5) * 65537;
        costs = ones(1,5) * 65537;
        costs(3) = cost;
       
        doneFlag = 0;   
        while doneFlag == 0
            for k = 1:5
                refBlkVer = y + SDSP(k,2);   % row/Vert co-ordinate for ref block
                refBlkHor = x + SDSP(k,1);   % col/Horizontal co-ordinate
                if ( refBlkVer < 1 || refBlkVer+mbSize-1 > row ...
                      || refBlkHor < 1 || refBlkHor+mbSize-1 > col)
                      continue;
                end

                if k == 3
                    continue;
                elseif refBlkHor < j-p || refBlkHor > j+p || refBlkVer < i-p ...
                            || refBlkVer > i+p
                    continue;
                elseif checkMatrix(y-i+SDSP(k,2)+p+1 , x-j+SDSP(k,1)+p+1) == 1
                    continue;
                end
            
                costs(k) = costFuncMAD(imgP(i:i+mbSize-1,j:j+mbSize-1), ...
                             imgI(refBlkVer:refBlkVer+mbSize-1, ...
                                 refBlkHor:refBlkHor+mbSize-1), mbSize);
                checkMatrix(y-i+SDSP(k,2)+p+1, x-j+SDSP(k,1)+p+1) = 1;
                computations =  computations + 1;
            end
            
            [cost, point] = min(costs);
           
            if point == 3
                doneFlag = 1;
            else
                x = x + SDSP(point, 1);
                y = y + SDSP(point, 2);
                costs = ones(1,5) * 65537;
                costs(3) = cost;
            end
        end  % while loop ends here
        
        MVy(round(i/mbSize)+1,round(j/mbSize)+1) = y - i; 
        MVx(round(i/mbSize)+1,round(j/mbSize)+1) = x - j; 
        vectors(1,mbCount) = y - i;    % row co-ordinate for the vector
        vectors(2,mbCount) = x - j;    % col co-ordinate for the vector            
        mbCount = mbCount + 1;
%         costs = ones(1,7) * 65537;
        costs = ones(1,11) * 65537;
        checkMatrix = zeros(2*p+1,2*p+1);
    end
end
    
motionVect = vectors;
Chencomputations = computations/(mbCount-1) ; 