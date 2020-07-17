

function [motionVect, FDGDScomputations] = MEFDGDS(imgP, imgI, mbSize, p)

[row,col] = size(imgI);
vectors = zeros(2,row*col/mbSize^2);

% 9 points 
SDSP(1,:) = [ 0 -1];
SDSP(2,:) = [-1  0];
SDSP(3,:) = [ 0  0];
SDSP(4,:) = [ 1  0];
SDSP(5,:) = [ 0  1];

SDSP(6,:) = [ 1  1];
SDSP(7,:) = [ 1  -1];
SDSP(8,:) = [ -1  -1];
SDSP(9,:) = [ -1  1];

costs = ones(1,9) * 65537; 
checkMatrix = zeros(2*p+1,2*p+1);
computations = 0;
mbCount = 1;
for i = 1 : mbSize : row-mbSize+1
    for j = 1 : mbSize : col-mbSize+1
        
        % the Adapive Rood Pattern search starts
        % we are scanning in raster order
        x = j;
        y = i;
        
        costs(3) = costFuncMAD(imgP(i:i+mbSize-1,j:j+mbSize-1), ...
                                    imgI(i:i+mbSize-1,j:j+mbSize-1),mbSize);
                                
%         if costs(3) <= 512 % Zero Motion Prejudgment
%             vectors(1,mbCount) = 0;    % row co-ordinate for the vector
%             vectors(2,mbCount) = 0;    % col co-ordinate for the vector   
%             mbCount = mbCount + 1;
%             computations =  computations + 1;
% %             costs = ones(1,7) * 65537;
%             costs = ones(1,9) * 65537;
%             checkMatrix = zeros(2*p+1,2*p+1);
%             continue;
            
%         else
            checkMatrix(p+1,p+1) = 1;
            computations =  computations + 1; 
            
            doneFlag = 0;   
            while doneFlag == 0
                for k = 1:size(SDSP,1)
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
                    
                    if costs(k) <= costs(3)/2
                        cost = costs(k);
                        costs = ones(1,9) * 65537;
                        costs(3) = cost;
                        x = x + SDSP(k, 1);
                        y = y + SDSP(k, 2);
                        break;
                    end
                end
                [cost, point] = min(costs);
           
                if point == 3
                    doneFlag = 1;
                else
                    x = x + SDSP(point, 1);
                    y = y + SDSP(point, 2);
                    costs = ones(1,9) * 65537;
                    costs(3) = cost;
                end
            end  % while loop ends here
        
        vectors(1,mbCount) = y - i;    % row co-ordinate for the vector
        vectors(2,mbCount) = x - j;    % col co-ordinate for the vector            
        mbCount = mbCount + 1;

        costs = ones(1,9) * 65537;
        checkMatrix = zeros(2*p+1,2*p+1);      
%         end
    end
end

motionVect = vectors;
FDGDScomputations = computations/(mbCount-1) ; 
end