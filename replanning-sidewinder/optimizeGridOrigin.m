function [gamma] = optimizeGridOrigin(gamma0, fp0, olapx, olapy, targetArea, dir)
%% Paula
w = fp0.sizex;
h = fp0.sizey;

deltax = 0.2*w/2;
deltay = 0.2*h/2;

%%
counter = 0;
gamma = gamma0;
opt = false;
delta = [0 0];
while ~opt && abs(gamma(1) - gamma0(1)) <= w/2 && abs(gamma(2) - gamma0(2)) <= h/2
    gamma = gamma + delta;
    grid = grid2D(fp0.sizex, fp0.sizey, olapx, olapy, gamma, targetArea);

    flag = 0;
    for i=1:size(grid,1)
        for j=1:size(grid,2)
            if ~isempty(grid{i,j})
                if norm(grid{i,j} - gamma) < 1e-3
                    ind = [i j];
                    flag = 1;
                    break;
                end
            end
        end
        if flag
            break;
        end
    end

    if exist('ind','var')
        if xor(dir(1), dir(2))
            if dir(1) > 0
                if ind(2) > 1
                    delta = deltax*[-1 0];
                elseif ind(1) == 2
                    delta = deltay*[0 1];
                elseif ind(1) == (size(grid,1) -  1)
                    delta = deltay*[0 -1];
                else
                    opt = true;
                end
            elseif dir(2) < 0
                if ind(1) > 1
                    delta = deltay*[0 1];
                end
            end
        else
            if dir(1) > 0 && dir(2) < 0
                if ind(1) > 1
                    delta = deltay*[0 1];
                elseif ind(2) == 2
                    delta = deltax*[-1 0];
                %elseif ind(1) == (size(grid,1) - 1)
                %    delta = delta + deltay*[0 -1];
                else
                    opt = true;
                end
            end
%             if dir(2) < 0 && abs(dir(1)) <= w
%                 if ind(1) > 1
%                     delta = delta + deltay*[0 1];
%                 else opt = true;
%                 end
%             elseif dir(1) < 0
%                 if ind(1) < sum(~isnan(grid(ind(1),:)))
%                     delta = delta + deltax*[1 0];
%                 else opt = true;
%                 end
%             elseif dir(1) > 0 && abs(dir(2)) <= h
%                 if ind(1) > 1
%                     delta = delta + deltay*[0 1];
%                 elseif ind(2) > 1
%                     delta = delta + deltax*[-1 0];
%                 else
%                     opt = true;
%                 end
%             end
        end
    end
end

%% Diego
%     % NO TILES IN PREVIOUS COLUMNS OR ROWS
%     ok = 0;
%     gamma = gamma0;
% 
%     %DIFFERENT CASES DEPENDING ON THE TOUR DIRECTION
%     %THE NEXT TILE IS TO THE RIGHT
%     if dir(2) == 0 && dir(1) > 0
%         sit = 1;
%     %THE NEXT TILE IS TO THE LEFT    
%     elseif dir(2) == 0 && dir(1) < 0
%         sit = 2;
%     %THE NEXT TILE IS TO THE BOTTOM 1    
%     elseif dir(2) < 0 && dir(1) >= 0  
%         sit = 3;
%     %THE NEXT TILE IS TO THE BOTTOM 2    
%     elseif dir(2) < 0 && dir(1) < 0   
%         sit = 4;
%     end
% 
%     %%%
% %     figure(2)
% %     hold off
% %     plot(targetArea(:,1),targetArea(:,2),'linewidth',3)
% %     hold on 
% %     xlim([0 360]);
% %     ylim([-90 90]);
% %     xlabel('longitude [°]','FontSize', 12)
% %     ylabel('latitude [°]','FontSize', 12)
%     %%%
% 
% 
%     while ok == 0 
% 
%         grid = grid_dicretization(w,h,gamma,targetArea);
%         dim1 = size(grid,1); 
%         dim2 = size(grid,2); 
%         flag = 0;
% 
% 
%         for i = 1:dim1
%             for j = 1:dim2
%                 if ~isempty(grid{i,j})
%                     if mean(grid{i,j} == gamma) == 1
%                         indexes = [i,j];
%                         flag = 1;
%                         break
%                     end    
%                 end
%             end 
%             if flag == 1
%                 break
% 
%             end    
%         end
% 
% 
%         if sit == 1
%             if indexes(1) == 1 && indexes(2) == 1
%                 ok = 1;
%             elseif indexes(1) == 1  &&  indexes(2) ~= 1 
%                 gamma = gamma + [-0.1*(h*0.5), 0];
%             elseif indexes(1) ~= 1 
%                 gamma = gamma + [0, 0.1*(h*0.5)];
%             end    
%         elseif sit == 2
%             if indexes(1) == 1 && indexes(2) == dim2
%                 ok = 1;
%             elseif indexes(1) == 1 && isempty(grid{indexes(1),indexes(2)+1})
%                 ok = 1;
%             elseif indexes(1) == 1 
%                 gamma = gamma + [0.5*(h*0.5), 0];
%             else
%                 gamma = gamma + [0, 0.1*(h*0.5)];
%             end 
%         elseif sit == 3
%             if indexes(1) == 1 && indexes(2) == dim2
%                 ok = 1;
%             elseif indexes(1) == 1 && indexes(2) == 1    
%                 ok = 1;
%             elseif indexes(1) == 1 && isempty(grid{indexes(1),indexes(2)+1})
%                 ok = 1;    
%             elseif indexes(1) ~= 1
%                 gamma = gamma + [0, 0.1*(h*0.5)];
%             end 
%         else
%             ok = 1;
%         end    
%     end    

end