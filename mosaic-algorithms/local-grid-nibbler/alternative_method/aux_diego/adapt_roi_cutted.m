function roi = adapt_roi_cutted(roi_old)
  
        k = 1;

        for i = 1:size(roi_old,1) 
            if i == 1
                roi(k,:) = roi_old(i,:);
                k = k + 1;
            elseif abs(roi_old(i,1)-roi_old(i-1,1))>200
                if roi_old(i-1,1)>0
                    midp1 = 0.5*[roi_old(i-1,1)+roi_old(i,1)+360, roi_old(i-1,2)+roi_old(i,2)]; 
                    midp2 = 0.5*[roi_old(i-1,1)+roi_old(i,1)-360, roi_old(i-1,2)+roi_old(i,2)];
                else
                    midp1 = 0.5*[roi_old(i-1,1)+roi_old(i,1)-360, roi_old(i-1,2)+roi_old(i,2)];
                    midp2 = 0.5*[roi_old(i-1,1)+roi_old(i,1)+360, roi_old(i-1,2)+roi_old(i,2)];
                end
                roi(k,:) = midp1;
                k = k + 1;
                roi(k,:) = [NaN, NaN];
                k = k + 1;
                roi(k,:) = midp2;
                k = k + 1;
                roi(k,:) = roi_old(i,:);
                k = k + 1;
            else
                roi(k,:) = roi_old(i,:);
                k = k + 1;
            end       
        end

        if sum(isnan(roi(:,1))) == 2
           roi_old_2 = roi;
           roi = [];
           indexes = find(isnan(roi_old_2(:,1))); 
           roi = [roi_old_2(indexes(1)+1:end,1) roi_old_2(indexes(1)+1:end,2);
                         roi_old_2(1:indexes(1)-1,1) roi_old_2(1:indexes(1)-1,2)];
        else
           if roi(end,1)>0
               roi(end+1,:) = 0.5*[roi(end,1)+roi(1,1)+360, roi(end,2)+roi(1,2)];
               roi = [0.5*(roi(end-1,1)+roi(1,1)-360), 0.5*(roi(end-1,2)+roi(1,2));
                         roi]; 
           else
               roi(end+1,:) = 0.5*[roi(end,1)+roi(1,1)-360, roi(end,2)+roi(1,2)];
               roi = [0.5*(roi(end-1,1)+roi(1,1)+360), 0.5*(roi(end-1,2)+roi(1,2));
                         roi];    
           end    
        end    

end