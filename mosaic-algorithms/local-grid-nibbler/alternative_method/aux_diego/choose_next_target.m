function [next_target] = choose_next_target(prev_target,prev_prev_target,foot_par)

        if prev_target == 0
            if foot_par(1)>foot_par(5)
               next_target = 1;
            else 
               next_target = 5; 
            end    
        elseif prev_target == 1
            if foot_par(1)>0
               next_target = 1;
            else
                if foot_par(8)>0
                    next_target = 8;
                elseif foot_par(7)>0  
                    next_target = 7;
                elseif foot_par(6)>0    
                    next_target = 6;
                elseif foot_par(2)>0  
                    next_target = 2;      
                elseif foot_par(3)>0  
                    next_target = 3;  
                elseif foot_par(4)>0  
                    next_target = 4;
                elseif foot_par(5)>0  
                    next_target = 5;    
                else
                    next_target = 1; 
                end    
            end    
        elseif prev_target > 5
            if foot_par(5)>0
               next_target = 5;
            elseif foot_par(1)>0  
               next_target = 1; 
            elseif prev_prev_target == 5
               next_target = 1;  
            elseif prev_prev_target == 1
               next_target = 5;   
            else
               next_target = 1;    
            end  
        elseif prev_target == 5
            if foot_par(5)>0
               next_target = 5;
            else
                if foot_par(6)>0
                    next_target = 6;
                elseif foot_par(7)>0  
                    next_target = 7;
                elseif foot_par(8)>0    
                    next_target = 8;
                elseif foot_par(4)>0    
                    next_target = 4;
                elseif foot_par(3)>0    
                    next_target = 3;
                elseif foot_par(2)>0    
                    next_target = 2;
                elseif foot_par(1)>0    
                    next_target = 1;  
                else
                    next_target = prev_target;     
                end    
            end
        elseif prev_target > 1 && prev_target < 5
            if prev_prev_target == 5 
                if foot_par(1)>0 
                    next_target = 1;
                elseif foot_par(4)>0 
                    next_target = 4;
                elseif foot_par(3)>0 
                    next_target = 3;
                elseif foot_par(2)>0     
                    next_target = 2;
                else
                    next_target = prev_target;
                end    
            elseif prev_prev_target == 1 
                if foot_par(5)>0   
                    next_target = 5;
                elseif foot_par(2)>0 
                    next_target = 2;
                elseif foot_par(3)>0 
                    next_target = 3;
                elseif foot_par(4)>0     
                    next_target = 4;
                else
                    next_target = prev_target;     
                end 
            else
                    next_target = prev_target;    
            end        
        end   

end