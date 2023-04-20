function [coordinate, intersection] = footprint_coverage_final(et,obs,scinst,rotation,target,targetframe,steps)

%Obtains the instrument code
[code_inst, ~] = cspice_bodn2c(scinst); 
%PRECISION REQUIRED FOR THE VALUES
room = 10; 
%Obtaining the FOV frame and the bounds vectors.
[~, FOV_frame, ~, bounds] = cspice_getfov(code_inst, room); 
%Obtaining body radii
bodyradii = cspice_bodvrd(target,'RADII',3);
averageradii = sum(bodyradii)/length(bodyradii);
%Obtaining observer position in the target body-fixed frame
[dobs, ~] = cspice_spkpos(obs,et,targetframe,'NONE',target);
error = 0;

if rotation == 0 %if no rotation matrix has been introduced the program will obtain that one to the input time from the C-kernels
   try  
    rotation = cspice_pxform(FOV_frame,targetframe,et);
   catch exception
    disp('Pointing data is not available');
    error = 1;
   end    
end    


%Set of boundary x and y values to define the FOV boundary vectors.
%steps = 200; %Number of boudary vectors chosen
ux = linspace(bounds(1,1),bounds(1,2),steps);
uy = linspace(bounds(2,1),bounds(2,3),steps);

top = zeros(3,steps);
bottom = zeros(3,steps);
left = zeros(3,steps);
right = zeros(3,steps);

top(1,:) = ux(:);
top(2,:) = uy(1);
top(3,:) = bounds(3,1);
top = rotation*top;

left(1,:) = ux(end);
left(2,:) = uy(:);
left(3,:) = bounds(3,1);
left = rotation*left;

bottom(1,:) = ux(end:-1:1);
bottom(2,:) = uy(end);
bottom(3,:) = bounds(3,1);
bottom = rotation*bottom;

right(1,:) = ux(1);
right(2,:) = uy(end:-1:1);
right(3,:) = bounds(3,1);
right = rotation*right;

L(:,:,1) = top./vecnorm(top);
L(:,:,2) = left./vecnorm(left);
L(:,:,3) = bottom./vecnorm(bottom);
L(:,:,4) = right./vecnorm(right);

n = 1;

if error == 0 %IF an exception appeared previously means that the data is not available and the program will show a message.
            %%alternative
        A = [1/((1000*bodyradii(1))^2) 0 0;
             0 1/((1000*bodyradii(2))^2) 0;
             0      0     1/((1000*bodyradii(3))^2)];

        S = 1000*dobs; 
       
        V = [0; 0; 0];    %body fixed centered in the target cente
        
        for t = 1:4 % 4 sides

            firstterm = diag(L(:,:,t).'*A*L(:,:,t));            
            secondterm = 2*S.'*A*L(:,:,t); %- 2*L.'*A*V;           
            indep = -2*S.'*A*V +S.'*A*S - 1; %V.'*A*V - 2*S.'*A*V + S.'*A*S - 1;          
            det = secondterm.^2 - 4*firstterm.'*indep;

            for i = 1:steps
                if det(i) == 0
                    k = -secondterm(i)/(2*firstterm(i));
                    intersection(:,n) = (S + k*L(:,i,t))/1000;
                    [~,intlon, intlat] = cspice_reclat(intersection(:,n));
                    coordinate(1,n) = intlon*(180/pi);
                    coordinate(2,n) = intlat*(180/pi);
                    if coordinate(1,n)<0
                       coordinate(1,n) = coordinate(1,n) + 360; 
                    end
                    found = 1;
                    n = n + 1;
                elseif det(i) > 0
                    k1 = (-secondterm(i) - sqrt(det(i)))/(2*firstterm(i));
                    k2 = (-secondterm(i) + sqrt(det(i)))/(2*firstterm(i));
                    k = min(k1,k2);
                    intersection(:,n) = (S + k*L(:,i,t))/1000;
                    [~,intlon, intlat] = cspice_reclat(intersection(:,n));
                    coordinate(1,n) = intlon*(180/pi);
                    coordinate(2,n) = intlat*(180/pi);
                    if coordinate(1,n)<0
                       coordinate(1,n) = coordinate(1,n) + 360; 
                    end
                    found = 1;
%                 else 
%                     coordinate(1,n) = NaN;
%                     coordinate(2,n) = NaN;
%                     found = 0;
                    n = n + 1;
                else 
                   intersection(:,n) = [NaN; NaN; NaN];
                   coordinate(:,n) = [NaN; NaN]; 
                   n = n + 1;
                end 
            end
        end  
    
    coordinate(1,n) = coordinate(1,1);
    coordinate(2,n) = coordinate(2,1);
    intersection(:,n) = intersection(:,1);
else
    coordinate(1,n) = NaN;
    coordinate(2,n) = NaN;
end


end