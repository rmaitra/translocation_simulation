function dxy = du_euler(u, parameters, dt, pore1, pore2)
    %parameters = [f_v gamma_p gamma_o omega];
    %molecule physical parameters
    f_v = parameters(1);
    gamma = parameters(2);
    omega = parameters(3); 
    k = parameters(4); 
    R_l = parameters(5);
    
    pore_outer_y = pore1(1,:);
    pore_outer_x = pore1(2,:);
    pore_inner_y = pore1(3,:);
    
    pore2_outer_y = pore2(1,:);
    pore2_outer_x = pore2(2,:);
    pore2_inner_y = pore2(3,:);
    

    for i = 1:length(u)   
        %initialize force voltage
        fx = 0;
        fy = 0;       
        
        %if particle is within the electric field, add force
        %inner pore
        if u(i,1) < max(pore_outer_x) && u(i,1) >= min(pore_outer_x) && u(i,2) > min(pore_inner_y) && u(i,2) < max(pore_inner_y) 
            fx = f_v; %convert to pN
        end

        
        %if just a point particle
        if length(u(:,1)) == 1
            dxy(i,1) = (fx*dt)/gamma;
            dxy(i,2) = (fy*dt)/gamma;
            break
        end
        
        %spring length and angle
        if i < length(u)
            L2 = sqrt((u(i,1)-u(i+1,1))^2+(u(i,2)-u(i+1,2))^2);
            c_theta2 = (u(i,2)-u(i+1,2))/L2;
            s_theta2 = (u(i,1)-u(i+1,1))/L2;
            K2 = (k)*(L2-R_l);
        end
        
        %if particle past the first, then add the pulling force of the
        %previos particle
        if i > 1
            L1 = sqrt((u(i-1,1)-u(i,1))^2+(u(i-1,2)-u(i,2))^2);
            c_theta1 = (u(i-1,2)-u(i,2))/L1;
            s_theta1 = (u(i-1,1)-u(i,1))/L1;
            K1 = (k)*(L1-R_l);
        end


        %deterimine dxy
        if i == 1
            %particle i == 1
            dxy(i,1) = (fx*dt - K2*s_theta2*dt + omega*(min(max(randn, -4), 4)))/gamma;
            dxy(i,2) = (fy*dt - K2*c_theta2*dt + omega*(min(max(randn, -4), 4)))/gamma;
        elseif i == length(u)
            %particle i == end
            dxy(i,1) = (fx*dt + K1*s_theta1*dt + omega*(min(max(randn, -4), 4)))/gamma;
            dxy(i,2) = (fy*dt + K1*c_theta1*dt + omega*(min(max(randn, -4), 4)))/gamma;
        else
            %particle 1 < i < end
            dxy(i,1) = (fx*dt*100 + K1*s_theta1*dt - K2*s_theta2*dt + omega*(min(max(randn, -4), 4)))/gamma;
            dxy(i,2) = (fy*dt*100 + K1*c_theta1*dt - K2*c_theta2*dt + omega*(min(max(randn, -4), 4)))/gamma;
        end
    end
end
%% misc
%     diam = diff(pore_outer_x);
%     l = diff(pore_outer_y);
%     
%     external_field = 1;

        %if outside pore
%         elseif u(i,1) < max(pore2_outer_x) && u(i,1) >= min(pore2_outer_x) && u(i,2) > min(pore2_inner_y) && u(i,2) < max(pore2_inner_y) 
%             fx = f_v;
%         elseif external_field
%             if u(i,1) < min(pore_outer_x)
%                 z = sqrt((u(i,1)-min(pore_outer_x))^2+(u(i,2))^2);
%                 F = (diam^2/(8*l*z))*f_v;
%                 fx = -F*(u(i,1)-min(pore_outer_x))/z;
%                 fy = -F*u(i,2)/z;
%                 if F > f_v
%                     fx = f_v;
%                     fy = 0;
%                 end
%             elseif u(i,1) > max(pore_outer_x) && u(i,1) < min(pore2_outer_x)
%                 %pore 1
%                 z = sqrt((u(i,1)-max(pore_outer_x))^2+(u(i,2))^2);
%                 F = (diam^2/(8*l*z))*f_v;
%                 fx = F*(u(i,1)-max(pore_outer_x))/z;
%                 fy = F*u(i,2)/z;
%                 
%                 %pore 2
%                 z2 = sqrt((u(i,1)-min(pore2_outer_x))^2+(u(i,2))^2);
%                 F2 = (diam^2/(8*l*z2))*f_v;
%                 fx2 = -F2*(u(i,1)-min(pore2_outer_x))/z2;
%                 fy2 = -F2*u(i,2)/z2;
% 
%                 if F > f_v
%                     fx = f_v;
%                     fy = 0;
%                 end
%                 if F2 > f_v
%                     fx2 = f_v;
%                     fy2 = 0;
%                 end
%                      
%                 fx = fx+fx2;
%                 fy = fy+fy2;
%             
%             elseif u(i,1) > max(pore2_outer_x)
%                 z2 = sqrt((u(i,1)-max(pore2_outer_x))^2+(u(i,2))^2);
%                 F2 = (diam^2/(8*l*z2))*f_v;
%                 fx2 = F2*(u(i,1)-max(pore2_outer_x))/z2;
%                 fy2 = F2*u(i,2)/z2;
%                 if F2 > f_v
%                     fx2 = f_v;
%                     fy2 = 0;
%                 end
%                 fx = fx2;
%                 fy = fy2;
%             end