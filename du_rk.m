function dxy = du_rk(u, parameters, h, pore1)
    %parameters = [f_v gamma_p gamma_o omega];
    %molecule physical parameters
    f_v = parameters(1);
    gamma = parameters(2);
    omega = parameters(3); 
    k = parameters(4); 
    R_l = parameters(5);
    dt = h;
    
    pore_outer_y = pore1(1,:);
    pore_outer_x = pore1(2,:);
    pore_inner_y = pore1(3,:);
    
    diam = diff(pore_outer_x);
    l = diff(pore_outer_y);
    
    external_field = 0;

    for i = 1:length(u)     
        %initialize force voltage
        fx = 0;
        fy = 0;       
        
        %if particle is within the electric field, add force
        %inner pore
        if u(i,1) < max(pore_outer_x) && u(i,1) >= min(pore_outer_x) && u(i,2) > min(pore_inner_y) && u(i,2) < max(pore_inner_y) 
            fx = f_v; %convert to pN
        %if outside pore
%         elseif external_field
%             if u(i,1) < min(pore_outer_x)
%                 z = sqrt(diff([u(i,1) min(pore_outer_x)])^2+(u(i,2))^2);
%                 F = (diam^2/(8*l*z))*f_v;
%                 fx = F*diff([u(i,1) min(pore_outer_x)])/z;
%                 fy = -F*u(i,2)/z;
%             elseif u(i,1) > max(pore_outer_x)
%                 z = sqrt(diff([u(i,1) max(pore_outer_x)])^2+(u(i,2))^2);
%                 F = (diam^2/(8*l*z))*f_v;
%                 fx = -F*diff([u(i,1) max(pore_outer_x)])/z;
%                 fy = F*u(i,2)/z;
%             end
        end

        
        %if just a point particle
        if length(u(:,1)) == 1
            dxy(i,1) = fx - omega*(randn)-gamma*v(i,1);
            dxy(i,2) = fy - omega*(randn)-gamma*v(i,2);
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
            %velocity particle i == 1
            dxy(i,1) = (fx - K2*s_theta2 + omega/sqrt(dt)*(randn))/gamma;
            dxy(i,2) = (fy - K2*c_theta2 + omega/sqrt(dt)*(randn))/gamma;
        elseif i == length(u)
            %velocity particle i == end
            dxy(i,1) = (fx + K1*s_theta1 + omega/sqrt(dt)*(randn))/gamma;
            dxy(i,2) = (fy + K1*c_theta1 + omega/sqrt(dt)*(randn))/gamma;
        else
            %velocity particle 1 < i < end
            dxy(i,1) = (fx + K1*s_theta1 - K2*s_theta2 + omega/sqrt(dt)*(randn))/gamma;
            dxy(i,2) = (fy + K1*c_theta1 - K2*c_theta2 + omega/sqrt(dt)*(randn))/gamma;
        end
    end
end
