%% Numerical Methods
% Runge Kutta
% Euler
function [positionf, velocityf] = numerical_solve(position, velocity, dt, boundary, pore1, pore2, parameters, order)
    if order == 2
        %step 1
        u = position;
        v = velocity;
        du1 = v;
        dv1 = dv_rk(u, v, parameters, pore1, pore2);

        %step 2
        u = u + ( dt / 2.0 ) * du1;
        v = v + ( dt / 2.0 ) * dv1;
        du2 = v;
        dv2 = dv_rk(u, v, parameters, pore1, pore2);

        %step 3
        u = u + ( dt / 2.0 ) * du2;
        v = v + ( dt / 2.0 ) * dv2;
        du3 = v;
        dv3 = dv_rk(u, v, parameters, pore1, pore2);

        %step 4
        u = u + dt * du3;
        v = v + dt * dv3;
        du4 = v;
        dv4 = dv_rk(u, v, parameters, pore1, pore2);

        %now combine the derivative estimates and compute new state
        positionf = position + ( dt / 6 ) * ( du1 + du2*2 + du3*2 + du4 );
        velocityf = velocity + ( dt / 6 ) * ( dv1 + dv2*2 + dv3*2 + dv4 );
    
    
    elseif order == 1
        runge = 0;
        euler = 1;
        if runge
            %step 1
            u = position;
            du1 = du_rk(u, parameters, dt, pore1);

            %step 2
            u = u + ( dt / 2.0 ) * du1;
            du2 = du_rk(u, parameters, dt, pore1);

            %step 3
            u = u + ( dt / 2.0 ) * du2;
            du3 = du_rk(u, parameters, dt, pore1);

            %step 4
            u = u + dt * du3;
            du4 = du_rk(u, parameters, dt, pore1);

            dxy = ( dt / 6 )*( du1 + du2*2 + du3*2 + du4 );
        end

        if euler
             u = position;
             dxy = du_euler(u, parameters, dt, pore1, pore2);
        end

        %now combine the derivative estimates and compute new state
        positionf = position + dxy;
        velocityf = velocity;
    end
    
    [positionf velocityf] = collision_detect(position, positionf, velocityf, boundary, pore1);
    [positionf velocityf] = collision_detect(position, positionf, velocityf, boundary, pore2);
end





