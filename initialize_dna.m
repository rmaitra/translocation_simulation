function [dna] = initialize_dna(N,kuhn)
%Generates initial position of dna
  
  x=zeros(1,N);   % place to store x locations
  y=zeros(1,N);   % place to store y locations

  x(1)=0.0;           % initial x location
  y(1)=0.0;           % initial y location
  bound = 200;                     % boundary for where DNA can't pass
                       
  for i=1:N          % take N steps
    theta = 360*rand;
    %First three joint
    if i*kuhn < bound
        x(i+1)=x(i)+kuhn*cosd(180);
        y(i+1)=y(i)+kuhn*sind(180);
        continue
    end
    %[x,y] = add_angle(theta,x,y,kuhn,i);
    if theta > 90 && theta < 180
        theta = 180-theta;
        %setting the angle to any number within the range
            x(i+1)=x(i)-kuhn*cosd(theta);
            y(i+1)=y(i)+kuhn*sind(theta);
    elseif theta > 180 && theta < 270
        theta = 270-theta;
        %setting the angle to any number within the range
            x(i+1)=x(i)-kuhn*cosd(theta);
            y(i+1)=y(i)-kuhn*sind(theta);
    elseif theta > 270 && theta < 360
        theta = 360-theta;
        %setting the angle to any number within the range
            x(i+1)=x(i)+kuhn*cosd(theta);
            y(i+1)=y(i)-kuhn*sind(theta);
    else
        %setting the angle to any number within the range
            x(i+1)=x(i)+kuhn*cosd(theta);
            y(i+1)=y(i)+kuhn*sind(theta);
    end
    if sqrt((x(i+1)+bound)^2+(y(i+1))^2) > 100 %if the strand goes past make the angle 180
        if x(i+1) < -bound
            if y(i+1) > 0
                x(i+1)= x(i)+kuhn*cosd(theta);
                y(i+1)= y(i)-kuhn*sind(theta);
            else
                x(i+1)= x(i)+kuhn*cosd(theta);
                y(i+1)= y(i)+kuhn*sind(theta);
            end
        else
            if y(i+1) > 0
                x(i+1)= x(i)-kuhn*cosd(theta);
                y(i+1)= y(i)-kuhn*sind(theta);
            else
                x(i+1)= x(i)-kuhn*cosd(theta);
                y(i+1)= y(i)+kuhn*sind(theta);
            end
        end
    end
  end
 
  
  dna = [x;y];
  plot(x,y);
  hold on
  plot(x,y,'.r');
  xlabel('position (nm)')
  ylabel('position (nm)')
  
  %draw the pore
  line([-100 -100],[10 100],'Color','k','LineWidth',2)
  line([0 0],[10 100],'Color','k','LineWidth',2)
  line([-80 -20],[5 5],'Color','k','LineWidth',2)
  line([-80 -100],[5 10],'Color','k','LineWidth',2)
  line([-20 0],[5 10],'Color','k','LineWidth',2)
  
  line([-100 -100],[-10 -100],'Color','k','LineWidth',2)
  line([0 0],[-10 -100],'Color','k','LineWidth',2)
  line([-80 -20],[-5 -5],'Color','k','LineWidth',2)
  line([-80 -100],[-5 -10],'Color','k','LineWidth',2)
  line([-20 0],[-5 -10],'Color','k','LineWidth',2)
  grid on
 
%% misc

% Angle Rotation     
%         display('in');
%         theta
%         if theta > 270 %if it is angled below 0% then rotate 1 degree per loop until it isn't past -100 nm
%             for j = theta+1:1:270
%                 j
%                 x(i+1)=x(i)+kuhn*cosd(j)
%                 y(i+1)=y(i)-kuhn*sind(j);
%                 if x(i+1) < -100 %once it passes, break from the loop
%                     break
%                 end
%             end
%         elseif theta < 90 %if it is angled below 0% then rotate 1 degree per loop until it isn't past -100 nm
%             for j = theta+1:1:90
%                 j
%                 x(i+1)=x(i)+kuhn*cosd(j);
%                 y(i+1)=y(i)+kuhn*sind(j);
%                 if x(i+1) < -100 %once it passes, break from the loop
%                     break
%                 end
%             end
%         end 

