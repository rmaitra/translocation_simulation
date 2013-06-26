function  animate_translocate(position,boundary,pore1,pore2)
    %initialize drawing window
    figure(1)
    clf
    
      %subplot(3,2,[1 3]);
        %[pvolumex pvolumey] = circle(position,r);
        plot(position(:,1),position(:,2), 'r.')
        hold on
        plot(position(:,1),position(:,2),'b')
        %plot(position(:,1),position(:,2),'b','LineWidth',2)
                
        %hold on
        %plot(pvolumex,pvolumey,'b.')
        grid on
        axis([min(boundary(:,1)) max(boundary(:,1)) min(boundary(:,2)) max(boundary(:,2))])
        %axis([-125 15 -25 20]);
        %axis([position(1,1)-20 position(1,1)+20 position(1,2)+20 position(1,2)+20]);
        xlabel('Position (nm)')
        ylabel('Position (nm)')
        title('Nanopore Translocation Simulation')
        axis([-300 50 -200 200])
        pore_color = 'k';

        %draw the pore1
        pore_outer_y = pore1(1,:);
        pore_outer_x = pore1(2,:);
        pore_inner_y = pore1(3,:);
        pore_inner_x = pore1(4,:);
        width =1;
        line([pore_outer_x(1,1) pore_outer_x(1,1)],[pore_outer_y(1,2) max(boundary(1:2,1))],'Color',pore_color,'LineWidth',width)
        line([pore_outer_x(1,2) pore_outer_x(1,2)],[pore_outer_y(1,2) max(boundary(1:2,1))],'Color',pore_color,'LineWidth',width)
        line(pore_inner_x,[pore_inner_y(1,2) pore_inner_y(1,2)],'Color',pore_color,'LineWidth',width)
        line([pore_inner_x(1,1) pore_outer_x(1,1)],[pore_inner_y(1,2) pore_outer_y(1,2)],'Color',pore_color,'LineWidth',width)
        line([pore_inner_x(1,2) pore_outer_x(1,2)],[pore_inner_y(1,2) pore_outer_y(1,2)],'Color',pore_color,'LineWidth',width)

        line([pore_outer_x(1,1) pore_outer_x(1,1)],[pore_outer_y(1,1) min(boundary(1:2,1))],'Color',pore_color,'LineWidth',width)
        line([pore_outer_x(1,2) pore_outer_x(1,2)],[pore_outer_y(1,1) min(boundary(1:2,1))],'Color',pore_color,'LineWidth',width)
        line(pore_inner_x,[pore_inner_y(1,1) pore_inner_y(1,1)],'Color',pore_color,'LineWidth',width)
        line([pore_inner_x(1,1) pore_outer_x(1,1)],[pore_inner_y(1,1) pore_outer_y(1,1)],'Color',pore_color,'LineWidth',width)
        line([pore_inner_x(1,2) pore_outer_x(1,2)],[pore_inner_y(1,1) pore_outer_y(1,1)],'Color',pore_color,'LineWidth',width)
    
        %draw the pore2
        pore_outer_y = pore2(1,:);
        pore_outer_x = pore2(2,:);
        pore_inner_y = pore2(3,:);
        pore_inner_x = pore2(4,:);
        width =1;
        line([pore_outer_x(1,1) pore_outer_x(1,1)],[pore_outer_y(1,2) max(boundary(1:2,1))],'Color',pore_color,'LineWidth',width)
        line([pore_outer_x(1,2) pore_outer_x(1,2)],[pore_outer_y(1,2) max(boundary(1:2,1))],'Color',pore_color,'LineWidth',width)
        line(pore_inner_x,[pore_inner_y(1,2) pore_inner_y(1,2)],'Color',pore_color,'LineWidth',width)
        line([pore_inner_x(1,1) pore_outer_x(1,1)],[pore_inner_y(1,2) pore_outer_y(1,2)],'Color',pore_color,'LineWidth',width)
        line([pore_inner_x(1,2) pore_outer_x(1,2)],[pore_inner_y(1,2) pore_outer_y(1,2)],'Color',pore_color,'LineWidth',width)

        line([pore_outer_x(1,1) pore_outer_x(1,1)],[pore_outer_y(1,1) min(boundary(1:2,1))],'Color',pore_color,'LineWidth',width)
        line([pore_outer_x(1,2) pore_outer_x(1,2)],[pore_outer_y(1,1) min(boundary(1:2,1))],'Color',pore_color,'LineWidth',width)
        line(pore_inner_x,[pore_inner_y(1,1) pore_inner_y(1,1)],'Color',pore_color,'LineWidth',width)
        line([pore_inner_x(1,1) pore_outer_x(1,1)],[pore_inner_y(1,1) pore_outer_y(1,1)],'Color',pore_color,'LineWidth',width)
        line([pore_inner_x(1,2) pore_outer_x(1,2)],[pore_inner_y(1,1) pore_outer_y(1,1)],'Color',pore_color,'LineWidth',width)
    
    
        
%     subplot(3,2,5);
%         plot([1:length(vel)]*h,position(:,1))
%         grid on
%         xlabel('Time (\mus)')
%         ylabel('Average Velocity (nm/\mus)')
%         title('Average Translocation Velocity in Left Pore')
%         
%     subplot(3,2,2)
%         plot([1:length(vel)]*h,vel)
%         grid on
%         xlabel('Time (\mus)')
%         ylabel('Average Velocity (nm/\mus)')
%         title('Average Translocation Velocity in Left Pore')
%     
%     subplot(3,2,4);
%          plot(([1:length(current(:,1))]*4), current(:,1),'Color','y','LineWidth',2)
%          hold on
%          plot(([1:length(current(:,1))]*4), current(:,2),'Color','m','LineWidth',2)
%          grid on
%          xlabel('Time (\mus)')
%          ylabel('Current (pA)')
%          title('Ionic Current in Pore 1 (left-blue) and Pore 2 (right-red)')
%     
%     subplot(3,2,6);
%          plot(([1:length(current(:,1))]*4), voltage(:,1),'Color','y','LineWidth',2)
%          hold on
%          plot(([1:length(current(:,1))]*4), voltage(:,2),'Color','m','LineWidth',2)
%          grid on
%          xlabel('Time (\mus)')
%          ylabel('Command Voltage (mV)')
%          title('Control Voltage in Pore 1 (left-blue) and Pore 2 (right-red)')
%          
%     pause(0.005)
%     
%     whitebg('k')
end