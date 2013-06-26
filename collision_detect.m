function   [positionf velocityf] = collision_detect(position, positionf, velocityf, boundary, pore)

    pore_outer_y = pore(1,:);
    pore_outer_x = pore(2,:);
    pore_inner_y = pore(3,:);
    pore_inner_x = pore(4,:);
    
    
    for i = 1:length(position(:,1))

        %if collided with x boundary\
        collide = find(positionf(:,1) < min(boundary(1:2,1)));
        if collide
            collide = positionf(:,1) < min(boundary(1:2,1));
            new_position = min(boundary(1:2,1))-(positionf(collide,1) - min(boundary(1:2,1)));
            positionf(collide,1) = new_position;
            velocityf(collide,1) = -velocityf(collide,1);
        end

        collide = find(positionf(:,1) > max(boundary(1:2,1)));
        if collide
            collide = positionf(:,1) > max(boundary(1:2,1));
            new_position = max(boundary(1:2,1))-(positionf(collide,1) - max(boundary(1:2,1)));
            positionf(collide,1) = new_position;
            velocityf(collide,1) = -velocityf(collide,1);
        end

        %if collided with y boundary
        collide = find(positionf(:,2) < min(boundary(1:2,2)));
        if collide
            collide = positionf(:,2) < min(boundary(1:2,2));
            new_position = min(boundary(1:2,2))-(positionf(collide,2) - min(boundary(1:2,2)));
            positionf(collide,2) = new_position;
            velocityf(collide,2) = -velocityf(collide,2);
        end

        collide = find(positionf(:,2) > max(boundary(1:2,2)));
        if collide
            collide = positionf(:,2) > max(boundary(1:2,2));
            new_position = max(boundary(1:2,2))-(positionf(collide,2) - max(boundary(1:2,2)));
            positionf(collide,2) = new_position;
            velocityf(collide,2) = -velocityf(collide,2);
        end

        %if collided with nanopore
        %if past left side of pore
        for i = 1:length(position(:,1))
            if positionf(i,1) > min(pore_outer_x) && positionf(i,1) < max(pore_outer_x) && positionf(i,2) < min(pore_inner_y)
                %we know it collided with the lower pore part
                %if collided with left side of pore 
                if positionf(i,2) < min(pore_outer_y)
                    if position(i,1) < min(pore_outer_x)
                        new_position = min(pore_outer_x)-(positionf(i,1) - min(pore_outer_x));
                        positionf(i,1) = new_position;
                        velocityf(i,1) = -velocityf(i,1);
                    %if collided with right side of pore 
                    elseif position(i,1) > max(pore_outer_x)
                        new_position = max(pore_outer_x)-(positionf(i,1) - max(pore_outer_x));
                        positionf(i,1) = new_position;
                        velocityf(i,1) = -velocityf(i,1);
                    end

                %if collided with angled bottom left side of pore
                elseif positionf(i,1) < min(pore_inner_x) && positionf(i,2) < -(pore_outer_y(1,2)-(1/4)*(positionf(i,1)-pore_outer_x(1,1)))
                    dy = abs(positionf(i,2)-((1/4)*(positionf(i,1)-pore_outer_x(1,1))-pore_outer_y(1,2)));
                    dx = abs(positionf(i,1)-((4)*(positionf(i,2)+pore_outer_y(1,2))+pore_outer_x(1,1)));
                    theta = 2*(180+26.565);
                    new_positionY = ((1/4)*(positionf(i,1)-pore_outer_x(1,1))-pore_outer_y(1,2)) + dy;
                    new_positionX = ((4)*(positionf(i,2)+pore_outer_y(1,2))+pore_outer_x(1,1)) + dx;
                    positionf(i,2) = new_positionY; %update y
                    positionf(i,1) = new_positionX; %update x
                    newvelocityX = velocityf(i,1)*cosd(theta)+velocityf(i,2)*sind(theta);
                    newvelocityY = velocityf(i,1)*sind(theta)-velocityf(i,2)*cosd(theta);
                    velocityf(i,1) = newvelocityX;
                    velocityf(i,2) = newvelocityY;


                %if collided with angled bottom right side of pore 
                elseif positionf(i,1) > max(pore_inner_x) && positionf(i,2) < (pore_inner_y(1,1)-(1/4)*(positionf(i,1)-pore_inner_x(1,2)))
                    dy = abs(positionf(i,2)-(-pore_inner_y(1,2)-(1/4)*(positionf(i,1)-pore_inner_x(1,2))));
                    dx = abs(positionf(i,1)+(-pore_inner_x(1,2)+(4)*(positionf(i,2)+pore_inner_y(1,2))));
                    theta = 2*(180-26.565);
                    new_positionY = (-pore_inner_y(1,2)-(1/4)*(positionf(i,1)-pore_inner_x(1,2))) + dy;
                    new_positionX = -(-pore_inner_x(1,2)+(4)*(positionf(i,2)+pore_inner_y(1,2))) - dx;
                    positionf(i,2) = new_positionY; %update y
                    positionf(i,1) = new_positionX; %update x
                    newvelocityX = velocityf(i,1)*cosd(theta)+velocityf(i,2)*sind(theta);
                    newvelocityY = velocityf(i,1)*sind(theta)-velocityf(i,2)*cosd(theta);
                    velocityf(i,1) = newvelocityX;
                    velocityf(i,2) = newvelocityY;
                %if collided with inside bottom side of pore 
                elseif position(i,2) > min(pore_inner_y) && position(i,1) < pore_inner_x(1,2) && position(i,1) > pore_inner_x(1,1)
                    new_position = min(pore_inner_y)-(positionf(i,2) - min(pore_inner_y));
                    positionf(i,2) = new_position;
                    velocityf(i,2) = -velocityf(i,2);
                end
            elseif positionf(i,1) > min(pore_outer_x) && positionf(i,1) < max(pore_outer_x) && positionf(i,2) > max(pore_inner_y)
                %we know it collided with the top pore part
                %if collided with left side of pore 
                if positionf(i,2) > max(pore_outer_y)
                    if position(i,1) < min(pore_outer_x)
                        new_position = min(pore_outer_x)-(positionf(i,1) - min(pore_outer_x));
                        positionf(i,1) = new_position;
                        velocityf(i,1) = -velocityf(i,1);
                    %if collided with right side of pore 
                    elseif position(i,1) > max(pore_outer_x)
                        new_position = max(pore_outer_x)-(positionf(i,1) - max(pore_outer_x));
                        positionf(i,1) = new_position;
                        velocityf(i,1) = -velocityf(i,1);
                    end



                %if collided with angled top left side of pore
                elseif positionf(i,1) < min(pore_inner_x) && positionf(i,2) > (pore_outer_y(1,2)-(1/4)*(positionf(i,1)-pore_outer_x(1,1)))
                    dy = abs(positionf(i,2)-(pore_outer_y(1,2)-(1/4)*(positionf(i,1)-pore_outer_x(1,1))));
                    dx = abs(positionf(i,1)+(-pore_outer_x(1,1)+(4)*(positionf(i,2)-pore_outer_y(1,2))));
                    theta = 2*(180-26.565);
                    new_positionY = (pore_outer_y(1,2)-(1/4)*(positionf(i,1)-pore_outer_x(1,1))) - dy;
                    new_positionX = -(-pore_outer_x(1,1)+(4)*(positionf(i,2)-pore_outer_y(1,2))) + dx;
                    positionf(i,2) = new_positionY; %update y
                    positionf(i,1) = new_positionX; %update x
                    newvelocityX = velocityf(i,1)*cosd(theta)+velocityf(i,2)*sind(theta);
                    newvelocityY = velocityf(i,1)*sind(theta)-velocityf(i,2)*cosd(theta);
                    velocityf(i,1) = newvelocityX;
                    velocityf(i,2) = newvelocityY;



                %if collided with angled top right side of pore 
                elseif positionf(i,1) > max(pore_inner_x) && positionf(i,2) > (pore_inner_y(1,2)+(1/4)*(positionf(i,1)-pore_inner_x(1,2)))
                    dx = abs(positionf(i,1)-(pore_inner_x(1,2)+(4)*(positionf(i,2)-pore_inner_y(1,2))));
                    dy = abs(positionf(i,2)-(pore_inner_y(1,2)+(1/4)*(positionf(i,1)-pore_inner_x(1,2))));
                    theta = 2*(180+26.565);
                    new_positionY = (pore_inner_y(1,2)+(1/4)*(positionf(i,1)-pore_inner_x(1,2))) - dy;
                    new_positionX = (pore_inner_x(1,2)+(4)*(positionf(i,2)-pore_inner_y(1,2))) - dx;
                    positionf(i,2) = new_positionY; %update y
                    positionf(i,1) = new_positionX; %update x
                    newvelocityX = velocityf(i,1)*cosd(theta)+velocityf(i,2)*sind(theta);
                    newvelocityY = velocityf(i,1)*sind(theta)-velocityf(i,2)*cosd(theta);
                    velocityf(i,1) = newvelocityX;
                    velocityf(i,2) = newvelocityY;
                %if collided with inside top side of pore 
                elseif position(i,2) < max(pore_inner_y) && position(i,1) < pore_inner_x(1,2) && position(i,1) > pore_inner_x(1,1)
                    new_position = max(pore_inner_y)-(positionf(i,2) - max(pore_inner_y));
                    positionf(i,2) = new_position;
                    velocityf(i,2) = -velocityf(i,2);
                end
            end
        end
    end

    
    %% SPHERICAL DETECTION
    
    %[pvolumex pvolumey] = circle(position,2.5);
    %[pfvolumex pfvolumey] = circle(positionf,2.5)
            %if collided with other sphere
%         for j = 1:length(position(:,1))
%             if j == i
%                 continue
%             end
%             L = sqrt(diff([positionf(i,1),positionf(j,1)])^2+diff([positionf(i,2),positionf(j,2)])^2);
%             if L < 5
%                z = 5-L;
%                x = abs(diff([positionf(i,1),positionf(j,1)]));
%                y = abs(diff([positionf(i,2),positionf(j,2)]));
%                theta = tand(y/x);
%                dx = z*cosd(theta);
%                dy = z*sind(theta);
%                if positionf(i,2) > positionf(j,2)
%                   positionf(i,2) = positionf(i,2)+dy;
%                   positionf(j,2) = positionf(j,2)-dy;
%                else
%                   positionf(i,2) = positionf(i,2)-dy;
%                   positionf(j,2) = positionf(j,2)+dy;
%                end
%                if positionf(i,1) > positionf(j,1)
%                   positionf(i,1) = positionf(i,1)+dx;
%                   positionf(j,1) = positionf(j,1)-dx;
%                else
%                   positionf(i,1) = positionf(i,1)-dx;
%                   positionf(j,1) = positionf(j,1)+dx;
%                end
%                newvelocity2 = velocityf(i,:);
%                newvelocity1 = velocityf(j,:);
%                velocityf(i,:) = newvelocity1;
%                velocityf(j,:) = newvelocity2;
%             end
%         end
%         
%         
%         %if collided with x boundary
%         if sum(pfvolumex(i,:) < min(boundary(1:2,1)))
%             dx = abs(mean(pfvolumex(i,pfvolumex(i,:) < min(boundary(1:2,1)))) - min(boundary(1:2,1)));
%             positionf(i,1) = positionf(i,1)+2*dx;
%             velocityf(i,1) = -velocityf(i,1);
%         end
%         if sum(pfvolumex(i,:) > max(boundary(1:2,1)))
%             dx = abs(mean(pfvolumex(i,pfvolumex(i,:) > max(boundary(1:2,1)))) - max(boundary(1:2,1)));
%             positionf(i,1) = positionf(i,1)-2*dx;
%             velocityf(i,1) = -velocityf(i,1);
%         end
% 
%         %if collided with y boundary
%         if sum(pfvolumey(i,:) < min(boundary(1:2,2)))
%             dy = abs(mean(pfvolumey(i,pfvolumey(i,:) < min(boundary(1:2,2)))) - min(boundary(1:2,2)));
%             positionf(i,:);
%             mean(pfvolumey(i,pfvolumey(i,:) < min(boundary(1:2,2))));
%             positionf(i,2) = positionf(i,2)+2*dy;
%             velocityf(i,2) = -velocityf(i,2);
%         end
%         if sum(pfvolumey(i,:) > max(boundary(1:2,2)))
%             dy = abs(mean(pfvolumey(i,pfvolumey(i,:) > max(boundary(1:2,2)))) - max(boundary(1:2,2)));
%             positionf(i,:);
%             mean(pfvolumey(i,pfvolumey(i,:) > max(boundary(1:2,2))));
%             positionf(i,2) = positionf(i,2)-2*dy;
%             velocityf(i,2) = -velocityf(i,2);
%         end
%         if sum(pfvolumex(i,:) > min(pore_outer_x)) && sum(pfvolumex(i,:) < max(pore_outer_x)) && sum(pfvolumey(i,:) < min(pore_inner_y))
%             %we know it collided with the lower pore part
%             %if collided with left side of pore 
%             if sum(pfvolumey(i,:) < min(pore_outer_y))
%                 if sum(pvolumex(i,:) < min(pore_outer_x))
%                     new_position = min(pore_outer_x)-(mean(pfvolumex(i,pfvolumex(i,:) > min(pore_outer_x))) - min(pore_outer_x))-2.5;
%                     positionf(i,1) = new_position;
%                     velocityf(i,1) = -velocityf(i,1);
%                 %if collided with right side of pore 
%                 elseif sum(pvolumex(i,:) > max(pore_outer_x))
%                     new_position = max(pore_outer_x)-(mean(pfvolumex(i,pfvolumex(i,:) < max(pore_outer_x))) - max(pore_outer_x))+2.5;
%                     positionf(i,1) = new_position;
%                     velocityf(i,1) = -velocityf(i,1);
%                 end
%         
%                 
%                 
%             %if collided with angled bottom left side of pore
%             elseif sum(pfvolumex(i,:) < min(pore_inner_x)) && sum(pfvolumey(i,:) < -(pore_outer_y(1,2)-(1/4)*(pfvolumex(i,:)-pore_outer_x(1,1))))
%                 theta = 2*(180+26.565);
%                 maxx = mean(pfvolumex(i,pfvolumex(i,:) < min(pore_inner_x)));
%                 maxy = mean(pfvolumey(i,pfvolumey(i,:) < -(pore_outer_y(1,2)-(1/4)*(pfvolumex(i,:)-pore_outer_x(1,1)))));
%                 dy = abs(maxy-((1/4)*(maxx-pore_outer_x(1,1))-pore_outer_y(1,2)));
%                 dx = abs(maxx-((4)*(maxy+pore_outer_y(1,2))+pore_outer_x(1,1)));
%                 positionf(i,2) = positionf(i,2)+dy; %update y
%                 positionf(i,1) = positionf(i,1); %update x
%                 newvelocityX = velocityf(i,1)*cosd(theta)+velocityf(i,2)*sind(theta);
%                 newvelocityY = velocityf(i,1)*sind(theta)-velocityf(i,2)*cosd(theta);
%                 velocityf(i,1) = newvelocityX;
%                 velocityf(i,2) = newvelocityY/2;
%             %if collided with angled bottom right side of pore 
%             elseif sum(pfvolumex(i,:) > max(pore_inner_x)) && sum(pfvolumey(i,:) < (pore_inner_y(1,1)-(1/4)*(pfvolumex(i,:)-pore_inner_x(1,2))))
%                 theta = 2*(180-26.565);
%                 maxx = mean(pfvolumex(i,pfvolumex(i,:) > max(pore_inner_x)));
%                 maxy = mean(pfvolumey(i,pfvolumey(i,:) < (pore_inner_y(1,1)-(1/4)*(pfvolumex(i,:)-pore_inner_x(1,2)))));
%                 dy = abs(maxy-(-pore_inner_y(1,2)-(1/4)*(maxx-pore_inner_x(1,2))));
%                 dx = abs(maxx+(-pore_inner_x(1,2)+(4)*(maxy+pore_inner_y(1,2))));
%                 positionf(i,2) = positionf(i,2)+dy; %update y
%                 positionf(i,1) = positionf(i,1); %update x
%                 newvelocityX = velocityf(i,1)*cosd(theta)+velocityf(i,2)*sind(theta);
%                 newvelocityY = velocityf(i,1)*sind(theta)-velocityf(i,2)*cosd(theta);
%                 velocityf(i,1) = newvelocityX;
%                 velocityf(i,2) = newvelocityY/2;
%                 
%             
%                 
%                 
%             %if collided with inside bottom side of pore 
%             elseif sum(pvolumey(i,:) > min(pore_inner_y))&& sum(pvolumex(i,:) < pore_inner_x(1,2)) && sum(pvolumex(i,:) > pore_inner_x(1,1))
%                 dy = abs(mean(pfvolumey(i,pfvolumey(i,:) < min(pore_inner_y))))+min(pore_inner_y);
%                 positionf(i,2) = positionf(i,2)+dy;
%                 velocityf(i,2) = -velocityf(i,2);
%             end
%         elseif sum(pfvolumex(i,:) > min(pore_outer_x)) && sum(pfvolumex(i,:) < max(pore_outer_x)) && sum(pfvolumey(i,:) > max(pore_inner_y))
%             %we know it collided with the top pore part
%             %if collided with left side of pore 
%             if sum(pfvolumey(i,:) > max(pore_outer_y))
%                 if sum(pvolumex(i,:) < min(pore_outer_x))
%                     new_position = min(pore_outer_x)-(mean(pfvolumex(i,pfvolumex(i,:) > min(pore_outer_x))) - min(pore_outer_x))-2.5;
%                     positionf(i,1) = new_position;
%                     velocityf(i,1) = -velocityf(i,1);
%                 %if collided with right side of pore 
%                 elseif sum(pvolumex(i,:) > max(pore_outer_x))
%                     new_position = max(pore_outer_x)-(mean(pfvolumex(i,pfvolumex(i,:) < max(pore_outer_x))) - max(pore_outer_x))+2.5;
%                     positionf(i,1) = new_position;
%                     velocityf(i,1) = -velocityf(i,1);
%                 end
%    
%                 
%                 
%             %if collided with angled top left side of pore
%             elseif sum(pfvolumex(i,:) < min(pore_inner_x)) && sum(pfvolumey(i,:) > (pore_outer_y(1,2)-(1/4)*(pfvolumex(i,:)-pore_outer_x(1,1))))                  
%                 theta = 2*(180-26.565);
%                 maxx = mean(pfvolumex(i,pfvolumex(i,:) < min(pore_inner_x)));
%                 maxy = mean(pfvolumey(i,pfvolumey(i,:) > (pore_outer_y(1,2)-(1/4)*(pfvolumex(i,:)-pore_outer_x(1,1)))));
%                 dy = abs(maxy-(pore_outer_y(1,2)-(1/4)*(maxx-pore_outer_x(1,1))));
%                 dx = abs(maxx+(-pore_outer_x(1,1)+(4)*(maxy-pore_outer_y(1,2))));
%                 positionf(i,2) = positionf(i,2)-dy; %update y
%                 positionf(i,1) = positionf(i,1); %update x
%                 newvelocityX = velocityf(i,1)*cosd(theta)+velocityf(i,2)*sind(theta);
%                 newvelocityY = velocityf(i,1)*sind(theta)-velocityf(i,2)*cosd(theta);
%                 velocityf(i,1) = newvelocityX;
%                 velocityf(i,2) = newvelocityY/2;
%             %if collided with angled top right side of pore 
%             elseif sum(pfvolumex(i,:) > max(pore_inner_x)) && sum(pfvolumey(i,:) > (pore_inner_y(1,2)+(1/4)*(pfvolumex(i,:)-pore_inner_x(1,2))))
%                 theta = 2*(180+26.565);
%                 maxx = mean(pfvolumex(i,pfvolumex(i,:) > max(pore_inner_x)));
%                 maxy = mean(pfvolumey(i,pfvolumey(i,:) > (pore_inner_y(1,2)+(1/4)*(pfvolumex(i,:)-pore_inner_x(1,2)))));
%                 dx = abs(maxx-(pore_inner_x(1,2)+(4)*(maxy-pore_inner_y(1,2))));
%                 dy = abs(maxy-(pore_inner_y(1,2)+(1/4)*(maxx-pore_inner_x(1,2))));
%                 positionf(i,2) = positionf(i,2)-dy; %update y
%                 positionf(i,1) = positionf(i,1); %update x
%                 newvelocityX = velocityf(i,1)*cosd(theta)+velocityf(i,2)*sind(theta);
%                 newvelocityY = velocityf(i,1)*sind(theta)-velocityf(i,2)*cosd(theta);
%                 velocityf(i,1) = newvelocityX;
%                 velocityf(i,2) = newvelocityY/2;
%             
%                 
%                 
%             %if collided with inside top side of pore 
%             elseif sum(pvolumey(i,:) < max(pore_inner_y)) && sum(pvolumex(i,:) < pore_inner_x(1,2)) && sum(pvolumex(i,:) > pore_inner_x(1,1))
%                 dy = abs(mean(pfvolumey(i,pfvolumey(i,:) > max(pore_inner_y))))-max(pore_inner_y);
%                 positionf(i,2) = positionf(i,2)-dy;
%                 velocityf(i,2) = -velocityf(i,2);
%             end
%         end
%     end