%     .................................................
%             ____  _       _   ____  _____   _        
%            |  _ \| |     |_| |  _ \|  ___| |_|       
%            | |_) | |___   _  | |_) | |___   _        
%            |  _ /|  _  | | | |  _ /|___  | | |       
%            | |   | | | | | | | |    ___| | | |       
%            |_|   |_| |_| |_| |_|   |_____| |_|       
%     .................................................
%     PhiPsi:     a general-purpose computational      
%                 mechanics program written in Fortran.
%     Website:    http://phipsi.top                    
%     Author:     Fang Shi  
%     Contact me: shifang@ustc.edu.cn     

function Plot_Shaped_Cracks(Shaped_Crack_Points)
% This function plots the shaped cracks.

	for i = 1:length(Shaped_Crack_Points)
	        c_Last_Shaped_Points = Shaped_Crack_Points{i};
			
			% Plot the division points.
			% plot(c_Last_Shaped_Points(:,1),c_Last_Shaped_Points(:,2),'o','Color','r','MarkerSize',5)

			% Patch the shaped crack.
			patch(c_Last_Shaped_Points(:,1),c_Last_Shaped_Points(:,2),'white','edgecolor','black','LineWidth',1.0)	
			
			% Patch the shaped crack.
			% patch(c_Last_Shaped_Points(:,1),c_Last_Shaped_Points(:,2),'white','edgecolor','white','LineWidth',1.0)				
	end 
end	
	
