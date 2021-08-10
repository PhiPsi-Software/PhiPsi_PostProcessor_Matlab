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

function Plot_Node_Disp_3D(isub,Crack_X,Crack_Y,Crack_Z,Post_Enriched_Nodes,POS,Elem_Type)
% This function plots the 3D nodal dispalcement contour.

global Node_Coor Elem_Node Key_POST_HF
global Num_Node Num_Elem
global Min_X_Coor Max_X_Coor Min_Y_Coor Max_Y_Coor Min_Z_Coor Max_Z_Coor
global Key_PLOT aveg_area_ele
global Size_Font Elem_Fontcolor Elem_Fontsize Node_Fontcolor Node_Fontsize
global num_Crack num_of_Material
global Color_Crack Width_Crack Full_Pathname
global Color_Backgro_Mesh_1 Color_Backgro_Mesh_2 Color_Backgro_Mesh_3 Color_Backgro_Mesh_4
global Color_Backgro_Mesh_5 Color_Backgro_Mesh_6 Color_Backgro_Mesh_7
global Color_Backgro_Mesh_8 Color_Backgro_Mesh_9 Color_Backgro_Mesh_10
global Elem_Material Num_Step_to_Plot DISP Num_Foc_z
global Num_Foc_x Num_Foc_y Foc_x Foc_y Num_Foc_z Foc_z FORCE_Matrix
global Key_Data_Format
global Crack_node_X Crack_node_Y Crack_node_Z
global Crack_Ele_1 Crack_Ele_2 Crack_Ele_3
global Model_Outline Model_OutArea
global Num_Bou_x Num_Bou_y Num_Bou_z Bou_x Bou_y Bou_z
global Crack_node_local_X Crack_node_local_Y Crack_node_local_Z
global Crack_Node_in_Ele
global Crack3D_Vector_S1_X Crack3D_Vector_S1_Y Crack3D_Vector_S1_Z
global Crack3D_CalP_X Crack3D_CalP_Y Crack3D_CalP_Z
global Crack3D_Vertex_Vector_X_X 
global Crack3D_Vertex_Vector_X_Y  
global Crack3D_Vertex_Vector_X_Z 
global Crack3D_Vertex_Vector_Y_X  
global Crack3D_Vertex_Vector_Y_Y  
global Crack3D_Vertex_Vector_Y_Z  
global Crack3D_Vertex_Vector_Z_X  
global Crack3D_Vertex_Vector_Z_Y  
global Crack3D_Vertex_Vector_Z_Z 	
global Crack3D_Meshed_Outline
global Num_Accuracy_Contour Key_Contour_Metd
global volume_ele max_volume_ele min_volume_ele aveg_volume_ele 
global Num_Contourlevel Key_Flipped_Gray Color_Contourlevel
global Color_Mesh
global Cracks_CalP_UpDis_3D_X Cracks_CalP_UpDis_3D_Y Cracks_CalP_UpDis_3D_Z
global Cracks_CalP_LowDis_3D_X Cracks_CalP_LowDis_3D_Y Cracks_CalP_LowDis_3D_Z
global Title_Font Key_Figure_Control_Widget

scale         = Key_PLOT(4, 6);
% Get the new coordinates of all nodes.
New_Node_Coor(1:Num_Node,1) = Node_Coor(1:Num_Node,1) + scale*DISP(1:Num_Node,2);
New_Node_Coor(1:Num_Node,2) = Node_Coor(1:Num_Node,2) + scale*DISP(1:Num_Node,3);
New_Node_Coor(1:Num_Node,3) = Node_Coor(1:Num_Node,3) + scale*DISP(1:Num_Node,4);


% Get the maximum and minimum value of the new coordinates of all nodes.
Min_X_Coor_New = min(min(New_Node_Coor(1:Num_Node,1)));
Max_X_Coor_New = max(max(New_Node_Coor(1:Num_Node,1)));
Min_Y_Coor_New = min(min(New_Node_Coor(1:Num_Node,2)));
Max_Y_Coor_New = max(max(New_Node_Coor(1:Num_Node,2)));
Min_Z_Coor_New = min(min(New_Node_Coor(1:Num_Node,3)));
Max_Z_Coor_New = max(max(New_Node_Coor(1:Num_Node,3)));

c_X_Length = Max_X_Coor_New-Min_X_Coor_New;
c_Y_Length = Max_Y_Coor_New-Min_Y_Coor_New;
c_Z_Length = Max_Z_Coor_New-Min_Z_Coor_New;

% Get resample coors.
delta = (aveg_volume_ele^(1/3))/Num_Accuracy_Contour;
gx    = Min_X_Coor:delta:Max_X_Coor; 
gy    = Min_Y_Coor:delta:Max_Y_Coor;
gz    = Min_Z_Coor:delta:Max_Z_Coor;

% 绘图设置.
Start_Dis_Type =1;     
End_Dis_Type   =3; 

if Key_PLOT(4,2)==1
    Start_Dis_Type =1; 
	End_Dis_Type   =1;
elseif Key_PLOT(4,2)==2
    Start_Dis_Type =2; 
	End_Dis_Type   =2;
elseif Key_PLOT(4,2)==3
    Start_Dis_Type =3; 
	End_Dis_Type   =3;
end

nel = Num_Elem;
nnel = 8;                % number of nodes per element

% 仅绘制非增强单元
if Key_PLOT(4,5)==3
	num_plot_el = 0;
	plot_eles =[];
	if isempty(Elem_Type)==0
		for iElem=1:nel	
			if sum(Elem_Type(iElem,1:size(Elem_Type,2)))==0  %非增强单元
				num_plot_el = num_plot_el +1;
				plot_eles(num_plot_el) = iElem;
			end
		end
	end
% 绘制全部单元
elseif Key_PLOT(4,5)~=3
	num_plot_el = nel;
	plot_eles =[];
	for iElem=1:nel	
		plot_eles(iElem) = iElem;
	end	
end
			


% Initialization of the required matrices
X = zeros(nnel,num_plot_el);
Y = zeros(nnel,num_plot_el);
Z = zeros(nnel,num_plot_el);
profile = zeros(nnel,num_plot_el);



% Plot the contours.
for i_Plot = Start_Dis_Type:End_Dis_Type
    % New figure.
    Tools_New_Figure
	hold on;
	axis off; axis equal;
	delta = sqrt(aveg_area_ele);
	axis([Min_X_Coor_New-delta Max_X_Coor_New+delta ...
		  Min_Y_Coor_New-delta Max_Y_Coor_New+delta ...
		  Min_Z_Coor_New-delta Max_Z_Coor_New+delta]);
	view(-37.5,30)  %视角	
	
	if i_Plot == 1
	    disp('    > Plotting nodal displacement contour in x direction....') 
	    title('Displacement contour in x direction','FontName',Title_Font,'FontSize',Size_Font)
        fm = [1 2 6 5; 2 3 7 6; 3 4 8 7; 4 1 5 8; 1 2 3 4; 5 6 7 8];
        XYZ = cell(1,num_plot_el) ;
        profile = XYZ;
        for iElem=1:num_plot_el	
            nd=Elem_Node(plot_eles(iElem),1:8);
			
            X = New_Node_Coor(nd,1);
            Y = New_Node_Coor(nd,2);
            Z = New_Node_Coor(nd,3);
            XYZ{iElem} = [X Y Z] ;
            profile{iElem} = DISP(nd,2);
        end
        cellfun(@patch,repmat({'Vertices'},1,num_plot_el),XYZ,.......
            repmat({'Faces'},1,num_plot_el),repmat({fm},1,num_plot_el),......
            repmat({'FaceVertexCdata'},1,num_plot_el),profile,......
            repmat({'FaceColor'},1,num_plot_el),repmat({'interp'},1,num_plot_el));	
		set(gca,'XTick',[]) ; set(gca,'YTick',[]); set(gca,'ZTick',[]) ;	
	elseif i_Plot == 2
	    disp('    > Plotting nodal displacement contour in y direction....') 
	    title('Displacement contour in y direction','FontName',Title_Font,'FontSize',Size_Font)
        fm = [1 2 6 5; 2 3 7 6; 3 4 8 7; 4 1 5 8; 1 2 3 4; 5 6 7 8];
        XYZ = cell(1,num_plot_el) ;
        profile = XYZ ;
		
        for iElem=1:num_plot_el	
            nd=Elem_Node(plot_eles(iElem),1:8);
            X = New_Node_Coor(nd,1);
            Y = New_Node_Coor(nd,2);
            Z = New_Node_Coor(nd,3);
            XYZ{iElem} = [X Y Z] ;
            profile{iElem} = DISP(nd,3);
        end
        cellfun(@patch,repmat({'Vertices'},1,num_plot_el),XYZ,.......
            repmat({'Faces'},1,num_plot_el),repmat({fm},1,num_plot_el),......
            repmat({'FaceVertexCdata'},1,num_plot_el),profile,......
            repmat({'FaceColor'},1,num_plot_el),repmat({'interp'},1,num_plot_el));	
		set(gca,'XTick',[]) ; set(gca,'YTick',[]); set(gca,'ZTick',[]) ;	
	elseif i_Plot == 3
	    disp('    > Plotting nodal displacement contour in z direction....') 
	    title('Displacement contour in z direction','FontName',Title_Font,'FontSize',Size_Font)
        fm = [1 2 6 5; 2 3 7 6; 3 4 8 7; 4 1 5 8; 1 2 3 4; 5 6 7 8];
        XYZ = cell(1,num_plot_el) ;
        profile = XYZ ;
		
        for iElem=1:num_plot_el	
            nd=Elem_Node(plot_eles(iElem),1:8);
            X = New_Node_Coor(nd,1);
            Y = New_Node_Coor(nd,2);
            Z = New_Node_Coor(nd,3);
            XYZ{iElem} = [X Y Z] ;
            profile{iElem} = DISP(nd,4);
        end
        cellfun(@patch,repmat({'Vertices'},1,num_plot_el),XYZ,.......
            repmat({'Faces'},1,num_plot_el),repmat({fm},1,num_plot_el),......
            repmat({'FaceVertexCdata'},1,num_plot_el),profile,......
            repmat({'FaceColor'},1,num_plot_el),repmat({'interp'},1,num_plot_el));	
		set(gca,'XTick',[]) ; set(gca,'YTick',[]); set(gca,'ZTick',[]) ;			
	end
	
	% Set colormap.
	if Key_Flipped_Gray==0
		colormap(Color_Contourlevel)
	elseif Key_Flipped_Gray==1
		colormap(flipud(gray))
	end
	
	

	%绘制坐标轴
	Arrow_length = (c_X_Length+c_Y_Length+c_Z_Length)/15;
	h = Tools_mArrow3([0 0 0],[Arrow_length 0 0],'color','red','stemWidth',Arrow_length/25.0,'tipWidth',Arrow_length/10.0,'facealpha',1.0);
	ts = text((c_X_Length+c_Y_Length+c_Z_Length)/14, 0, 0,"x",'Color','red','FontSize',15,'FontName','Consolas','FontAngle','italic');
	h = Tools_mArrow3([0 0 0],[0 Arrow_length 0],'color','green','stemWidth',Arrow_length/25.0,'tipWidth',Arrow_length/10.0,'facealpha',1.0);
	ts = text(0,(c_X_Length+c_Y_Length+c_Z_Length)/14,0,"y",'Color','green','FontSize',15,'FontName','Consolas','FontAngle','italic');
	h = Tools_mArrow3([0 0 0],[0 0 Arrow_length],'color','blue','stemWidth',Arrow_length/25.0,'tipWidth',Arrow_length/10.0,'facealpha',1.0);
	ts = text(0,0,(c_X_Length+c_Y_Length+c_Z_Length)/14,"z",'Color','blue','FontSize',15,'FontName','Consolas','FontAngle','italic');
	


	
    % Set colorbar.
	cbar = colorbar;
	brighten(0.5); 
	set(get(cbar,'title'),'string','m');
	% get the color limits
	clim = caxis;
	ylim(cbar,[clim(1) clim(2)]);
	numpts = 10 ;    % Number of points to be displayed on colorbar
	kssv = linspace(clim(1),clim(2),numpts);
	set(cbar,'YtickMode','manual','YTick',kssv); % Set the tickmode to manual
	for i = 1:numpts
		imep = num2str(kssv(i),'%+3.2E');
		vasu(i) = {imep} ;
	end
	set(cbar,'YTickLabel',vasu(1:numpts),'FontName',Title_Font,'FontSize',Size_Font);
	

	% Active Figure control widget (2021-08-01)
	% Ref: https://ww2.mathworks.cn/matlabcentral/fileexchange/38019-figure-control-widget
	% Press q to exit.
	% Press r (or double-click) to reset to the initial.
	if Key_Figure_Control_Widget==1
		fcw(gca);
	end
	
	% Save pictures.
    if i_Plot == 1
	    Save_Picture(c_figure,Full_Pathname,'Dpxn')
	elseif i_Plot ==2
	    Save_Picture(c_figure,Full_Pathname,'Dpyn')
	elseif i_Plot ==3
	    Save_Picture(c_figure,Full_Pathname,'Dpzn')		
	end	

end


