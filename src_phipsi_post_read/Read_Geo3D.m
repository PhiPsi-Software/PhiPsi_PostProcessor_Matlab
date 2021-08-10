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

function Read_Geo3D
% This function read information about nodes and elements from *node and *elem files which can be obtained by ansys mac file, Ansys2PhiPsi_3D.mac.
% After the reading, the extreme values will be calculated and stored as global values.

global Full_Pathname Node_Coor Elem_Node Bou_x Bou_y Bou_z Foc_x Foc_y Foc_z
global Num_Node Num_Elem Num_Bou_x Num_Bou_y Num_Bou_z Num_Foc_x Num_Foc_y Num_Foc_z
global Min_X_Coor Max_X_Coor Min_Y_Coor Max_Y_Coor Min_Z_Coor Max_Z_Coor Outline Inline
global max_area_ele min_area_ele aveg_area_ele Centroid_elem
global Elemet_Info Node_Info 
global Yes_Checked Crack Elem_Material num_of_Material
global G_NN G_X_NODES G_Y_NODES G_Z_NODES G_X_Min G_X_Max G_Y_Min G_Y_Max G_Z_Min G_Z_Max
global volume_ele max_volume_ele min_volume_ele aveg_volume_ele 

% Read files
disp('    Reading node file....') 

% Check if Full_Pathname exist or not.
if size(Full_Pathname,2) ==1
    disp(['    Error :: Filename is not defined.'])
	Error_Message
end

% Check if *.node file exists or not.
Check_exist = [Full_Pathname,'.node'];
if exist(Check_exist,'file') ==0
    disp(['    Error :: Can not find input files.'])
	Error_Message
end

Node_Coor = load([Full_Pathname,'.node']);

disp('    Reading elem file....') 
Elem_Info = load([Full_Pathname,'.elem']);
Elem_Node = Elem_Info(:,1:8);
Elem_Material = Elem_Info(:,9);       % Material number of elements
num_of_Material = max(Elem_Material); % Number of materials

disp('    Reading boux file....') 
Bou_x = load([Full_Pathname,'.boux']);
disp('    Reading bouy file....') 
Bou_y = load([Full_Pathname,'.bouy']);
disp('    Reading bouz file....') 
Bou_z = load([Full_Pathname,'.bouz']);
disp('    Reading focx file....') 
Foc_x = load([Full_Pathname,'.focx']);
disp('    Reading focy file....') 
Foc_y = load([Full_Pathname,'.focy']);
disp('    Reading focz file....') 
Foc_z = load([Full_Pathname,'.focz']);
disp('    Geometry files reading done.') 

% Get the numbers of nodes, elements, boundaries and forces.
Num_Node = size(Node_Coor,1);
Num_Elem = size(Elem_Node,1);
Num_Bou_x= size(Bou_x,1);
Num_Bou_y= size(Bou_y,1);
Num_Bou_z= size(Bou_z,1);
Num_Foc_x= size(Foc_x,1);
Num_Foc_y= size(Foc_y,1);
Num_Foc_z= size(Foc_z,1);

% Initial Element_Info Matrix.
Elemet_Info = zeros(Num_Elem,10);
for i=1:Num_Elem
    Elemet_Info(i,1) = i;
end
% Initial Node_Info Matrix.
Node_Info   = zeros(Num_Node,10);
for i=1:Num_Node
    Node_Info(i,1) = i;
end

% Initial "Yes_Checked" Matrix, flag for initiation check.
Yes_Checked = zeros(Num_Elem,1);

% Get the max and min value of node coordinates.
Max_X_Coor = max(max(Node_Coor(:,1)));
Min_X_Coor = min(min(Node_Coor(:,1)));
Max_Y_Coor = max(max(Node_Coor(:,2)));
Min_Y_Coor = min(min(Node_Coor(:,2)));
Max_Z_Coor = max(max(Node_Coor(:,3)));
Min_Z_Coor = min(min(Node_Coor(:,3)));
% Check is initial crack out of range or not.
% for i=1:length(Crack)
% end

% Find outline of quadrilateral mesh.
% tem = sort([Elem_Node(:,1) Elem_Node(:,2); Elem_Node(:,2) Elem_Node(:,3); 
          % Elem_Node(:,3) Elem_Node(:,4); Elem_Node(:,4) Elem_Node(:,1) ]')';
% [u,m,n] = unique(tem,'rows');      % determine uniqueness of edges
% counts  = accumarray(n(:), 1);     % determine counts for each unique edge
% Outline = u(counts==1,:);          % extract edges that only occurred once and store as a public value 

% Sort Outline by end to end.
% [New_Outline] = Tools_Srot_by_End_to_End(Outline);

% Check if hole exists or not, if exists, then output "Inline".
% if size(New_Outline,1) < size(Outline,1)
    % Hole exists.
	% Inline=Outline;
	% for ii=1:size(New_Outline,1)
	    % case1 = [New_Outline(ii,1) New_Outline(ii,2)];
		% case2 = [New_Outline(ii,2) New_Outline(ii,1)];
		% for jj =1:size(Outline,1)
		    % if case1 == Outline(jj,:) | case2 == Outline(jj,:)
			    % Inline(jj,:)=[0 0];
			% end
		% end
	% end
	% Inline(all(Inline==0,2),:)=[]; 
	% Sort Inline by end to end.
	% Tem_Inline = Inline(1,:);
	% c_Inline   =1;
	% for i = 1:size(Inline,1);
		% if i==1
			% c_Inline = 1;
			% c_node_num_2 = Inline(1,2);
		% end
		% for j = 1:size(Inline,1)
			% if j ~= c_Inline
				% if c_node_num_2 == Inline(j,1)
					% if ismember([Inline(j,1) Inline(j,2)],Tem_Inline,'rows') ==0
						% Tem_Inline = [Tem_Inline; [Inline(j,1) Inline(j,2)]];
						% c_Inline = j;
						% c_node_num_2 = Inline(j,2);
						% continue
					% end
				% elseif c_node_num_2 == Inline(j,2)
					% if ismember([Inline(j,2) Inline(j,1)],Tem_Inline,'rows') ==0
						% Tem_Inline = [Tem_Inline; [Inline(j,2) Inline(j,1)]];
						% c_Inline = j;
						% c_node_num_2 = Inline(j,1);
						% continue
					% end
				% end
			% end
		% end
	% end
	% Outline = New_Outline;
	% Inline  = Tem_Inline;
% else
    % Hole does not exist.
    % Outline = New_Outline;
	% Inline  = [];
% end

% Get the maximum & minimum & average volume.
for i=1:Num_Elem
    N1  = Elem_Node(i,1);                                                  % Node 1 for current element
    N2  = Elem_Node(i,2);                                                  % Node 2 for current element
    N3  = Elem_Node(i,3);                                                  % Node 3 for current element
    N4  = Elem_Node(i,4);                                                  % Node 4 for current element
    N5  = Elem_Node(i,5);                                                  % Node 1 for current element
    N6  = Elem_Node(i,6);                                                  % Node 2 for current element
    N7  = Elem_Node(i,7);                                                  % Node 3 for current element
    N8  = Elem_Node(i,8);                                                  % Node 4 for current element	
	NN = [N1 N2 N3 N4 N5 N6 N7 N8];                                        % 8 nodes of the current element  
	X_NODES = Node_Coor(NN,1);                                             % Coordinates of the four nodes of the current element
	Y_NODES = Node_Coor(NN,2);
	Z_NODES = Node_Coor(NN,3);
	%--------------------------------
	%Get the volume of each element.
	%--------------------------------
	shp = alphaShape(X_NODES,Y_NODES,Z_NODES); %Matlab内置函数,根据坐标生成几何体
	volume_ele(i) = volume(shp);               %Matlab内置函数,计算生成的几何体的体积
end

max_volume_ele = max(volume_ele);
min_volume_ele = min(volume_ele);
aveg_volume_ele = mean(volume_ele);

% Get G_NN,G_X_NODES,G_Y_NODES,G_X_Min,G_X_Max,G_Y_Min,G_Y_Max,G means Global.
G_NN      = zeros(8,Num_Elem);
G_X_NODES = zeros(8,Num_Elem);
G_Y_NODES = zeros(8,Num_Elem);
G_Z_NODES = zeros(8,Num_Elem);
G_X_Min   = zeros(1,Num_Elem);
G_X_Max   = zeros(1,Num_Elem);
G_Y_Min   = zeros(1,Num_Elem);
G_Y_Max   = zeros(1,Num_Elem);
G_Z_Min   = zeros(1,Num_Elem);
G_Z_Max   = zeros(1,Num_Elem);

G_NN = [Elem_Node(:,1)';Elem_Node(:,2)';Elem_Node(:,3)';Elem_Node(:,4)';
        Elem_Node(:,5)';Elem_Node(:,6)';Elem_Node(:,7)';Elem_Node(:,8)';];

G_X_NODES = [Node_Coor(G_NN(1,:),1)';Node_Coor(G_NN(2,:),1)';Node_Coor(G_NN(3,:),1)';Node_Coor(G_NN(4,:),1)';
             Node_Coor(G_NN(5,:),1)';Node_Coor(G_NN(6,:),1)';Node_Coor(G_NN(7,:),1)';Node_Coor(G_NN(8,:),1)'];
G_Y_NODES = [Node_Coor(G_NN(1,:),2)';Node_Coor(G_NN(2,:),2)';Node_Coor(G_NN(3,:),2)';Node_Coor(G_NN(4,:),2)';
             Node_Coor(G_NN(5,:),2)';Node_Coor(G_NN(6,:),2)';Node_Coor(G_NN(7,:),2)';Node_Coor(G_NN(8,:),2)'];
G_Z_NODES = [Node_Coor(G_NN(1,:),3)';Node_Coor(G_NN(2,:),3)';Node_Coor(G_NN(3,:),3)';Node_Coor(G_NN(4,:),3)';
             Node_Coor(G_NN(5,:),3)';Node_Coor(G_NN(6,:),3)';Node_Coor(G_NN(7,:),3)';Node_Coor(G_NN(8,:),3)'];
G_X_Min   = min(G_X_NODES,[],1);
G_X_Max   = max(G_X_NODES,[],1);
G_Y_Min   = min(G_Y_NODES,[],1);
G_Y_Max   = max(G_Y_NODES,[],1);
G_Z_Min   = min(G_Z_NODES,[],1);
G_Z_Max   = max(G_Z_NODES,[],1);