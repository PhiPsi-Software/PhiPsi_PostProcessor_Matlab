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

function Read_Geo
% This function read information about nodes and elements from *node and *elem files which can be obtained by ansys mac file, Ansys2PhiPsi_2D.mac.
% After the reading, the extreme values will be calculated and stored as global values.

global Full_Pathname Node_Coor Elem_Node Bou_x Bou_y Foc_x Foc_y
global Num_Node Num_Elem Num_Bou_x Num_Bou_y Num_Foc_x Num_Foc_y
global Min_X_Coor Max_X_Coor Min_Y_Coor Max_Y_Coor Outline Inline
global max_area_ele min_area_ele aveg_area_ele Centroid_elem
global Elemet_Info Node_Info 
global Yes_Checked Crack Elem_Material num_of_Material
global G_NN G_X_NODES G_Y_NODES G_X_Min G_X_Max G_Y_Min G_Y_Max
global Yes_has_FZ frac_zone_min_x frac_zone_max_x frac_zone_min_y frac_zone_max_y

% Read files
disp('    Reading node file....') 

% Check if Full_Pathname exist or not.
if size(Full_Pathname,2) ==1
    disp(['    Error :: Filename is not defined.'])
	%Error_Message
	warndlg('Filename is not defined!','ERROR')
end

% Check if *.node file exists or not.
Check_exist = [Full_Pathname,'.node'];
if exist(Check_exist,'file') ==0
    disp(['    Error :: Can not find input files.'])
	%Error_Message
	warndlg('Can not find input files!','ERROR')
end

Node_Coor = load([Full_Pathname,'.node']);

disp('    Reading elem file....') 
Elem_Info = load([Full_Pathname,'.elem']);
Elem_Node = Elem_Info(:,1:4);
Elem_Material = Elem_Info(:,5);       % Material number of elements
num_of_Material = max(Elem_Material); % Number of materials

if exist([Full_Pathname,'.boux'], 'file') ==2  
    disp('    Reading boux file....') 
    Bou_x = load([Full_Pathname,'.boux']);
end
if exist([Full_Pathname,'.bouy'], 'file') ==2  
    disp('    Reading bouy file....') 
    Bou_y = load([Full_Pathname,'.bouy']);
end
if exist([Full_Pathname,'.focx'], 'file') ==2  
    disp('    Reading focx file....') 
    Foc_x = load([Full_Pathname,'.focx']);
end
if exist([Full_Pathname,'.focy'], 'file') ==2  
    disp('    Reading focy file....') 
    Foc_y = load([Full_Pathname,'.focy']);
end
disp('    Geometry files reading done.') 

% Get the numbers of nodes, elements, boundaries and forces.
Num_Node = size(Node_Coor,1);
Num_Elem = size(Elem_Node,1);
Num_Bou_x= size(Bou_x,1);
Num_Bou_y= size(Bou_y,1);
Num_Foc_x= size(Foc_x,1);
Num_Foc_y= size(Foc_y,1);

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

% Check is initial crack out of range or not.
% for i=1:length(Crack)
% end

% Find outline of quadrilateral mesh.
tem = sort([Elem_Node(:,1) Elem_Node(:,2); Elem_Node(:,2) Elem_Node(:,3); 
          Elem_Node(:,3) Elem_Node(:,4); Elem_Node(:,4) Elem_Node(:,1) ]')';
[u,m,n] = unique(tem,'rows');      % determine uniqueness of edges
counts  = accumarray(n(:), 1);     % determine counts for each unique edge
Outline = u(counts==1,:);          % extract edges that only occurred once and store as a public value 

% Sort Outline by end to end.
[New_Outline] = Tools_Srot_by_End_to_End(Outline);

% Check if hole exists or not, if exists, then output "Inline".
if size(New_Outline,1) < size(Outline,1)
    % Hole exists.
	Inline=Outline;
	for ii=1:size(New_Outline,1)
	    case1 = [New_Outline(ii,1) New_Outline(ii,2)];
		case2 = [New_Outline(ii,2) New_Outline(ii,1)];
		for jj =1:size(Outline,1)
		    if case1 == Outline(jj,:) | case2 == Outline(jj,:)
			    Inline(jj,:)=[0 0];
			end
		end
	end
	Inline(all(Inline==0,2),:)=[]; 
	% Sort Inline by end to end.
	Tem_Inline = Inline(1,:);
	c_Inline   =1;
	for i = 1:size(Inline,1);
		if i==1
			c_Inline = 1;
			c_node_num_2 = Inline(1,2);
		end
		for j = 1:size(Inline,1)
			if j ~= c_Inline
				if c_node_num_2 == Inline(j,1)
					if ismember([Inline(j,1) Inline(j,2)],Tem_Inline,'rows') ==0
						Tem_Inline = [Tem_Inline; [Inline(j,1) Inline(j,2)]];
						c_Inline = j;
						c_node_num_2 = Inline(j,2);
						continue
					end
				elseif c_node_num_2 == Inline(j,2)
					if ismember([Inline(j,2) Inline(j,1)],Tem_Inline,'rows') ==0
						Tem_Inline = [Tem_Inline; [Inline(j,2) Inline(j,1)]];
						c_Inline = j;
						c_node_num_2 = Inline(j,1);
						continue
					end
				end
			end
		end
	end
	Outline = New_Outline;
	Inline  = Tem_Inline;
else
    % Hole does not exist.
    Outline = New_Outline;
	Inline  = [];
end

% Get the maximum & minimum & average area, centroid of all elements.
% Filter edge elements whose nodes belongs to "outline". 
Centroid_elem    = zeros(Num_Elem,2);
for i=1:Num_Elem
    N1  = Elem_Node(i,1);                                                  % Node 1 for current element
    N2  = Elem_Node(i,2);                                                  % Node 2 for current element
    N3  = Elem_Node(i,3);                                                  % Node 3 for current element
    N4  = Elem_Node(i,4);                                                  % Node 4 for current element
	NN = [N1 N2 N3 N4];                                                    % Four nodes of the current element  
	X_NODES = Node_Coor(NN,1);                                             % Coordinates of the four nodes of the current element
	Y_NODES = Node_Coor(NN,2);
	
	% Get the area of each element.
	area(i) = polyarea(X_NODES,Y_NODES);
	
	% Get the coordinates of centroid of each element.
	[center_x,center_y] = Cal_Centroid_of_Polygon(X_NODES,Y_NODES);
	Centroid_elem(i,1) = center_x;
	Centroid_elem(i,2) = center_y;
	
	% Filter edge elements whose nodes belongs to "outline". 
	Yes(1) = ismember([N1 N2],Outline,'rows');
	Yes(2) = ismember([N2 N3],Outline,'rows');
	Yes(3) = ismember([N3 N4],Outline,'rows');
	Yes(4) = ismember([N4 N1],Outline,'rows');
	Yes(5) = ismember([N2 N1],Outline,'rows');
	Yes(6) = ismember([N3 N2],Outline,'rows');
	Yes(7) = ismember([N4 N3],Outline,'rows');
	Yes(8) = ismember([N1 N4],Outline,'rows');

	if max(Yes)==1
	    Elemet_Info(i,2) = 1;                                              % Edge element
	end
end

max_area_ele = max(area);
min_area_ele = min(area);
aveg_area_ele = mean(area);

% Get G_NN,G_X_NODES,G_Y_NODES,G_X_Min,G_X_Max,G_Y_Min,G_Y_Max,G means Global.
G_NN      = zeros(4,Num_Elem);
G_X_NODES = zeros(4,Num_Elem);
G_Y_NODES = zeros(4,Num_Elem);
G_X_Min   = zeros(1,Num_Elem);
G_X_Max   = zeros(1,Num_Elem);
G_Y_Min   = zeros(1,Num_Elem);
G_Y_Max   = zeros(1,Num_Elem);

G_NN = [Elem_Node(:,1)' ;Elem_Node(:,2)' ;Elem_Node(:,3)' ;Elem_Node(:,4)'];

G_X_NODES = [Node_Coor(G_NN(1,:),1)';Node_Coor(G_NN(2,:),1)';Node_Coor(G_NN(3,:),1)';Node_Coor(G_NN(4,:),1)'];
G_Y_NODES = [Node_Coor(G_NN(1,:),2)';Node_Coor(G_NN(2,:),2)';Node_Coor(G_NN(3,:),2)';Node_Coor(G_NN(4,:),2)'];

G_X_Min   = min(G_X_NODES,[],1);
G_X_Max   = max(G_X_NODES,[],1);
G_Y_Min   = min(G_Y_NODES,[],1);
G_Y_Max   = max(G_Y_NODES,[],1);

% 读取破裂区文件
if exist([Full_Pathname,'.fraz'], 'file') ==2
    Yes_has_FZ = 1;
	disp('    > Reading fraz file....') 
	namefile= [Full_Pathname,'.fraz'];
	data=fopen(namefile,'r'); 
	lineNum =0;
	while ~feof(data)
		lineNum = lineNum+1;
		TemData = fgetl(data);    
		if lineNum>=2   %第一行是文件标识行,不予读取   
			ttt_DATA(1:4)  = str2num(TemData);
		end
	end
	fclose(data); 
	frac_zone_min_x = ttt_DATA(1); 
	frac_zone_max_x = ttt_DATA(2); 
	frac_zone_min_y = ttt_DATA(3); 
	frac_zone_max_y = ttt_DATA(4); 
end
