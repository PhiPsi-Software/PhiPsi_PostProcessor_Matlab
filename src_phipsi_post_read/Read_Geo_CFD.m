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

function Read_Geo_CFD
% Read MESH FILE for CFD Analysis.

global Full_Pathname Node_Coor Mesh_Connects Bou_x Bou_y Foc_x Foc_y
global Num_Node Num_Elem Num_Bou_x Num_Bou_y Num_Foc_x Num_Foc_y
global Min_X_Coor Max_X_Coor Min_Y_Coor Max_Y_Coor Outline Inline
global max_area_ele min_area_ele aveg_area_ele Centroid_elem
global Elemet_Info Node_Info
global Yes_Checked Crack Elem_Material num_of_Material
global G_NN G_X_NODES G_Y_NODES G_X_Min G_X_Max G_Y_Min G_Y_Max
global Yes_has_FZ frac_zone_min_x frac_zone_max_x frac_zone_min_y frac_zone_max_y
global Mesh_Points number_mesh_points number_mesh_connect Mesh_Connects

% Read files
disp('    Reading mesh point file....')

% Check if Full_Pathname exist or not.
if size(Full_Pathname,2) ==1
    disp(['    Error :: Filename is not defined.'])
	%Error_Message
	warndlg('Filename is not defined!','ERROR')
end

% Check if *.mesh file exists or not.
Check_exist = [Full_Pathname,'.mesh'];
if exist(Check_exist,'file') ==0
    disp(['    Error :: Can not find input mesh files.'])
	%Error_Message
	warndlg('Can not find input mesh files!','ERROR')
else
    %第一步:读取第一行的数据,分别为网格点数和连接单元数
    fpn=fopen(Check_exist,'rt');
	c_file_line =fgetl(fpn);         %读取第一行
	c = textscan(c_file_line,'%6.1f');
	number_mesh_points  = c{1,1}(1);
	number_mesh_connect = c{1,1}(2);
	fclose(fpn);
	
	%然后先把第一行去掉,分别存入新文件mesp,mesc,
    new_tem_file1 =[Full_Pathname,'.mesp'];  %点的坐标和类型文件
	new_tem_file2 =[Full_Pathname,'.mesc'];  %点连接文件
    fid1=fopen(new_tem_file1,'wt');
	fid2=fopen(new_tem_file2,'wt');
    fpn=fopen(Check_exist,'rt');
	cc_count =0;
    while feof(fpn)~=1
	   cc_count = cc_count +1;
       c_file_line =fgetl(fpn); 
       if cc_count>=2 && cc_count <= number_mesh_points+1
           fprintf(fid1,'%s\n',c_file_line);
       end	   
       if cc_count> number_mesh_points+1
           fprintf(fid2,'%s\n',c_file_line);
       end
    end
    fclose(fid1);
	fclose(fid2);
	fclose(fpn);

    Mesh_Points   = load([Full_Pathname,'.mesp']);
    Mesh_Connects = load([Full_Pathname,'.mesc']);
end


% Get the max and min value of points coordinates.
Max_X_Coor = max(max(Mesh_Points(:,2)));
Min_X_Coor = min(min(Mesh_Points(:,2)));
Max_Y_Coor = max(max(Mesh_Points(:,3)));
Min_Y_Coor = min(min(Mesh_Points(:,3)));

% -------------------------------------
% Find outline of quadrilateral mesh.
% -------------------------------------
tem = sort([Mesh_Connects(:,2) Mesh_Connects(:,3); Mesh_Connects(:,3) Mesh_Connects(:,4); Mesh_Connects(:,4) Mesh_Connects(:,2)]')';
[u,m,n] = unique(tem,'rows');      % determine uniqueness of edges
counts  = accumarray(n(:), 1);     % determine counts for each unique edge
Outline = u(counts==1,:);          % extract edges that only occurred once and store as a public value 
% Sort Outline by end to end.
[New_Outline] = Tools_Srot_by_End_to_End(Outline);
% 检查是否存在孔洞.
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


