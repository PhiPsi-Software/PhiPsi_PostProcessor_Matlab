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

% This function deletes all the result files in specified folder.

% Clear and close.
clear all; close all; clc; format compact;  format long;
disp(['  ======================================================================'])
disp(['  ****                                                              ****'])
disp(['  ****                FraxFem2D Results File Deleter                ****']) 
disp(['  ****                                                              ****'])
disp(['  ======================================================================']) 
disp([' > Features:                                                            ']) 
disp(['   This tool is used to delete results files produced by FraxFem2D.     ']) 
disp(['   This program is written by Matlab.                                   '])
disp(['   -------------------------------------------------------------------- ']) 
disp([' > Author:  SHI Fang, University of Science & Technology of China       ']) 
disp([' > Website: http://www.FraxFem2D.com                                    ']) 
disp([' > Email:   FraxFem2D@sina.com                                          ']) 
disp(['  ======================================================================']) 
disp(['  '])
% Get folder name.
% ----------------------
%      Option 1
% ----------------------
% name_folder ='C:\Matlab work\FraxFEM work\';
% if exist(name_folder,'dir') ==7
	% str = sprintf('    %s is found!', name_folder);
	% disp(str);
% end

% ----------------------
%      Option 2
% ----------------------
% Read folder name form user's input. 
for i=1:20     % Try 20 times to read the right the file name
	name_folder = input(' >> Please input the folder name:\n', 's');
	
	if name_folder(length(name_folder)) ~='\'
	   name_folder = [name_folder '\'];
	end
	
	if exist(name_folder,'dir') ==0
	    disp('    Error :: Input folder name is not found!')
		if i ~=20
	       continue
		elseif i==20
		   disp('    Error :: Input folder name is not found after 20 attempts!')
		   Error_Message
		end
	elseif exist(name_folder,'dir') ==7
		str = sprintf('    %s is found!', name_folder);
		disp(str);
		break
	end
end

disp(['  '])
% Check again.
disp([' >> Warning: please make sure the folder name is right!'])
Make_sure = input(' >> Press Enter if the folder name is right, type NO to cancel:\n', 's');
if isempty(Make_sure)==1
    Yes_Delete = 1;
elseif lower(Make_sure) =='no'
    Yes_Delete = 0;
end

% Delete now.
if Yes_Delete==1
	disp(['  '])
	disp([' >> Now begin to delete all results files......'])

	% Read all files in the folder.
	FileList=dir(name_folder); 

	i_folder=0;
	% Loop through each file.
	for i=1:length(FileList)
		% Check if the file is a folder.
		if(FileList(i).isdir==1&&~strcmp(FileList(i).name,'.')&&~strcmp(FileList(i).name,'..'))
			i_folder = i_folder+1;
			 % Store folder name.
			fileFolder{i_folder}=[FileList(i).name];  
		end
	end

	% Loop through each folder.
	Size_of_deleted_files =0;
	count_f =0;
	for j=1:length(fileFolder)
		c_folder = [name_folder fileFolder{j}];
		c_FileList=dir(c_folder);
		% Loop through each file.
		for k=1:length(c_FileList)
			filename = c_FileList(k).name;
			% extension name.
			if length(filename) >2
				full_path_name = [c_folder '\' filename];
				ext_name = filename(strfind(filename,'.'):end);
				
				if length(ext_name) <=4   % *.png and *.gif
					delete(full_path_name);
					Size_of_deleted_files =Size_of_deleted_files+c_FileList(k).bytes;
					disp(['    ',filename,' deleted.']);
					count_f =count_f + 1;
				elseif length(ext_name) >4 % FraxFEM files.
					l_ext_name = lower(ext_name(2:5));
					FraxFEM_ext = {'apdl';'boux';'bouy';'elem';'node';'focx';'focy';'ivex';'ivey'};	
					if ismember(l_ext_name,FraxFEM_ext)==0
						Size_of_deleted_files =Size_of_deleted_files+c_FileList(k).bytes;
						delete(full_path_name);
						disp(['    ',filename,' deleted.']);
						count_f =count_f + 1;
					end				
				end
				
			end
		end
	end

	Size_of_deleted_files=Size_of_deleted_files/1024/1024;
	disp([' >> ',num2str(count_f ),' files were deleted.']);
	disp([' >> Size of all deleted result files is ',num2str(Size_of_deleted_files),' MB']);
end

Exit = input(' >> Press any key to exit!\n', 's');