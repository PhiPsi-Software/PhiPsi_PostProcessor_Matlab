## PhiPsi Post-Processor v1.18.2 (2021-08-07) new features:
1. Added support for 3D XFEM refinement of enriched elements.
   
## PhiPsi Post-Processor v1.18.1 (2021-08-03) new features:
1. Function enhancement of Plot_Crack_Contour_3D.m.

## PhiPsi Post-Processor v1.18.0 (2021-08-01) new features:
1. Added a figure control widget fcw.m and camview.m. Controlled by Key_Figure_Control_Widget in PhiPsi_Color_and_Font_Settings.m.
   Ref: https://ww2.mathworks.cn/matlabcentral/fileexchange/38019-figure-control-widget.
   Matlab code to active the figure control widget:
		% Active Figure control widget (2021-08-01)
		% Press q to exit.
		% Press r (or double-click) to reset to the initial.
		fcw(gca); 
2. Matlab codes are placed into individual folder. 
		% Add path of source files.
		addpath('src_fcw')
		addpath('src_geom3d')
		addpath('src_meshes3d')
		addpath('src_phipsi_post_animate')
		addpath('src_phipsi_post_cal')
		addpath('src_phipsi_post_main')
		addpath('src_phipsi_post_plot')
		addpath('src_phipsi_post_read')
		addpath('src_phipsi_post_tool')		
3. Performance improvement for the 3D drawing. fill3 is replaced by patch.
		fill3(c_x_1(1:4,1:count_1),c_y_1(1:4,1:count_1),c_z_1(1:4,1:count_1),Color_3D_ele_face,'FaceAlpha',FaceAlpha_3D_ele_face) 
   has been changed to:
        patch(c_x_1(1:4,1:count_1),c_y_1(1:4,1:count_1),c_z_1(1:4,1:count_1),Color_3D_ele_face,'FaceAlpha',FaceAlpha_3D_ele_face) 
4. Performance improvement for the 3D drawing. Reduce the number of objects to be drawn when using plot3. For example, the following codes:
		%绘制单元面
		nele = size(Crack_Ele_1{i},2);
		for j=1:nele
			c_x = [Crack_node_X_new{i}(Crack_Ele_1{i}(j)),Crack_node_X_new{i}(Crack_Ele_2{i}(j)),Crack_node_X_new{i}(Crack_Ele_3{i}(j))];
			c_y = [Crack_node_Y_new{i}(Crack_Ele_1{i}(j)),Crack_node_Y_new{i}(Crack_Ele_2{i}(j)),Crack_node_Y_new{i}(Crack_Ele_3{i}(j))];
			c_z = [Crack_node_Z_new{i}(Crack_Ele_1{i}(j)),Crack_node_Z_new{i}(Crack_Ele_2{i}(j)),Crack_node_Z_new{i}(Crack_Ele_3{i}(j))];
			patch(c_x,c_y,c_z,[1,1,0],'FaceAlpha',0.5,'FaceLighting','gouraud')					
		end	
   have been changed to:
		%绘制单元面
		nele = size(Crack_Ele_1{i},2);
		c_x =[];c_y =[];c_z =[];
		for j=1:nele
			c_x(j,1:3) = [Crack_node_X_new{i}(Crack_Ele_1{i}(j)),Crack_node_X_new{i}(Crack_Ele_2{i}(j)),Crack_node_X_new{i}(Crack_Ele_3{i}(j))];
			c_y(j,1:3) = [Crack_node_Y_new{i}(Crack_Ele_1{i}(j)),Crack_node_Y_new{i}(Crack_Ele_2{i}(j)),Crack_node_Y_new{i}(Crack_Ele_3{i}(j))];
			c_z(j,1:3) = [Crack_node_Z_new{i}(Crack_Ele_1{i}(j)),Crack_node_Z_new{i}(Crack_Ele_2{i}(j)),Crack_node_Z_new{i}(Crack_Ele_3{i}(j))];					
		end		
		patch(c_x',c_y',c_z',[1,1,0],'FaceAlpha',0.5,'FaceLighting','gouraud')		

## PhiPsi Post-Processor v1.17.0 (2021-07-27) new features:
1. 解决Colorbar无法显示极值的问题:
		 t1=caxis;
		 t1=linspace(t1(1),t1(2),8);
		 my_handle=colorbar('ytick',t1);
2. PhiPsi_Color_and_Font_Settings.m中增加全局变量Num_Colorbar_Level.
     Num_Colorbar_Level  = 8; % Number of colorbar level of values, default: 8, 2021-07-27		 
3. PhiPsi_Color_and_Font_Settings.m中增加全局变量Colorbar_Font、Title_Font.

## PhiPsi Post-Processor v1.16.0 (2021-06-30) new features:
1. 2D新增局部加密支持.

## PhiPsi Post-Processor v1.15.8 (2021-02-11) new features:
1. 3D绘图功能增加：增加了Plot_Deformation3D.m绘制裂尖增强节点对应的单元, Key_PLOT(2,10)==4时可用.
2. 3D绘图功能增加：增加了Plot_Deformation3D.m绘制三维裂缝体的功能, Key_PLOT(2,5)==2时可用.

## PhiPsi Post-Processor v1.15.7 (2021-02-07) new features:
1. 3D绘图功能增加：增加了Plot_Deformation3D.m绘制裂尖加权应力球的功能, Key_PLOT(2,12)==1或Key_PLOT(2,12)==2时可用.

## PhiPsi Post-Processor v1.15.6 (2021-01-21) new features:
1. 3D绘图bug修复：修复了无法绘制3D载荷箭头的问题,Plot_Deformation3D.m.

## PhiPsi Post-Processor v1.15.5 (2020-12-27) new features:
1. 3D功能增强：位移云图、应力云图不绘制增强单元，从而显示出裂缝面位置.
    Key_PLOT(3,5)=3
	  Key_PLOT(4,5)=3

## PhiPsi Post-Processor v1.15.4 (2020-08-10) new features:
1. 增加椭圆空缺支持.
2. 增加压应力为正云图绘制控制：Key_Posivite_Compressive = 1.

## PhiPsi Post-Processor v1.15.3 (2020-03-12) new features:
1. 新增3D裂缝面云图绘制功能(新增Key_PLOT(5,:)控制变量).

## PhiPsi Post-Processor v1.15.2 (2020-03-11) new features:
1. 修复了调用Tools_mArrow3.m绘制3D箭头时的错误.
2. 新增3D节点位移云图绘制和3D节点应力云图绘制.

## PhiPsi Post-Processor v1.15.1 (2020-02-26) new features:
1. 修复了变形很小的情况下2D位移云图、应力云图显示空白的Bug.

## PhiPsi Post-Processor v1.15.0 (2020-01-05) new features:
1. 3D功能添加.

## PhiPsi Post-Processor v1.14.3 (2019-08-10) new features:
1. 3D开度云图初步.
2. 通过删除重复部分大幅提高了3D网格和单元的绘制速度.
3. 新增3D裂缝开度矢量绘制(含Slice切片中的开度矢量).
4. 新增3D切片(Slice)绘制功能.
5. 支持生死单元后处理.
6. 新增流体单元绘制.

## PhiPsi Post-Processor v1.13.1 (2019-08-01) new features:
1. Post-Process of excavation.

## PhiPsi Post-Processor v1.10.0 (2019-04-16) new features:
1. Post-Process of 3D crack propagation.

## PhiPsi Post-Processor v1.9.1 (2018-09-26) new features:
1. Bug fixes.

## PhiPsi Post-Processor v1.9.0 (2018-09-10) new features:
1. Support field problems with cracks and holes.

## PhiPsi Post-Processor v1.8.4 (2018-08-28) new features:
1. Bug fixes.

## PhiPsi Post-Processor v1.8.3 (2018-06-19) new features:
1. Bug fixes.

## PhiPsi Post-Processor v1.8.1 (2018-01-11) new features:
1. Speed improvement of binary files reading.

## PhiPsi Post-Processor v1.8.0 (2017-10-31) new features:
1. Add supports for the reading of binary result files generated by PhiPsi.
2. Add support for molecular dynamics analysis.
3. Replace PhiPsi_Input_Control.m with PhiPsi_Color_and_Font_Settings.m.
4. Plot contour of displacement, i.e., sqrt(disx^2+disy^2).
5. Add support for the curved cracks.
6. Bug fixes.

## PhiPsi Post-Processor v1.4.2 (2017-05-21) new features:
1. Add support for the intersection of crack and void.
2. Plot the load-displacement curves.
3. Add support for the crossing cracks.
4. Add support for the cohesive cracks.
5. Bug fixes.

## PhiPsi Post-Processor v1.2.3 (2017-04-16) new features:
1. Add support for the post-processing of static field problems.
2. Add support for the post-processing of molecular dynamics analysis.
3. Plot contour of the equivalent plastic strain.
4. Bug fixes.

## PhiPsi Post-Processor v1.2.2 (2017-02-04) new features:
1. 绘图默认打开openGL渲染器.
2. 增加等效塑性应变的绘制.
3. 增加分子动力学模拟后处理.

## PhiPsi Post-Processor v1.1.8 (2016-11-21) new features:
1. 新增某点压力曲线的绘制.
2. 新增页岩气产量曲线绘制.
3. 场问题流量矢量动画绘制.
4. 场问题动画绘制.
5. 场问题边界条件的绘制.

## PhiPsi Post-Processor v1.1.7 (2016-11-01) new features:
1. Plot 3D mesh.
2. Plot inclusions.
3. Bug fixes.

## PhiPsi Post-Processor v1.1.5.3 (2016-10-12) new features:
1. 新增多边形夹杂绘制.
2. 优化.

## PhiPsi Post-Processor v1.1.4.3 (2016-10-3) new features:
1. 圆形夹杂相关优化.
2. 新增夹杂之增强节点的绘制.
3. 新增圆形夹杂的绘制.

## PhiPsi Post-Processor v1.1.3.5 (2016-09-28) new features:
1. 修复了变形图和非变形图一起画时出错的bug.
2. -999表示智能处理最后一个计算步的结果.
3. 增加了场问题场值的云图绘制.

## PhiPsi Post-Processor v1.1.0.0 (2016-07-11) new features:
1. Junction增强修复.

## PhiPsi Post-Processor v1.0.0.68 (2016-06-16) new features:
1. 新增裂缝切向开度绘制功能.
2. bug修复.
3. 优化.

## PhiPsi Post-Processor v1.0.0.63 (2016-05-15) new features:
1. 支持XFEM孔洞的后处理.
2. 绘制破裂区边界线.

## PhiPsi Post-Processor v1.0.0.48 (2016-03-13) new features:
1. 增加位移动画和应力动画的绘制功能.
2. 新增变形动画绘制功能.
3. 新增刚体小球的绘制.
4. 绘制注水点水压变化曲线.
5. 绘制应力强度因子各破裂步变化曲线.
6. 新增裂缝编号绘制功能.
7. 网格变形图支持绘制闭合单元(粉色).

## PhiPsi Post-Processor v1.0.0.30 (2016-02-05) new features:
1. 新增水压线绘制功能.
2. 增加天然裂缝的绘制功能.
3. 增加了KGD模型理论解.
4. 绘制支撑剂浓度曲线.
5. 绘制流量曲线.
6. 新增绘制流速曲线功能.
7. 新增支撑剂颗粒的绘制功能.

## PhiPsi Post-Processor v1.0.0.18 (2015-10-14) new features
1. 裂纹面接触状态分两种,接触和支撑剂支撑,用两种不同颜色标记接触单元.
2. 修复了Gauss点位移云图绘制的一个小bug.
3. 新增3D网格绘制.
4. 计算裂纹变形形态时包含裂纹片段端点.
5. 绘制各个单元的接触状态.
6. 支持主应力云图和方向的绘制.
7. 支持Gauss点位移云图的绘制.
8. 支持Gauss点应力云图的绘制.
9. 区分不同的裂尖增强方案.
10. 新增计算点编号绘制功能(Key_PLOT(1,6)=1).
11. 新增HF计算点的绘制(Key_PLOT(1,5)=1)
