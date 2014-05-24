function create_retro_flow_image_collagen(rf_path,im,bin_im,index,x_flow,y_flow,x_shift,y_shift,mue2pix_ratio,fps,max_r_flow,flow_field_arrow_distance, files, varargin)

nVarargs = length(varargin)
if nVarargs > 0 
    lunit_str = varargin{1};
    tunit_str = varargin{2};
else
    lunit_str = 'µm';
    tunit_str = 'minute';
end

switch lunit_str
    case 'm'
        length_scale_unit = 1e6;
    case 'mm'
        length_scale_unit = 1e3;
    case 'µm'
        length_scale_unit = 1;
    case 'nm'
        length_scale_unit = 1e-3;    
    otherwise
        assert(false,'unknown lenghtscale')
end




switch tunit_str
    case 'day'
        time_scale_unit = 1440;
    case 'hour'
        time_scale_unit = 60;
    case 'minute'
        time_scale_unit = 1;
    case 'second'
        time_scale_unit = 1/60;
    otherwise
        assert(false,['unknown timescale:' tunit_str])
end




hold off;
tic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Define parameters (are saved in analysis_parameters.mat):
%flow_field_arrow_distance=25;       %Flow field arrow distance
%max_r_flow=5;                       %max retrograde flow amplitude for colorcoding
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%first run: save the parameters defined above
if index==1
    save([rf_path,'/analysis_parameters.mat'],'flow_field_arrow_distance','max_r_flow','-append');
end

%create a new figure to display the flow fiels in.
%h=fig_handle;
h=figure;
%axes(h);

%get the images
im1=im;
[y_max,x_max]=size(im1);
%now take care of the x_shift, y_shift;
im_int(1:y_shift+y_max,x_shift+x_max)=0;
im_int(y_shift+1:y_shift+y_max,x_shift+1:x_shift+x_max)=im1;
im1=im_int;

%This is already prepared for Dannys Labview version, with the binary shapes stored in rf_path/binary_shape.
%If the folder exists, read from there, otherwise do it on your own
% if isdir([rf_path,filesep,'binary_shape']
%     files_bin=dir([rf_path,filesep,'binary_shape',filesep,'corr*']);
%     r_im_erode=double(imread([rf_path,filesep,'binary_shape',filesep,files_bin(index).name]));
%     r_im_erode(find(r_im_erode==255))=1;
%     r_im_erode_int(1:y_shift+y_max,x_shift+x_max)=0;
%     r_im_erode_int(y_shift+1:y_shift+y_max,x_shift+1:x_shift+x_max)=r_im_erode;
%     r_im_erode=r_im_erode_int;
% else
%     r_im_erode=better_edge(im1);    %edge detection function written by Danny.
% end

%again, take care of the image shift
im_int(1:y_shift+y_max,x_shift+x_max)=0;
im_int(y_shift+1:y_shift+y_max,x_shift+1:x_shift+x_max)=bin_im;
r_im_erode=im_int;

[y_max,x_max]=size(im1);
max_size=max(size(im1));

x_flow=double(x_flow);
y_flow=double(y_flow);
[y_max_field,x_max_field]=size(x_flow);

%Make sure im1, and the flow data have the same size, if not reshape the flow data accordingly
if (~(size(x_flow)==size(im1)))
    intermediate(1:y_max,1:x_max)=NaN;
    intermediate(1:min([y_max,y_max_field]),1:min([x_max,x_max_field]))=x_flow(1:min([y_max,y_max_field]),1:min([x_max,x_max_field]));
    x_flow=intermediate;
    intermediate(1:y_max,1:x_max)=NaN;
    intermediate(1:min([y_max,y_max_field]),1:min([x_max,x_max_field]))=y_flow(1:min([y_max,y_max_field]),1:min([x_max,x_max_field]));
    y_flow=intermediate;
end
clear('intermediate');

%invert_im1
im1_inv=double(im1);
im1_inv=255-im1_inv;
im1_inv=uint8(im1_inv);

%use the shape, to cut of whatever of the retro flow maps is outside the real shape
x_flow=r_im_erode.*x_flow;
y_flow=r_im_erode.*y_flow;
x_flow(isnan(x_flow))=0;
y_flow(isnan(y_flow))=0;
x_flow_i=x_flow(x_shift+1:end,y_shift+1:end);
y_flow_i=y_flow(x_shift+1:end,y_shift+1:end);
x_flow=x_flow_i;y_flow=y_flow_i;
[th_flow,r_flow]=cart2pol(x_flow,y_flow);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Now we should have everything together to start the drawing the image
%First make the color coded growth retrograde flow
r_flow(isnan(r_flow))=0;
%r_flow=r_flow(x_shift+1:end,y_shift+1:end);
imshow(r_flow/length_scale_unit*time_scale_unit)
colormap('gray')
set(gca,'CLim',[0,max_r_flow/length_scale_unit*time_scale_unit])
display(max_r_flow/length_scale_unit*time_scale_unit)


mkdir([rf_path,filesep,'gray-scale']);
mkdir([rf_path,filesep,'color_code']);
mkdir([rf_path,filesep,'color_code_and_arrows']);


out0=[rf_path,filesep,'gray-scale',filesep, files(index).name(1:end-4),'_','grayscale.png']
print(h,'-dpng','-r150',out0);

%here I will write the colorbar
if index==1
    fig1=figure;
    left=100; bottom=100 ; width=20 ; height=500;
    pos=[left bottom width height];
    axis off
    set(gca,'CLim',[0,max_r_flow/length_scale_unit*time_scale_unit])
    cb=colorbar
    ylabel(cb,['speed ' lunit_str '/' tunit_str])
    colormap gray
    set(fig1,'OuterPosition',pos)
    color_bar=[rf_path,filesep,'gray-scale',filesep, 'colorbar_grayscale.pdf']
    print(fig1,'-dpdf','-r150',color_bar);
    close(fig1)
end
    


colormap('jet')

set(gca,'CLim',[0,max_r_flow/length_scale_unit*time_scale_unit])
display(max_r_flow/length_scale_unit*time_scale_unit)
hold on

%now get and draw the flow field arrows
plot_data=generate_plot_normalized_collagen(x_flow,y_flow,flow_field_arrow_distance);
hold on


out=[rf_path,filesep,'color_code',filesep, files(index).name(1:end-4),'_','colorcode.png'];
print(h,'-dpng','-r150',out);

%here I will write the colorbar
if index==1
    fig1=figure;
    left=100; bottom=100 ; width=20 ; height=500;
    pos=[left bottom width height];
    axis off
    set(gca,'CLim',[0,max_r_flow/length_scale_unit*time_scale_unit])
    cb=colorbar
    ylabel(cb,['speed ' lunit_str '/' tunit_str])
    colormap jet
    set(fig1,'OuterPosition',pos)
    color_bar=[rf_path,filesep,'color_code',filesep, 'colorbar.pdf']
    print(fig1,'-dpdf','-r150',color_bar);
    close(fig1)
end
quiver(plot_data(:,1),plot_data(:,2),plot_data(:,3),plot_data(:,4),0.3,'k');
%cb = colorbar();


out2=[rf_path,'/color_code_and_arrows',filesep, files(index).name(1:end-4),'_','colorcode_arrows.png']
print(h,'-dpng','-r150',out2);
if index==1
    fig1=figure;
    left=100; bottom=100 ; width=20 ; height=500;
    pos=[left bottom width height];
    axis off
    set(gca,'CLim',[0,max_r_flow/length_scale_unit*time_scale_unit])
    colorbar
    colormap jet
    set(fig1,'OuterPosition',pos)
    color_bar=[rf_path,filesep,'color_code_and_arrows',filesep, 'colorbar.pdf']
    print(fig1,'-dpdf','-r150',color_bar);
    close(fig1)
end

%Draw a 10µm scale bar
fill([x_max-x_shift-5-10/mue2pix_ratio,x_max-x_shift-5-10/mue2pix_ratio,x_max-5-x_shift,x_max-5],[5,15,15,5],'w')
%text(x_max-5-10/mue2pix_ratio,20,'10 µm','Color','w')


%Save the images, the display_edge_array, and of course the final_results_out

%now correct for the stupid matlab exporting and resizeing
for cout = {out2 out out0}
    im_corr = imread(cout{1});
    im_bw   = mean(im_corr,3);
    x_mean=mean(im_bw,2);
    y_mean=mean(im_bw,1);
    x_im_f=find(~(x_mean==255));
    y_im_f=find(~(y_mean==255));
    im_corr_cut = im_corr(x_im_f(1):x_im_f(length(x_im_f)),y_im_f(1):y_im_f(length(y_im_f)),:);
    [y_im_corr_cut,x_im_corr_cut,depth]=size(im_corr_cut);
    im1_resized=imresize(im1_inv,[y_im_corr_cut,x_im_corr_cut]);
    im_corr_resized=imresize(im_corr_cut,size(im));
    imwrite(im_corr_resized,cout{1});
end
%finally get an overlay with the original image, to give the structural information!!
%im1_inv_d(:,:,1)=double(im1_inv);im1_inv_d(:,:,2)=double(im1_inv);im1_inv_d(:,:,3)=double(im1_inv);
%im1_inv_resize=im1_inv_d(x_shift+1:end,y_shift+1:end,:);
%im_corr_cut_im1=uint8(ceil(double(im_corr_resized).*im1_inv_resize./255));
%imwrite(im_corr_cut_im1,[rf_path,filesep,'interpolation_image',filesep,'plus_actin',filesep,'all_overlay',num2str(index+1000),'all_filt.png']);
data_out=[rf_path,'/interpolation_data',filesep, files(index).name(1:end-4),'_','colorcode_arrows.png']
save([rf_path,filesep,'interpolation_data',filesep, files(index).name(1:end-4),'_','interpolated_flow.mat'],'x_flow','y_flow','mue2pix_ratio','fps');

close(h)
toc
pause(0.1)



