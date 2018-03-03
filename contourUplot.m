% April 12, 2011: single_system_simulation_v2:  Add APD vs. next DI
% diagnostic.

% April 19, 2011: simple_plot_v2_0: Just plot the membrane potential.

% Nov 16, 2015: simple_plot_v2_0s: Simplify a little.

% Parameters:

% The following should be copied from the corresponding parameters in the simulation:
Nx = 144; Ny = 144; %72% Number of gridpoints in the x and y directions
Dt = 0.0016; % Timestep size (in ms)
Napd = 240; Ndi = 300; Nvabs = 100; % Dimensions used to make phase space plots
Dx = 0.1667; Dy = 0.1667; %*2 % Grid spacing (in cm)
ntplot = 50; % Plots are made every ntplot timesteps in the simulation
x0 =0; y0=-7;obs2_radius = 0.99;
% User may choose the following:

% Folder containing the data files:
data_folder_name = 'data_1580_1790';

% Files to process in the data folder:
% Files uxxxxx through uyyyyy where xxxxx = timestep_start,
% yyyyy = timestep_end, with increment timestep_inc:
timestep_start =0 ; timestep_end =20000; timestep_inc = 50;
uhist = zeros(1,timestep_end);
vhist = zeros(1,timestep_end);

% To make a movie file set equal to 1; to see on the screen (faster) set
% equal to 0.  Needs to be updated to use videowriter when equal to 1:
make_movie = 0;

xx = ((0:(Nx-1))-0.5*Nx)*Dx;
yy = ((0:(Ny-1))-0.5*Ny)*Dy;
if (make_movie)
    vidobj = VideoWriter('G2.avi');
    open(vidobj);

    
    fig = figure('Position',[37 824 542 477],'visible','off'); clf;
else
    fig=figure(1);
    axis equal;
    %fig = figure('Position',[37 824 542 477],'visible','on'); clf;
end
set(gca,'FontSize',20);


obs = zeros(Nx,Ny);
for i = 1:Nx
    for j = 1:Ny
        x = (i-1-Nx/2)*Dx;
        y = (j-1-Ny/2)*Dy;
        if (x-x0)*(x-x0)+(y-y0)*(y-y0)< obs2_radius*obs2_radius
            obs(i,j) = 1;
        end
    end
end


% Loop over plotting times:
for it = timestep_start:timestep_inc:timestep_end
    
    % Membrane potential snapshot:;
    fp = fopen(sprintf('%s/u%07i',data_folder_name,it),'rb');
    if (fp==-1)
        fprintf('Can''t find file %s/u%07i\n',subfolder,it); break;
    end;
    u = fread(fp,[Nx,Ny],'double');
    if (it>0)
        uhist(it)=u(1,1);
    end
    
    fclose(fp);
    
    fp = fopen(sprintf('%s/v%07i',data_folder_name,it),'rb');
    if (fp==-1)
        fprintf('Can''t find file %s/v%07i\n',subfolder,it); break;
    end;
    v = fread(fp,[Nx,Ny],'double');
    
    
    fclose(fp);
    
    clf;
    whitebg('w');
    
     umin = 0.0;
     umax = 1.0;
     umid = 0.5;
     vmin = -0.075;
     vmid = 0.5*(0.6)-0.075;
     vmax = 1.0-5*0.075;
     
     red= (u-umin)/(umax-umin);
     blue= (v-vmin)/(vmax-vmin);
     green= 0.5*(1+tanh(10*(vmid-v))).*exp(-20*(u-umid).*(u-umid));
%     alpha = 0.1;
%     a = 0.6;
%     b = 0.075; % a,b need to be consistent with the corresproning values
%     red = (u > alpha);
%     green = (v > (a/2-b))*0.5;
%     blue = 0;
    

    C = zeros(Nx,Ny,3);
    C(:,:,1) = red';
    C(:,:,2) = green';
    C(:,:,3) = blue';
    
    C(:,:,1) = max(C(:,:,1),obs');
    C(:,:,2) = max(C(:,:,2),obs');
    C(:,:,3) = max(C(:,:,3),obs');
    
    
    image(xx,yy,C);
    
    %imagesc(x,y,v');
    axis xy; %caxis([0 1.4]); 
    %colormap(jet); shading interp;
    xlabel('x'); ylabel('y');
    axis square;
    %h=colorbar;
    %caxis([0 0.3]);
    %hold on;
    %contour(x,y,u',[0.2,0.2],'w');
  
    
    %set(h,'fontsize',20);
    title(sprintf('Membrane potential at timestep = %i',it));
    
    drawnow;
    %pause;
    
    
    hold off;
    
    if (make_movie)
        writeVideo(vidobj,getframe(gcf));
    end
    
end


close(vidobj);
%if (make_movie)
%    mov = close(mov);
%end
%figure(2);
%plot(uhist,'o');
