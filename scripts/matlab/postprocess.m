function postprocess(foldername, xmin, xmax, ymin, ymax, omgmax, domg, pmax, dp)

% default values --------------------------------------------------------------
if nargin == 7
	xmin = str2num(xmin);
	xmax = str2num(xmax);
	ymin = str2num(ymin);
	ymax = str2num(ymax);
	omgmax = str2num(omgmax);
	domg = str2num(domg);
	pmax = 1.0;
	dp = 0.1;
elseif nargin == 5
	xmin = str2num(xmin);
	xmax = str2num(xmax);
	ymin = str2num(ymin);
	ymax = str2num(ymax);
	omgmax = 3.0;
	domg = 0.4;
	pmax = 1.0;
	dp = 0.1;
elseif nargin == 1
	xmin = -3;
	xmax = 3;
	ymin = -3;
	ymax = 3;
	omgmax = 3.0;
	domg = 0.4;
	pmax = 1.0;
	dp = 0.1;
elseif nargin == 0
	foldername = 'new';
	xmin = -3;
	xmax = 3;
	ymin = -3;
	ymax = 3;
	omgmax = 3.0;
	domg = 0.4;
	pmax = 1.0;
	dp = 0.1;
end

% obtain nt and nsave ---------------------------------------------------------

info_file = strcat('../../', foldername, '/run.info');
file = fopen(info_file, 'r');
count=1;
nt=0;
nsave=0;
start_step=0;
for i=1:6
	[opt, count] = fscanf(file, '%s', 1);
	[value, count] = fscanf(file, '%f', 1);
	if (strcmpi(opt,'--nt')==1)
		nt = value;
	end
	if (strcmpi(opt,'--nsave')==1)
		nsave = value;
	end
	if (strcmpi(opt,'--start_step')==1)
		start_step = value;
	end
end

% obtain grid data ------------------------------------------------------------

grid_file = strcat('../../', foldername, '/grid');
file = fopen(grid_file, 'r');
nx = fscanf(file, '%d', 1);
x = fscanf(file, '%f', nx+1);
ny = fscanf(file, '%d', 1);
y = fscanf(file, '%f', ny+1);

% pressure
xp = zeros(1,nx);
yp = zeros(1,ny);
xp(1:nx) = 0.5*(x(2:nx+1)+x(1:nx));
yp(1:ny) = 0.5*(y(2:ny+1)+y(1:ny));
[Xp, Yp] = meshgrid(xp, yp);

xo = zeros(1,nx-1);
yo = zeros(1,ny-1);
xo(1:nx-1) = 0.5*(xp(1:nx-1) + xp(2:nx));
yo(1:ny-1) = 0.5*(yp(1:ny-1) + yp(2:ny));
%
%
[Xo, Yo] = meshgrid(xo, yo);

[Xu, Yu] = meshgrid(xo, yp);
[Xv, Yv] = meshgrid(xp, yo);

% step through different save points ------------------------------------------

N_u = (nx-1)*ny;
N_v = nx*(ny-1);
N_p = nx*ny;
N_omg=(nx-1)*(ny-1);
f = [];
for n=start_step+nsave:nsave:nt
	subfolder = num2str(n, '%07d');
	
	u = [];
	v = [];
	q_file = strcat('../../', foldername, '/', subfolder, '/q');
	file = fopen(q_file, 'r');
	nq = fscanf(file, '%d', 1);
	
	% populate u
	for I=0:N_u-1
		q = fscanf(file, '%f', 1);
		j = floor(I/(nx-1)) + 1;
		u = [u, q/(y(j+1)-y(j))];
	end
	
	% populate v
	for I=0:N_v-1
		q = fscanf(file, '%f', 1);
		i = rem(I, nx) + 1;
		v = [v, q/(x(i+1)-x(i))];
	end
	u = reshape(u, nx-1, ny);
	v = reshape(v, nx, ny-1);
	
	omg = [];
	for I=0:N_omg-1
		i = rem(I, nx-1);
		j = floor(I/(nx-1));
		value = (v(i+2,j+1)-v(i+1,j+1))/(xp(i+2)-xp(i+1)) - (u(i+1,j+2)-u(i+1,j+1))/(yp(j+2)-yp(j+1));
		omg = [omg, value];
	end
	
	omg = transpose(reshape(omg, nx-1, ny-1));
	contourf(Xo, Yo, omg, [-omgmax:domg:omgmax]);
	%	shading flat
	%	shading interp
	%	h = pcolor(Xo, Yo, omg);
	%	set(h,'edgecolor','none') 
	set(gca, 'FontName', 'Arial', 'FontSize', 14);
	set(gcf, 'PaperPositionMode', 'auto', 'PaperSize', [6 4.5]);
	axis equal;
	axis([xmin xmax ymin ymax]);
	caxis([-omgmax omgmax])
	out_file = strcat('../../', foldername, '/o', subfolder, '.png');
	print('-dpng', '-r300', out_file);
	
	lambda_file = strcat('../../', foldername, '/', subfolder, '/lambda');
	file = fopen(lambda_file, 'r');
	nlambda = fscanf(file, '%d', 1);
	% populate p	
	p = fscanf(file, '%f', N_p);
	
	p = transpose(reshape(p, nx, ny));
	contourf(Xp, Yp, p, [-pmax:dp:pmax]);
	set(gca, 'FontName', 'Arial', 'FontSize', 14);
	set(gcf, 'PaperPositionMode', 'auto', 'PaperSize', [5 4.5]);
	axis equal;
	axis([xmin xmax ymin ymax]);
	caxis([-pmax pmax]);
	out_file = strcat('../../', foldername, '/p', subfolder, '.png');
	print('-dpng', '-r300', out_file);
	
	% plot u
	%mesh(Xu, Yu, transpose(u), 'EdgeColor', 'black')
	%out_file = strcat('../../', foldername, '/umesh', subfolder, '.png');
	%print('-dpng', '-r300', out_file);
	
	contourf(Xu, Yu, transpose(u), [-1:0.1:1])
	colorbar;
	axis equal;
	axis([xmin xmax ymin ymax]);
	caxis([-1 1]);
	out_file = strcat('../../', foldername, '/u', subfolder, '.png');
	print('-dpng', '-r300', out_file);
	
	% plot v
	%mesh(Xv, Yv, transpose(v), 'EdgeColor', 'black')
	%out_file = strcat('../../', foldername, '/vmesh', subfolder, '.png');
	%print('-dpng', '-r300', out_file);
	
	contourf(Xv, Yv, transpose(v), [-1:0.1:1])
	colorbar;
	axis equal;
	axis([xmin xmax ymin ymax]);
	caxis([-1 1]);
	out_file = strcat('../../', foldername, '/v', subfolder, '.png');
	print('-dpng', '-r300', out_file); 
end
