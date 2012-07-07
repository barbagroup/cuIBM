%{
  Copyright (C) 2011 by Anush Krishnan, Simon Layton, Lorena Barba

  Permission is hereby granted, free of charge, to any person obtaining a copy
  of this software and associated documentation files (the "Software"), to deal
  in the Software without restriction, including without limitation the rights
  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
  copies of the Software, and to permit persons to whom the Software is
  furnished to do so, subject to the following conditions:

  The above copyright notice and this permission notice shall be included in
  all copies or substantial portions of the Software.

  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
  THE SOFTWARE.	
%}

function cavity(foldername)

if nargin == 0
	foldername = 'new';
end

yCen = [0.0000 0.0547 0.0625 0.0703 0.1016 0.1719 0.2813 0.4531 0.5000 0.6172 0.7344 0.8516 0.9531 0.9609 0.9688 0.9766 1.0000];
xCen = [0.0000 0.0625 0.0703 0.0781 0.0938 0.1563 0.2266 0.2344 0.5000 0.8047 0.8594 0.9063 0.9453 0.9531 0.9609 0.9688 1.0000];
uCavCen100 = [0.00000 -0.03717 -0.04192 -0.04775 -0.06434 -0.10150 -0.15662 -0.21090 -0.20581 -0.13641 0.00332 0.23151 0.68717 0.73722 0.78871 0.84123 1.00000];
uCavCen1000= [0.00000 -0.18109 -0.20196 -0.22220 -0.29730 -0.38289 -0.27805 -0.10648 -0.06080  0.05702 0.18719 0.33304 0.46604 0.51117 0.57492 0.65928 1.00000];
vCavCen100 = [0.00000 0.09233 0.10091 0.10890 0.12317 0.16077 0.17507 0.17527 0.05454 -0.24533 -0.22445 -0.16914 -0.10313 -0.08864 -0.07391 -0.05906 0.00000];
vCavCen1000= [0.00000 0.27485 0.29012 0.30353 0.32627 0.37095 0.33075 0.32235 0.02426 -0.31966 -0.42665 -0.51550 -0.39188 -0.33714 -0.27669 -0.21388 0.00000];

% default values --------------------------------------------------------------
xmin = 0;
xmax = 1;
ymin = -0.7;
ymax = 1.3;

% obtain nt and nsave ---------------------------------------------------------

info_file = strcat('../../', foldername, '/run.info');
file = fopen(info_file, 'r');
count=1;
nt=0;
nsave=0;
for i=1:5
	[opt, count] = fscanf(file, '%s', 1);
	[value, count] = fscanf(file, '%f', 1);
	if (strcmpi(opt,'--nt')==1)
		nt = value;
	end
	if (strcmpi(opt,'--nsave')==1)
		nsave = value;
	end
end

% obtain grid data ------------------------------------------------------------

grid_file = strcat('../../', foldername, '/grid');
file = fopen(grid_file, 'r');
nx = fscanf(file, '%d', 1);
x = [];
for i=1:nx+1
	x = [x, fscanf(file, '%f', 1)];
end
	ny = fscanf(file, '%d', 1);
y = [];
for i=1:ny+1
	y = [y, fscanf(file, '%f', 1)];
end

% pressure
xp = [];
yp = [];
for i=1:nx
	xp = [xp, 0.5*(x(i+1) + x(i))];
end
for i=1:ny
	yp = [yp, 0.5*(y(i+1) + y(i))];
end
[Xp, Yp] = meshgrid(xp, yp);

% step through different save points ------------------------------------------

N_u = (nx-1)*ny;
N_v = nx*(ny-1);

subfolder = num2str(nt, '%07d');

U = [];
V = [];
q_file = strcat('../../', foldername, '/', subfolder, '/q');
file = fopen(q_file, 'r');
nq = fscanf(file, '%d', 1);
for I=0:N_u-1
	q = fscanf(file, '%f', 1);
	j = floor(I/(nx-1)) + 1;
	U = [U, q/(y(j+1)-y(j))];
end
for I=0:N_v-1
	q = fscanf(file, '%f', 1);
	i = rem(I, nx) + 1;
	V = [V, q/(x(i+1)-x(i))];
end

u=[];
v=[];
for j=0:ny-1
	u = [u, U( j*(nx-1) + nx/2 )];
end

for i=1:nx
	v = [v, V( (ny/2-1)*nx + i)];
end

plot(yp, u, 'r', yCen, uCavCen1000, 'o');
%axis equal;
axis([xmin xmax ymin ymax]);
out_file = strcat('../../', foldername, '/u.png');
print('-dpng', '-r600', out_file);

plot(xp, v, 'r', xCen, vCavCen1000, 'o');
%axis equal;
axis([xmin xmax ymin ymax]);
out_file = strcat('../../', foldername, '/v.png');
print('-dpng', '-r600', out_file);