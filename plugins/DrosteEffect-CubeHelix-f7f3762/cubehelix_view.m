function [map,prm] = cubehelix_view(N,start,rots,sat,gamma,irange,domain)
% An interactive figure for Cubehelix colormap parameter selection. With demo!
%
% (c) 2015 Stephen Cobeldick
%
% View Dave Green's Cubehelix colorschemes in a figure.
%
% * Two colorbars give the colorscheme in color and grayscale.
% * A button toggles between 3D-cube and 2D-lineplot of the RGB values.
% * A button toggles an endless demo of randomly generated colorschemes.
% * Nine sliders allow real-time adjustment of the Cubehelix parameters.
% * Warning text if any RGB values are clipped.
% * Warning text if the grayscale is not monotonic increasing/decreasing.
%
% Syntax:
%  cubehelix_view
%  cubehelix_view(N)
%  cubehelix_view(N,start,rots,sat,gamma)
%  cubehelix_view(N,start,rots,sat,gamma,irange)
%  cubehelix_view(N,start,rots,sat,gamma,irange,domain)
%  cubehelix_view(N,[start,rots,sat,gamma],...)
%  cubehelix_view([],...)
%  cubehelix_view({axes/figure handles},...) % see "Adjust External Colormaps"
%  [map,prm] = cubehelix_view(...)
%
% Calling the function with an output argument blocks MATLAB execution until
% the figure is deleted: the final colormap and parameters are then returned.
%
% Cubehelix is defined here: http://astron-soc.in/bulletin/11June/289392011.pdf
% For more information and examples: http://www.mrao.cam.ac.uk/~dag/CUBEHELIX/
%
% See also CUBEHELIX BREWERMAP RGBPLOT COLORMAP COLORMAPEDITOR COLORBAR UICONTROL ADDLISTENER
%
% ### Adjust External Colormaps ###
%
% % Example:
% load spine
% image(X)
% cubehelix_view({gca})
%
% Very useful! Simply provide a cell array of axes or figure handles when
% calling this function, and their colormaps will be updated in real-time.
% Note that MATLAB versions <=2010 only support axes handles for this option.
%
% ### Input and Output Arguments ###
%
% Input Arguments (*==default):
%  N     = NumericScalar, an integer to define the colormap length.
%        = *[], colormap length of one hundred and twenty-eight (128).
%        = {axes/figure handles}, their colormaps will be updated by this function.
%  start = NumericScalar, the helix's start color, with R=1, G=2, B=3 (modulus 3).
%  rots  = NumericScalar, the number of R->G->B rotations over the scheme length.
%  sat   = NumericScalar, controls how saturated the colors are (saturation).
%  gamma = NumericScalar, can be used to emphasize low or high intensity values.
%  irange = NumericVector, range of brightness levels of the colormap's endnodes. Size 1x2.
%  domain = NumericVector, domain of the cubehelix calculation (endnode positions). Size 1x2.
%
% Output Arguments (block execution until the figure is deleted!):
%  map = NumericMatrix, the cubehelix colormap defined when the figure is closed.
%  prm = NumericVector, the parameters of <map>: [start,rots,sat,gamma,irange,domain].
%
% [map,prm] = cubehelix_view(N, start, rots, sat, gamma, *irange, *domain)
% OR
% [map,prm] = cubehelix_view(N, [start,rots,sat,gamma], *irange, *domain)

% ### Input Wrangling ###
%
if nargin==0 || isnumeric(N)&&isempty(N)
	N = 128;
elseif iscell(N)&&numel(N)
	tmp = all(1==cellfun('prodofsize',N)&cellfun(@ishghandle,N));
	assert(tmp,'Input <N> may be a cell array of axes or figure handles.')
else
	assert(isnumeric(N)&&isscalar(N),'Input <N> must be a scalar numeric.')
	assert(isreal(N)&&fix(N)==N&&N>0,'Input <N> must be positive real integer: %g+%gi',N,imag(N))
	N = double(N);
end
%
% Random Cubehelix parameters based on the current time:
T = sum(clock*100);
foo = @(n)min(3,max(0,sqrt(-log(rem(T/pow2(n),1))*2)));
bar = @(n)min(3,max(0,rem(T/pow2(n),1)^36));
dfa(5:8,1) = [bar(3);1-bar(2);bar(1);1-bar(0)]; % irange and domain
dfa(3:4,1) = [foo(5);foo(4)]; % sat and gamma
dfa(2,1) = min(3,max(-3,log10(rem(T,1)/(1-rem(T,1))))); % rots
dfa(1,1) = rem(T,3); % start
%
%      [sta; rots; sat; gam; yrng; domn]
%dfa = [0.5; -1.5;   1;   1; 0; 1; 0; 1];
%
stp = '%s input can be a vector of the four Cubehelix parameters.';
str = '%s input can be a vector of the endnode brightness levels (range).';
std = '%s input can be a vector of the endnode relative positions (domain).';
%
switch nargin
	case {0,1}
		[~,~,h] = chvUpDt(N, dfa(1:8));
	case 2
		start = chvChk(4,start,stp,'Second');
		[~,~,h] = chvUpDt(N, [start;dfa(5:8)]);
	case 3
		start = chvChk(4,start,stp,'Second');
		rots  = chvChk(2,rots, str,'Third');
		[~,~,h] = chvUpDt(N, [start;rots;dfa(7:8)]);
	case 4
		start = chvChk(4,start,stp,'Second');
		rots  = chvChk(2,rots, str,'Third');
		sat   = chvChk(2,sat,  std,'Fourth');
		[~,~,h] = chvUpDt(N, [start;rots;sat]);
	case 5
		[~,~,h] = chvUpDt(N, [chvC2V(start,rots,sat,gamma);dfa(5:8)]);
	case 6
		irange = chvChk(2,irange,str,'Sixth');
		[~,~,h] = chvUpDt(N, [chvC2V(start,rots,sat,gamma);irange;dfa(7:8)]);
	case 7
		irange = chvChk(2,irange,str,'Sixth');
		domain = chvChk(2,domain,std,'Seventh');
		[~,~,h] = chvUpDt(N, [chvC2V(start,rots,sat,gamma);irange;domain]);
	otherwise
		error('Wrong number of inputs. Enter parameters individually or in a vector.')
end
%
if nargout
	waitfor(h);
	[~,prm,map] = chvUpDt();
end
%
end
%----------------------------------------------------------------------END:cubehelix_view
function V = chvC2V(varargin)
% Check that all of the input variables are real scalar numerics.
str = 'Input Cubehelix parameters must be %s values.';
assert(all(cellfun(@isnumeric,varargin)),str,'numeric')
assert(all(cellfun(@isscalar,varargin)),str,'scalar')
assert(all(cellfun(@isfinite,varargin)),str,'finite')
assert(all(cellfun(@isreal,varargin)),str,'real')
V = cellfun(@double,varargin(:));
end
%----------------------------------------------------------------------END:chvC2V
function x = chvChk(n,x,msg,ord)
% Check that the input variable <x> is real numeric vector with <n> elements.
assert(isnumeric(x)&&isreal(x)&&all(isfinite(x))&&isvector(x)&&numel(x)==n,msg,ord)
x = double(x(:));
end
%----------------------------------------------------------------------END:chvChk
function [N,V,ghnd] = chvUpDt(N,V)
% Draw a new figure or update an existing figure. Callback for sliders & demo.
%
persistent H prvV prvN extH
%
% LHS and RHS slider bounds/limits, and slider step sizes:
lbd = [  1, 0,-3, 0, 0, 0, 0, 0, 0].';
rbd = [128, 3, 3, 3, 3, 1, 1, 1, 1].';
stp = [  1,[5, 5, 5, 5, 1, 1, 1, 1]/100;... % minor
        10,[5, 5, 5, 5, 1, 1, 1, 1]/10].';  % major
%       [N,st,ro,sa,ga,i1,i2,d1,d2]
%
switch nargin
	case 0 % Demo initialize OR return colormap
		N = prvN;
		V = prvV;
		if nargout==3 % Return colormap
			ghnd = cubehelix(N, V(1:4), V(5:6), V(7:8));
			return
		end
	case 1 % Slider callback
		if get(H.demo,'Value') % Demo
			return
		elseif N % cubehelix parameters
			V = prvV;
			V(N) = get(H.vSld(N+1),'Value');
			N = prvN;
		else % parameter N only
			V = prvV;
			N = round(get(H.vSld(1),'Value'));
		end
	case 2 % Demo update OR function call
		if iscell(N) % External axes/figure
			extH = N;
			N = size(colormap(extH{1}),1);
		end
		if isempty(H) || ~ishghandle(H.fig)
			if nargout<3 % Demo
				return
			end
			% Create a new figure:
			H = chvPlot(lbd,rbd,stp,[N;V]);
		else % Update slider positions:
			set(H.vSld, {'Value'},num2cell(max(lbd,min(rbd,[N;V]))));
		end
		%
	otherwise
		error('How did this happen? This should not be possible.')
end
%
prvN = N;
prvV = V;
ghnd = H.fig;
%
% Update parameter value text:
set(H.vTxt(1), 'String',sprintf('%.0f',N));
set(H.vTxt(2:end), {'String'},sprintfc('%.2f',V));
%
% Get Cubehelix colormap and grayscale equivalent:
[map,lo,hi]  = cubehelix(N, V(1:4), V(5:6), V(7:8));
mag = sum(map*[0.298936;0.587043;0.114021],2);
%
% Update colorbar values:
set(H.cbAx, 'YLim',[0,abs(N)+(N==0)]+0.5);
set(H.cbIm(1), 'CData',reshape(map,[],1,3))
set(H.cbIm(2), 'CData',repmat(mag,[1,1,3]))
%
% Update 2D line / 3D patch values:
if get(H.D2D3,'Value')
	set(H.ln2D, 'XData',linspace(0,1,abs(N)));
	set(H.ln2D, {'YData'},num2cell([map,mag],1).');
else
	set(H.pt3D, 'XData',map(:,1), 'YData',map(:,2), 'ZData',map(:,3), 'FaceVertexCData',map)
end
%
% Update warning text:
mad = diff(mag);
str = {'Not Monotonic';'Clipped'};
set(H.warn,'String',str([any(mad<=0)&&any(0<=mad);any(lo(:))||any(hi(:))]));
%
% Update external axes/figure:
if ~isempty(extH)
	for k = find(cellfun(@ishghandle,extH))
		colormap(extH{k},map);
	end
end
%
end
%----------------------------------------------------------------------END:chvUpDt
function H = chvPlot(lbd,rbd,stp,V)
% Draw a new figure with RGBplot axes, ColorBar axes, and uicontrol sliders.
%
% Parameter names for each slider:
names = {'N';'start';'rotations';'saturation';'gamma';'irange(1)';'irange(2)';'domain(1)';'domain(2)'};
%
M = numel(names); % number of sliders
gap = 0.01; % gaps
bth = 0.04; % demo height
btw = 0.09; % demo width
uih = 0.40; % height of UI control group
cbw = 0.24; % width of both colorbars
axh = 1-uih-2*gap; % axes height
axw = 1-cbw-2*gap; % axes width
slh = uih/M - gap; % slider height
%
H.fig = figure('HandleVisibility','callback', 'Color','white',...
	'IntegerHandle','off', 'NumberTitle','off',...
	'Name','Cubehelix Interactive Parameter Selector');
%
% Add 2D lineplot:
H.ax2D = axes('Parent',H.fig, 'Position',[gap, uih+gap, axw, axh],...
	'ColorOrder',[1,0,0; 0,1,0; 0,0,1; 0.6,0.6,0.6], 'HitTest','off',...
	'Visible','off', 'XLim',[0,1], 'YLim',[0,1], 'XTick',[], 'YTick',[]);
H.ln2D = line([0,0,0,0;1,1,1,1],[0,0,0,0;1,1,1,1], 'Parent',H.ax2D, 'Visible','off');
%
% Add 3D scatterplot:
H.ax3D = axes('Parent',H.fig, 'OuterPosition',[0, uih, axw+2*gap, 1-uih],...
	'Visible','on', 'XLim',[0,1], 'YLim',[0,1], 'ZLim',[0,1], 'HitTest','on');
H.pt3D = patch('Parent',H.ax3D, 'XData',[0;1], 'YData',[0;1], 'ZData',[0;1],...
	'Visible','on', 'LineStyle','none', 'FaceColor','none', 'MarkerEdgeColor','none',...
	'Marker','o', 'MarkerFaceColor','flat', 'MarkerSize',10, 'FaceVertexCData',[1,1,0;1,0,1]);
view(H.ax3D,3);
grid(H.ax3D,'on')
xlabel(H.ax3D,'Red')
ylabel(H.ax3D,'Green')
zlabel(H.ax3D,'Blue')
%
% Add warning text:
H.warn = text('Parent',H.ax2D, 'Units','normalized', 'Position',[0,1],...
	'HorizontalAlignment','left', 'VerticalAlignment','top', 'Color','r');
%
% Add demo button:
H.demo = uicontrol(H.fig, 'Style','togglebutton', 'Units','normalized',...
	'Position',[axw-btw+gap,uih+gap+0*bth,btw,bth], 'String','Demo',...
	'Max',1, 'Min',0, 'Callback',@chvDemo);
% Add 2D/3D button:
H.D2D3 = uicontrol(H.fig, 'Style','togglebutton', 'Units','normalized',...
	'Position',[axw-btw+gap,uih+gap+1*bth,btw,bth], 'String','2D / 3D',...
	'Max',1, 'Min',0, 'Callback',@(h,~)chv2D3D(H,h));
%
% Add colorbars:
C = reshape([1,1,1],1,[],3);
H.cbAx(2) = axes('Parent',H.fig, 'Visible','off', 'Units','normalized',...
	'Position',[1-cbw/2,gap,cbw/2-gap,1-2*gap], 'YLim',[0.5,1.5], 'HitTest','off');
H.cbAx(1) = axes('Parent',H.fig, 'Visible','off', 'Units','normalized',...
	'Position',[1-cbw/1,gap,cbw/2-gap,1-2*gap], 'YLim',[0.5,1.5], 'HitTest','off');
H.cbIm(2) = image('Parent',H.cbAx(2), 'CData',C);
H.cbIm(1) = image('Parent',H.cbAx(1), 'CData',C);
%
% Add parameter sliders, listeners, and corresponding text:
V = max(lbd,min(rbd,V));
H.vLab(M) = 0;
H.vTxt(M) = 0;
H.vSld(M) = 0;
for m = 1:M
	Y = gap+(M-m)*(slh+gap);
	H.vLab(m) = uicontrol(H.fig,'Style','text', 'Units','normalized',...
		'Position',[gap,Y,btw,slh], 'String',names{m});
	H.vTxt(m) = uicontrol(H.fig,'Style','text', 'Units','normalized',...
		'Position',[gap+btw,Y,btw,slh], 'String','X');
	H.vSld(m) = uicontrol(H.fig,'Style','slider', 'Units','normalized',...
		'Position',[2*btw+gap,Y,axw-2*btw,slh], 'Min',lbd(m), 'Max',rbd(m),...
		'SliderStep',stp(m,:)/(rbd(m)-lbd(m)), 'Value',V(m));
	addlistener(H.vSld(m), 'Value', 'PostSet',@(a,b)chvUpDt(m-1));
end
%
end
%----------------------------------------------------------------------END:chvPlot
function chv2D3D(H,tgh)
% Switch between 2D-line and 3D-cube: swap visibility and hittest, then update.
%
if get(tgh,'Value')% 2D
	set(H.ax3D, 'HitTest','off', 'Visible','off')
	set(H.ax2D, 'HitTest','on')
	set(H.pt3D, 'Visible','off')
	set(H.ln2D, 'Visible','on')
else % 3D
	set(H.ax2D, 'HitTest','off')
	set(H.ax3D, 'HitTest','on', 'Visible','on')
	set(H.ln2D, 'Visible','off')
	set(H.pt3D, 'Visible','on')
end
chvUpDt();
%
end
%----------------------------------------------------------------------END:chv2D3D
function chvDemo(tgh,~)
% While the toggle button is depressed run a loop showing random Cubehelix schemes.
%
% Initial values and goal values:
[N,V] = chvUpDt();
Ng = N;
Vg = V;
%
% Step length between updates:
stpL = 0.03;
%
% Functions to randomly find new values:
randfn(7:8) = {@()rand(1,1).^42,@()1-rand(1,1).^42};
randfn(3:4) = {@()sqrt(-log(rand(1,1))*2)};
randfn(1:2) = {@()3*rand(1,1),@()randn(1,1)};
randfn(5:6) = randfn(7:8);
%
% While the toggle button is down, step values:
while ishghandle(tgh)&&get(tgh,'Value')
	%
	isnw = V==Vg; % create new goal
	Vg(isnw) = round(100*cellfun(@(fn)fn(),randfn(isnw)))/100;
	%
	isnx = abs(Vg-V)<=stpL; % close to goal
	V(isnx) = Vg(isnx);
	%
	isgo = ~(isnw|isnx); % on its way to goal
	V(isgo) = V(isgo) + stpL*sign(Vg(isgo)-V(isgo));
	%
	if N==Ng % new goal
		Ng = randi(128);
	elseif abs(N-Ng)<=1 % close to goal
		N = Ng;
	else % on its way to goal
		N = N + sign(Ng-N);
	end
	%
	% Update figure:
	[N,V] = chvUpDt(N,V);
	%
	% Frame-rate faster/slower:
	pause(0.09);
end
%
end
%----------------------------------------------------------------------END:chvDemo
% Copyright (c) 2015 Stephen Cobeldick
%
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
%
% http://www.apache.org/licenses/LICENSE-2.0
%
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and limitations under the License.
%----------------------------------------------------------------------END:license
