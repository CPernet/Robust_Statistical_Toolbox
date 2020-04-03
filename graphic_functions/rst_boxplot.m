function rst_boxplot(varargin)

% routine to make box plots specifying several parameters for display
%
% FORMAT rst_boxplot(data,'LineWidth',l,'BoxColors',rgb,'BoxFillColors',rgb,'NewFig','on')
%
% INPUT data a n*m matrix of data to plot using box plot
%       optional:
%       LineWidth the thickness of box and wisker lines (default = 3)
%       BoxColors a RGB matrix (n*3) of box colors (default from rst_colour_map)
%       BoxFillColors a RGB matrix (n*3) of box colors or [] (fefault same of BoxColors)
%       NewFig 'on' (default) or 'off' on plot in current figure
%
% Cyril Pernet RS toolbox
% -----------------------

%% deal with inputs
data          = varargin{1};
if ~isnumeric(data)
    error('1st variable must be numerical')
end
linewidth     = 3;
boxcolors     = rst_colour_maps(size(data,2)+2);
boxfillcolors = flipud(boxcolors);
newfig        = 'on';

for v=2:nargin
    if strcmpi(varargin{v},'LineWidth')
        linewidth = varargin{v+1}; 
        if any(size(linewidth) == [1 1])
            warn('LineWidth should be a single value, swithching to default');
            linewidth = 3;
        end
        
    elseif strcmpi(varargin{v},'BoxColor') || strcmpi(varargin{v},'BoxColour')
            boxcolors = varargin{v+1}; 
            if size(boxcolors,1) == 3 && size(boxcolors,2) ~= 3
                boxcolors = boxcolors';
            end
            if size(boxcolors,2) ~= 3
                warn('BoxColor matrix not recognized, swithing to defaults')
                boxcolors     = rst_colour_maps(size(data,2)+2);
            end
                
    elseif strcmpi(varargin{v},'BoxColor') || strcmpi(varargin{v},'BoxColour')
            boxfillcolors = flipud(varargin{v+1});
            if size(boxfillcolors,1) == 3 && size(boxfillcolors,2) ~= 3
                boxfillcolors = boxfillcolors';
            end
            if size(boxfillcolors,2) ~= size(boxcolors)
                warn('BoxFillColor matrix not recognized, swithing to defaults')
                boxfillcolors = flipud(boxcolors);
            end
    end
end

%% plot

if strcmpi(newfig,'on')
    figure('Name','rst boxplot')
end

boxplot(data, ...
    'boxstyle','outline', ...
    'whisker',1.5, ...
    'positions',1:1:size(data,2), ...
    'widths',0.2,...
    'colors',boxcolors);
    
% lines
a = findall(gca,'type','line');
for l=1:size(a,1)
    if strcmp(a(l).LineStyle,'-') || strcmp(a(l).LineStyle,'--')
        set(a(l), 'linewidth',linewidth);
    end
end

% box color
if ~isempty(boxfillcolors)
    h = findobj(gca,'Tag','Box');
    for j=length(h):-1:1
        patch(get(h(j),'XData'),get(h(j),'YData'),boxfillcolors(j+2,:),'FaceAlpha',.4);
    end
end 

% outline
grid on; box on
set(gca,'FontSize',12,'Layer','Top')



