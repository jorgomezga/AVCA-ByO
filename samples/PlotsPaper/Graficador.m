function Graficador( name, xLab, yLab, titulo, axisIn, lineW, fontS, bColor, str )

if nargin < 9, str = []; end
if nargin < 8, bColor = 1; end
if nargin < 7, fontS = 30; end
if nargin < 6, lineW = 3; end
if nargin < 5, axisIn = []; end
if nargin < 4, titulo = ''; end
if nargin < 3, yLab   = ''; end
if nargin < 2, xLab   = ''; end

nameeps = [name, '.eps'];

if ~isempty( titulo )
    title ( titulo, 'Interpreter', 'LaTex')
end
if ~isempty( yLab )
    ylabel( yLab, 'Interpreter', 'LaTex' )
end
if ~isempty( xLab )
    xlabel( xLab, 'Interpreter', 'LaTex' )
end
pretty_xyplot

% Graphics
p = gca;
set( gcf, 'color', 'w' );
% set( gca, 'fontsize', fontS )
set(findall(gcf,'-property','FontSize'),'FontSize',fontS)

try
    set( p.Children, 'LineWidth', lineW ) % gcf
catch e
    warning( e.message )
end
    
if ~isempty( axisIn )
    axis( axisIn )
end

if ~isempty( str )
    TextLocation( str,'Interpreter','latex', 'Location', 'Best' )
end

set(findall(gcf,'-property','XColor'),'XColor',[.3 .3 .3])
set(findall(gcf,'-property','YColor'),'YColor',[.3 .3 .3])

set( gcf,'WindowState','maximized')

savefig( gcf, name )
if bColor
    saveas( gcf, nameeps,'epsc' )
else
    saveas( gcf, nameeps )
end
close(gcf)

system( ['epstopdf ', nameeps] )