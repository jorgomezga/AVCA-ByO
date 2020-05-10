pd=pwd;

output = pwd;%dirPlus(pwd, 'ReturnDirs', true);

m2html('mfiles','./src','htmldir',[ 'doc' ], 'recursive', 'on',...
    'graph', 'on', 'source', 'off');