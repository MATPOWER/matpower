% modify the path temporarily and silently to install MATPOWER
cwd = pwd;              % current working directory
imp = which('install_matpower');    % install_matpower path
if isempty(imp)
    error('matpower_project_startup: install_matpower() not found');
end
[p, n, e] = fileparts(imp);
cd(p);
install_matpower(1, 0, 0);
cd(cwd);
