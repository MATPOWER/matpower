function t_ok(cond, msg)

global t_counter;
global t_ok_cnt;
global t_not_ok_cnt;

if nargin < 2 | strcmp(msg, '')
	msg = '';
else
	msg = [' - ', msg];
end
if cond
	t_ok_cnt = t_ok_cnt + 1;
else
	t_not_ok_cnt = t_not_ok_cnt + 1;
	fprintf('not ');
end
fprintf('ok %d%s\n', t_counter, msg);
t_counter = t_counter + 1;

return;
