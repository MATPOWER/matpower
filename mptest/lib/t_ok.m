function t_ok(cond, msg)
	global t_counter;
	if nargin < 2 | strcmp(msg, '')
		msg = '';
	else
		msg = [' - ', msg];
	end
	if ~cond
		fprintf('not ');
	end
	fprintf('ok %d%s\n', t_counter, msg);
	t_counter = t_counter + 1;
return;
