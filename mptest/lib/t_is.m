function t_is(got, expected, prec, msg)

if nargin < 4
	msg = '';
end
if nargin < 3 | isempty(prec)
	prec = 5;
end

cond = ( max(max(abs(got - expected))) < 10^(-prec) );

t_ok(cond, msg);
if ~cond
	got
	expected
	got_minus_expected = got - expected
end

return;
