define_constants;
mpopt = mpoption('out.all', 0, 'verbose', 2);
mpopt = mpoption(mpopt, 'cpf.stop_at', 'FULL', 'cpf.step', 0.2);
mpopt = mpoption(mpopt, 'cpf.plot.level', 2);
mpcb = loadcase('t_case9_pfv2');                    % load base case
mpct = mpcb;                                        % set up target case with
mpct.gen(:, [PG QG]) = mpcb.gen(:, [PG QG]) * 2.5;  % increased generation
mpct.bus(:, [PD QD]) = mpcb.bus(:, [PD QD]) * 2.5;  % and increased load
results = runcpf(mpcb, mpct, mpopt);
