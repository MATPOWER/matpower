pw = [1 0 0];
mpopt1 = mpoption('exp.sys_wide_zip_loads.pw',pw);
mpopt2 = mpoption(mpopt1,'pf.alg','PQSUM','pf.radial.vcorr',1);
mpc = case4_dist;
% mpc = case33bw_radial;
% mpc = scale_load(3,mpc);
mpc1 = runpf(mpc,mpopt1);
mpc2 = runpf(mpc,mpopt2);
compare_case(mpc1,mpc2);