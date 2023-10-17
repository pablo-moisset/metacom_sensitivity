%unstable #9
%r = load_serialized('tseries_15664.24_tp4489817858643590_P04_p100le200S45c5x5u_var0.5_T1_D-0.333_px2_CON0.2.hlp.gz') ;

%stable #3
r=load_serialized('MOD_landscape/tseries_15538.20_tp6360292031686362_P00_p100le200S45c5x5u_var0.5_T0.2_D-0.333_px2_CON0.2.hlp.gz')

f=fopen('tmp_matrices.txt', 'w');
inputForm(r.M{1,3}.A,f)
inputForm(r.M{1,3}.CM,f)
inputForm(r.M{1,3}.R,f)
fclose(f)
