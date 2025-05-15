function [A36, qoi] = run_forward_model(parameters, qoi_option)

% Runs the coare 3.6 algorithm given the parameters [a1,a2,Ad,Bd]
y=load('VP_test_data_1.txt');
jdy=y(:,1)-y(1,1)+datenum(2000,1,8)-datenum(1999,12,31);%julian day
%COLUMN A: DATE & TIME (UTC)
P=y(:,2);%COLUMN B: Atmospheric Pressure mbs
ta=y(:,3);%COLUMN C: Air Temperature (2m)
ts=y(:,4);%COLUMN D: Sea temperature (deprth 3m)
u=y(:,5);%COLUMN E: Wind speed (m/sec, 2m)
qa=y(:,6);%COLUMN F: Specific Humidity (gr/Kgr, 10m)
cloudf=y(:,7);%COLUMN G: Total Cloud Coverage (8 times the same daily value)

%**************  set height of the input data  ********
zu=2;
zt=2;
zq=10;
%******************************   Bogus in downward solar, IR flux, and BL
%height
Rs=200*ones(size(jdy));
Rl=400*ones(size(jdy));
zi=600*ones(size(jdy));
lat=35;%latitude of the site
%**************  set reference height for output of mean variables (e.g., 10-m)  ********
zref_u=10;
zref_t=10;
zref_q=10;
%%************************

rh=relhum5([ta qa P]);
%A=coare30vn_ref(u,zu,ta,zt,rh,zq,P,ts,Rs,Rl,lat,zi,zref_u,zref_t,zref_q);%%% original 3.0 veresion
ss=35*ones(length(u),1);
cp=NaN;
sigH=NaN;
Rainrate=0;
A36=coare36vn_zrf(u,zu,ta,zt,rh,zq,P,ts,Rs,Rl,lat,zi,Rainrate,ss,cp,sigH,zref_u,zref_t,zref_q, parameters);

% Bring in data
yy=load('VP_test_results_1.txt');
%COLUMN A: DATE & TIME
Rsday=yy(:,2);%COLUMN B: DAILY INSOLATION (8 same values)
Rl=yy(:,3);%COLUMN C: LONG-WAVE RADIATION
Hs=-yy(:,4);%COLUMN D: SENSIBLE HEAT
Hl=-yy(:,5);%COLUMN E: LATENT HEAT
Tau=yy(:,6);%COLUMN F: WIND STRESS


% quantities of interest
if qoi_option == 1
    % RMSE of stress, sensible heat, latent heat against data
    qoi = [rmse(A36(:,2),Tau), rmse(A36(:,3), Hs), rmse(A36(:,4), Hl)];  
elseif qoi_option == 2
    % mean values of Cd, Ch, Ce
    qoi = [mean(A36(:,13)), mean(A36(:,14)), mean(A36(:,15))];  
end

end