function A=coare36vn_zrf(u,zu,t,zt,rh,zq,P,ts,Rs,Rl,lat,zi,rain,Ss,cp,sigH,zrf_u,zrf_t,zrf_q, parameters)
%
% Vectorized version of COARE 3.6 code (Fairall et al, 2003) with 
% modification based on the CLIMODE, MBL and CBLAST experiments 
% (Edson et al., 2012). The cool skin option is retained but warm layer 
% and surface wave options removed. 
%
% Latest update: added parameters [a1,a2,gamma] as input (EYS: 04292025)
%
% This version include parameterizations using wave height and wave
% slope using cp and sigH.  If these are set to NaN, then the wind
% speed dependent formulation is used.  The parameterizations are based
%on fits to the Banner-Norison wave model and the Fairall-Edson flux
%database.  It also allows salinity as a input.  Open ocean Ss=35; Great
%Lakes Ss=0;
%
%********************************************************************
% An important component of this code is whether the inputed ts 
% represents the skin temperature of a near surface temperature.  
% How this variable is treated is determined by the jcool parameter:
% set jcool=1 if Ts is bulk ocean temperature (default),
%     jcool=0 if Ts is true ocean skin temperature. 
%********************************************************************

jcoolx=1;%set to 1.0 if you want the cool skin correction  computed for ocean surface temperature

% The code assumes u,t,rh,ts are vectors; 
% sensor heights zu,zt,zl, latitude lat, and PBL height zi may be constants;
% air pressure P and radiation Rs,Rl may be vectors or constants. 
% Default values are assigned for P,Rs,Rl,lat,and zi if these data are not 
% available.  Input NaNs to indicate no data. Defaults should be set to 
% representative regional values if possible.
%
% Input:  
%
%     u = water-relative mean wind vector magnitude (m/s) at height zu(m)
%     t = bulk air temperature (degC) at height zt(m)
%    rh = relative humidity (%) at height zq(m)
%     P = surface air pressure (mb) (default = 1015)
%    ts = water temperature (degC) see jcool below
%    Rs = downward shortwave radiation (W/m^2) (default = 150) 
%    Rl = downward longwave radiation (W/m^2) (default = 370)
%   lat = latitude (default = +45 N)
%    zi = PBL height (m) (default = 600m)
%  rain = rain rate (mm/hr)
%    Ss = sea surface salinity (PSU)
%    cp = phase speed of dominant waves (m/s)  
%  sigH =  significant wave height (m)
%  zu, zt, zq heights of the observations (m)
%  zrf_u, zrf_t, zrf_q  reference height for profile.  Use this to compare observations at different heights  
%
% The user controls the output.  This is currently set as:
%
%Output:   A=[usr tau hsb hlb hbb hsbb hlwebb tsr qsr zo  zot zoq  Cd  Ch  Ce  L zet dterx dqerx tkt Urf Trf Qrf RHrf UrfN Rnl  Le rhoa UN U10 U10N Cdn_10 Chn_10 Cen_10 RF Qs Evap T10 Q10 RH10 ug Whf Edis];
%              1   2   3   4   5    6     7    8   9  10  11   12  13  14  15 16  17   18   19    20  21  22  23  24   25   26  27  28  29  30  31    32     33     34   35 36  37   38  39  40  41  42  43 
%where
%
%   usr = friction velocity that includes gustiness (m/s)
%   tau = wind stress (N/m^2)
%   hsb = sensible heat flux into ocean (W/m^2)
%   hlb = latent heat flux into ocean (W/m^2)
%   hbb = buoyany flux into ocean (W/m^2)
%   hsbb = "sonic" buoyancy flux measured directly by sonic anemometer 
%   tsr = temperature scaling parameter (K)
%   qsr = specific humidity scaling parameter (g/Kg)
%   zot = thermal roughness length (m)
%   zoq = moisture roughness length (m)
%   Cd = wind stress transfer (drag) coefficient at height zu   
%   Ch = sensible heat transfer coefficient (Stanton number) at height zu   
%   Ce = latent heat transfer coefficient (Dalton number) at height zu
%    L = Obukhov length scale (m) 
%  zet = Monin-Obukhov stability parameter zu/L 
% dter = cool-skin temperature depression (degC)
% dqer = cool-skin humidity depression (degC)
%  tkt = cool-skin thickness (m)
%  Urf = wind speed at reference height (user can select height below)
%  Tfr = temperature at reference height
%  Qfr = specific humidity at reference height
% RHfr = relative humidity at reference height
% UrfN = neutral value of wind speed at reference height
%  Rnl = Upwelling IR radiation computed by COARE
%   Le = latent heat of vaporization
% rhoa = density of air
%   UN = neutral value of wind speed at zu
%  U10 = wind speed adjusted to 10 m
% UN10 = neutral value of wind speed at 10m
%Cdn_10 = neutral value of drag coefficient at 10m    
%Chn_10 = neutral value of Stanton number at 10m    
%Cen_10 = neutral value of Dalton number at 10m    
%Rf     = Rain heat flux (W/m^2)
%Qs     = surface specific humidity (g/g)
%Evap   = evaporation rate (mm/h)
%T10    = air temperature at 10m
%Q10    = air specific humidity at 10m
%RH10   = air relative humidity at 10m
%uq     = gustiness velocity (m/s)
%Whf    = whitecap fraction
%Edis   = energy dissipated by wave breaking (W/m^2)

% Notes: 1) u is the surface-relative wind speed, i.e., the magnitude of the
%           difference between the wind (at zu) and ocean surface current 
%           vectors.
%        2) Set jcool=0 in code if ts is true surface skin temperature,
%           otherwise ts is assumed the bulk temperature and jcool=1.
%        3) Set P=NaN to assign default value if no air pressure data 
%           available. 
%        4) Set Rs=NaN, Rl=NaN if no radiation data available.  This assigns 
%           default values to Rs, Rl so that cool skin option can be applied. 
%        5) Set lat=NaN and/or zi=NaN to assign default values if latitude
%           and/or PBL height not given. 
%        6) The code to compute the heat flux caused by precipitation is 
%           included if rain data are available (default is no rain).
%        7) Code updates the cool-skin temperature depression dter and thickness
%           tkt during iteration loop for consistency.
%        8) Number of iterations set to nits = 6.

% Reference:
%
%  Fairall, C.W., E.F. Bradley, J.E. Hare, A.A. Grachev, and J.B. Edson (2003),
%  Bulk parameterization of air sea fluxes: updates and verification for the 
%  COARE algorithm, J. Climate, 16, 571-590.
%
%Edson, J.B., J. V. S. Raju, R.A. Weller, S. Bigorre, A. Plueddemann, C.W. Fairall, 
%S. Miller, L. Mahrt, Dean Vickers, and Hans Hersbach, 2013: On the Exchange of momentum
%over the open ocean. J. Phys. Oceanogr., 43, 1589–1610. doi: http://dx.doi.org/10.1175/JPO-D-12-0173.1 

% Code history:
% 
% 1. 12/14/05 - created based on scalar version coare26sn.m with input
%    on vectorization from C. Moffat.  
% 2. 12/21/05 - sign error in psiu_26 corrected, and code added to use variable
%    values from the first pass through the iteration loop for the stable case
%    with very thin M-O length relative to zu (zetu>50) (as is done in the 
%    scalar coare26sn and COARE3 codes).
% 3. 7/26/11 - S = dt was corrected to read S = ut.
% 4. 7/28/11 - modification to roughness length parameterizations based 
%    on the CLIMODE, MBL, Gasex and CBLAST experiments are incorporated
% 5. New wave parameterization added 9/20/2017  based on fits to wave model
%-----------------------------------------------------------------------
%display('flockewe');
% convert input to column vectors
u=u(:);t=t(:);rh=rh(:);P=P(:);ts=ts(:);
Rs=Rs(:);Rl=Rl(:);lat=lat(:);zi=zi(:);
zu=zu(:);zt=zt(:);zq=zq(:);
zrf_u=zrf_u(:);zrf_t=zrf_t(:);zrf_q=zrf_q(:);
rain=rain(:);  
Ss=Ss(:);cp=cp(:);sigH=sigH(:);
N=length(u);
jcool=jcoolx*ones(N,1);jcool=jcool(:);

% set local variables to default values if input is NaN
ii=find(isnan(P)); P(ii)=1013;               % pressure
ii=find(isnan(Rs)); Rs=200;                  % incident shortwave radiation
ii=find(isnan(lat)); lat(ii)=45;             % latitude
ii=find(isnan(Rl)); Rl(ii)=400-1.6*abs(lat(ii)); % incident longwave radiation
ii=find(isnan(zi)); zi(ii)=600;              % PBL height
ii=find(isnan(Ss)); Ss(ii)=35;               % Salinity


% input variable u is assumed relative wind speed (magnitude of difference
% between wind and surface current vectors). to follow orginal Fairall code, set
% surface current speed us=0. if us data are available, construct u prior to
% using this code.
us = 0*u;

% convert rh to specific humidity
Tf=-0.0575*Ss+1.71052E-3*Ss.^1.5-2.154996E-4*Ss.*Ss;%freezing point of seawater
Qs = qsat26sea(ts,P,Ss,Tf)./1000;    % surface water specific humidity (g/kg)
[Q,Pv]  = qsat26air(t,P,rh);  % specific humidity of air (g/kg).  Assumes rh relative to ice T<0
Q=Q./1000;

ice=zeros(size(u));
iice=find(ts<Tf);ice(iice)=1;jcool(iice)=0;
zos=5E-4;
%***********  set constants **********************************************
zref=10;
Beta = 1.2;
von  = 0.4;
fdg  = 1.00; % Turbulent Prandtl number
tdk  = 273.16;
grav = grv(lat);

%***********  air constants **********************************************
Rgas = 287.1;
Le   = (2.501-.00237*ts)*1e6;
cpa  = 1004.67;
cpv  = cpa*(1+0.84*Q);
rhoa = P*100./(Rgas*(t+tdk).*(1+0.61*Q));
rhodry = (P-Pv)*100./(Rgas*(t+tdk));
visa = 1.326e-5*(1+6.542e-3.*t+8.301e-6*t.^2-4.84e-9*t.^3);

%***********  cool skin constants  ***************************************
%%%%%%%%%%%%%  salinity dependent thermal expansion coeff for water
tsw=ts;ii=find(ts<Tf);tsw(ii)=Tf(ii);
Al35   = 2.1e-5*(tsw+3.2).^0.79;
Al0   =(2.2*real((tsw-1).^0.82)-5)*1e-5;
Al=Al0+(Al35-Al0).*Ss/35;
%%%%%%%%%%%%%%%%%%%
bets=7.5e-4;%salintity expansion coeff; assumes beta independent of temperature
be   = bets*Ss;%be is beta*Salinity
%%%%  see Computing the seater expansion coefficients directly from the
%%%%  1980 equation of state.  J. Lillibridge, J.Atmos.Oceanic.Tech, 1980.
cpw  = 4000;
rhow = 1022;
visw = 1e-6;
tcw  = 0.6;
bigc = 16*grav*cpw*(rhow*visw)^3./(tcw.^2*rhoa.^2);
wetc = 0.622*Le.*Qs./(Rgas*(ts+tdk).^2);

%***********  net radiation fluxes ***************************************
Rns = 0.945.*Rs; % albedo correction
% IRup = eps*sigma*T^4 + (1-eps)*IR
% Rnl = IRup - IR
% Rnl = eps*sigma*T^4 - eps*IR  as below

Rnl = 0.97*(5.67e-8*(ts-0.3*jcool+tdk).^4-Rl); % initial value

% IRup = Rnl + IR

%****************  begin bulk loop ********************************************

gamma = parameters(3);    % EYS 04292025

%***********  first guess ************************************************
du = u-us;
dt = ts-t-.0098.*zt;
dq = Qs-Q;
ta = t+tdk;
ug = 0.5;
dter  = 0.3;
ut    = sqrt(du.^2+ug.^2);
u10   = ut.*log(10/1e-4)./log(zu/1e-4);
usr   = 0.035*u10;
zo10  = 0.011*usr.^2./grav + gamma*visa./usr;
Cd10  = (von./log(10./zo10)).^2;
Ch10  = 0.00115;
Ct10  = Ch10./sqrt(Cd10);
zot10 = 10./exp(von./Ct10);
Cd    = (von./log(zu./zo10)).^2;
Ct    = von./log(zt./zot10);
CC    = von*Ct./Cd;
Ribcu = -zu./zi./.004/Beta^3;
Ribu  = -grav.*zu./ta.*((dt-dter.*jcool)+.61*ta.*dq)./ut.^2;
zetu = CC.*Ribu.*(1+27/9*Ribu./CC);
k50=find(zetu>50); % stable with very thin M-O length relative to zu
k=find(Ribu<0); 
if length(Ribcu)==1
    zetu(k)=CC(k).*Ribu(k)./(1+Ribu(k)./Ribcu); clear k;
else
    zetu(k)=CC(k).*Ribu(k)./(1+Ribu(k)./Ribcu(k)); clear k;
end
L10 = zu./zetu;
gf=ut./du;
usr = ut.*von./(log(zu./zo10)-psiu_40(zu./L10));
tsr = -(dt-dter.*jcool).*von*fdg./(log(zt./zot10)-psit_26(zt./L10));
qsr = -(dq-wetc.*dter.*jcool)*von*fdg./(log(zq./zot10)-psit_26(zq./L10));
tkt = 0.001*ones(N,1);

%**********************************************************
%  The following gives the new formulation for the
%  Charnock variable
%**********************************************************
%%%%%%%%%%%%%   COARE 3.5 wind speed dependent charnock
charnC = 0.011*ones(N,1);
umax=19;
% a1=0.0017;
% a2=-0.0050;
a1 = parameters(1);   % EYS 04292025
a2 = parameters(2);   % EYS 04292025
charnC=a1*u10+a2;
k=find(u10>umax);
charnC(k)=a1*umax+a2;


%%%%%%%%%   if wave age is given but not wave height, use parameterized
%%%%%%%%%   wave height
    hsig=(0.02*(cp./u10).^1.1-0.0025).*u10.^2;
    hsig=max(hsig,.25);
    ii=find(~isnan(cp) & isnan(sigH));
    sigH(ii)=hsig(ii);
    
Ad=0.2;  %Sea-state/wave-age dependent coefficients from wave model
%Ad=0.73./sqrt(u10);
Bd=2.2;
zoS=sigH.*Ad.*(usr./cp).^Bd;
charnS=zoS.*grav./usr./usr;


nits=10; % number of iterations
charn=charnC;
ii=find(~isnan(cp));charn(ii)=charnS(ii);
%**************  bulk loop **************************************************

for i=1:nits
    zet=von.*grav.*zu./ta.*(tsr +.61*ta.*qsr)./(usr.^2);
    L=zu./zet;
    zo=charn.*usr.^2./grav+gamma*visa./usr; % surface roughness
    zo(iice)=zos;
    rr=zo.*usr./visa;
    zoq=min(1.6e-4,5.8e-5./rr.^.72);       % These thermal roughness lengths give Stanton and
     
     %%%%%%%%%%%%%  %Andreas 1987 for snow/ice
    ik=find(rr(iice)<=.135);
   		rt(iice(ik))=rr(iice(ik))*exp(1.250);
     	rq(iice(ik))=rr(iice(ik))*exp(1.610);
     ik=find(rr(iice)>.135 & rr(iice)<=2.5);
        rt(iice(ik))=rr(iice(ik)).*exp(0.149-.55*log(rr(iice(ik))));
     	rq(iice(ik))=rr(iice(ik)).*exp(0.351-0.628*log(rr(iice(ik))));
     ik=find(rr(iice)>2.5 & rr(iice)<=1000);
     	rt(iice(ik))=rr(iice(ik)).*exp(0.317-0.565*log(rr(iice(ik)))-0.183*log(rr(iice(ik))).*log(rr(iice(ik))));
      	rq(iice(ik))=rr(iice(ik)).*exp(0.396-0.512*log(rr(iice(ik)))-0.180*log(rr(iice(ik))).*log(rr(iice(ik))));
 
        zot=zoq;                               % Dalton numbers that closely approximate COARE 3.0
    cdhf=von./(log(zu./zo)-psiu_26(zu./L));
    cqhf=von.*fdg./(log(zq./zoq)-psit_26(zq./L));
    cthf=von.*fdg./(log(zt./zot)-psit_26(zt./L));
    usr=ut.*cdhf;
    qsr=-(dq-wetc.*dter.*jcool).*cqhf;
    tsr=-(dt-dter.*jcool).*cthf;
    tvsr=tsr+0.61*ta.*qsr;
    tssr=tsr+0.51*ta.*qsr;
    Bf=-grav./ta.*usr.*tvsr;
    ug=0.2*ones(N,1);
    k=find(Bf>0); 
    if length(zi)==1;
        ug(k)=Beta*(Bf(k).*zi).^.333; clear k;
    else
        ug(k)=Beta*(Bf(k).*zi(k)).^.333; clear k;
    end
    ut=sqrt(du.^2+ug.^2);
    gf=ut./du;
    hsb=-rhoa*cpa.*usr.*tsr;
    hlb=-rhoa.*Le.*usr.*qsr;
    qout=Rnl+hsb+hlb;
    dels=Rns.*(0.065+11*tkt-6.6e-5./tkt.*(1-exp(-tkt/8.0e-4)));
    qcol=qout-dels;
    alq=Al.*qcol+be.*hlb.*cpw./Le;
    xlamx=6.0*ones(N,1);
    tkt=min(0.01, xlamx.*visw./(sqrt(rhoa./rhow).*usr));
    k=find(alq>0);
    xlamx(k)=6./(1+(bigc(k).*alq(k)./usr(k).^4).^0.75).^0.333;
    tkt(k)=xlamx(k).*visw./(sqrt(rhoa(k)./rhow).*usr(k)); clear k;
    dter=qcol.*tkt./tcw;
    dqer=wetc.*dter;
    Rnl=0.97*(5.67e-8*(ts-dter.*jcool+tdk).^4-Rl); % update dter
    if i==1; % save first iteration solution for case of zetu>50;
        usr50=usr(k50);tsr50=tsr(k50);qsr50=qsr(k50);L50=L(k50);
        zet50=zet(k50);dter50=dter(k50);dqer50=dqer(k50);tkt50=tkt(k50);
    end
    u10N = usr./von./gf.*log(10./zo);
    charnC=a1*u10N+a2;
    k=find(u10N>umax);
    charnC(k)=a1*umax+a2;
    charn=charnC;
     zoS=sigH.*Ad.*(usr./cp).^Bd;%-0.11*visa./usr;
    charnS=zoS.*grav./usr./usr;
     ii=find(~isnan(cp));charn(ii)=charnS(ii);
end

% insert first iteration solution for case with zetu>50
usr(k50)=usr50;tsr(k50)=tsr50;qsr(k50)=qsr50;L(k50)=L50;
zet(k50)=zet50;dter(k50)=dter50;dqer(k50)=dqer50;tkt(k50)=tkt50;

%****************  compute fluxes  ********************************************
tau=rhoa.*usr.*usr./gf;      % wind stress
hsb=-rhoa.*cpa.*usr.*tsr;     % sensible heat flux
hlb=-rhoa.*Le.*usr.*qsr;      % latent heat flux
hbb=-rhoa.*cpa.*usr.*tvsr;    % buoyancy flux
hsbb=-rhoa.*cpa.*usr.*tssr;   % sonic heat flux
wbar=1.61*hlb./Le./(1+1.61*Q)./rhoa+hsb./rhoa./cpa./ta;
hlwebb=rhoa.*wbar.*Q.*Le;
Evap=1000*hlb./Le./1000*3600;   %mm/hour

%*****  compute transfer coeffs relative to ut @ meas. ht  ********************
Cd=tau./rhoa./ut./max(.1,du);
Ch=-usr.*tsr./ut./(dt-dter.*jcool);
Ce=-usr.*qsr./(dq-dqer.*jcool)./ut;

%***  compute 10-m neutral coeff relative to ut (output if needed) ************
Cdn_10=1000*von.^2./log(10./zo).^2;
Chn_10=1000*von.^2.*fdg./log(10./zo)./log(10./zot);
Cen_10=1000*von.^2.*fdg./log(10./zo)./log(10./zoq);

%***  compute 10-m neutral coeff relative to ut (output if needed) ************
%  Find the stability functions
%*********************************
%zrf_u=10;             %User defined reference heights
%zrf_t=10;
%zrf_q=10;
psi=psiu_26(zu./L);
psi10=psiu_26(10./L);
psirf=psiu_26(zrf_u./L);
psiT=psit_26(zt./L);
psi10T=psit_26(10./L);
psirfT=psit_26(zrf_t./L);
psirfQ=psit_26(zrf_q./L);
gf=ut./du;

%*********************************************************
%  Determine the wind speeds relative to ocean surface
%  Note that usr is the friction velocity that includes 
%  gustiness usr = sqrt(Cd) S, which is equation (18) in
%  Fairall et al. (1996)
%*********************************************************
S = ut;
U = du;
S10 = S + usr./von.*(log(10./zu)-psi10+psi);
U10 = S10./gf;
% or U10 = U + usr./von./gf.*(log(10/zu)-psi10+psi);
Urf = U + usr./von./gf.*(log(zrf_u./zu)-psirf+psi);
UN = U + psi.*usr/von./gf;
U10N = U10 + psi10.*usr/von./gf;
UrfN = Urf + psirf.*usr/von./gf;

UN2 = usr/von./gf.*log(zu./zo);
U10N2 = usr./von./gf.*log(10./zo);
UrfN2  = usr./von./gf.*log(zrf_u./zo);

%******** rain heat flux (save to use if desired) *****************************
if isnan(rain(1));
    RF=zeros(size(usr));
else
    dwat=2.11e-5*((t+tdk)./tdk).^1.94; %! water vapour diffusivity
    dtmp=(1. + 3.309e-3*t - 1.44e-6.*t.*t).*0.02411./(rhoa.*cpa); %! heat diffusivity
    dqs_dt=Q.*Le./(Rgas.*(t+tdk).^2); %! Clausius-Clapeyron
    alfac= 1./(1+0.622*(dqs_dt.*Le.*dwat)./(cpa.*dtmp)); %! wet bulb factor
    RF= rain.*alfac.*cpw.*((ts-t-dter.*jcool)+(Qs-Q-dqer.*jcool).*Le./cpa)./3600;
end

lapse=grav/cpa;
SST=ts-dter.*jcool;

T = t;
%[size(-psi10T+psiT + lapse*(zt-10))]
T10 = T + tsr./von.*(log(10./zt)-psi10T+psiT) + lapse.*(zt-10);
Trf = T + tsr./von.*(log(zrf_t./zt)-psirfT+psiT) + lapse.*(zt-zrf_t);
TN = T + psiT.*tsr/von;
T10N = T10 + psi10T.*tsr/von;
TrfN = Trf + psirfT.*tsr/von;

TN2 = SST + tsr/von.*log(zt./zot)-lapse.*zt;
T10N2 = SST + tsr/von.*log(10./zot)-lapse.*10;
TrfN2 = SST + tsr/von.*log(zrf_t./zot)-lapse.*zrf_t;

dqer=wetc.*dter.*jcool;
SSQ=Qs-dqer;
SSQ=SSQ*1000;
Q=Q*1000;
qsr=qsr*1000;
Q10 = Q + qsr./von.*(log(10./zq)-psi10T+psiT);
Qrf = Q + qsr./von.*(log(zrf_q./zq)-psirfQ+psiT);
QN = Q + psiT.*qsr/von./sqrt(gf);
Q10N = Q10 + psi10T.*qsr/von;
QrfN = Qrf + psirfQ.*qsr/von;

QN2 = SSQ + qsr/von.*log(zq./zoq);
Q10N2 = SSQ + qsr/von.*log(10./zoq);
QrfN2 = SSQ + qsr/von.*log(zrf_q./zoq);
RHrf=RHcalc(Trf,P,Qrf/1000,Tf);
RH10=RHcalc(T10,P,Q10/1000,Tf);

%%%%%%%%%%%%%%%%%  Other wave breaking statistics fromn Banner-Morison wave
%%%%%%%%%%%%%%%%%  model

Whf=1.6e-3*U10N.^1.1./sqrt(cp./U10N);      %  whitecap fraction
ii=find(isnan(cp));Whf(ii)=7.3E-4*(U10N(ii)-2).^1.43;jj=find(U10N(ii)<2.1);Whf(ii(jj))=1e-5;
Edis=0.095*rhoa.*U10N.*usr.^2;             %  energy dissipation rate from breaking waves W/m^2
Whf(iice)=0;Edis(iice)=0;

%****************  output  ****************************************************
dterx=dter.*jcool;dqerx=dqer.*jcool;
A=[usr tau hsb hlb hbb hsbb hlwebb tsr qsr zo  zot zoq Cd Ch Ce  L zet dterx dqerx tkt Urf Trf Qrf RHrf UrfN Rnl Le rhoa UN U10 U10N Cdn_10 Chn_10 Cen_10 RF Qs Evap T10 Q10 RH10 ug Whf Edis];
%   1   2   3   4   5   6    7      8   9  10  11  12  13 14 15  16  17   18   19   20  21  22  23   24   25 26  27  28  29  30  31     32     33   34    35 36  37  38   39  40  41  42  43
end
%------------------------------------------------------------------------------
function psi=psit_26(zet)
% computes temperature structure function
dzet=min(50,0.35*zet); % stable
psi=-((1+0.6667*zet).^1.5+0.6667*(zet-14.28).*exp(-dzet)+8.525);
k=find(zet<0); % unstable
x=(1-15*zet(k)).^0.5;
psik=2*log((1+x)./2);
x=(1-34.15*zet(k)).^0.3333;
psic=1.5*log((1+x+x.^2)./3)-sqrt(3)*atan((1+2*x)./sqrt(3))+4*atan(1)./sqrt(3);
f=zet(k).^2./(1+zet(k).^2);
psi(k)=(1-f).*psik+f.*psic;
end
%------------------------------------------------------------------------------
function psi=psiu_26(zet)
% computes velocity structure function
dzet=min(50,0.35*zet); % stable
a=0.7;
b=3/4;
c=5;
d=0.35;
psi=-(a*zet+b*(zet-c/d).*exp(-dzet)+b*c/d);
k=find(zet<0); % unstable
x=(1-15*zet(k)).^0.25;
psik=2*log((1+x)/2)+log((1+x.*x)/2)-2*atan(x)+2*atan(1);
x=(1-10.15*zet(k)).^0.3333;
psic=1.5*log((1+x+x.^2)/3)-sqrt(3)*atan((1+2*x)./sqrt(3))+4*atan(1)./sqrt(3);
f=zet(k).^2./(1+zet(k).^2);
psi(k)=(1-f).*psik+f.*psic;
end
%------------------------------------------------------------------------------
function psi=psiu_40(zet)
% computes velocity structure function
dzet=min(50,0.35*zet); % stable
a=1;
b=3/4;
c=5;
d=0.35;
psi=-(a*zet+b*(zet-c/d).*exp(-dzet)+b*c/d);
k=find(zet<0); % unstable
x=(1-18*zet(k)).^0.25;
psik=2*log((1+x)/2)+log((1+x.*x)/2)-2*atan(x)+2*atan(1);
x=(1-10*zet(k)).^0.3333;
psic=1.5*log((1+x+x.^2)/3)-sqrt(3)*atan((1+2*x)./sqrt(3))+4*atan(1)./sqrt(3);
f=zet(k).^2./(1+zet(k).^2);
psi(k)=(1-f).*psik+f.*psic;
end
%------------------------------------------------------------------------------
function exx=bucksat(T,P,Tf)
% computes saturation vapor pressure [mb]
% given T [degC] and P [mb] Tf is freezing pt 
exx=6.1121.*exp(17.502.*T./(T+240.97)).*(1.0007+3.46e-6.*P);
ii=find(T<Tf);
exx(ii)=(1.0003+4.18e-6*P(ii)).*6.1115.*exp(22.452.*T(ii)./(T(ii)+272.55));%vapor pressure ice
end
%------------------------------------------------------------------------------
function qs=qsat26sea(T,P,Ss,Tf)
% computes surface saturation specific humidity [g/kg]
% given T [degC] and P [mb]
ex=bucksat(T,P,Tf);
fs=1-0.02*Ss/35;% reduction sea surface vapor pressure by salinity
es=fs.*ex; 
qs=622*es./(P-0.378*es);
end
%------------------------------------------------------------------------------
function [q,em]=qsat26air(T,P,rh)
% computes saturation specific humidity [g/kg]
% given T [degC] and P [mb]
Tf=0;%assumes relative humidity for pure water
es=bucksat(T,P,Tf);
em=0.01*rh.*es;
q=622*em./(P-0.378*em);
end
%------------------------------------------------------------------------------
function g=grv(lat)
% computes g [m/sec^2] given lat in deg
gamma=9.7803267715;
c1=0.0052790414;
c2=0.0000232718;
c3=0.0000001262;
c4=0.0000000007;
phi=lat*pi/180;
x=sin(phi);
g=gamma*(1+c1*x.^2+c2*x.^4+c3*x.^6+c4*x.^8);
end
%------------------------------------------------------------------------------
function RHrf=RHcalc(T,P,Q,Tf)
% computes relative humidity given T,P, & Q
es=6.1121.*exp(17.502.*T./(T+240.97)).*(1.0007+3.46e-6.*P);
ii=find(T<Tf);%ice case
es(ii)=6.1115.*exp(22.452.*T(ii)./(T(ii)+272.55)).*(1.0003+4.18e-6*P(ii));
em=Q.*P./(0.378.*Q+0.622);
RHrf=100*em./es;
end


