%% gait.m -- Starter File for Gait Data Analysis
%% B17 Biomechanics -- Hilary Term 2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all

% Load data into arrays Xmk (markers) and Xfp (forceplate) for Subject X
% You must replace 'X' below with the code for your Subject (A, B, C, D)
mk = csvread('Subject-B-markers.csv');
fp = csvread('Subject-B-forceplate.csv');

% Find the number of rows and columns in the input files
[rmk,cmk] = size(mk);    [rfp,cfp] = size(fp);

% initialize arrays of datapt number in marker file & assign to an array
% note that ':' means 'all the rows (or columns)'
datapt = zeros(rmk,1);
datapt(:,1) = mk(:,1);

% Initialize arrays for marker data (right leg only)
sac = zeros(rmk,3);      % sacrum (SACR)
rasi = zeros(rmk,3);      % anterior superior iliac spine (ASIS)
rthi = zeros(rmk,3);      % thigh wand marker (THI) 
rkne = zeros(rmk,3);      % lateral femoral condyle (KNE)
rtib = zeros(rmk,3);      % tibia wand marker (TIB)
rank = zeros(rmk,3);      % later malleolus (LMA)
rhee = zeros(rmk,3);      % heel (HEE)
rtoe = zeros(rmk,3);      % 2nd metatarsal head (TOE)

% Assign xyz coordinates of markers to right side sacrum, asis, thigh, knee,
% ankle, heel, and toe arrays for right leg
sac(:,1:3) = mk(:,2:4);
rasi(:,1:3) = mk(:,8:10);             
rthi(:,1:3) = mk(:,29:31);      
rkne(:,1:3) = mk(:,32:34); 
rtib(:,1:3) = mk(:,35:37);
rank(:,1:3) = mk(:,38:40);      
rhee(:,1:3) = mk(:,41:43);      
rtoe(:,1:3) = mk(:,44:46);    

%% question 1

%load left leg 
%initialize
lasi = zeros(rmk,3);      % anterior superior iliac spine (ASIS)
lthi = zeros(rmk,3);      % thigh wand marker (THI) 
lkne = zeros(rmk,3);      % lateral femoral condyle (KNE)
ltib = zeros(rmk,3);      % tibia wand marker (TIB)
lank = zeros(rmk,3);      % later malleolus (LMA)
lhee = zeros(rmk,3);      % heel (HEE)
ltoe = zeros(rmk,3);      % 2nd metatarsal head (TOE)

%assign coordinates
lasi(:,1:3) = mk(:,5:7);             
lthi(:,1:3) = mk(:,11:13);      
lkne(:,1:3) = mk(:,14:16); 
ltib(:,1:3) = mk(:,17:19);
lank(:,1:3) = mk(:,20:22);      
lhee(:,1:3) = mk(:,23:25);      
ltoe(:,1:3) = mk(:,26:28); 


% Plot yz trajectories
     figure(1)
     plot(sac(:,2),sac(:,3))
     hold on
     text(sac(rmk,2),sac(rmk,3),'SACRUM')
     plot(rasi(:,2),rasi(:,3))
     text(rasi(rmk,2),rasi(rmk,3),'ILIAC')
% continue with rest of markers     
     plot(rthi(:,2),rthi(:,3))
     text(rthi(rmk,2),rthi(rmk,3),'THIGH')
     plot(rkne(:,2),rkne(:,3))
     text(rkne(rmk,2),rkne(rmk,3),'KNEE')
     plot(rtib(:,2),rtib(:,3))
     text(rtib(rmk,2),rtib(rmk,3),'TIBIA')
     plot(rank(:,2),rank(:,3))
     text(rank(rmk,2),rank(rmk,3),'ANKLE')
     plot(rhee(:,2),rhee(:,3))
     text(rhee(rmk,2),rhee(rmk,3),'HEEL')
     plot(rtoe(:,2),rtoe(:,3))
     text(rtoe(rmk,2),rtoe(rmk,3),'TOE')

     %AXIS
     xlabel('Y (Posterior - Anterior)')
     ylabel('Z (Inferior - Superior)')
     title('Subject B Marker Trajectories')
     %axis('equal')
     print -depsc epsFigq1
     hold off

 %% question 2
 %find COM - midpoint between LASI and RASI
 COM = mean(0.5*(rasi+lasi));
 
%find and plot velocity
%use SACR point 

%time as 1/f
t=0.01;
%y values for sac
[m,n]=size(sac);
y=sac(2:m,2);
yn=sac(1:m-1,2);
%velocity
v=(y-yn)/t;
%datapoint number
[b,~]=size(datapt);
d=datapt(2:b,:);

%find mean velocity
V=mean(v);
sz=size(d);
p=ones(sz);
Vp=V*p;

%plot
figure(2)
plot(d,v)
ylim([1250,1800])
xlabel('Data Point Number')
ylabel('Velocity (mm/s)')
title('Velocity Against Data Point Number')
hold on 
plot(d,Vp)
legend('Instantaneous Velocity','Average Velocity')

print -depsc epsFigq2


%% question 3
%calculate lengths of shank
%right leg in 2d
rank2=rank(:,2:3);
rkne2=rkne(:,2:3);
rminus2=rank2-rkne2;
rsquare2=rminus2.^2;
rsum2=sum(rsquare2,2);
rlength2=rsum2.^0.5;

%left leg in 2d
lank2=lank(:,2:3);
lkne2=lkne(:,2:3);
lminus2=lank2-lkne2;
lsquare2=lminus2.^2;
lsum2=sum(lsquare2,2);
llength2=lsum2.^0.5;

%right leg in 3d
rminus=rank-rkne;
rsquare=rminus.^2;
rsum=sum(rsquare,2);
rlength=rsum.^0.5;

%left leg in 3d
lminus=lank-lkne;
lsquare=lminus.^2;
lsum=sum(lsquare,2);
llength=lsum.^0.5;

%plot
figure(3)
plot(datapt,rlength)
ylim([400,445])
xlabel('Data Point Number')
ylabel('Length (mm/s)')
title('Length of Shank')
hold on 
plot(datapt,llength)
plot(datapt,rlength2)
plot(datapt,llength2)
legend({'Right Leg in 3D','Left Leg in 3D','Right Leg in 2D','Left leg in 2D'},'Location','southwest')

print -depsc epsFigq3

%% question 4
%a) vertical marker of ank, hee & toe
rankz=rank(:,3);
rheez=rhee(:,3);
rtoez=rtoe(:,3);

%plot
figure(4)
plot(datapt,rankz)
hold on
plot(datapt,rheez)
plot(datapt,rtoez)
xlabel('Data Point Number')
ylabel('Vertical Position/mm')
title('Vertical Position against Data Point Number')


%plot labels
xline(220,'--r')
xline(232,'--g')
xline(258,'--c')
xline(278,'--')

legend('Ankle','Heel','Toe','HS','FF','HO','TO')

print -depsc epsFigq4a

%b) vertical component of Fz
Fz=fp(:,4);

%data points for force plate
Fdatapt = zeros(rfp,1);
Fdatapt(:,1) = fp(:,1);

%plot
figure(5)
plot(Fdatapt,Fz)
xlabel('Data Point Number')
ylabel('Force/N')
title('Vertical Component of Ground Reaction Force')

print -depsc epsFigq4b



%% question 6
%find angle of knee 
%get yz
yzthi = rthi(:,2:3);
yzkne = rkne(:,2:3);
yztib = rtib(:,2:3);

%contstruct vectors
BA = yzkne-yzthi;
BC = yztib-yzkne;

%find dot product
dot=(BA(:,1).*BC(:,1)+BA(:,2).*BC(:,2));

%find magnitudes
BAsquare=BA.^2;
BAsum=sum(BAsquare,2);
BAmag=BAsum.^0.5;

BCsquare=BC.^2;
BCsum=sum(BCsquare,2);
BCmag=BCsum.^0.5;

%combine 
mag=BAmag.*BCmag;
comb=dot./mag;

%find angle
theta=acosd(comb);

%plot
figure(6)
plot(datapt,theta)
hold on
ylim([-15, 70])
xlabel('Data Point Number')
ylabel('Flexion/degrees')
title('2D Angle of Knee Joint')

%plot heelstrikes
plot(220,19.6699,'d')
plot(311,22.4091,'d')


%plot toe offs
plot(184,32.1078,'o')
plot(278,29.4181,'o')


print -depsc epsFigq6





