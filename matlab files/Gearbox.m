
%     "GEARBOX"
 clc;  
fprintf('\nWelcome to GEARBOX design script!\n\n')
FID=fopen('Notes.txt');
while ~feof(FID)
    TL=fgetl(FID);
    disp(TL)
end
fclose(FID);
fprintf('\n\nAre you ready? Press "Enter" key')
pause;
clc
clear
close all

% overal data of project
H_nom=1360;         % nominal power of electromotor     (W)
omega_m=1450;       % nominal RPM of electromotor       (rpm)
phi_n=20;           % normal pressure angle of teeth    (degree)
n_d=3;              % overall safety factor
HB=325;             % Brinell hardness of core of the gears

% data from outside of script
EL_tot=load('EstimatedLength.txt');  % estimated length of shafts
EL1=EL_tot(1,:);
EL2=EL_tot(2,:);

EL3=EL_tot(3,:);
shaft_disp_total=load('shaft_displacement.txt');    % component displacement

Fid=fopen('ComponentNames.txt');    % a file encompassing all component names in the gearbox
i=0;
while ~feof(Fid)
    i=i+1;
    TL=fgetl(Fid);
    COM(i)=string(TL);
end
fclose(Fid);
COMs=split(COM');

% % ----- 1st phase: Belt-Pulley subsystem
% fprintf('*** You got into "BELT-PULLEY" design scope ***\n')
% [P_width,F1,F2,alpha,TR_BP]=BeltPulley(H_nom,omega_m,n_d);

% ----- 2nd phase: HelicalGears subsystem
fprintf('\n\n*** You got into "HELICAL GEAR" design scope ***\n\n')
[Teeth,mG,b,C_G,m_n,Dia,Gear]=HelicalGear(H_nom,omega_m,n_d);



% ----- 4th phase: Shaft-Layouts
fprintf('\n\n*** You got into "SHAFT" design scope ***\n\n')
[Reactions1,W_a1,M_tot1]=FirstShaft(H_nom,omega_m,Gear,EL1,phi_n);
[Reactions2,W_a2,M_tot2]=SecondShaft(H_nom,omega_m,Gear,EL2,phi_n,mG);
[Reactions3,W_a3,M_tot3]=ThirdShaft(H_nom,omega_m,Gear,EL3,phi_n,mG);


T_maximum=(H_nom/omega_m)*(30/pi)*[1 mG(1) prod(mG)];   % max. torque in each shaft
M={M_tot1;M_tot2;M_tot3};

D_shaft_total=zeros(3,5);       % a container for whole diameters of the three shafts
bearing_seat_dia=zeros(3,2);    % diameters of all six bearing seats
for v=1:3   % iteration on three shafts
    shaft_disp=shaft_disp_total(8*(v-1)+1:8*v,:);     % characteristic table containing shafts dimensions
    D_shaft_total(v,:)=ShaftDiameters(shaft_disp,M{v},T_maximum(v),EL_tot(v,:),n_d,COMs(:,v)');
    bearing_seat_dia(v,:)=D_shaft_total(v,[max([3-v,1]),min([7-v,5])]); % assignment of bearing seat diameter
end
D_G=D_shaft_total([5 6 10 11]); % diameters of all four gear seats
DIA=Dia';
for i=1:4   % rim-thickness analysis
    mB(i)=(DIA(i)-D_G(i))/(2.25*m_n(ceil(i/2))); % backup ratio (add+ded=2.25*m_n)
    if mB<1.2
        fprintf('\nIn the gear "G%s" there is rim-thickness!\nrim-thicknees factor is K_b = %.2f\n',64+i,1.6*log(2.242/mB(i)))
    else
        disp('All gears in the gearbox are free from rim-thickness')
    end
end

%% % ----- 5th phase: bearing-selection
fprintf('\n\n*** You got into "BEARING" design scope ***\n\n')
omega=omega_m./[1,mG(1),prod(mG)];    % concatenation of shafts RPM
reaction=[Reactions1;Reactions2;Reactions3]/1000;       % reactoin force on bearings in kN
Flag_axial=[1 0;1 0;0 1];                               % guide for indication of axial-supporting bearings
F_a=abs([W_a1, W_a2(1)-W_a2(2), W_a3])/1000;            % concatenation of axial force of shafts in kN

[BearingType,R]=BearingSelection(omega,reaction,bearing_seat_dia,Flag_axial,F_a);

fprintf('\n\n*** Design process of all five phase finished ***\n\n now, you can refine your gearbox-designing\n\n')

