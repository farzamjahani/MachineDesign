function [Reactions2,W_a2,M_tot2]=SecondShaft(H_nom,omega_m,Gear,EL2,phi_n,mG)

    phi_t=atand(tand(phi_n)*secd(Gear(:,3)));
    % Gears Force Analysis (GB & GC)
    W_t=(60000*H_nom/(Gear(1,2)*pi*omega_m/(1*mG(1))))*[1 (Gear(1,2)/Gear(2,1))];
    W_r=W_t.*tand(phi_t);
    W_a2=W_t.*tand(Gear(:,3));
    % Second Bearing Force Analysis (B22)
    F_zB22=(W_a2(1)*(Gear(1,2)/2)-W_r(1)*(EL2(2)-EL2(1))-W_a2(2)*(Gear(2,1)/2)-W_r(2)*(EL2(3)-EL2(1)))/(EL2(4)-EL2(1));
    F_yB22=(W_t(1)*(EL2(2)-EL2(1))-W_t(2)*(EL2(3)-EL2(1)))/(EL2(4)-EL2(1));
    F_B22=sqrt(F_zB22^2+F_yB22^2);
    % First Bearing Force Analysis (B12)
    F_zB12=-W_r(1)-W_r(2)-F_zB22;
    F_yB12=W_t(1)-W_t(2)-F_yB22;
    F_B12=sqrt(F_zB12^2+F_yB12^2);
    Reactions2=[F_B12 F_B22];   % concatenating of reaction force in 2nd shaft
    % domains along second shaft
    L1=0:0.01:EL2(1);
    L2=EL2(1):0.01:EL2(2);
    L3=EL2(2):0.01:EL2(3);
    L4=EL2(3):0.01:EL2(4);
    L5=EL2(4):0.01:EL2(5);
    L=[L1 L2 L3 L4 L5];
    % bending moment in Y-direction along 2nd shaft
    M1_y=zeros(size(L1));
    M2_y=-F_zB12*(L2-EL2(1));
    M3_y=-F_zB12*(L3-EL2(1))-W_r(1)*(L3-EL2(2))-W_a2(1)*(Gear(1,2)/2);
    M4_y=-F_zB12*(L4-EL2(1))-W_r(1)*(L4-EL2(2))-W_a2(1)*(Gear(1,2)/2)-W_r(2)*(L4-EL2(3))+W_a2(2)*(Gear(2,1)/2);
    M5_y=zeros(size(L5));
    M_y2=[M1_y M2_y M3_y M4_y M5_y]/1000;
    % bending moment in Z-direction along 2nd shaft
    M1_z=zeros(size(L1));
    M2_z=F_yB12*(L2-EL2(1));
    M3_z=F_yB12*(L3-EL2(1))-W_t(1)*(L3-EL2(2));
    M4_z=F_yB12*(L4-EL2(1))-W_t(1)*(L4-EL2(2))+W_t(2)*(L4-EL2(3));
    M5_z=zeros(size(L5));
    M_z2=[M1_z M2_z M3_z M4_z M5_z]/1000;
    
    M_tot2=sqrt(M_y2.^2+M_z2.^2);	% magnitude of bending moment
    
    figure(4);
    plot(L,M_y2,'k','LineWidth',1.5)
    title('Bending Moment in Second Shaft')
    xlabel('x (mm)')
    ylabel('M_y2 (N.m)')
    grid on
    
    figure(5);
    plot(L,M_z2,'k','LineWidth',1.5)
    title('Bending Moment in Second Shaft')
    xlabel('x (mm)')
    ylabel('M_z2 (N.m)')
    grid on
    
    figure(6);
    plot(L,M_tot2,'k','LineWidth',1.5)
    title('Bending Moment in Second Shaft')
    xlabel('x (mm)')
    ylabel('M_{tot2} (N.m)')
    grid on

end