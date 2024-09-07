function [Reactions1,W_a1,M_tot1]=FirstShaft(H_nom,omega_m,Gear,EL1,phi_n)

    phi_t=atand(tand(phi_n)*secd(Gear(1,3)));	% transvers pressure angle
    % Gear Force Analysis (GA)
    W_t=60000*H_nom/(Gear(1,1)*pi*omega_m);
    W_r=W_t*tand(phi_t);
    W_a1=W_t*tand(Gear(1,3));
    % reaction force of bearing B21 (second bearing of first shaft)
    F_zB21=(W_r*(EL1(2)-EL1(1))+W_a1*(Gear(1,1)/2))/(EL1(3)-EL1(1));
    F_yB21=(-W_t*(EL1(2)-EL1(1)))/(EL1(3)-EL1(1));
    F_B21=sqrt(F_zB21^2+F_yB21^2);
    % reaction force of bearing B11 (first bearing of first shaft)
    F_zB11=W_r-F_zB21;
    F_yB11=-W_t-F_yB21;
    F_B11=sqrt(F_zB11^2+F_yB11^2); 
    Reactions1=[F_B11 F_B21];   % concatenating of reaction force in 1st shaft
    % domains along first shaft
    L1=0:0.01:EL1(1);
    L2=EL1(1):0.01:EL1(2);
    L3=EL1(2):0.01:EL1(3);
    L4=EL1(3):0.01:EL1(4);
    L=[L1 L2 L3 L4];
    % bending moment in Y-direction along 1st shaft
    M1_y=zeros(size(L1));
    M2_y=-F_zB11*(L2-EL1(1));
    M3_y=-F_zB11*(L3-EL1(1))-W_a1*(Gear(1,1)/2)+W_r*(L3-EL1(2));
    M4_y=zeros(size(L4));
    M_y1=[M1_y M2_y M3_y M4_y]/1000;
    % bending moment in Z-direction along 1st shaft
    M1_z=zeros(size(L1));
    M2_z=F_yB11*(L2-EL1(1));
    M3_z=F_yB11*(L3-EL1(1))+W_t*(L3-EL1(2));
    M4_z=zeros(size(L4));
    M_z1=[M1_z M2_z M3_z M4_z]/1000;    % N.mm to N.m
    
    
    M_tot1=sqrt(M_y1.^2+M_z1.^2);	% magnitude of bending moment
    
    figure(1);
    plot(L,M_y1,'k','LineWidth',1.5)
    title('Bending Moment in First Shaft')
    xlabel('x (mm)')
    ylabel('M_y1 (N.m)')
    grid on
    
    figure(2);
    plot(L,M_z1,'k','LineWidth',1.5)
    title('Bending Moment in First Shaft')
    xlabel('x (mm)')
    ylabel('M_z1 (N.m)')
    grid on
    
    
    figure(3);
    plot(L,M_tot1,'k','LineWidth',1.5)
    title('Bending Moment in First Shaft')
    xlabel('x (mm)')
    ylabel('M_{tot1} (N.m)')
    grid on
    
end