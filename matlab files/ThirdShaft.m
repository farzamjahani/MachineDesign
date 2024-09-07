function [Reactions3,W_a3,M_tot3]=ThirdShaft(H_nom,omega_m,Gear,EL3,phi_n,mG)

    phi_t=atand(tand(phi_n)*secd(Gear(2,3)));
    % Gear Force Analysis
    W_t=60000*H_nom/(Gear(2,2)*pi*omega_m/(prod(mG)));
    W_r=W_t*tand(phi_t);
    W_a3=W_t*tand(Gear(2,3));
    % reaction force of bearing B23 (second bearing of third shaft)
    F_zB23=((W_r*(EL3(2)-EL3(1))-W_a3*(Gear(2,2)/2)))/(EL3(3)-EL3(1));
    F_yB23=(W_t*(EL3(2)-EL3(1)))/(EL3(3)-EL3(1));
    F_B23=sqrt(F_zB23^2+F_yB23^2);
    % reaction force of bearing B13 (first bearing of third shaft)
    F_zB13=W_r-F_zB23;
    F_yB13=W_t-F_yB23;
    F_B13=sqrt(F_zB13^2+F_yB13^2);
    Reactions3=[F_B13 F_B23];   % concatenating of reaction force in 3rd shaft
    % domain along third shaft
    L1=0:0.01:EL3(1);
    L2=EL3(1):0.01:EL3(2);
    L3=EL3(2):0.01:EL3(3);
    L4=EL3(3):0.01:EL3(4);
    L=[L1 L2 L3 L4];
    % bending moment in Y-direction along 3rd shaft
    M1_y=zeros(size(L1));
    M2_y=-F_zB13*(L2-EL3(1));
    M3_y=-F_zB13*(L3-EL3(1))+W_r*(L3-EL3(2))+W_a3*(Gear(2,2)/2);
    M4_y=zeros(size(L4));
    M_y3=[M1_y M2_y M3_y M4_y]/1000;
    % bending moment in Z-direction along 3rd shaft
    M1_z=zeros(size(L1));
    M2_z=F_yB13*(L2-EL3(1));
    M3_z=F_yB13*(L3-EL3(1))-W_t*(L3-EL3(2));
    M4_z=zeros(size(L4));
    M_z3=[M1_z M2_z M3_z M4_z]/1000;
    
    M_tot3=sqrt(M_y3.^2+M_z3.^2);     % magnitude of bending moment
    
    figure(7);
    plot(L,M_y3,'k','LineWidth',1.5)
    title('Bending Moment in Third Shaft')
    xlabel('x (mm)')
    ylabel('M_y3 (N.m)')
    grid on
    
    figure(8);
    plot(L,M_z3,'k','LineWidth',1.5)
    title('Bending Moment in Third Shaft')
    xlabel('x (mm)')
    ylabel('M_z3 (N.m)')
    grid on
    
    figure(9);
    plot(L,M_tot3,'k','LineWidth',1.5)
    title('Bending Moment in Third Shaft')
    xlabel('x (mm)')
    ylabel('M_{tot3} (N.m)')
    grid on
    
    

end