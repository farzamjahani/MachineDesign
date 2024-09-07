function [Teeth,mG,b,C_G,m_n,Dia,Gear]=HelicalGear(H_nom,omega_m,n_d)
 
 phi_n=20; % normal pressure angle
 k=1; % full depth teeth assumption
 e=11.9; % total gearbox transmission ratio
 HB=325; % Brinell hardness of core of the gears
 mG_raw=sqrt(e); % initial transmission ratio as an equal-reduction assumption
 Omega=omega_m; % TR_BP is transmission ratio of belt-pulley subsystem
 
 Stage=['1st';'2nd'];
 GM=load('Geo&ModFactors.txt'); % table that contains geometry and modifying factor of helical gears
 
 Flag=1;
 for i=[1,2] % for-loop for both two reduction-stage designing
    N=10*300*8*60*Omega; % desired life for pinion
    % ------ "first principle": intercutting avoidance
    while true
        psi=input(['Enter helix angle of "' Stage(i,:) '" stage of gears subsystem in degree \n']);
        phi_t=atand(tand(phi_n)*secd(psi)); % transverse pressure angle of teeth
        Np_min=(2*k*cosd(psi)*(mG_raw+sqrt(mG_raw^2+(1+2*mG_raw)*(sind(phi_t))^2)))/((1+2*mG_raw)*sind(phi_t)^2);
        disp(['minimum No. teeth of pinion in the "' Stage(i,:) '" stage is Np_min = ' num2str(Np_min)])
        Np=input('Choose an integer Np that must be as small as possible :\n');
        Ng_raw=Np*mG_raw; % initial No. gear teeth
        disp(['Number of teeth of gear in the "' Stage(i,:) '" stage is Ng = ' num2str(Ng_raw)])
        Ng=input('Please round Ng to an integer number :\n');
        mG=Ng/Np; % final trans
 
        if gcd(Np,Ng)~=1 % hunter Analysis 
            warning(['The hunter phenomenon will occur in the "' Stage(i,:) '" of gear subsystem!'])
            ANS1=input(['Do you accept hunter phenomenon in the "' Stage(i,:) '" of gears subsystem? 1.Yes 2.No \n']);
            if ANS1==2
                continue
            end
        else
            disp(['The "' Stage(i,:) '" stage of gear subsystem is free from hunter phenomenon'])
        end
        if (mG^2<(e-0.02*e) || mG^2>(e+0.02*e)) % Local checking of transmission ratio
            warning(['The transmission ratio in the "' Stage(i,:) ...
            '" of gears subsystems out of employer desired (with equal-stage-reduction consideration)!'])
            ANS2=input('Is it ok? 1.Yes, it will checked later and maybe there is no problem! 2.No');
            if ANS2==2
                continue
            else
                break
            end
        else
            disp(['The transmission ratio in the "' Stage(i,:) '" of gears subsystem is in allowable range that employer desired'])
            break
        end
    end

    % ------ "second principle": bending fatigue analysis
    S_t=(0.703*HB+113)*10^6; % allowable bending stress number
    S_F=1; % applying safety factor on load, S_F is one;
    Y_theta=1; % assuming temperature in the system is under 120 degree Celsius
    R=input('Please Enter the "percentage" of Reliability: \n');
    R=R/100;
    if R>0.5 && R<0.99
        Y_z=0.658-0.0759*log(1-R);
    elseif R>=0.99 && R<=0.9999
        Y_z=0.5-0.109*log(1-R);
    end
 
    Y_N=(1.6831*N^(-0.0323)+1.3558*N^(-0.0178))/2; % cycle factor with midrange consideration!
    sigma_b_all=(S_t/S_F)*(Y_N/(Y_z*Y_theta)); % allowable bending stress in teeth
    K_o=1.5; % overload factor
    K_s=1; % not knowing any available data, assume "size factor" equal to one
    K_B=1; % not knowing shaft dia. at gear seat, assume "rim-thickness factor" is one!
    IND=find(GM(:,1)==psi);
    for j=2:7
        if GM(1,j)<=Np && GM(1,j+1)>Np
            J_prime=GM(IND,j)+(Np-GM(1,j))*(GM(IND,j+1)-GM(IND,j))/(GM(1,j+1)-GM(1,j)); % raw geometry facto
        end
        if GM(1,j+5)<=Ng && GM(1,j+6)>Ng && j>=4
            K_J=GM(IND,j+5)+(Ng-GM(1,j+5))*(GM(IND,j+6)-GM(IND,j+5))/(GM(1,j+6)-GM(1,j+5)); % modifying factor
        end
    end
    Y_J=K_J*J_prime; % finalized geometry factor
    x=input(['Enter width ratio (x) of the "' Stage(i,:) '" of gears subsystem: \n']); % width ratio of helical gear
    T=(60*H_nom)/(2*pi*Omega); % torque in shaft containing pinion (N.m)
 
    FLAG=0;
    while true % check for variation of load-distribution factor with width of gears
        if FLAG==0
            K_H=1.6; % initial load-distribution factor
        else
            if b>=50
                K_H_iterate=1.7;
            else
                K_H_iterate=1.6;
            end
            if K_H_iterate==K_H
                break
            else
                K_H=K_H_iterate;
            end
        end
        m_n_guess=input(['Enter initial guess for normal module in mm in the "' Stage(i,:) '" of gears subsystem: \n']);

        flag=1;
        while flag~=0 % check for variation of dynamic factor with module of gears
            V=pi*m_n_guess*Np*secd(psi)*Omega/60000; % velocity of pitch point (m/s)
            K_v=(3.56+sqrt(V))/3.56; % dynamic factor for hob-manufactured teeth
            m_n_minb=(cosd(psi))*(2*n_d*T*K_o*K_B*K_H*K_v*K_s/(x*Np*Y_J*sigma_b_all))^(1/3)*10^3; % min. normal module (mm)
            disp(['"BENDING", min. value for normal module from bending analysis is ' num2str(m_n_minb) ' mm'])
            m_nb=input('Please correct this module from Tab. 13-2 [1] and select a preferred one : \n'); % correcting the normal module
            if m_nb==m_n_guess
                flag=0;
            else
                m_n_guess=m_nb;
            end
        end
        % ------ "third principle": surf. contact fatigue analysis
        S_c=(2.41*HB+237)*10^6; % allowable contact stress number
        S_H=1; % applying safety factor on load, S_H is unity;
        Z_N=(2.466*N^(-0.056)+1.4488*N^(-0.023))/2; % cycle factor with midrange consideration!
        Z_R=1; % surface condition factor
        sigma_c_all=(S_c/S_H)*(Z_N*Z_R/(Y_z*Y_theta)); % allowable contact stress
        Z_E=187*10^3; % elastic factor
        m_n_it=m_nb; % m_n_it is a type of m_n_guess for iteration on surf. contact analysis
        while flag==0
            a=k*m_n_it; % estimated addendum (mm) 
            m_t=m_n_it*secd(psi); % estimated transverse module (mm)
            d=Np*m_t; % estimated diameter of pitch circle for pinion (mm) 
            pn=pi*m_n_it; % normal Circular Pitch (mm) 
            pN=pn*cosd(phi_n); % (mm) 
            rP=d/2; % radius of circular pitch for Pinion (mm) 
            rbP=rP*cosd(phi_t); % radius of base circle for Pinion (mm)
            rG=rP*mG; % radius of circular pitch for Gear (mm) 
            rbG=rG*cosd(phi_t); % radius of base circle for Gear (mm) 
            Z1=((rP+a)^2-rbP^2)^0.5;
            Z2=((rG+a)^2-rbG^2)^0.5;
            Z3=(rP+rG)*sind(phi_t); 
            Z=min(Z1,Z3)+min(Z2,Z3)-Z3; % length of the line of action
            mN=pN/(0.95*Z); % load-sharing ratio (mm) 
            Z_I=(cosd(phi_t)*sind(phi_t)*(mG))/(2*mN*(mG+1)); % geometry factor
            m_n_minc=(cosd(psi))*((2*n_d*T*Z_R*K_o*K_H*K_v*K_s*(Z_E/sigma_c_all)^2)/(Z_I*x*(Np^2)))^(1/3)*10^3;
            disp(['"SURF. CONTACT", min. value for normal module from surf. contact analysis is ' num2str(m_n_minc)])
            m_nc=input('"SURF. CONTACT", please correct this module from Tab. 13-2 and select a preferred one \n');
            if m_nc>m_nb
                if m_nc==m_n_it
                    flag=1;
                end
                m_n_it=m_nc;
                m_n_guess=m_nc;
                V=pi*m_n_guess*Np*secd(psi)*Omega/60000;
                K_v=(3.56+sqrt(V))/3.56;
            else
                flag=1;
            end
        end
        m_n=max(m_nb,m_nc); % finalizing the normal module
        C=(Np+Ng)*m_n*secd(psi)/2; % center Distance of gears
        b=x*m_n; % Final Width of Gears (mm) 
        disp(['Final width of the "' Stage(i,:) '" of gears subsystem is b = ' num2str(b) ' mm'])
        m_t=m_n*secd(psi); %Final transverse Module (mm) 
        d=Np*m_t; % Final Diameter of Pitch Circle for Pinion (mm) 
        D=Ng*m_t; % Diameter of Pitch Circle For Gear (mm) 
        ded=1.25*m_n; % Deddendum (mm) 
        dr=d-2*ded; % Diameter of Root Circle for Pinion (mm)
        Dr=D-2*ded; % Diameter of Root Circle for Gear (mm)
        FLAG=1;
    end
    OutPuts(i,:)={psi,Np,Ng,mG,b,C,d,D,dr,Dr,m_n}; % initial container for outputs
    if Flag~=0
        fprintf('\n"First stage" designing of gears successfully finished!\n"Second stage" designing of gears started\n')
        Omega=Omega/(mG(i));
    end
    Flag=0;
 end
 for i=[1,2] % for-loop of output data assignment
    [psi(i),Np(i),Ng(i),mG(i),b(i),C_G(i),d(i),D(i),dr(i),Dr(i),m_n(i)]=OutPuts{i,:};
 end
 
 if (prod(mG)<(e-0.02*e) || prod(mG)>(e+0.02*e)) % Final Check of total Transmission Ratio of gearbox
    warning('Total transmission ratio of "GEARBOX" e = %.3f is out of employer desired range!\n',prod(mG))
    warning('You must redesign one or both stage of reduction!')
 else
    fprintf('Total transmission ratio of "GEARBOX" e = %.3f is in allowable range that employer needs and there is no problem \n',prod(mG))
 end
 
 Teeth=[Np(1) Ng(1);Np(2) Ng(2)];
 Dia=[dr(1) Dr(1);dr(2) Dr(2)];
 Gear=[d(1) D(1) psi(1);d(2) D(2) psi(2)];
 
end
