function D_shaft=ShaftDiameters(shaft_disp,M_a,T_maximum,EL,n_d,COM)

    Su=660*10^6;                % ultimate tension stress of shaft material
    Sy=470*10^6;                % yield strength of shaft material
    Mm=0;                       % Midrange Moment   
    T_minimum=0;                % Minimum Torque
    n=numel(shaft_disp(:,1));   % No. initial stress concentration point on each shaft
    T_min=zeros(1,n);           % container for lower bound of torque
    T_max=zeros(1,n);           % container for upper bound of torque
    NUM=importdata('Numbers.txt');  % number file
    SNO=1+(shaft_disp(4,2)==3)+(shaft_disp(7,2)==2);  % detection which shaft is in progress 
    SHAFT={'FIRST','SECOND','THIRD'};   % shaft ordering
    
    % ----- torque assignment   
    for i=1:n
            T_min(i)=0;   
    end        
    for i=1:n 
        if shaft_disp(i,4)==0        
            T_max(i)=0;
        elseif shaft_disp(i,4)==1        
            T_max(i)=T_maximum;        
            T_min(i)=T_minimum;     
        end
    end
    T_a=(T_max+T_min)./2;	% amplitude of Torque (N.m) 
    T_m=(T_max-T_min)./2;	% midrange of Torque (N.m) 
    
    % ----- endurance limit calculations
    Se_prime=min([0.5*Su,700*10^6]); 
    a=4.51;     % a-Factor for Marin surf. mod. factor (machined)
    b=-0.265;   % b-exponent for Marin surf. mod. factor (machined)
    ka=a*(Su/(10^6))^(b);	% surface modification factor   
    kc=1;       % loading factor (combination of loading)
    TABLE1=load('table_6-5.txt');   % data of tab. 6-5 [1]
    OPT=num2cell(sort([50,95,100-logspace(1,-4,6)]));   % available percentage of reliability
    R_num=menu(['Please choose the reliability of "' SHAFT(SNO) '" shaft from the MENU'],OPT);   
    ke=TABLE1(R_num,2);	% reliability factor 
    
   
    kd=1;
    kf=1;	% not knowing any available information, We assume miscellaneous-effects factor is 1
    
    % stress concentration factor for "Steps" 
    kf_steps=2.05;
    kfs_steps=1.85;
    d_guess=input(['"Shaft", Estimate the initial diameter for "' SHAFT{SNO} '" shaft in mm: \n']);     %initial estimation for Shaft Diameter
    
    while true  % iteration for endurance limit estimation
        if d_guess>=2.79 && d_guess<=51     
            kb=(d_guess/7.62)^(-0.107);
        elseif d_guess>51 && d_guess<=254    
            kb=1.51*d_guess^(-0.157);
        end
        i=5;    % fifth point on the shaft
        ki=[ka kc kd ke kf];	% concatenating of Marin Coefficients except size factor
        Se_it=Se_prime*(prod([ki kb]));
        D=((16*n_d/pi*((1/Se_it)*(4*(kf_steps*M_a(round(shaft_disp(i,3)*100)))^2 + ...
          3*(kfs_steps*T_a(i))^2)^0.5+(1/Sy)*(4*(kf_steps*Mm)^2+3*(kfs_steps*T_m(i))^2)^0.5))^(1/3))*1000;
        if (abs(d_guess-D)/D)<=0.01
            Se=Se_it; % final endurance limit
            break
        else
            d_guess=D;
        end
    end
    Se=Se*ones(1,8);
    
    % ----- Stress Concentration Analysis
    % initial stress concentration factor for "Steps" 
    kf_steps=2.05*ones(1,4);
    kfs_steps=1.85*ones(1,4);
    % initial stress concentration factor for "Retaining Rings" 
    kf_rings=2*ones(1,4);
    kfs_rings=1.5*ones(1,4);
    % initial stress concentration factor for "Keyseats" 
    kf_key=1.8*ones(1,2);
    kfs_key=1.7*ones(1,2);
    
    c=0;	% counter for steps
    p=0;	% counter for keyseats
    k=0;	% counter for retaining rings
    d=zeros(1,n);   % a container for shaft section diameters
    
    flag=0;
    NS=load('NotchSensitivity.txt');            % a table containing data from figs. 6-21,22 [1]
    SSCS=load('StaticStressConStep.txt');       % a table containing data from figs. A-15-8,9 [1]
    SSCG=load('StaticStressConGroove.txt');     % a table containing data from figs. A-15-16,17 [1]
    d_shaft_holder=zeros(n,1);
    
    while true
        I=[];   % container for keyseats
        J=[];   % container for retaining rings
        % ***** 1st set of critical points: stress concentration points
        for i=1:n   % iteration for "The Ten" stress concentration points    
            if flag~=0  % final endurance limit calculation
                if d_shaft(i)>=2.79 && d_shaft(i)<=51     
                    kb(i)=(d_shaft(i)/7.62)^(-0.107);
                elseif d_shaft(i)>51 && d_shaft(i)<=254    
                    kb(i)=1.51*(d_shaft(i))^(-0.157);
                end
                Se(i)=Se_prime*(prod([ki kb(i)]));  % final endurance limit
            end
            % diameter assignment + update of stress concentration factors
             switch shaft_disp(i,2) 
                 % *** STEPS
                 case 1                 % Steps on shaft 
                     c=c+1;
                     if flag~=0
                         r1(c)=0.02*d_shaft(i);      % racord fillet of steps
                         for j=1:11
                             if (NS(j,1)<=r1(c) && NS(j+1,1)>r1(c))
                                 qt1(c)=NS(j,2)+(r1(c)-NS(j,1))*(NS(j+1,1)-NS(j,1))/(NS(j+1,2)-NS(j,2));      % notch sensitivity q
                                 qts1(c)=NS(j,3)+(r1(c)-NS(j,1))*(NS(j+1,1)-NS(j,1))/(NS(j+1,3)-NS(j,3));     % notch sensitivity qs
                             end
                         end
                         switch c
                             case 1         % first step on shaft in Left-to-Right direction
                                 Diaratio(c)=d_shaft(i+1)/d_shaft(i);
                             case {2,3}    % second and third step on shaft in Left-to-Right direction
                                 Diaratio(c)=d_mean/d_shaft(i);
                             case 4         % fourth step on shaft in Left-to-Right direction
                                 Diaratio(c)=d_shaft(i-1)/d_shaft(i);
                         end
                         for j=1:4  % static stress concentration factors (Kt & Kts)
%                              if (SSCS(j,1)<=Diaratio(c) && SSCS(j+1,1)>Diaratio(c))
                                 Kt1(c)=SSCS(j,2)+(Diaratio(c)-SSCS(j,1))*(SSCS(j+1,1)-SSCS(j,1))/(SSCS(j+1,2)-SSCS(j,2));
                                 Kts1(c)=SSCS(j,3)+(Diaratio(c)-SSCS(j,1))*(SSCS(j+1,1)-SSCS(j,1))/(SSCS(j+1,3)-SSCS(j,3));
%                              end
                         end
                         kf_steps(c)=1+qt1(c)*(Kt1(c)-1);     % dynamic stress concentration factor Kf for steps
                         kfs_steps(c)=1+qts1(c)*(Kts1(c)-1);  % dynamic stress concentration factor Kfs for steps
                     end

                      d(i)=(16*n_d/pi*((1/Se(i))*(4*(kf_steps(c)*M_a(round(shaft_disp(i,3)*100)))^2 +3*(kfs_steps(c)*T_a(i))^2)^0.5+(1/Sy)*(4*(kf_steps(c)*Mm)^2+3*(kfs_steps(c)*T_m(i))^2)^0.5))^(1/3);

                 % *** RETAINING RINGS     
%                  case 2         % Retaining Rings on shaft 
%                      p=p+1; 
%                      J=[J i];   % saving No. retaining ring on each shaft
%                      if flag~=0
%                         r2(p)=0.1*m(p);     % internal racord fillet of retaining ring seat
%                         for j=1:10
%                              if (NS(j,1)<=r2(p) && NS(j+1,1)>r2(p))
%                                  qt2(p)=NS(j,2)+(r2(p)-NS(j,1))*(NS(j+1,1)-NS(j,1))/(NS(j+1,2)-NS(j,2));
%                                  qts2(p)=NS(j,3)+(r2(p)-NS(j,1))*(NS(j+1,1)-NS(j,1))/(NS(j+1,3)-NS(j,3));
%                              end
%                         end
%                         r_to_t(p)=r2(p)/(0.5*(d_shaft(i)-d_ret(p)));    % r/t ratio for retaining rings
%                         for j=1:5   % static stress concentration factors (Kt & Kts)
% %                             if (SSCG(j,1)<=r_to_t(p) && SSCG(j+1,1)>r_to_t(p))
%                                  Kt2(p)=SSCG(j,2)+(r_to_t(p)-SSCG(j,1))*(SSCG(j+1,1)-SSCG(j,1))/(SSCG(j+1,2)-SSCG(j,2));
%                                  Kts2(p)=SSCG(j,3)+(r_to_t(p)-SSCG(j,1))*(SSCG(j+1,1)-SSCG(j,1))/(SSCG(j+1,3)-SSCG(j,3));
% %                             end
%                         end
%                         kf_rings(p)=1+qt2(p)*(Kt2(p)-1);
%                         kfs_rings(p)=1+qts2(p)*(Kts2(p)-1);
%                      end
%                      d(i)=(16*n_d/pi*((1/Se(i))*(4*(kf_rings(p)*M_a(round(shaft_disp(i,3)*100)))^2+ ...
%                          3*(kfs_rings(p)*T_a(i))^2)^0.5+(1/Sy)*(4*(kf_rings(p)*Mm)^2+3*(kfs_rings(p)*T_m(i))^2)^0.5))^(1/3);    

                 % *** KEYSEATS
                 case 3          % Keyseats on shaft (in center of keyseats, torque is half!)
                     I=[I i];    % saving No. keyseats on each shaft
                     k=k+1;
                     if flag~=0
                         r3(k)=0.02*d_shaft(i);
                         for j=1:11
                             if (NS(j,1)<=r3(k) && NS(j+1,1)>r3(k))
                                 qt3(k)=NS(j,2)+(r3(k)-NS(j,1))*(NS(j+1,1)-NS(j,1))/(NS(j+1,2)-NS(j,2));      % notch sensitivity q
                                 qts3(k)=NS(j,3)+(r3(k)-NS(j,1))*(NS(j+1,1)-NS(j,1))/(NS(j+1,3)-NS(j,3));     % notch sensitivity qs
                             end
                         end
                         Kt3=2.2;
                         Kts3=2;
                         kf_key(k)=1+qt3(k)*(Kt3-1);
                         kfs_key(k)=1+qts3(k)*(Kts3-1);
                     end
                     
                     d(i)=(16*n_d/pi*((1/Se(i))*(4*(kf_key(k)*max(M_a(shaft_disp(i,3)*100+(-10:10))))^2+ ...
                         3*(kfs_key(k)*T_a(i)/2)^2)^0.5+(1/Sy)*(4*(kf_key(k)*Mm)^2+3*(kfs_key(k)*T_m(i)/2)^2)^0.5))^(1/3); 
             end  
        end

        % ***** 2nd critical point: maximum bending moment point along the shaft
        d_int_bearing=0;            % diameter of internal bearing seat (only 4 shafts number 1 & 3)
        Ind=find(M_a==max(M_a));
          if any(shaft_disp(:,3)==round(Ind*0.01))
              fprintf('\nlocation of max. bending moment in the "%s" is a stress concentration point and it will check in above! \n\n',SHAFT{SNO})
          else 
              fprintf('\nlocation of max. bending moment in the "%s" is the internal bearing and free from stress concentration! \n\n',SHAFT{SNO})
              d_int_bearing=((16*n_d/pi*((1/Se(i))*(4*(max(M_a))^2+ ...
                  3*(T_maximum/2)^2)^0.5+(1/Sy)*(4*(Mm)^2+3*(T_maximum/2)^2)^0.5))^(1/3))*10^3;
              if flag==1
                  n_D0=((((16/(pi*(d_shaft(I(1)+2))^3))*((1/Se(i))*(4*(max(M_a))^2 + ...
                              3*(T_maximum/2)^2)^0.5+(1/Sy)*(4*(Mm)^2+3*(T_maximum/2)^2)^0.5)))^(-1))/10^9;
              end
          end
          d=d*10^3;   % (mm)

       %  Show the diameter of Stress Concentration Points to the User and correction of them 
        if d_int_bearing~=0
            S1=[(shaft_disp(1:I(1)+2,2))' "none" (shaft_disp(I(1)+3:end,2))'; ...
                round([d(1:I(1)+2),d_int_bearing,d(I(1)+3:end)],1);string([COM(1:I(1)+2),COM(I(1)+2),COM(I(1)+3:end)])];
            Table_dia=array2table(S1,'VariableNames',[NUM(1:I(1)+2)','max_bending',NUM(I(1)+3:end)'], ...
                'RowNames',{'StrCon Guide','Diameter (mm)','Component'});
            disp(Table_dia)
        else
            S2=[(shaft_disp(:,2))';round(d(:)',2);COM];
            Table_dia=array2table(S2,'VariableNames',NUM,'RowNames',{'StrCon Guide','Diameter (mm)','Component'});
            disp(Table_dia)
        end
        disp('"Try to select same diameter for numbers that are locating in a same section, e.g. numbers on gears seat"')
        d_shaft=input(['Correct "' SHAFT{SNO} '" shaft diameter that are showed above in format of a "1-by-10" horizontal vector: \n']);

        d_mean=1.2*max(d_shaft([n/2,(n/2)+1]));     % diameter of the biggest section on shaft
        % final outputs
        if (flag==1 && all(d_shaft_holder==d_shaft))   % "if" condition for check the convergance of shaft-diameter-finding process
           D_raw_shaft=[unique(d_shaft(1:5)) unique(d_shaft(6:8),'stable')];
           D_shaft=[D_raw_shaft(1:2) ,d_mean ,D_raw_shaft(3:4)];  
           break
        end
        d_shaft_holder=d_shaft;      % a container for saving diameters of shafts
        flag=1;
        
        % ----- Key and Keyseat calculation:
        % 1.rectangular key
        Flag1=0;
        Flag2=0;
        k=0;    % reset counter
        RecKey=load('RectangularKey.txt');  % data of required retaining rings from std. tables
        for i=I     % iteration for rectangular keys selection
            k=k+1;
            w=RecKey(RecKey(:,1)==d_shaft(i),2);
            h=RecKey(RecKey(:,1)==d_shaft(i),3);
            L_min=RecKey(RecKey(:,1)==d_shaft(i),4);

            S_y_key=Sy;
            S_y_shaft=Sy;
            S_y_com=shaft_disp(i,6);

            S_y_min=min([S_y_key,S_y_com,S_y_shaft]);
            L_s=(4*T_maximum*n_d*1000)/(w*d_shaft(i)*S_y_key);
            L_p=(4*T_maximum*n_d*1000)/(h*d_shaft(i)*S_y_min);

            L_key=max([L_s,L_p,L_min]);
            L_keyseat=L_key+w;

            if L_keyseat+15>=shaft_disp(i,5)    % check for hub-needing (at least be 15 mm difference)
                warning('\n\nOops! You need hub in the component %s of %s shaft! ',COM(i),SHAFT{SNO})
                if L_keyseat+15<(shaft_disp(i+1,3)-shaft_disp(i-1,3))
                    fprintf('but width that considered for this component seat is sufficient\n')
                else
                   fprintf('and width that considered for this component seat is "NOT" sufficient\n') 
                end
            end

            shaft_disp(i,7)=w;
            shaft_disp(i,8)=h;
            shaft_disp(i,9)=L_keyseat;

            switch i    % check for "The Two" critical point that appear after rec. key calculation  
                case I(1)
                    d1=(16*n_d/pi*((1/Se(i))*(4*(kf_key(k)*M_a(round((shaft_disp(i,3)+L_keyseat/2)*100)))^2+ ...
                      3*(kfs_key(k)*T_a(i))^2)^0.5+(1/Sy)*(4*(kf_key(k)*Mm)^2+3*(kfs_key(k)*T_m(i))^2)^0.5))^(1/3);   
                    if d1>d_shaft(i)
                        warning('diameter of shaft at the end of first keyseat is dangerous')
                        disp('You should repeat the shaft-diameter-finding process!')
                    else
                        n1_D=((((16/(pi*(d_shaft(i))^3))*((1/Se(i))*(4*(kf_key(k)*max(M_a(round((shaft_disp(i,3)+L_keyseat/2)*100))))^2+ ...
                         3*(kfs_key(k)*T_a(i)/2)^2)^0.5+(1/Sy)*(4*(kf_key(k)*Mm)^2+3*(kfs_key(k)*T_m(i)/2)^2)^0.5)))^(-1))/10^9; 
                        Flag1=1;
                    end
                case I(2)
                    d2=(16*n_d/pi*(1/Se(i)*(4*(kf_key(k)*M_a(round((shaft_disp(i,3)-L_keyseat/2)*100)))^2+ ...
                     3*(kfs_key(k)*T_a(i))^2)^0.5+1/Sy*(4*(kf_key(k)*Mm)^2+3*(kfs_key(k)*T_m(i))^2)^0.5))^(1/3);   
                    if d2>d_shaft(i)
                        warning('diameter of shaft at the beginning of second keyseat is dangerous') 
                        disp('You should repeat the shaft-diameter-finding process!')
                    else
                        n2_D=((((16/(pi*(d_shaft(i))^3))*((1/Se(i))*(4*(kf_key(k)*max(M_a(round((shaft_disp(i,3)-L_keyseat/2)*100))))^2+ ...
                         3*(kfs_key(k)*T_a(i)/2)^2)^0.5+(1/Sy)*(4*(kf_key(k)*Mm)^2+3*(kfs_key(k)*T_m(i)/2)^2)^0.5)))^(-1))/10^9; 
                        Flag2=1;
                    end
            end       
        end

        % 2. Retaining ring
        RetRingKey=load('RetainingRingKey.txt');
        c=0;
        for i=J
            c=c+1;
            m(c)=RetRingKey(RetRingKey(:,1)==d_shaft(i),2);
            n_min(c)=RetRingKey(RetRingKey(:,1)==d_shaft(i),3);
            d_ret(c)=RetRingKey(RetRingKey(:,1)==d_shaft(i),4);
            chamfer=2;
            switch i
                case J(1)
                    if (m(c)+n_min(c)+chamfer)>(shaft_disp(i,3))
                       fprintf('distance from "%s" retaining ring and edge is not sufficient!\n','first') 
                    end
                case J(2)
                    if (m(c)+n_min(c)+chamfer)>(shaft_disp(i,3)-shaft_disp(i-1,3))
                       fprintf('distance from "%s" retaining ring and edge is not sufficient!\n','second') 
                    end
                case J(3)
                    if (m(c)+n_min(c)+chamfer)>(shaft_disp(i,3)-shaft_disp(i-1,3))
                       fprintf('distance from "%s" retaining ring and edge is not sufficient!\n','third') 
                    end
                case J(4)
                    if (m(c)+n_min(c)+chamfer)>(EL(end)-shaft_disp(i,3))
                       fprintf('distance from "%s" retaining ring and edge is not sufficient!','fourth') 
                    end
            end
            
            shaft_disp(i,10)=m(c);
            shaft_disp(i,11)=d_ret(c);
        end
           c=0;     % reset counter
           p=0;     % reset counter
           k=0;     % reset counter
    end
    
    % ----- safety factor of each section
    c=0;
    p=0;
    k=0;
    n_D=zeros(1,n);
    if (Flag1==1 && Flag2==1)
        for i=1:n      
             switch shaft_disp(i,2)              
                 case 1     % Steps 
                     c=c+1;
                     n_D(i)=((((16/(pi*(d_shaft(i))^3))*((1/Se(i))*(4*(kf_steps(c)*M_a(round(shaft_disp(i,3)*100)))^2 + ...
                          3*(kfs_steps(c)*T_a(i))^2)^0.5+(1/Sy)*(4*(kf_steps(c)*Mm)^2+3*(kfs_steps(c)*T_m(i))^2)^0.5)))^(-1))/10^9;
                 case 2     % Retaining Rings 
                     p=p+1;
                     n_D(i)=((((16/(pi*(d_shaft(i))^3))*((1/Se(i))*(4*(kf_rings(p)*M_a(round(shaft_disp(i,3)*100)))^2+ ...
                         3*(kfs_rings(p)*T_a(i))^2)^0.5+(1/Sy)*(4*(kf_rings(p)*Mm)^2+3*(kfs_rings(p)*T_m(i))^2)^0.5)))^(-1))/10^9;     
                 case 3     % Keyseats
                     k=k+1;
                     n_D(i)=((((16/(pi*(d_shaft(i))^3))*((1/Se(i))*(4*(kf_key(k)*max(M_a(shaft_disp(i,3)*100+(-10:10))))^2+ ...
                         3*(kfs_key(k)*T_a(i)/2)^2)^0.5+(1/Sy)*(4*(kf_key(k)*Mm)^2+3*(kfs_key(k)*T_m(i)/2)^2)^0.5)))^(-1))/10^9; 
             end
        end
        % Show the safety factors versus stress concentration point indicator on shaft
        if d_int_bearing~=0     % max. bending is on internal bearing
            N_D=[n_D(1:I(1)),n1_D,n_D(I(1)+1:I(1)+2),n_D0,n_D(I(1)+3:I(2)-1),n2_D,n_D(I(2):end)];
            S4=array2table(round(N_D,2),'VariableNames',[NUM(1:I(1))','rightkey1',NUM(I(1)+1:I(1)+2)',...
                'max_bending',NUM(I(1)+3:I(2)-1)','leftkey2',NUM(I(2):end)'],'RowNames',{'SafetyFactor'});
            disp(S4)
        else    % max. bending is a stress concentration point
            N_D=[n_D(1:I(1)) n1_D n_D(I(1)+1:I(2)-1) n2_D n_D(I(2):end)];
            S3=array2table(round(N_D,2),'VariableNames',[NUM(1:I(1))','rightkey1',NUM(I(1)+1:I(2)-1)',...
                'leftkey2',NUM(I(2):end)'],'RowNames',{'SafetyFactor'});
            disp(S3)
        end
        if all(N_D>=n_d)
            disp('Shaft-diameter-finding process successfully finished!')
        else
            Ind2=find(n_D<n_d);
            x=num2str(Ind2(2:end)');
            y=[repmat(' & ',[numel(x) 1]) x]';
            
%             warning(['The point(s) number ' num2str(Ind2(1)) y(:)' ...
%                 ' in shaft is(are) not safe because of low safety factor and you should repeat or refine shaft-diameter-finding process'])
        end
    end
end