function [BearingType,R]=BearingSelection(omega,reaction,bearing_seat_dia,Flag_axial,F_a)

    % SKF parameters;
    L10=10^6;       % rating life
    x0=0.02;
    theta=4.459;
    b=1.483;
    a=3;
    af=1;           % application factor
    BearingType=repmat("DGBB",[3 2]);   % default assumption
    
    R_set=input('Please input "percentage" of reliability of bearing set \n');
    R_set=R_set/100;
    R0=(R_set^(1/6))*ones(3,2);      % initial assumption
    
    DGBB_SKF_tab=load('DGBB_GuideTable.txt');           % guide table for deep groove ball bearings (table8-SKF) 
    DGBB_SKF=load('DeepGrooveBallBearing.txt');         % table encompassing some required DGBB data
    CRB_SKF=load('CylindricalRollerBearing.txt');       % table encompassing some required CRB data
    TAB_SKF=DGBB_SKF;                                   % working table in calculation flow
    FLAG=[];
    flag=0;                        % "flag" specifies type of bearing calculation: 1=CRB & 0=DGBB
    t=0;
    while true
        for i=1:3           % iteration on No. shafts
            LD=10*300*8*60*omega(i);   % absolute desired life
            xD=LD/L10;                  % rational desired life
            for j=[1 2]     % iteration on No. bearings
                if flag~=0 && FLAG(i,j)==0
                    continue
                end
                Fr=reaction(i,j);                       % radial force on bearing
                d_b=bearing_seat_dia(i,j);              % bearing seat diameter on shaft 
                C10=TAB_SKF(TAB_SKF(:,1)==d_b,2);       % rating load of bearing
                
                if Flag_axial(i,j)==1 && flag==0
                    Fa=F_a(i);                          % axial load of bearing that supports axial load                 
                    C0=TAB_SKF(TAB_SKF(:,1)==d_b,3);    % static rating load    
                    f0=TAB_SKF(TAB_SKF(:,1)==d_b,4);    % f0-factor from SKF catalogue
                    DGB_Ref=(f0*Fa)/(C0);               % reference for comparison with threshold "e"
                    
                    for k=1:1:8   % interpolation
                        if (DGB_Ref>DGBB_SKF_tab(k,1)) 
                            if(DGB_Ref<DGBB_SKF_tab(k+1,1))
                                e=DGBB_SKF_tab(k+1,2)-(DGBB_SKF_tab(k+1,1)-DGB_Ref)*(DGBB_SKF_tab(k+1,2)-DGBB_SKF_tab(k,2))/(DGBB_SKF_tab(k+1,1)-DGBB_SKF_tab(k,1));
                                X=DGBB_SKF_tab(k+1,3)-(DGBB_SKF_tab(k+1,1)-DGB_Ref)*(DGBB_SKF_tab(k+1,3)-DGBB_SKF_tab(k,3))/(DGBB_SKF_tab(k+1,1)-DGBB_SKF_tab(k,1));
                                Y=DGBB_SKF_tab(k+1,4)-(DGBB_SKF_tab(k+1,1)-DGB_Ref)*(DGBB_SKF_tab(k+1,4)-DGBB_SKF_tab(k,4))/(DGBB_SKF_tab(k+1,1)-DGBB_SKF_tab(k,1));
                            end
                        elseif DGB_Ref<0.172
                            e=0.19;
                            X=0.56;
                            Y=2.3;
                        elseif DGB_Ref>6.89
                            e=0.44;
                            X=0.56;
                            Y=1;
                        end
                    end
                    P=Fr;                   
                else
                    P=Fr;   % equivalent force of non-axial-supporting bearings
                end                
                C10_min=af*P*(xD/(x0+(theta-x0)*(1-R0(i,j))^(1/b)))^(1/a);
                if C10_min<C10
                    fprintf('Ball Bearing "B%.0f%.0f" with respect to its own estimated reliability is Ok! \n',i,j)
                end
                R(i,j)=0.99;   % collection of whole reliability of bearing set
            end
        end
        if (prod(R(:))<R_set && flag~=0)
            t=t+1;
            fprintf('In the try "%.0f" you can not reach to set reliability of bearings!! R_SET = %.3f\n',t,prod(R(:)))
            disp(R)
            ANS1=input('Do you want to continue? 1.Yes 2.No \n');
            if ANS1==2
                warning('You should refine shaft diameters and repeat all stages!')
                break
            end 
        end        
        if (prod(R(:))>=R_set && flag==0)   % check the total reliability required for set
            fprintf('\nMinimum set reliability satisfied! R_SET = %.3f\nAll bearings are deep groove ball bearings\n',prod(R(:)))
            break
        elseif (prod(R(:))>=R_set && flag~=0)
            fprintf('\nMinimum set reliability satisfied! R_SET = %.3f \n',prod(R(:)))
            [X, Y]=ind2sub([3 2],find(FLAG==1));
            BearingType(X,Y)="CRB";
            break   
        end
            if t==0
                ANS2=input('We offer you to select roller bearing instead of ball bearing, Do you accept? 1.Yes 2.No \n');
            end
            if ANS2==1
                a=10/3;
                TAB_SKF=CRB_SKF;    % change the working table for CR bearing calculation
                flag=1;
                disp('location guide map : [B11 B21; B12 B22; B13 B23]')
                disp(R)
                FLAG=input('Bearings that you want to substitute by roller bearing, please indicate them by "1" and others by "0"\n');
            else
                warning('You should refine shaft diameters and repeat all stages!')
                break
            end
    end
end