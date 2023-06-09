classdef corrugated_composite
    %"corrugated_compoiste" is a object a class used to analyze bending and
    %free vibration response of corrugated core sandwitch panel based on
    %model of equivalent elastic constants and first order shear deformable
    %plate theory (FSDT). The model includes functionality to model cross
    %ply laminates
    
    properties
        geom    struct      % Corrugation/Plate Geometry 
        mat_props   struct  % Material Properties
        equ_props           % Equivalent laminate material properties
        D           % Plate bending stiffness matrix
        M = 1       % Mode shape factor in x-direction (default to 1)
        N = 1       % Mode shape factor in y-direction (default to 1)
        L_eff_x = 1;      % Effective Length in x-direction (default 1 for simply supported on both ends)
        L_eff_y = 1;      % Effecitve Length in y-direction (default 1 for simply supported on both ends)
        w_n                 % Plate global natural frequency
        bending_rslt        % Structure with results from bending analysis
        wn_strip            % Local behavior natural frequency
    end
    
    methods
        function obj = def_geom(obj,geom)
            % def geometry method used to define the corrugation geometry 
            % and stores them as property of object
            obj.geom = geom;
            obj.geom.f = (geom.p-geom.h_c/tand(geom.theta));
        end
        function obj = def_mat_props(obj,mat_props)
            % define material properties, caulcate equivalent material
            % properties and store as property of object
            obj.mat_props   = mat_props;
            obj.equ_props.E = mean([mat_props.E1,mat_props.E2]);
            obj.equ_props.v = sqrt(mat_props.v12*mat_props.v21);
            obj.equ_props.G = mat_props.G12;
            obj.equ_props.rho = mat_props.rho;
        end
        function obj = calc_stiffness(obj)
            % Call to comp_equ_stiff_fun to calculate the equivalent
            % stiffness of a composite panel with the equivalent material
            % properties. Assumed to be cross ply layup with equal number
            % of 0 and 90 plys
            obj.D = Comp_Equ_Stiff_Fun(obj.geom,obj.equ_props);
        end
        function obj = Modal_Analysis(obj)
            % the modal analysis method calculates the global and local 
            % natural frequincies based on FSDT
            
            % Calculate equivalent panel properties
            obj = calc_stiffness(obj);  
            A_c = (obj.geom.f + obj.geom.h_c/sind(obj.geom.theta))*obj.geom.t_c/obj.geom.p; % Area of corrugation per pitch
            A_f = 2*obj.geom.t_f;  % Area of face sheets per pitch
            % Equivalent density assuming the plate is the thickness of corrugation]
            rho_equ = obj.mat_props.rho*(A_c+A_f)/obj.geom.h_c;
            rho_equ = rho_equ/(obj.L_eff_x*obj.L_eff_y);    % Scale by effective length
            I_0 = rho_equ*obj.geom.h_c;
            I_2 = rho_equ*obj.geom.h_c^3/12;
            
            % Scale by effective length for boundary condition
            a_eff = obj.geom.a*obj.L_eff_x;
            b_eff = obj.geom.b*obj.L_eff_y;

            alpha = obj.M*pi/a_eff;
            beta = obj.N*pi/b_eff;

            % calculate the S matrix 
            s = zeros(3,3);
            s(1,1) = obj.D.D_qx*alpha^2 + obj.D.D_qy*beta^2;
            s(1,2) = obj.D.D_qx*alpha;
            s(1,3) = obj.D.D_qy*beta;
            s(2,2) = (obj.D.D_11*alpha^2 + obj.D.D_66*beta^2 + obj.D.D_qx);
            s(2,3) = (obj.D.D_12 + obj.D.D_66)*alpha*beta;
            s(3,3) = obj.D.D_66*alpha^2 + obj.D.D_22*beta^2 + obj.D.D_qy;
            s(2,1) = s(1,2);
            s(3,1) = s(1,3);
            s(3,2) = s(2,3);

            mass = [I_0 0 0; 0 I_2 0; 0 0 I_2]; % define the mass matrix
            [A,B] = eig(s,mass);  % cauculate eigen values of mass and stiffness
            FREQ = sqrt(min(B(B>0)))/(2*pi); % post process to find natural frequency
            obj.w_n = FREQ(1);      % store natural freqiency in object properties

            % locic to eliminate results from physcally unobtainable
            % geometry
            if obj.geom.f < 0       
                obj.w_n = NaN;
                obj.wn_strip = NaN;
            else
                % for physically obtainable goemetry call to strio_vibe
                % subroutine to find local natural freqiency
                [obj.wn_strip] = Strip_Vibe(obj.geom,obj.equ_props);
            end
        end
        
        function obj = Bending_Analysis(obj,q_0)
            % sub-routine to calculate the bending response to a uniformly
            % distributed load 'q_0'
                obj.bending_rslt.q_0 = 0;
                a = obj.geom.a;
                b = obj.geom.b;
                q_mn = zeros(15,15);
                W_mn = zeros(15,15);
                x = linspace(0,a,50);
                y = linspace(0,b,50);
                w = zeros(50,50);
            for m = 1:30
                for n = 1:30       
                    q_mn(m,n) = 16*q_0/(pi^2*m*n);
                    alpha = m*pi/a;
                    beta = n*pi/b;
                    s = zeros(3,3);
                    s(1,1) = obj.D.D_qx*alpha^2 + obj.D.D_qy*beta^2;
                    s(1,2) = obj.D.D_qx*alpha;
                    s(1,3) = obj.D.D_qy*beta;
                    s(2,2) = (obj.D.D_11*alpha^2 + obj.D.D_66*beta^2 + obj.D.D_qx);
                    s(2,3) = (obj.D.D_12 + obj.D.D_66)*alpha*beta;
                    s(3,3) = obj.D.D_66*alpha^2 + obj.D.D_22*beta^2 + obj.D.D_qy;
                    s(2,1) = s(1,2);
                    s(3,1) = s(1,3);
                    s(3,2) = s(2,3);       
                    Unk = (s)\[q_mn(m,n);0;0];
                    W_mn(m,n) = Unk(1);
                end
            end
            for i = 1:50
                for j = 1:50
                    for m = 1:15
                        for n = 1:15
                            temp = W_mn(m,n)*sin(m*pi*x(i)/a)*sin(n*pi*y(j)/b);
                            w(i,j) = w(i,j) + temp;
                        end
                    end
                end
                
            end
            obj.bending_rslt.w = w;
            obj.bending_rslt.x = x;
            obj.bending_rslt.y = y;
        end  
    end
        
    
end
% Functions Referenced Within The Object
function[w_n] = Strip_Vibe(geom,equ_props)

    % sub-routine to calcualte the local natural freqiencies based on
    % assumed SSC plate between corrugation cells
    
    % Define corrugation geometry variables
    
    
    p     = geom.p;
    t_f   = geom.t_f;
    t_c   = geom.t_c;
    E   = equ_props.E;
    v   = equ_props.v;
    rho   = equ_props.rho/0.75;
    G = equ_props.G;
    f = geom.f;
    t_eff = t_f;

    D_11 = E*t_eff^3/(12*(1-v^2));
    D_12 = v*D_11;
    D_22 = D_11;
    D_66 = G*t_eff^3/12;
    A_44 = 5/6*G*t_eff;
    A_55 = A_44;
    
    a = geom.a;
    b = (2*p-f)*0.75;
        
    I_0 = rho*t_eff;
    I_2 = rho*t_eff^3/12;

    alpha = 1*pi/(a);
    beta = 1*pi/(b);

    s = zeros(3,3);
    s(1,1) = A_44*alpha^2 + A_55*beta^2;
    s(1,2) = A_44*alpha;
    s(1,3) = A_55*beta;
    s(2,2) = (D_11*alpha^2 + D_66*beta^2 + A_44);
    s(2,3) = (D_12 + D_66)*alpha*beta;
    s(3,3) = D_66*alpha^2 + D_22*beta^2 + A_44;
    s(2,1) = s(1,2);
    s(3,1) = s(1,3);
    s(3,2) = s(2,3);
    if b > 0
        mass = [I_0 0 0; 0 I_2 0; 0 0 I_2];
        [A,B] = eig(s,mass);
        FREQ = sqrt(min(B(B>0)))/(2*pi);
        w_n = FREQ(1); 
    else
        w_n = NaN;
    end
    
end



function S = S_Fun(h_c,h,t_c,t_f,theta,p,v_c,v_f,E_c,E_f)   
% S_fun subroutine used to calculate the shear correction factor for a
% symmetrical sandwitch panel as defined in appendix D of Libove and Hubka
R_c1 = 0;
k_y = 1;
k_z = 1;
a_1 = (1-k_z/2)*h_c - R_c1;
e_1 = R_c1*cosd(theta);
g_1 = R_c1*sind(theta);
j_1 = a_1 + e_1;
k_1 = j_1*cotd(theta);
d_1 = j_1*cscd(theta);
b_1 = k_1 + g_1;
f_1 = 2*((1-k_y/2)*p-b_1);
   
KI_z = 2/3*(k_1/h_c)^2*(d_1/h_c) + 2/3*(1/8*(p/h_c)^3-(b_1/h_c)^3); 
KI_yz = 2/3*(j_1/h_c)*b_1/h_c*d_1/h_c + 1/2*(1/4*(p/h_c)^2 - (b_1/h_c)^2);
KI_y = 2/3*(j_1/h_c)^2*d_1/h_c + f_1/h_c/4;
KA_z = 0;
KA_y = 0;
K_L = 2*d_1/h_c+f_1/h_c;
KL_y = f_1/h_c + 2*d_1/h_c*cosd(theta)^2;
KL_yz = 2*d_1/h_c*sind(theta)*cosd(theta);
KL_z = 2*d_1/h_c*sind(theta)^2;

B_3 = KI_z + 1/12*(t_c/h_c)^2*KL_z;
B_4 = KI_yz - 1/12*(t_c/h_c)^2*KL_yz;
B_6 = KI_y + 1/12*(t_c/h_c)^2*KL_y;
B_7 = E_f/E_c*(1-v_c^2)/(1-v_f^2)*(t_f/t_c)^3;

NUM = 6*h_c/p*B_3*B_7+(p/h_c)^2;
T1 = -2*(p/h_c)^2*B_4;
T2 = h_c/h*((6*B_7*(B_3*B_6-B_4^2)+(p/h_c)^3*B_6));
T3 = h/h_c*p/h_c*B_3;

S = NUM/(12*(T1+T2+T3));
end
% Functions Referenced Within The Object
function [D] = Comp_Equ_Stiff_Fun(geom,equ_props)

    theta = geom.theta;
    h_c   = geom.h_c;
    p     = geom.p;
    t_c   = geom.t_c;
    t_f   = geom.t_f;
    
    E     = equ_props.E;
    v     = equ_props.v;  
    G     = equ_props.G;
    
    h   = h_c + (t_c+t_f);
    f = (p-h_c/tand(theta));
%Find the moment of inertia per unit width of corrugation
    I_c = ((f/2*t_c*h_c^2/2)+2*t_c*h_c^3/(12*sind(theta)))/(2*p);
    EI_y = E/2*t_f*h^2;
%Find the flexural rigidity in the x direction
    D_x = E*I_c + E/2*t_f*h^2;
    D_y = EI_y/(1-v^2*(1-EI_y/D_x));
    D_xy = G*t_f*h^2;
    v_x = v;
    v_y = v*D_y/D_x;
    A_c = (f + h_c/sind(theta))*t_c/p;

    S = S_Fun(h_c,h,t_c,t_f,theta,p,v,v,E,E);
%Find transverse shear stiffness
    D_qy = S*h*(E/(1-v^2))*(t_c/h_c)^3;
    D_qx = G*t_c^2/A_c*(h/p)^2;

    D_11 = D_x/(1-v_x*v_y);
    D_22 = D_y/(1-v_x*v_y);
    D_66 = D_xy/2;
    D_12 = D_11*v_y;

    D = struct;
    D.D_11 = D_11;
    D.D_22 = D_22;
    D.D_66 = D_66;
    D.D_12 = D_12;
    D.D_qy = D_qy;
    D.D_qx = D_qx;
end


