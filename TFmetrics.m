%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computes various metrics of a transfer function
%
% SYNTAX: M = TFmetrics(sys,"type")
% INPUTS:
% G - A transfer function of Continous LTI type.
% type - A string indicating the type of system. "UD" for underdamped, "CD"
%          for critically damped, "OD" for overdamped.
% 
%
% OUTPUTS: (Exhaustive list of outputs, some of these metrics are neglected
%           if type = "OD" or "CD")
% "M", which is a struct containing the metrics below:

% UNDERDAMPED CASE:
% M.Peak - Peak value
% M.Fv - Final value 
% M.Tau - Time constant (s)
% M.Pos - Percentage overshoot (%)
% M.Gos - Goal overshoot (%)
% M.Zeta - Damping coefficient (scalar between [0,1])
% M.Beta - Beta
% M.Tr - Rise time (s)
% M.Wnr - Rise time natural frequency (rad/s)
% M.Ts - Settle time (s)
% M.Wns - Settle time natural frequency (rad/s)
% M.Tp - Peak time (s)
% M.Wnp - Peak time natural frequency (rad/s)
%----------------------------------------------
% OVERDAMPED/CRITICALLY DAMPED CASE:
% M.Taue - Effective time constant (s)
% M.Tau - Time constant (s)
% M.Tr1 - Modified rise time (s)
% M.Zeta - Damping coefficient (scalar >= 1)
% M.Wn - Natural frequency (rad/s)
%----------------------------------------------
% Author: John Song
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function M = TFmetrics(G, type)
    t = 0:0.001:10; % long enough time vector for step to stabilize
    stepG = step(G,t);
    FV = stepG(end);
    peak = max(stepG);

    % Tau
    [~,tau_idx]=min(abs(stepG-FV*0.63)); 
    Tau = t(tau_idx);
    
    if (type == "UD")
        % Pos
        Pos = (peak - FV)*100/FV;
    
        % Ess
        Ess = (1-FV)*100;
    
        % Gos
        %Gos = (peak - 1)*100;
        Gos = (Pos*(FV/100) - Ess/100)*100;
    
        % Zeta & Beta
        Zeta = sqrt((log(Pos/100)^2)/(pi^2+(log(Pos/100))^2));
        Beta = sqrt(1-Zeta^2);

            
        % Check if Zeta is between 0 and 1
        if (Zeta < 0 || Zeta > 1) 
            error('Zeta out of bounds, double check your system!');
        end

        % Peak time
        Tp_idx = find(stepG==peak);
        Tp = t(Tp_idx);
        Wnp = pi/(Tp*Beta);
    
        % Rise time
        [~,Tr_idx] = min(abs(stepG(1:Tp_idx) - FV));
        Tr = t(Tr_idx);
        Wnr = (pi-atan(Beta/Zeta))/(Tr*Beta);
    
        % Settle Time
        
        min_idx = find(stepG(Tp_idx:end)==min(stepG(Tp_idx:end)));
        %max_val = arr(max_index);
        min_idx = Tp_idx + min_idx;
        %{
        [~,Ts_idx] = min(abs(stepG(Tp_idx:min_idx) - FV));
        Ts = t(Tp_idx + Ts_idx);
        %}
        [~,Ts_idx] = min(abs(stepG(min_idx:end) - FV*1.02));
        Ts = t(min_idx + Ts_idx);
        Wns = log(50/Beta)/(Zeta*Ts);

        % Assign outputs to struct
        M.Peak = peak;
        M.Fv = FV;
        M.Tau = Tau;
        M.Pos = Pos;
        M.Gos = Gos;
        M.Ess = Ess;
        M.Zeta = Zeta;
        M.Beta = Beta;
        M.Tr = Tr;
        M.Wnr = Wnr;
        M.Ts = Ts;
        M.Wns = Wns;
        M.Tp = Tp;
        M.Wnp = Wnp;
        return
    end
    
    if (type == "CD" || type == "OD") 
        % Taue
        P = (pole(G)).'; % two poles
        tau1 = -1/P(1);
        tau2 = -1/P(2);
        Taue = 1.05*(tau1 + tau2);

        % Modified rise time
        [~, t1_idx] = min(abs(stepG - FV * 0.1));
        t1 = t(t1_idx);
        [~, t2_idx] = min(abs(stepG - FV * 0.9));
        t2 = t(t2_idx);
        Tr1 = t2 - t1;

        % Zeta
        if (type == "CD") 
            Zeta = 1;
        end
        if (type == "OD")
            Zeta = Taue/(3.86*Taue - 1.83*Tr1);
        end

        % Wn
        Wn = (2.1*Zeta)/Taue;

        % Assign outputs to struct
        M.Tau = Tau;
        M.Taue = Taue;
        M.Tr1 = Tr1;
        M.Zeta = Zeta;
        M.Wn = Wn;
        return
    end
end
