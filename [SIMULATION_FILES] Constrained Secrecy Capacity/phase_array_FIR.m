function [FIR_rx1,FIR_rx2] = phase_array_FIR(angle_rx1,angle_rx2)

    N = 16;
    dist_el = 5.3e-3;
    
    TD_el_rx1 = (dist_el * sin(angle_rx1))/(3e8);
    TD_el_rx2 = (dist_el * sin(angle_rx2))/(3e8);
    
    time_rx1 = 0:TD_el_rx1:TD_el_rx1 * (N-1);
    time_rx2 = 0:TD_el_rx2:TD_el_rx2 * (N-1);
    
    signal_rx1 = exp(-2*1i*pi*28e9*time_rx1);
    signal_rx2 = exp(-2*1i*pi*28e9*time_rx2);
    
    symb_time = 1/10e9;


    FIR_rx1(1) = sum(signal_rx1(time_rx1 < 0.5 * symb_time));
    shftDex = 0;
    while sum(signal_rx1( ( (shftDex+0.5) * symb_time <= time_rx1) & (time_rx1 < (shftDex+1.5) * symb_time) )) ~= 0
        FIR_rx1(shftDex+2) = sum(signal_rx1( ( (shftDex+0.5) * symb_time <= time_rx1) & (time_rx1 < (shftDex+1.5) * symb_time) ));
        shftDex=shftDex+1;
    end

    FIR_rx2(1) = sum(signal_rx2(time_rx2 < 0.5 * symb_time));
    shftDex = 0;
    while sum(signal_rx2( ( (shftDex+0.5) * symb_time <= time_rx2) & (time_rx2 < (shftDex+1.5) * symb_time) )) ~= 0
        FIR_rx2(shftDex+2) = sum(signal_rx2( ( (shftDex+0.5) * symb_time <= time_rx2) & (time_rx2 < (shftDex+1.5) * symb_time) ));
        shftDex=shftDex+1;
    end
    
    % FIR_rx1(1) = sum(signal_rx1(time_rx1 < 0.5 * symb_time));
    % FIR_rx1(2) = sum(signal_rx1( (0.5 * symb_time <= time_rx1) & (time_rx1 < 1.5 * symb_time) ));
    % FIR_rx1(3) = sum(signal_rx1( (1.5 * symb_time <= time_rx1) & (time_rx1 < 2.5 * symb_time) ))
    % 
    % FIR_rx2(1) = sum(signal_rx2(time_rx2 < 0.5 * symb_time));
    % FIR_rx2(2) = sum(signal_rx2( (0.5 * symb_time <= time_rx2) & (time_rx2 < 1.5 * symb_time) ));
    % FIR_rx2(3) = sum(signal_rx2( (1.5 * symb_time <= time_rx2) & (time_rx2 < 2.5 * symb_time) ))

end