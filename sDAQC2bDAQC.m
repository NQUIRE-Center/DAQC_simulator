function [gate_type,time,angles] = sDAQC2bDAQC(gate_type,time,angles,rxgatetime,rzgatetime)
%SDAQC2BDAQC converts a stepwise schedule to a banged schedule
% We ask for the rxgatetime and rzgatetime for simplicity of the code

% TODO: SIMPLIFY THE SQG BLOCKS, BUT THIS IS NOT 100% NECCESARY

% We count the number of analog blocks
K = sum(gate_type=="zz");

% Counter of the time at which the "current analog blocks stars"
time_current = 0;

% First digital block
% We implement this full block in the begining of the first analog block (for this we actually dont have to change anything)
deltaT = 0;
i = 1;
while gate_type(i) ~= "zz"
    deltaT = time(i);
    i = i+1;
end
% Obtain the time for the analog block
T = time(i+1)-time(i);
% Initilize the output variable for the time
timeOUT = time(1:i);

i = i+1;

% Run over all analog blocks
for k = 1:K-1
    % Obtain the total time for the digital block
    deltaT = 0;
    auxTimes = 0;
    while gate_type ~= "zz"
        deltaT = time(i);
        if gate_type(i) == "rx"
            auxTimes = [auxTimes auxTimes(end)+rxgatetime];
        elseif gate_type(i) == "rz"
            auxTimes = [auxTimes auxTimes(end)+rzgatetime];
        else
            error("Unknown SQG: "+gate_type(i))
        end
        i = i+1;
    end
    
    % Advance the time
    time_current = time_current+T;

    % Write the SQGs in the new position
    timeOUT = [timeOUT time_current-deltaT/2+auxTimes(1:end-1)];

    % Obtain the time for the next analog block
    T = time(i+1)-time(i);

    i = i+1;
end


% Last analog block
% In this case we have to displace the SQGs deltaT instead of deltaT/2
% Obtain the total time for the digital block
deltaT = 0;
auxTimes = 0;
while gate_type ~= "zz"
    deltaT = time(i);
    if gate_type(i) == "rx"
        auxTimes = [auxTimes auxTimes(end)+rxgatetime];
    elseif gate_type(i) == "rz"
        auxTimes = [auxTimes auxTimes(end)+rzgatetime];
    else
        error("Unknown SQG: "+gate_type(i))
    end
    i = i+1;
end
% Write the SQGs in the new position
timeOUT = [timeOUT time_current-deltaT+auxTimes(1:end-1)];

% Delete all zz blocks from gate_type and angles
% We are assuming that any "good" stepwise schedule cannot end in a zz block
i = 1;
while i < length(gate_type)
    % Removes the "zz" blocks
    if gate_type(i) == "zz"
        gate_type(i) = [];
        angles(:,i) = [];
    else
        i = i+1;
    end
end

end