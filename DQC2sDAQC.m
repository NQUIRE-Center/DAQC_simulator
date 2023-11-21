%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% By: Mikel Garcia de Andoin, mikelgda@gmail.com
% Licensed under CC BY 4.0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [gate_type,time,angles] = DQC2sDAQC(gate_type,angles,CNOT_ctrl,CNOT_targ,Hs,minstep,rxgatetime,rzgatetime,timeCutThres)
%D2DA Function to conver a digital circuit to a stepwise digital-analog schedule
%(For the moment) It accepts just rx, rz and CNOT gates
%   gate_type            	Vector of strings that defines what type of block is the next. It can be CNOT, rz or rx
%   angles                  Matrix of the angles of rotations of the SQG. Each row corresponds to one qubit, and the columns match the ones in gate_type
%   CNOT_(ctrl/targ)        Matrix that labels the control/target qubits of the CNOT gates
%   Hs                      Coefficients of the source Hamiltonian. In units of 1/2pi*s.
%   minstep                 Discretization in time of the DA process. This can be the smallest time scale for which we can control the system. In units of s.
%   r(z/x)gatetime          Time for the application of the corresponding SQG. In units of s.
%   timeCutThres            Minimum time for considering the application of an analog block
%   
%   OUTPUT:
%    gate_type              Vector of strings that defines what type of block is the next. It can be zz, rz or rx
%    time                   Vector that defines the time stamps for changing each block type. It starts at 0, so it has one extra element compared to gate_type. In units of s/hbar
%    angles                 Matrix of the angles of rotations of the SQG. Each row corresponds to one qubit, and the columns match the ones in gate_type

% Size of the system
N = size(angles,1);


%---------------------------------------------------------------------------------------------------------------------------
% DECOMPOSE THE CNOT GATES INTO SQGs AND ZZ BLOCKS
%---------------------------------------------------------------------------------------------------------------------------
% We go through all the digital gates in order
NEWgate_type = strings(0);
NEWangles_SQG = [];
angles_2QG = [];
for i = 1:length(gate_type)
    switch gate_type(i)
        case {"rx","rz"}
            NEWgate_type(end+1) = gate_type(i);%#ok<*AGROW>
            NEWangles_SQG(:,end+1) = angles(:,i);
            if i == 1
                % Special initialization
                angles_2QG = zeros(N);
            else
                angles_2QG(:,:,end+1) = zeros(N);
            end
        case "CNOT"
            % We decompose the CNOT gate into a zz block and SQGs
            NEWgate_type(end+1:end+6) = ["rz" "rx" "rz" "zz" "rx" "rz"];
            NEWangles_SQG(:,end+1:end+6) = zeros(N,6);
            if i == 1
                angles_2QG = zeros(N,N,6);
            else
                angles_2QG(:,:,end+1:end+6) = zeros(N,N,6);
            end
            % We take into account that there might be CNOTs in parallel
            for j = nonzeros(CNOT_ctrl(:,i))'
                ctrlindex = find(CNOT_ctrl(:,i)==j);
                targindex = find(CNOT_targ(:,i)==j);
                NEWangles_SQG(ctrlindex,end-3) = -pi/2;
                NEWangles_SQG(targindex,end-5) = pi/2;
                NEWangles_SQG(targindex,end-4) = pi/2;
                NEWangles_SQG(targindex,end-3) = pi/2;
                NEWangles_SQG(targindex,end-1) = pi/2;
                NEWangles_SQG(targindex,end) = pi/2;
                angles_2QG(min([ctrlindex targindex]),max([ctrlindex targindex]),end-2) = pi/4;
            end
        otherwise
            error("Invalid gate_type - "+gate_type(i))
    end
end
gate_type = NEWgate_type;
angles = NEWangles_SQG;
clear NEWangles_SQG NEWgate_type


%---------------------------------------------------------------------------------------------------------------------------
% DECOMPOSE ALL ZZ BLOCKS IN TERMS OF DA BLOCKS USING ADRIÁN PR ALGORITHM
%---------------------------------------------------------------------------------------------------------------------------
% For this, we need to simulate arbitrary Hamiltonians with some source Hamiltonian
% Because of the way we did the decomposition of the problem Hamiltonian in vertical and horizantal terms, the topology of the zz Hamiltonians we need to simulate "coincide" with the hardware topology
% For the topologies to coincide, we will supose a continious NN Hamiltonian that takes the parallel rows/columns in order, so that the d coefficient for qubits belonging to different rows/columns is 0
% In order to not have infinitum problems, we need that the terms that are not connected have a b coefficient equal to 0 (thus, 0/0=0)
% In this section, first we solve the system of equations from ep.2 in Adrián PR
% The formula is slightly modified so that the source Hamiltonian can have arbitrary coefficients
% Then, we will introduce the x gates along the evolution
block_type = strings(0);
NEWangles_SQG = [];
analog_time = [];
% Since the M matrix is the same for all the blocks, we can generate it first
% This M matrix is generated so that it only accounts for the pairs of qubits that are present in the resource system
% Thus, we have to keep track of a new mapping, from i,j to alpha and again to alpha'. We have that alpha'=[1:number of couplings], so that the system of equations is solvable.
[Mmat,pairs_list] = M_alphabeta(Hs);
for i = 1:length(gate_type)
    if strncmp(gate_type(i),"zz",2)       
        b = b_beta(angles_2QG(:,:,i),Hs);
        t_alpha = linsolve(Mmat,b)';
        % We purge the times that are too small (hopefully this does not intruduce a lot of error)
        t_alpha_validIndex = [];
        for j = 1:length(t_alpha)
            if abs(t_alpha(j)) > timeCutThres
                t_alpha_validIndex(end+1) = j;
            end
        end
        % Now we have to introduce the X gates in between the analog blocks
        % We check if we have negative times, and if so, we use apendix A to work around it
        tmin = min(t_alpha(t_alpha_validIndex));
        if tmin < 0            
        	t_alpha = t_alpha+abs(tmin);
            % We purge again the times to eliminate the new ceros
         	t_alpha_validIndex = [];
            for j = 1:length(t_alpha)
                if abs(t_alpha(j)) > timeCutThres
                    t_alpha_validIndex(end+1) = j;
                end
            end
        end     
        if sum(t_alpha>=0) ~= numel(t_alpha)
            warning("We have negative times for the analog blocks!")
        end
        % We write the analog blocks with the x gate sandwiches
        for j = t_alpha_validIndex
            block_type(end+1:end+3) = ["rx" "zz" "rx"];
            auxAngles = zeros(N,1);
            if pairs_list(1,j) ~= 0
                n = pairs_list(1,j);
                m = pairs_list(2,j);
                auxAngles(n) = pi;
                auxAngles(m) = pi;
            end
            NEWangles_SQG(:,end+1:end+3) = [auxAngles zeros(N,1) auxAngles];
            analog_time(end+1:end+3) = [0 t_alpha(j) 0];
        end
        % If we need to simulate a negative time evolution, we need an extra block without any sandwiching
        if tmin < 0
            block_type(end+1) = "zz";
            NEWangles_SQG(:,end+1) = zeros(N,1);
            analog_time(end+1) = abs(((N*(N-9)/2)+8)*tmin);
        end
    else
        % This is a local rotation block, so we dont need to do anything special
        block_type(end+1) = gate_type(i); 
        NEWangles_SQG(:,end+1) = angles(:,i);
        analog_time(end+1) = 0;
    end        
end
angles_SQG = NEWangles_SQG;
clear Mmat b t_alpha tmin angles_2QG i j k t_alpha_validIndex n m pairs_list gate_type


%---------------------------------------------------------------------------------------------------------------------------
% MERGE ALL ADJACENT ROTATION GATES OF THE SAME TYPE
%---------------------------------------------------------------------------------------------------------------------------
% We parse the block_type vector and group together all the adjacet blocks
% To avoid any problem with the pulses, we write all the angles so that for each rx gate, the angle of rotation is in the range [0,2pi)
% At this step we remove all indications of vertical of horizontal blocks (this can be reverted for further optimization in the future)
NEWblock_type = strings(0);
NEWangles_SQG = [];
NEWanalog_time = [];
for i = 1:length(block_type)
    if strncmp(block_type(i),"rx",2) || strncmp(block_type(i),"rz",2)
        auxblock = split(block_type(i));
        auxblock = auxblock(1);
        % We search for the same block ahead, and sum it to the current block
        NEWblock_type(end+1) = auxblock;
        NEWangles_SQG(:,end+1) = angles_SQG(:,i);
        NEWanalog_time(end+1) = 0;
        auxi = 1;
        while i+auxi <= length(block_type) && strncmp(block_type(i+auxi),auxblock,2)
            % We mark the merged blocks, if there are no adjacent blocks we leave it as it was
            block_type(i+auxi) = "delete";
            NEWangles_SQG(:,end) = NEWangles_SQG(:,end)+angles_SQG(:,i+auxi);
            auxi = auxi+1;
        end
    elseif block_type(i) == "delete"
        % If we encounter a deleted block, we skip it
        continue
    else
        % If we find any zz block, we do nothing, as it is imposible to have 2 zz blocks one close to another
        NEWblock_type(end+1) = block_type(i);
        NEWangles_SQG(:,end+1) = angles_SQG(:,i);
        NEWanalog_time(end+1) = analog_time(i);
    end
end
block_type = NEWblock_type;
% We put all rotation angles in the range [0,2pi)
angles = mod(NEWangles_SQG,2*pi);
analog_time = NEWanalog_time;
clear auxblock NEWblock_type NEWangles_SQG NEWanalog_time i auxi
            

%---------------------------------------------------------------------------------------------------------------------------
% WRITE EVERYTHING AS A STEPWISE DIGITAL ANALOG BLOCKS SCHEDULE
%---------------------------------------------------------------------------------------------------------------------------
% For this we take into account the length of the SQGs and the timing for the zz blocks
% Since the rz blocks can be done virtualy, these gates will the smallest possible gate time (this is, minstep, that is the step used for the numerical integration or the application of pulses in qiskit)
% The time variable marks the start of a SQG block
% The duration of a SQG is minstep for rz gates and rxgatetime for rx gates
% The gates are timed for the rx gates, the adjacent rz gates will be acomodated ASAP. Since the rz gates commute with the evolution of the zz Hamiltonian, we can do this trick.
time = zeros(1,length(block_type)+1);
angles = angles+zeros(N,1);
gate_type = block_type;
% We have to parse the block type vector
% Each time we encounter a zz block we advance the current time
% Then, we time all the adjacent SQGs so that the middle pint of the total blocks time is centered in current time
for i = 1:length(block_type)
    % Depending of the block type, its duration is different
    time(i+1) = time(i) + round(analog_time(i)/minstep)*minstep*(block_type(i)=="zz")+round(rxgatetime/minstep)*minstep*(block_type(i)=="rx")+round(rzgatetime/minstep)*minstep*(block_type(i)=="rz");
end
warning("on")
end




%---------------------------------------------------------------------------------------------------------------------------
%---------------------------------------------------------------------------------------------------------------------------
% AUXILIARY FUNCTIONS
%---------------------------------------------------------------------------------------------------------------------------
%---------------------------------------------------------------------------------------------------------------------------
function res = b_beta(angles,Hs)
% Auxiliary function that takes the Hamiltonian coefficients as matrices and returns the b_beta vector, and the pairs that are present in the system
T = length(angles);
MatInc = Hs~=0;
res = zeros(T*(T-1)/2,1);                                                               
for i = 1:T-1
    for j = i+1:T
        if Hs(i,j) == 0
            res(T*(i-1)-i*(i+1)/2+j) = 0;
        else
            res(T*(i-1)-i*(i+1)/2+j) = angles(i,j)/Hs(i,j);
        end
    end
end
% We need to take out all interactions that are not in the hardware, and only run the algorithm in those we can control
delete_list = [];
for i = 1:T-1
    for j = i+1:T
        if ~MatInc(i,j)
            delete_list(end+1) = T*(i-1)-i*(i+1)/2+j;
        end
    end
end
for i = sort(delete_list,'descend')
    res(i) = [];
end
end
%---------------------------------------------------------------------------------------------------------------------------
function [res,pair_list] = M_alphabeta(Hs)
% Auxiliary function to write the M matrix from eq.3 in Adrian PR
T = length(Hs);
res = zeros(T*(T-1)/2);
% Case for N==4 ATA
if length(Hs) == 4 && sum(sum(triu(Hs~=0))) == 6
    res = [ 1 -1 -1 1 1 -1;...
            1 -1 1 -1 -1 1;...
            1 -1 1 1 -1 -1;...
            1 1 -1 -1 -1 -1;...
            1 1 -1 1 -1 1;...
            1 1 1 -1 1 -1]; 
    pair_list = [0 0; 1 1; 2 2; 3 3; 1 2; 2 3];
    return
end
% Generate the M matrix as in Adrian PR
for n = 1:T-1
    for m = n+1:T
        alpha = T*(n-1)-n*(n+1)/2+m;
        for i = 1:T-1
            for j = i+1:T
                beta = T*(i-1)-i*(i+1)/2+j;
                res(alpha,beta) = (-1)^((n==i)+(n==j)+(m==i)+(m==j));
            end
        end
    end
end
% We need to take out all interactions that are not in the hardware, and only run the algorithm in those we can control
delete_list = [];
pair_list = [];
MatInc = Hs~=0;
for i = 1:T-1
    for j = i+1:T
        if ~MatInc(i,j)
            delete_list(end+1) = T*(i-1)-i*(i+1)/2+j;
        else
            pair_list(:,end+1) = [i; j];
        end
    end
end
for i = sort(delete_list,'descend')
    res(i,:) = [];
    res(:,i) = [];
end
end
