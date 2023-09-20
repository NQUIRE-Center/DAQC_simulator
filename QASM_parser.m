function [gate_type,angles,CNOT_ctrl,CNOT_targ,measqubits] = QASM_parser(filename)
%QASM_PARSER parse QASM text file into an input adequate for D_simulator
%The program accepts QASM files with rx, rz and CNOT gates
%The parser is naive in the sense that it does not parallelize the gates (yet)
%   filename        Name of the text file
%   
%   OUTPUT:
%    gate_type              Vector of strings that defines what type of block is the next. It can be CNOT, rz or rx
%    angles                 Matrix of the angles of rotations of the SQG. Each row corresponds to one qubit, and the columns match the ones in gate_type
%    CNOT_(ctrl/targ)       Matrix that labels the control/target qubits of the CNOT gates
%    measqubits             List of qubits that are measured at the end of the circuit

% Open the file
fileID = fopen(filename,"r");

% Initialize the output
gate_type = strings(0);
angles = [];
CNOT_ctrl = [];
CNOT_targ = [];

measqubits = [];

% Read the text file line by line
while true
    line = fgetl(fileID);
    if line==-1
        % Stop at the end of file
        break
    elseif strncmp(line,"#",1)
        % We dont read the comments
        continue
    end
    % We split the line as it is intended to be read 
    elems = split(line,["(",")"," ","q[","],q[","];","] -> c[","], ",") "]);
    % We clean elems of empty strings
    elems = elems(elems~="");
    if isempty(elems)
        % Skip empty lines
        continue
    end
    % Gate parser
    switch elems{1}
        case "qreg"
            % Obtain the number of qubits of the system
            N = str2double(elems{2});
        case {"rx","rz"}
            gate_type(end+1) = elems{1}; %#ok<*AGROW>
            angles(:,end+1) = zeros(N,1);
            angles(str2double(elems{3})+1,end) = eval(elems{2});
            CNOT_ctrl(:,end+1) = zeros(N,1);
            CNOT_targ(:,end+1) = zeros(N,1);
        case "ry"
            gate_type(end+1:end+3) = ["rz" "rx" "rz"];
            angles(:,end+1:end+3) = zeros(N,3);
            angles(str2double(elems{3})+1,end-2:end) = [-pi/2 eval(elems{2}) pi/2];
            CNOT_ctrl(:,end+1:end+3) = zeros(N,3);
            CNOT_targ(:,end+1:end+3) = zeros(N,3);
        case {"x","z"}
            gate_type(end+1) = "r"+elems{1};
            angles(:,end+1) = zeros(N,1);
            angles(str2double(elems{2})+1,end) = pi;
            CNOT_ctrl(:,end+1) = zeros(N,1);
            CNOT_targ(:,end+1) = zeros(N,1);
        case "h"
            gate_type(end+1:end+3) = ["rz" "rx" "rz"];
            angles(:,end+1:end+3) = zeros(N,3);
            angles(str2double(elems{2})+1,end-2:end) = pi/2;
            CNOT_ctrl(:,end+1:end+3) = zeros(N,3);
            CNOT_targ(:,end+1:end+3) = zeros(N,3);
        case "cx"
            gate_type(end+1) = "CNOT";
            angles(:,end+1) = zeros(N,1);
            CNOT_ctrl(:,end+1) = zeros(N,1);
            CNOT_targ(:,end+1) = zeros(N,1);
            CNOT_ctrl(str2double(elems{2})+1,end) = 1;
            CNOT_targ(str2double(elems{3})+1,end) = 1;
        case "sx"
            gate_type(end+1) = "rx";
            angles(:,end+1) = zeros(N,1);
            angles(str2double(elems{2})+1,end) = pi/2;
            CNOT_ctrl(:,end+1) = zeros(N,1);
            CNOT_targ(:,end+1) = zeros(N,1);
        case "measure"
            measqubits(end+1) = str2double(elems{2})+1;            
        % Here we can add more gates to parse        
        
        % These are the instructions we dont care about
        case {"include","creg","OPENQASM"}
            continue
        otherwise
            warning(elems{1}+" instruction not implemented")
    end
end

%---------------------------------------------------------------------------------------------------------------------------
% MERGE ALL ADJACENT ROTATION GATES OF THE SAME TYPE
%---------------------------------------------------------------------------------------------------------------------------
% We parse the gate_type vector and group together all the adjacet blocks
% To avoid any problem with the pulses, we write all the angles so that for each rx gate, the angle of rotation is in the range [0,2pi)
% At this step we remove all indications of vertical of horizontal blocks (this can be reverted for further optimization in the future)
NEWgate_type = strings(0);
NEWangles = [];
NEWCNOT_ctrl = [];
NEWCNOT_targ = [];
for i = 1:length(gate_type)
    if strncmp(gate_type(i),"rx",2) || strncmp(gate_type(i),"rz",2)
        auxblock = split(gate_type(i));
        auxblock = auxblock(1);
        % We search for the same block ahead, and sum it to the current block
        NEWgate_type(end+1) = auxblock;
        NEWangles(:,end+1) = angles(:,i);
        NEWCNOT_ctrl(:,end+1) = zeros(N,1);
        NEWCNOT_targ(:,end+1) = zeros(N,1);
        auxi = 1;
        while i+auxi <= length(gate_type) && strncmp(gate_type(i+auxi),auxblock,2)
            % We mark the merged blocks, if there are no adjacent blocks we leave it as it was
            gate_type(i+auxi) = "delete";
            NEWangles(:,end) = NEWangles(:,end)+angles(:,i+auxi);
            auxi = auxi+1;
        end
    elseif gate_type(i) == "delete"
        % If we encounter a deleted block, we skip it
        continue
    else
        % If we find any zz block, we do nothing, as it is imposible to have 2 zz blocks one close to another
        NEWgate_type(end+1) = gate_type(i);
        NEWangles(:,end+1) = angles(:,i);
        NEWCNOT_ctrl(:,end+1) = CNOT_ctrl(:,i);
        NEWCNOT_targ(:,end+1) = CNOT_targ(:,i);
    end
end
% Reasing the variables
gate_type = NEWgate_type;
% We put all rotation angles in the range [0,2pi)
angles = mod(NEWangles,2*pi);
CNOT_ctrl = NEWCNOT_ctrl;
CNOT_targ = NEWCNOT_targ;
clear auxblock NEWgate_type NEWangles NEWanalog_time i auxi
end