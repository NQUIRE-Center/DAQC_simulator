function [rho,rhomeas,Uideal] = sDA_simulator(gate_type,time,angles,Hs,measqubits,rxgatetime,rzgatetime,rD,rB,pbf,pmeas,pth,T1,T2,digitalOnlyBF,rho)
%DA_SIMULATOR Function to simulate the application of a stepwise Digital-Analog circuit in a noisy platform.
%   gate_type               Type of SQG that is applied.
%   time                    Times at which the digital blocks starts. In units of s.
%   angles                  Angles for the SQGs.
%   Hs                      Coefficients of the source Hamiltonian. In units of 1/2pi*s.
%   measqubits              List of the indeces of the qubits that are measured at the end of the circuit.
%   (rz/rx)gatetime         Time for the application of the corresponding SQG. In units of s.
%   rD                      Noise: deviation ratio for the SQGs, comming from an uniform magnetic field noise.
%   rB                      Noise: deviation ratio for the source Hamiltonian, comming from normal fluctuations in the coupling strength.
%   pbf                     Noise: bit-flip error during the rz gate. Here we assume that rzgatetime<rxgatetime
%   pmeas                   Noise: measurement error modeled as bit-flip errors
%   pth                     Noise: thermal popilation of the ground state of the qubits
%   T1                      Noise: decoherence time. In units of s.
%   T2                      Noise: dephasing time. In units of s.
%   digitalOnlyBF           Optional: if true, the bit flip error are only applied to the digital blocks (default=false)
%   rho                     Optional: Initial state for the circuit, given as a density matrix
%
%   OUTPUT:
%    rho                    Final density matrix traced over measqubits.
%    rhomeas                Density matrix after simulating a noisy measurement traced over measqubits.
%    Uideal                 Unitary matrix of the ideal circuit
% Note: this program could be further optimized with gpuarrays

% Number of qubits
N = size(angles,1);

% Initialize density matrix (if not already given)
if ~exist('rho','var')
    rho = zeros(2^N);
    rho(1,1) = 1;
end

% Initilize digitalOnlyBf variable
if ~exist('digitalOnlyBF','var')
    digitalOnlyBF = false;
end

% Initialization of the unitary
Uideal = eye(2^N);

% Auxiliary parameter for the timing error
timingError = max([rxgatetime rzgatetime]);%Este parámetro igual es mejor tomarlo como un parámetro de entrada

% Generate the vector with the diagonal of the source Hamiltonian
Hs = genHs(Hs);

% We can precalculate the bit-flip channel renormalization matrix, exactly for the SQGs
% This is trivial because the renormalization matrix is just a constant
pbf_rx = 0;
for i = 1:2:ceil(rxgatetime/rzgatetime)
    pbf_rx = pbf_rx*nchoosek(ceil(rxgatetime/rzgatetime),i)*pbf^i*(1-pbf)^(ceil(rxgatetime/rzgatetime)-i);
end
if pbf_rx == 0
    pbf_rx = pbf;
end

% Check if there is any checkpointFile
% % % if exist('chkS.mat','file')==2
% % %     load('chkS.mat')
% % %     str = chk;
% % % else
str = 1;
% % % end

% We go through all the steps
for i = str:length(gate_type)
% % %     tic%kk
    if gate_type(i)=="rz" || gate_type(i)=="rx"
        %---------------------------------------------------------------------------------------------------------------------------
        % DIGITAL BLOCK
        %---------------------------------------------------------------------------------------------------------------------------
        % We calculate the digital block time
        deltaT = (gate_type(i)=="rz")*rzgatetime+(gate_type(i)=="rx")*rxgatetime;   
        % We calculate the unitary evolution with the corresponding unitary errors. In sDAQC the interaction is off while on a digital block
        %U = expm(-1i*(deltaT*normrnd(0,rB*deltaT)*Hs+gen_HSQG(angles(i),rD,deltaT,gate_type(i))));%Mirar a ver si esto es invent, o si tiene sentido simular estas fluctuaciones
        U = USQG(angles(:,i),rD,deltaT,gate_type(i));
        if gate_type(i)=="rz"
            % If the unitary matrix is diagonal, we wan make the operation much faster
            U = sparse(diag(U));
            rho = U*rho*U';
        else
            % We check if we can write this matrix as sparse to multiply it faster
            if nnz(U)/numel(U) < 1/7.9
                U = sparse(U);
            end
            % Update the new density matrix with the noisy unitary evolution and the ideal evolution operator
            rho = U*rho*U';
        end
        clear U
        % Update ideal evolution
        Uidealaux = Uideal_step(gate_type(i),angles(:,i),Hs,[]);
        % We check if we can write this matrix as sparse to multiply it faster
        if gate_type(i)=="rz"
            % If the unitary matrix is diagonal, we wan make the operation much faster
            Uidealaux = sparse(diag(Uidealaux));
            Uideal = Uidealaux*Uideal;
        else
            if nnz(Uidealaux)/numel(Uidealaux) < 1/7.9
                Uidealaux = sparse(Uidealaux);
            end
            Uideal = Uidealaux*Uideal;
        end
        clear Uidealaux         

        % QUANTUM CHANNELS
        if pbf > 0
            rho = bitFlipChannel(rho,pbf,pbf_rx,gate_type(i));   
        end

        if pth == 0 && T1 < Inf
            if pth == 0
                rho = amplitudeDampingChannel(rho,deltaT,T1);
            elseif pth > 0
                rho = generalizedAmplitudeDampingChannel(rho,deltaT,T1,pth);
            end
        end

        if T2>0
            rho = dephasingChannel(rho,deltaT,T2);
        end


    elseif gate_type(i)=="zz"
        %-----------------------------------------------------------------------------------------------------------------------
        % ANALOG BLOCK
        %-----------------------------------------------------------------------------------------------------------------------
        % We calculate the analog block time
        deltaT = time(i+1)-time(i);
        % We calculate the diagonal of the unitary evolution with the corresponding unitary error
        U = spdiags(exp(-1i*normrnd(deltaT,rB*timingError)*diag(Hs)),0,2^N,2^N);
        % We update the new density matrix with the noisy unitary evolution
        rho = U*rho*U';
        clear U
        % Update ideal evolution
        Uidealaux = sparse(Uideal_step(gate_type(i),angles(:,i),Hs,deltaT));
        % The unitary matrix is diagonal, we wan make the operation much faster
        Uideal = Uidealaux*Uideal;
        clear Uidealaux
        
        % QUANTUM CHANNELS
        if pbf > 0 && ~digitalOnlyBF
            rho = bitFlipChannel(rho,pbf,pbf_rx,gate_type(i));   
        end

        if pth == 0 && T1 < Inf
            if pth == 0
                rho = amplitudeDampingChannel(rho,deltaT,T1);
            elseif pth > 0
                rho = generalizedAmplitudeDampingChannel(rho,deltaT,T1,pth);
            end
        end        

        if T2 < Inf
            rho = dephasingChannel(rho,deltaT,T2);
        end
    end    


    if abs(1-sum(diag(rho)))>1e-4
        warning("Error in tr(rho) > 10^-4",'traza')
        warning("off",'last')
    end
    
% % %     aux = toc;%kk
% % %     fprintf("%.2f%% %.2fs\n",[100*i/length(gate_type) aux])%kk
% % %     % Checkpoint each 10 iterations
% % %     if mod(i,10)==0
% % %         chk = i+1;
% % %         save('chkS.mat','-v7.3')
% % %     end
end
% Delete checkpoint after finishing
% % % delete chkS.mat
clear rhoaux

%-----------------------------------------------------------------------------------------------------------------------
% MEASUREMENT
%-----------------------------------------------------------------------------------------------------------------------
% We have to trace the qubits that are not measured
% It is more efficient to trace out the qubits one by one, (2^(5+2N)-2^(5+2N-2m)/3 VS 2^(4m+2N+1)
auxN = N;
rhomeas = rho;
for delIndex = flip(setdiff(1:N,measqubits))
    % We trace the last possible qubit, so we don't mess the order of the qubits
    rhomeas = kron(eye(2^(delIndex-1)),kron([1 0],eye(2^(auxN-delIndex))))*rhomeas*kron(eye(2^(delIndex-1)),kron([1; 0],eye(2^(auxN-delIndex))))...
            + kron(eye(2^(delIndex-1)),kron([0 1],eye(2^(auxN-delIndex))))*rhomeas*kron(eye(2^(delIndex-1)),kron([0; 1],eye(2^(auxN-delIndex))));    
    auxN = auxN-1;
end

if pmeas > 0
    % MEASUREMENT ERROR (BIT-FLIP CHANNEL)
    rhomeas = bitFlipChannel(rhomeas,pmeas,[],"meas");
end
end




%---------------------------------------------------------------------------------------------------------------------------
%---------------------------------------------------------------------------------------------------------------------------
% AUXILIARY FUNCTIONS
%---------------------------------------------------------------------------------------------------------------------------
%---------------------------------------------------------------------------------------------------------------------------
function res = USQG(angles,rD,deltaT,gate_type)
% Auxiliary function to generate the noisy unitary evolution implementing the gates
N = length(angles);
if gate_type == "rx"
    res = 1;
    for i = 1:N
        if angles(i) == 0
            % If we don't apply any gate we assume we have exactly the identity
            res = kron(res,eye(2));
        else
            theta = angles(i) + (rand()-1/2)*rD;
            res = kron(res,[cos(theta/2) -1i*sin(theta/2);-1i*sin(theta/2) cos(theta/2)]);
        end
    end
elseif gate_type == "rz"
    res = 1;
    for i = 1:N
        if angles(i) == 0
            % If we don't apply any gate we assume we have exactly the identity
            res = kron(res,ones(1,2));
        else
            theta = angles(i) + (rand()-1/2)*rD;
            res = kron(res,[exp(-1i*theta/2) exp(1i*theta/2)]);
        end
    end
else
    warning("Invalid single qubit rotation gate")
end
end
%---------------------------------------------------------------------------------------------------------------------------
function res = genHs(H)
% Auxiliary function to write the full source Hamiltonian from the coefficients matrix
N = length(H);
res = zeros(1,2^N);
for i = 1:N-1
    for j = i+1:N
        if H(i,j)~=0
            res = res + H(i,j)*kron(ones(1,2^(i-1)),kron([1 -1],kron(ones(1,2^(j-i-1)),kron([1 -1],ones(1,2^(N-j))))));
        end
    end
end
res = diag(res);
end
%---------------------------------------------------------------------------------------------------------------------------
function res = Uideal_step(gate_type,angles,Hs,time)
% Auxiliary function to generate the unitary ideal evolution of the sDAQC circuit
if gate_type=="rx"
    res = 1;
    for theta = angles'
        res = kron(res,[cos(theta/2) -1i*sin(theta/2);-1i*sin(theta/2) cos(theta/2)]);
    end
elseif gate_type=="rz"
    res = 1;
    for theta = angles'
        res = kron(res,[exp(-1i*theta/2) exp(1i*theta/2)]);
    end
elseif gate_type=="zz"
    res = diag(exp(-1i*time*diag(Hs)));
else
    disp(gate_type)
    error("Invalid gate")
end
end
%---------------------------------------------------------------------------------------------------------------------------
function rho = bitFlipChannel(rho,pbf,pbf_rx,gate_type)
% BIT-FLIP CHANNEL
    N = round(log2(height(rho)));
    % We need first the bit flip probability on the current block (this has already been precalculated)
    if gate_type == "rx"
        PBF = pbf_rx;
    else
        PBF = pbf;
    end
    % We divide the full channel, into N channels for each qubit
    for k = 1:N
        Eaux = sqrt(pbf)*kron(speye(2^(k-1)),kron(sparse([0 1;1 0]),speye(2^(N-k))));
        % We apply the E0 and the E1 operators 
        rho = (1-pbf)*rho + Eaux*rho*Eaux;
    end
end
%---------------------------------------------------------------------------------------------------------------------------
function rho = amplitudeDampingChannel(rho,deltaT,T1)
% AMPLITUDE-DAMPING CHANNEL
    N = round(log2(height(rho)));
    % Calculate the gamma for the current block
    gamma = 1-exp(-deltaT/T1);
    % We calculate the Krauss operators
    E0 = sparse([1 0;0 sqrt(1-gamma)]);
    E1 = sparse([0 sqrt(gamma);0 0]);
    % We divide the full channel, into N channels for each qubit
    for k = 1:N
        Eaux0 = kron(speye(2^(k-1)),kron(E0,speye(2^(N-k))));
        Eaux1 = kron(speye(2^(k-1)),kron(E1,speye(2^(N-k))));
        rho = Eaux0*rho*Eaux0 + Eaux1*rho*Eaux1';
    end
end
%---------------------------------------------------------------------------------------------------------------------------
function rho = generalizedAmplitudeDampingChannel(rho,deltaT,T1,pth)
% GENERALIZED AMPLITUDE-DAMPING CHANNEL
    N = round(log2(height(rho)));
    % Calculate the gamma for the current block
    gamma = 1-exp(-deltaT/T1);
    % Calculate the Krauss operators
    E0 = sqrt(pth)*sparse([1 0;0 sqrt(1-gamma)]);
    E1 = sqrt(pth)*sparse([0 sqrt(gamma);0 0]);
    E2 = sqrt(1-pth)*sparse([sqrt(1-gamma) 0;0 1]);
    E3 = sqrt(1-pth)*sparse([0 0;sqrt(gamma) 0]);
    % Then, we apply the N*2^N combinations of a single sqrt(p) operator and the rest sqrt(1-p) operators
    for k = 1:N % Run over the qubits
        Eaux0 = kron(speye(2^(k-1)),kron(E0,speye(2^(N-k))));
        Eaux1 = kron(speye(2^(k-1)),kron(E1,speye(2^(N-k))));
        Eaux2 = kron(speye(2^(k-1)),kron(E2,speye(2^(N-k))));
        Eaux3 = kron(speye(2^(k-1)),kron(E3,speye(2^(N-k))));
        rho = Eaux0*rho*Eaux0' + Eaux1*rho*Eaux1' + Eaux2*rho*Eaux2' + Eaux3*rho*Eaux3';
    end
end
%---------------------------------------------------------------------------------------------------------------------------
function rho = dephasingChannel(rho,deltaT,T2)
% DEPHASING CHANNEL
    N = round(log2(height(rho)));
    % Calculate the gamma for the current block
    alpha = (1+exp(-deltaT/(2*T2)))/2;
    % We calculate the Krauss operators
    E0 = sqrt(alpha)*sparse([1 0;0 1]);
    E1 = sqrt(1-alpha)*sparse([1 0;0 -1]);
    % We divide the full channel, into N channels for each qubit
    for k = 1:N
        Eaux0 = kron(speye(2^(k-1)),kron(E0,speye(2^(N-k))));
        Eaux1 = kron(speye(2^(k-1)),kron(E1,speye(2^(N-k))));
        rho = Eaux0*rho*Eaux0 + Eaux1*rho*Eaux1';
    end
end







% %---------------------------------------------------------------------------------------------------------------------------
% % OLD FUNCTIONS FOR THE CHANNELS
% %---------------------------------------------------------------------------------------------------------------------------
% function rhoaux = bitFlipChannel(rho,pbf,pbf_rx,gate_type)
% % BIT-FLIP CHANNEL
%     N = round(log2(height(rho)));
%     % We need first the bit flip probability on the current block (this has already been precalculated)
%     if gate_type == "rx"
%         PBF = pbf_rx;
%     else
%         PBF = pbf;
%     end
%     % The renormalization matrix for the bit-flip channel is proportional to the identity matrix
%     OmegaBF = (1-PBF)^N+N*PBF*(1-PBF)^(N-1);
%     % We use this renormalization matrix to "renormalize" the density matrix beforehand
%     rho = rho/OmegaBF;
%     % We first apply the identity operator (no bit-flip), E0 on every qubit
%     rhoaux = rho*(1-pbf)^N;
%     % Then, the operators where only one bit flip happend
%     for k = 1:N
%         Eaux = sqrt(pbf)*sqrt(1-pbf)^(N-1)*kron(speye(2^(k-1)),kron(sparse([0 1;1 0]),speye(2^(N-k))));
%         rhoaux = rhoaux + Eaux*rho*Eaux;
%     end
% end
% %---------------------------------------------------------------------------------------------------------------------------
% function rhoaux = amplitudeDampingChannel(rho,deltaT,T1)
% % AMPLITUDE-DAMPING CHANNEL
%     N = round(log2(height(rho)));
%     % We can write write this as the combinations of all E_0 E_1 operations
%     % But since gamma is small, we can apply the same cutoff as for the bit-flip channel
%     % Calculate the gamma for the current block
%     gamma = 1-exp(-deltaT/T1);
%     % Obtain the renormalization matrix (this is a diagonal matrix)
%     aux = zeros(1,2^N);
%     E0 = [1 1-gamma];
%     E1 = [0 gamma];
%     auxx = 1;
%     for j = 1:N
%         auxx = kron(auxx,E0);
%     end
%     aux = aux + auxx;
%     for j = 1:N
%         auxx = 1;
%         for k = 1:N
%             if j == k
%                 auxx = kron(auxx,E1);
%             else
%                 auxx = kron(auxx,E0);
%             end
%         end
%         aux = aux + auxx;
%     end
%     % The renormalization matrix is the inverse of the square root of this sum
%     aux = sparse(diag(1./sqrt(aux)));
%     % Apply the renormalization
%     rho = aux*rho*aux;
%     clear aux auxx 
%     % We apply first the krauss operator with E0
%     E0 = sparse([1 0;0 sqrt(1-gamma)]);
%     E1 = sparse([0 sqrt(gamma);0 0]);
%     Eaux = 1;
%     for j = 1:N
%         Eaux = kron(Eaux,E0);
%     end
%     rhoaux = Eaux*rho*Eaux;
%     % Then, we apply the krauss operators with one instance of a sqrt(gamma) operator
%     for j = 1:N
%         Eaux = 1;
%         for k = 1:N
%             if k == j
%                 Eaux = kron(Eaux,E1);
%             else
%                 Eaux = kron(Eaux,E0);
%             end
%         end
%         rhoaux = rhoaux + Eaux*rho*Eaux';
%     end
% end
% %---------------------------------------------------------------------------------------------------------------------------
% function rhoaux = generalizedAmplitudeDampingChannel(rho,deltaT,T1,pth)
% % GENERALIZED AMPLITUDE-DAMPING CHANNEL
%     N = round(log2(height(rho)));
%     % Calculate the gamma for the current block
%     gamma = 1-exp(-deltaT/T1);
%     % Apply the renormalization matrix to the density matrix
%     aux = zeros(1,2^N);
%     % First the operators with probability sqrt(1-p)
%     for j = 0:2^N-1
%         auxx = 1;
%         for k = dec2bin(j,N)
%             if k == '0'
%                 auxx = kron(auxx,(1-pth)*[(1-gamma) 1]);
%             else
%                 auxx = kron(auxx,(1-pth)*[gamma 1]);
%             end
%         end
%         aux = aux+auxx;
%     end
%     % Then, we apply the N*2^N combinations of a single sqrt(p) operator and the rest sqrt(1-p) operators
%     for j = 1:N % Run over the qubits
%         for op = [0 1] % Run over the two possible sqrt(p) operators
%             for k = 0:2^(N-1) % Run over the combinations of the sqrt(1-p) operators
%                 % Obtain the corresponding combination of operators
%                 comb = dec2bin(k,N-1);
%                 auxx = 1;
%                 for l = 1:N
%                     if l < j
%                         if comb(l) == '0'
%                             auxx = kron(auxx,(1-pth)*[(1-gamma) 1]);
%                         else
%                             auxx = kron(auxx,(1-pth)*[gamma 1]);
%                         end
%                     elseif l == j
%                         if op
%                             auxx = kron(auxx,pth*[1 1-gamma]);
%                         else
%                             auxx = kron(auxx,pth*[0 gamma]);
%                         end
%                     else
%                         if comb(l-1) == '0'
%                             auxx = kron(auxx,(1-pth)*[(1-gamma) 1]);
%                         else
%                             auxx = kron(auxx,(1-pth)*[gamma 1]);
%                         end
%                     end
%                 end
%                 aux = aux + auxx;
%             end
%         end
%     end
%     % The renormalization matrix is the inverse of the square root of this sum
%     aux = sparse(diag(1./sqrt(aux)));
%     % Apply the renormalization
%     rho = aux*rho*aux;
%     clear aux auxx 
%     % We first apply all the 2^N combinations of operators with probability sqrt(1-p)
%     rhoaux = zeros(2^N);
%     for j = 0:2^N-1
%         % Obtain the corresponding combination of operators
%         comb = dec2bin(j,N);
%         % Construct the operator by running over the character variable
%         Eaux = sparse(1);
%         for k = comb
%             if k == '0'
%                 Eaux = kron(Eaux,sparse((1-pth)*[sqrt(1-gamma) 0;0 1]));
%             else
%                 Eaux = kron(Eaux,sparse((1-pth)*[0 0;sqrt(gamma) 0]));
%             end
%         end
%         % Add the action to the sum
%         rhoaux = rhoaux + Eaux*rho*Eaux';
%     end
%     % Then, we apply the N*2^N combinations of a single sqrt(p) operator and the rest sqrt(1-p) operators
%     for j = 1:N % Run over the qubits
%         for op = [0 1] % Run over the two possible sqrt(p) operators
%             for k = 0:2^(N-1) % Run over the combinations of the sqrt(1-p) operators
%                 % Obtain the corresponding combination of operators
%                 comb = dec2bin(k,N-1);
%                 Eaux = sparse(1);
%                 for l = 1:N
%                     if l < j
%                         if comb(l) == '0'
%                             Eaux = kron(Eaux,sparse((1-pth)*[sqrt(1-gamma) 0;0 1]));
%                         else
%                             Eaux = kron(Eaux,sparse((1-pth)*[0 0;sqrt(gamma) 0]));
%                         end
%                     elseif l == j
%                         if op
%                             Eaux = kron(Eaux,sparse(pth*[1 0;0 sqrt(1-gamma)]));
%                         else
%                             Eaux = kron(Eaux,sparse(pth*[0 sqrt(gamma);0 0]));
%                         end
%                     else
%                         if comb(l-1) == '0'
%                             Eaux = kron(Eaux,sparse((1-pth)*[sqrt(1-gamma) 0;0 1]));
%                         else
%                             Eaux = kron(Eaux,sparse((1-pth)*[0 0;sqrt(gamma) 0]));
%                         end
%                     end
%                 end
%                 rhoaux = rhoaux + Eaux*rho*Eaux';
%             end
%         end
%     end
% end