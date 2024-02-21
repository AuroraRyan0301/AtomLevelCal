oma=1; omb=1; % express omeg rabi freqs, in units of gamma
dels=0; % common detuning set to zero
N=3; % number of energy levels
R=401; % number of points to plot
%initialize and set dimensions for all matrices
del=zeros(1,R); %diff detuning array
M=zeros(N^2,N^2); %M-matrix
rho=zeros(N,N); %density matrix
Ham=zeros(N,N); %Hamiltonian with decay
Q=zeros(N,N); %matrix representing derivative of density matrix
W=zeros((N^2-1),(N^2-1)); %W-matrix
S=zeros((N^2-1),1); %S-vector
B=zeros((N^2-1),1); %B-vector
A=zeros(N^2,R); %A-vectors, for all detunings
for m=1:R %start the overall-loop
    del(1,m)=(m-(R+1)/2)/10; %define the detuning
    Ham=[del(1,m)/2 0 oma/2; 0 del(1,m)*(-1)/2 omb/2; ...
        oma/2 omb/2 (dels+0.5i)*(-1)];
    for n=1:N^2 %start the outer-loop for finding elements of M;
        for p=1:N^2 %start inner-loop for finding elements of M;
            %finding alpha and beta
            remain=rem(n,N);
            if remain==0
                beta=N;
            else beta=remain;
            end
            alpha=(1+(n-beta)/N);
            %finding epsilon and sigma
            remain=rem(p,N);
            if remain==0
            sigma=N;
            else sigma=remain;
            end
            epsilon=(1+(p-sigma)/N);
            rho=zeros(N,N); %reset rho to all zeros
            rho(epsilon,sigma)=1; %pick one element to unity
            Q=(Ham*rho-rho*conj(Ham))*(0-1i); %find first part of Q matrix
            disp(n);
            disp(p);
            Q(1,1)=Q(1,1)+rho(3,3)/2; %add pop source term to Q
            Q(2,2)=Q(2,2)+rho(3,3)/2; %add pop source term to Q
            %Modify as needed for general %systems
            M(n,p)=Q(alpha,beta);
        end %end the inner-loop for finding elements of M
    end %end of the outer-loop for finding elements of M
    S=M(1:(N^2-1),N^2:N^2); %find S-vector
    W=M(1:(N^2-1),1:(N^2-1)); %initialize W-matrix

    for d=1:(N-1)
    W(:,((d-1)*N+d))=W(:,((d-1)*N+d))-S; %update W by subtracting
    %from selected columns
    end
    B=(W\S)*(-1); %find B-vector: primary solution
    rhonn=1; %initialize pop of N-th state
    %determine pop of N-th state
    for f=1:(N-1)
    rhonn=rhonn-B(((f-1)*N+f), 1);
    end
    %determine elements of A vector
    A(1:(N^2-1),m)=B;
    A(N^2,m)=rhonn;
end %end of over-all loop
plot(del,real( A ( (N^2-0),: ) ) )