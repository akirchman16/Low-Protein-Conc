clearvars;
close all;

% This code is meant to track a single protein (RAD51 only; no
% competition) that binds to and unbinds from a ssDNA lattice at low
% concentrations. Rather than having a constant concentration over time,
% individual protein populations will be tracked. This utilizes the SSA for
% stochastic simulations. No cooperativity will be included in this model
% (a future model will include cooperativity).

N = 1000;    %length of ssDNA
n = 3;  %length of protein (RAD51)

MaxEvents = 100;   %maximum number of events to occur
InitialFree = 100;    %initial number of proteins that are around the lattice

k_on = 10;   %kinetic rate constant for protein binding
k_off = 1;  %kinetic rate constant for protein unbinding

% Memory Allocation
X_Locations = zeros(1,MaxEvents);   %number of free locations on the DNA lattice over time
X_Free = zeros(1,MaxEvents);    %number of free proteins with each event
X_Free(1) = InitialFree;    %records initial amount of free proteins
X_Bound = zeros(1,MaxEvents);   %number of bound proteins on the lattice with each event (lattice is initially empty)
Populations = zeros(3,MaxEvents);   %combined array tracking populations of bound/free proteins and free locations
a_P = zeros(2,MaxEvents);   %propensity functions
t = zeros(1,MaxEvents+1); %time tracker
dt = zeros(1,MaxEvents);    %records each time interval in the simulation
EventHistory = zeros(2,MaxEvents);  %records where each event occurs at each step
FracCover = zeros(1,MaxEvents+1);   %tracks saturation of ssDNA lattice

DNA = zeros(1,N);   %ssDNA lattice
BoundAtSpot = zeros(1,N);   %used to record where proteins are bound

Event = 0;  %counts number of events which have occured
for i = 1:MaxEvents
    Event = Event+1;    %increases event counter
    
    % Location Counter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    AvailableLocations = [];    %initializes an array to track available locations
    for Loc = 1:N-(n-1)   %currently checks each DNA location individually (Should be optimized in the future)
        if DNA(Loc:Loc+(n-1)) == 0 %if location is free...
            AvailableLocations = [AvailableLocations,Loc];  %...add location to the list of available ones
        end
    end
    X_Locations(Event) = numel(AvailableLocations); %number of free locations on the lattice
    X_Bound(Event) = numel(find(BoundAtSpot == 1)); %number of proteins bound to the lattice
    BoundLocations = find(BoundAtSpot == 1);    %list of locations where proteins are bound
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    Populations(:,Event) = [X_Free(Event);X_Bound(Event);X_Locations(Event)];   %tracks population of each species
    a_P(:,Event) = [k_on*Populations(1,Event)*Populations(3,Event);k_off*Populations(2,Event)]; %propensity functions of each reaction
        % 1 - Binding Reaction
        % 2 - Unbinding Reaction
    Rands = [rand,rand];    %random numbers for Monte Carlo step
    dt(Event) = (1/sum(a_P(:,Event)))*log(1/Rands(1));   %time step calculation
    
    if Rands(2)*sum(a_P(:,Event)) < a_P(1,Event)    %binding reactions (via SSA)
        BindLocation = AvailableLocations(randi(X_Locations(Event)));   %selects random location for binding
        DNA(BindLocation:BindLocation+(n-1)) = 1;   %binds protein
        X_Free(Event+1) = X_Free(Event)-1;  %updates free protein tracker (decrease by 1)
        BoundAtSpot(BindLocation) = 1;  %records where proteins are bound
        EventHistory(1,Event) = BindLocation;   %records where event happened
    else   %unbinding reaction
        UnbindLocation = BoundLocations(randi(X_Bound(Event)));     %selects random protein to unbind from lattice
        DNA(UnbindLocation:UnbindLocation+(n-1)) = 0;   %unbinds protein
        X_Free(Event+1) = X_Free(Event)+1;  %updates free protein tracker (increase by 1)
        BoundAtSpot(UnbindLocation) = 0;    %updates BoundAtSpot
        EventHistory(2,Event) = UnbindLocation;
    end
    
    t(Event+1) = t(Event)+dt(Event);    %updates timer by advancing time according to the calculated time step
    FracCover(Event+1) = numel(find(DNA == 1))/N;   %calcualtes saturation level of ssDNA lattice
    if mod(numel(find(DNA == 1)),n) ~= 0  %if there isn't an integer number of bound proteins...
        disp('BROKEN PROTEIN');
        break  %ERROR STOP
    end
end

figure(1);
subplot(2,2,[3,4]);
scatter(t,FracCover,5,'r','filled');
hold on;
yline(n*InitialFree/N,'--k','Maximum Saturation');
xlim([0 t(end)]);
xlabel('Time, t');
ylim([0 1]);
ylabel('Saturation Level');
title(['Low Protein Population (Initial X: ', num2str(InitialFree), ')']);
box on;
subplot(2,2,1);
scatter(t(1:end-1),X_Bound,5,'r','filled');
hold on;
scatter(t,X_Free,5,'b','filled');
xlim([0 t(end)]);   xlabel('Time, t');
ylim([0 InitialFree]); ylabel('Population, X');
title('Protein Populations');
box on;
legend('Bound','Free','Location','east');
subplot(2,2,2);
scatter(t(1:end-1),X_Locations,5,'g','filled');
hold on;
xlim([0 t(end)]); xlabel('Time, t');
ylim([0 N]);    ylabel('Free Locations');
title('Available Locations');
box on;