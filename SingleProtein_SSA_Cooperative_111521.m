clearvars;
close all;

% This code is meant to track a single protein (RAD51 only, no
% competition) that binds to and unbinds from a ssDNA lattice at low
% concentrations. Rather than having a constant concentration over time,
% individual protein populations will be tracked. This utilizes the SSA for
% stochastic simulations. Now cooperativity will be included, but still no
% competition.

N = 50;    %length of ssDNA
n = 3;  %length of protein (RAD51)

MaxEvents = 100;   %maximum number of events to occur
InitialFree = 10;    %initial number of proteins that are around the lattice

k_on = 10;   %kinetic rate constant for protein binding
k_off = 1;  %kinetic rate constant for protein unbinding
w = 1;  %cooperativity constant

% Memory Allocation
X_Locations = zeros(3,MaxEvents);   %number of free locations on the DNA lattice over time
X_Free = zeros(1,MaxEvents);    %number of free proteins with each event
X_Free(1) = InitialFree;    %records initial amount of free proteins
X_Bound = zeros(1,MaxEvents);   %number of bound proteins on the lattice with each event (lattice is initially empty)
Populations = zeros(3,MaxEvents);   %combined array tracking populations of bound/free proteins and free locations
a_P = zeros(4,MaxEvents);   %propensity functions
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
    clear LeftEdges RightEdges GapSize AvailablePerGap DC_AvailableLocations SC_AvailableLocations I_AvailableLocations
    AvailableLocations = [];    %initializes an array to track available locations
    LeftEdges = find(diff([1 DNA 1]) == -1);    %all left edges of any gaps
    RightEdges = find(diff([1 DNA 1]) == 1);    %all right edges of any gaps
    GapSize = RightEdges-LeftEdges; %gap size of each gap on the lattice
    LeftEdges(GapSize < n) = [];    %clears position of any gap that is too small for binding
    GapSize(GapSize < n) = [];      %clears any gap that is too small for binding
    AvailablePerGap = GapSize-(n-1);    %number of available spaces per gap
    for a = 1:numel(GapSize)
        AvailableLocations = [AvailableLocations,LeftEdges(a):1:LeftEdges(a)+AvailablePerGap(a)-1]; %all available binding locations
    end
    DC_AvailableLocations = LeftEdges(GapSize == n);    %list of all doubly contiguous binding locations
    DC_AvailableLocations(DC_AvailableLocations == 1 | DC_AvailableLocations == N-(n-1)) = [];  %clears first and last position of lattice from DC
    
    SC_AvailableLocations = sort([AvailableLocations(diff([0,diff(AvailableLocations) == 1,0])>0),AvailableLocations(diff([0,diff(AvailableLocations) == 1,0])<0)]);    %calcualtes positon of all SC locations
    if ismember(1,SC_AvailableLocations)==1 && DNA(n) == 0  %if first location isn't really SC
        SC_AvailableLocations(SC_AvailableLocations==1) = [];   %...clear it
    end
    if ismember(N-(n-1),SC_AvailableLocations)==1 && DNA(N-n) == 0  %do the same thing with the last location on the lattice
        SC_AvailableLocations(SC_AvailableLocations==N-(n-1)) = [];
    end
    
    I_AvailableLocations = AvailableLocations(~ismember(AvailableLocations,DC_AvailableLocations)); %removes DC locations from available list of isolated locations
    I_AvailableLocations = I_AvailableLocations(~ismember(I_AvailableLocations,SC_AvailableLocations)); %removes SC locations from available list
    
    X_Locations(:,Event) = [numel(I_AvailableLocations);numel(SC_AvailableLocations);numel(DC_AvailableLocations)]; %number of free locations on the lattice
        %1 - Isolated
        %2 - Singly Contiguous
        %3 - Doubly Contiguous
    X_Bound(Event) = numel(find(BoundAtSpot == 1)); %number of proteins bound to the lattice
    BoundLocations = find(BoundAtSpot == 1);    %list of locations where proteins are bound
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    a_P(:,Event) = [k_on*X_Locations(1,Event)*X_Free(Event);k_on*X_Locations(2,Event)*X_Free(Event)*w;k_on*X_Locations(3,Event)*X_Free(Event)*(w^2);k_off*X_Bound(Event)];
        % 1 - Binding Reaction (I)
        % 2 - Binding Reaction (SC)
        % 3 - Binding Reaction (DC)
        % 4 - Unbinding Reaction
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