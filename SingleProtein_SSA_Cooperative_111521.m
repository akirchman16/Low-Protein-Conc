clearvars;
close all;

% This code is meant to track a single protein (RAD51 only, no
% competition) that binds to and unbinds from a ssDNA lattice at low
% concentrations. Rather than having a constant concentration over time,
% individual protein populations will be tracked. This utilizes the SSA for
% stochastic simulations. Now cooperativity will be included, but still no
% competition.

N = 5000;    %length of ssDNA
n = 3;  %length of protein (RAD51)

MinEvents = 100;   %minimum number of events to occur (typically keep at 100)
InitialFree = 500;    %initial number of proteins that are around the lattice

k_on = 1;   %kinetic rate constant for protein binding
k_off = 100;  %kinetic rate constant for protein unbinding
w = 1;  %cooperativity constant

% Memory Allocation
X_Locations = zeros(3,MinEvents);   %number of free locations on the DNA lattice over time
X_Free = zeros(1,MinEvents);    %number of free proteins with each event
X_Free(1) = InitialFree;    %records initial amount of free proteins
X_Bound = zeros(1,MinEvents);   %number of bound proteins on the lattice with each event (lattice is initially empty)
Populations = zeros(3,MinEvents);   %combined array tracking populations of bound/free proteins and free locations
a_P = zeros(4,MinEvents);   %propensity functions
t = zeros(1,MinEvents+1); %time tracker
dt = zeros(1,MinEvents);    %records each time interval in the simulation
EventHistory = zeros(4,MinEvents);  %records where each event occurs at each step
    % 1 - Isolated Binding
    % 2 - Singly Contiguous Binding
    % 3 - Doubly Contiguous Binding
    % 4 - Unbinding
FracCover = zeros(1,MinEvents+1);   %tracks saturation of ssDNA lattice

DNA = zeros(1,N);   %ssDNA lattice
BoundAtSpot = zeros(1,N);   %used to record where proteins are bound

Event = 0;  %counts number of events which have occured
Equilibrium = 0;
while Equilibrium == 0 %runs until the system reaches equilibrium
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
    
    if Rands(2)*sum(a_P(:,Event)) < a_P(1,Event)    %isolated binding reactions (via SSA)
        I_BindLocation = I_AvailableLocations(randi(numel(I_AvailableLocations)));   %selects random location for binding
        DNA(I_BindLocation:I_BindLocation+(n-1)) = 1;   %binds protein
        X_Free(Event+1) = X_Free(Event)-1;  %updates free protein tracker (decrease by 1)
        BoundAtSpot(I_BindLocation) = 1;  %records where proteins are bound
        EventHistory(1,Event) = I_BindLocation;   %records where event happened
    elseif Rands(2)*sum(a_P(:,Event)) < sum(a_P(1:2,Event))    %singly contiguous binding reaction
        SC_BindLocation = SC_AvailableLocations(randi(numel(SC_AvailableLocations)));   %selects random location for binding
        DNA(SC_BindLocation:SC_BindLocation+(n-1)) = 1; %binds protein
        X_Free(Event+1) = X_Free(Event)-1;  %updates free protein tracker
        BoundAtSpot(SC_BindLocation) = 1;   %records where protein is now bound
        EventHistory(2,Event) = SC_BindLocation;    %records where the event happened
    elseif Rands(2)*sum(a_P(:,Event)) < sum(a_P(1:3,Event))    %doubly contiguous binding event
        DC_BindLocation = DC_AvailableLocations(randi(numel(DC_AvailableLocations)));   %selects random location for binding
        DNA(DC_BindLocation:DC_BindLocation+(n-1)) = 1; %binds protein
        X_Free(Event+1) = X_Free(Event)-1;  %updates free protein tracker
        BoundAtSpot(DC_BindLocation) = 1;   %records where protein is now bound
        EventHistory(3,Event) = DC_BindLocation;    %records where the event happened
    else   %unbinding reaction
        UnbindLocation = BoundLocations(randi(X_Bound(Event)));     %selects random protein to unbind from lattice
        DNA(UnbindLocation:UnbindLocation+(n-1)) = 0;   %unbinds protein
        X_Free(Event+1) = X_Free(Event)+1;  %updates free protein tracker (increase by 1)
        BoundAtSpot(UnbindLocation) = 0;    %updates BoundAtSpot
        EventHistory(4,Event) = UnbindLocation;
    end
    
    t(Event+1) = t(Event)+dt(Event);    %updates timer by advancing time according to the calculated time step
    FracCover(Event+1) = numel(find(DNA == 1))/N;   %calcualtes saturation level of ssDNA lattice
    if mod(numel(find(DNA == 1)),n) ~= 0  %if there isn't an integer number of bound proteins...
        disp('BROKEN PROTEIN');
        break  %ERROR STOP
    end
    if Event >= MinEvents
        QuarterTime = t(end)-(0.25*t(end));   %time that is 25% of the way back in the simulation
        TimeDiff = t-QuarterTime;
        TestPosition_t = find(abs(TimeDiff) == min(abs(TimeDiff)));  %position of time that is closest to 25% of the way back
        t_Equilibrium_Test = t(TestPosition_t:end);  %set of time we're looking for and testing for equilibrium
        Sat_Equilibrium_Test = FracCover(TestPosition_t:end);    %set of FracCover values to test for equilibrium in
        Avg_Saturation = sum(Sat_Equilibrium_Test)/numel(Sat_Equilibrium_Test); %average saturation level
        CurveFit = polyfit(t_Equilibrium_Test,Sat_Equilibrium_Test,1);  %linear fit to the last quarter of events
        Y_Int_Error(Event) = abs(Avg_Saturation-CurveFit(2))/Avg_Saturation;   %percent error in y-intercept fit
        if CurveFit(1) < (1/mean(dt)) & (Y_Int_Error(Event) < 0.05 | isnan(Y_Int_Error(Event)))  %if fitted slope is basically zero and the fitted y-int is basically the average saturation...
            Equilibrium = 1;    %...then the system is considered to be at equilibirum
        end
    end
end

figure(1);
subplot(2,2,[3,4]);
scatter(t,FracCover,5,'r','filled');
hold on;
yline(n*InitialFree/N,'--k',['Max. Saturation: ', num2str(round(n*InitialFree/N,2))],'LabelHorizontalAlignment','left');
yline(Avg_Saturation,'--r',['Eq. Saturation: ', num2str(round(Avg_Saturation,2))],'LabelHorizontalAlignment','right','LabelVerticalAlignment','bottom');
% xline(QuarterTime,'--r',['Eq. Time: ', num2str(round(QuarterTime,2))],'LabelVerticalAlignment','top','LabelHorizontalAlignment','left');
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
scatter(t(1:end-1),X_Locations(1,:),5,'r','filled'); hold on;    %Isolated locations
scatter(t(1:end-1),X_Locations(2,:),5,'b','filled');    %Singly Contiguous locations
scatter(t(1:end-1),X_Locations(3,:),5,'g','filled');    %Doubly Contiguous locations
scatter(t(1:end-1),sum(X_Locations(1:3,:)),5,'k','filled');   %total free locations
xlim([0 t(end)]); xlabel('Time, t');
ylim([0 N]);    ylabel('Free Locations');
legend('Isolated','Singly Contiguous','Doubly Contiguous','Total','Location','Northeast');
title('Available Locations');
box on;