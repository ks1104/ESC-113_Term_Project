function enzyme_gui()
    % Create the main GUI figure window
    f = figure('Name','Michaelis-Menten Simulator','Position',[200 200 800 400]);

    % --- Input section: labels and edit boxes for all user inputs ---

    uicontrol('Style','text','Position',[50 350 150 20],'String','Initial Substrate (S0):');
    s0_input = uicontrol('Style','edit','Position',[210 350 100 25],'String','10');

    uicontrol('Style','text','Position',[50 310 150 30],'String','Enzyme Concentration (E0):');
    e0_input = uicontrol('Style','edit','Position',[210 310 100 25],'String','1');

    uicontrol('Style','text','Position',[50 270 150 20],'String','Michaelis Constant (KM):');
    km_input = uicontrol('Style','edit','Position',[210 270 100 25],'String','5');

    uicontrol('Style','text','Position',[50 230 150 20],'String','Catalytic Constant (k2):');
    k2_input = uicontrol('Style','edit','Position',[210 230 100 25],'String','1');

    uicontrol('Style','text','Position',[50 190 150 20],'String','Simulation Time (t_max):');
    tmax_input = uicontrol('Style','edit','Position',[210 190 100 25],'String','50');

    % --- Button to trigger the simulation ---
    uicontrol('Style','pushbutton','String','Run Simulation',...
              'Position',[100 140 150 30],'Callback',@simulate);

    % --- Axes for plotting the output graph ---
    ax = axes('Units','pixels','Position',[360 95 420 280]);

    % --- Callback function triggered when user clicks 'Run Simulation' ---
    function simulate(~,~)
        % Retrieve and convert user input strings to numbers
        S0 = str2double(get(s0_input, 'String'));
        E0 = str2double(get(e0_input, 'String'));
        KM = str2double(get(km_input, 'String'));
        k2 = str2double(get(k2_input, 'String'));
        tmax = str2double(get(tmax_input, 'String'));

        % Compute Vmax from enzyme concentration and catalytic constant
        Vmax = k2 * E0;

        % Set time step and time vector for simulation
        h = 0.01;
        t = [0:h:tmax];

        % Initial condition: [S0; P0] (P0 = 0 initially)
        y0 = [S0; 0];

        % Solve using custom RK4 integrator
        y = rk4_solver(@(t,y) odes(t,y,Vmax,KM), y0, t);

        % --- Plot the simulation results ---
        plot(ax, t, y(1,:), 'b', 'DisplayName','[S] Substrate');  % Substrate curve
        hold(ax, 'on');
        plot(ax, t, y(2,:), 'r', 'DisplayName','[P] Product');    % Product curve
        legend(ax);                                               % Add legend
        xlabel(ax,'Time'); ylabel(ax,'Concentration');
        title(ax,'Michaelis-Menten Kinetics');
        hold(ax, 'off');
    end

    % --- Nested function defining the system of ODEs ---
    function dydt = odes(~, y, Vmax, KM)
        % y(1) = S, y(2) = P
        S = y(1);
        dSdt = -Vmax * S / (KM + S);   % Substrate consumption
        dPdt =  Vmax * S / (KM + S);   % Product formation
        dydt = [dSdt; dPdt];
    end
end
