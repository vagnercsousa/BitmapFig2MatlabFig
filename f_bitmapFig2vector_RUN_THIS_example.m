clc;
clear;


% moment of inertia

imgfile = 'Fig_Airfoil_MomIn3.jpg';
[x1, y1] = flBitmapFig2matlabFig(imgfile, [0 9.8], [-100 94], 1);
set(gca, 'FontSize', 15);
set(gcf, 'position', [582   332   420   322]);
set(gca, 'position', [0.2120    0.2050    0.7523    0.7613]);
xlabel('$t$ [s]', 'Interpreter', 'latex', 'FontSize', 18);
ylabel('$x$ [mm]', 'Interpreter', 'latex', 'FontSize', 18);
flSetMinorTicks();
xlim([0 10]);
ylim([-100 100]);



% pitch damping

imgfile = 'Fig_Airfoil_freeVib_pit.jpg';
[x2, y2] = flBitmapFig2matlabFig(imgfile, [0 2], [-62 31], 1);
set(gca, 'FontSize', 15);
xlabel('$t$ [s]', 'Interpreter', 'latex', 'FontSize', 18);
ylabel('$\theta$ [deg]', 'Interpreter', 'latex', 'FontSize', 18);
flSetMinorTicks();
xlim([0 2]);
ylim([-70 70]);



% plunge damping

imgfile = 'Fig_Airfoil_freeVib_plu.jpg';
[x3, y3] = flBitmapFig2matlabFig(imgfile, [0 7.8], [-38 34.5], 1);
set(gca, 'FontSize', 15);
xlabel('$t$ [s]', 'Interpreter', 'latex', 'FontSize', 18);
ylabel('$h$ [mm]', 'Interpreter', 'latex', 'FontSize', 18);
flSetMinorTicks();
xlim([0 8]);
ylim([-50 50]);



figure();

subplot(2, 1, 1); % pitch damping
hold('on');
plot(x2, y2, 'LineWidth', 1.2);
set(gca, 'FontSize', 15);
set(gcf, 'position', [403   234   560   420]);
set(gca, 'position', [0.1746    0.6452    0.7750    0.3250]);
ylabel('$\alpha$ [deg]', 'Interpreter', 'latex', 'FontSize', 18');
ylabel('$\theta$ [deg]', 'Interpreter', 'latex', 'FontSize', 18);
flSetMinorTicks();
xlim([0 2]);
ylim([-70 70]);

subplot(2, 1, 2); % plunge damping
hold('on');
plot(x3, y3, 'LineWidth', 1.2);
set(gca, 'FontSize', 15);
set(gca, 'position', [0.1750    0.1714    0.7746    0.3310]);
xlabel('$t$ [s]', 'Interpreter', 'latex', 'FontSize', 18);
ylabel('$h$ [mm]', 'Interpreter', 'latex', 'FontSize', 18);
flSetMinorTicks();
xlim([0 8]);
ylim([-50 50]);


