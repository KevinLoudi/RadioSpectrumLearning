clc
clear all
close all

run('Scripts\add_dstemPath.m');
while (1)
    disp(' ');
    disp('****************************************************************');
    pause(0.05);
    disp('* JSS paper - D-STEM: A Software for the Analysis and Mapping  *');
    pause(0.05);
    disp('*             of Environmental Space-Time Variables            *');
    pause(0.05);
    disp('****************************************************************');
    pause(0.05);
    disp(' ');
    disp('Choose which script to run to reproduce the paper results:');
    pause(0.05);
    disp(' ');
    disp('[1] - Section 4  : Univariate model (< 3 minutes)');
    pause(0.05);
    disp('[2] - Section 5  : Bivariate  model (< 3 minutes)');
    pause(0.05);
    disp('[3] - Section 6  : Downscaler model (> 10 minutes)');
    pause(0.05);
    disp('[4] - Section 7.1: Downscaler model of Section 6 with tapering enabled (> 10 minutes)');
    pause(0.05);
    disp('[5] - Section 7.2: Downscaler model of Section 6 with distributed computing (> 10 minutes)');
    pause(0.05);
    disp('[6] - Section 7.3: Univariate model of Section 4 with log-likelihood computing disabled (> 5 minutes)');
    pause(0.05);
    disp(' ');
    disp('[7] - Exit');
    disp(' ');
    choice  =  input('Insert a number and press Enter: ','s');
    switch choice
        case '1'
            run('Scripts\demo_section4');
        case '2'
            run('Scripts\demo_section5');
        case '3'
            run('Scripts\demo_section6');
        case '4'
            run('Scripts\demo_section7_1');
        case '5'
            %delete all the MATLAB files in the Distributed folder
            files=dir('../Distributed/*.mat');
            for i=1:length(files)
                delete(['../Distributed/',files(i).name]);
            end
            %run the slave script
            if ispc
                fullpath=['"',matlabroot,'\bin\matlab.exe','"',' -r run(''Scripts\demo_runslave'')'];
            else
                fullpath=['nohup "',matlabroot,'/bin/./matlab','"',' -r run(''/Scripts/demo_runslave'') &'];
            end
            [status, res] = system(fullpath,'-echo');
            %wait for the slave to open and run
            pause(25);
            %run the master script
            run('Scripts\demo_section7_2');
        case '6'
            run('Scripts\demo_section7_3');
        case '7'
            break;
        otherwise
            disp(' ');
            disp('Please insert a number between 1 and 7');
            pause(2.5);
    end
end