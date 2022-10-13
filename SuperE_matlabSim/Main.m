%% Explanations

% This program provides displacement and efforts in element for any 3D
% tubular structure.
% The behaviour is elastic - large displacement
%function [Screening] = Main (ShipAngle,ShipSpeed)
%% Modifications by Jonathan

%%%% All comments with (%%%%) were hidden from the original code
    %in the following functions
    %Main,solve,DataJacket,automaticContact,ComputeZoneAB,impactPointTim,draw
function [Screening] = Main(ShipAngle,ShipSpeed)
%% Cleaning

close all
%clear all
clearvars -except ShipAngle ShipSpeed %samples
clc
disp(['Ship Angle = ',num2str(ShipAngle),'    Ship Speed = ',num2str(ShipSpeed)])
%% Data (geometric, mechanical, boundary conditions, forces)

drawInit = 0; % draw initial structure : 1 -> yes ; 0 -> no

%%%%%%%%%% Choose of the right file %%%%%%%%%%
DataJacket3D
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

solveOption.az = -20;
solveOption.el = 35;

[properties,dofs] = autoCalculation(data,properties,dofs,solveOption,drawInit);

%% Boundary conditions

dofs.freeDofs = freeDofs(dofs);

%% Solve

solveOption.behaviour = 1; % 0 -> elastic ; 1 -> elasto-plastic
solveOption.deltaT = 0.001; % second
solveOption.deltaTFigures = 0.3;
solveOption.epsilon = 10^(-3);
solveOption.drawAllStep = 1;
solveOption.amplification = 40;

solveOption.fontsizeTitle = 24;
solveOption.fontsizeLabel = 20;
solveOption.fontsizeLegend = 14;
solveOption.fontsizeAxis = 16;

%JONATHAN (INPUT)---------------------------------------------------------
%variables taken from DataJacket3D
data.shipTrajectory = ShipAngle;
data.shipSpeed = ShipSpeed;
data.pointTrajectory = [0 0]; %In the Jacket frame of reference
%No bulb is considered 
%JONATHAN (INPUT)---------------------------------------------------------
%%%%tic
[properties,propertiesInit,contactElement,contactNode,dofs,dispTot,matrices,plasticity,cylinder,cylinderTot,ship,impact,punchProp,punchPropLeft,punchPropRight,basisJacket,Output] = solve(data,ship,properties,dofs,solveOption);  %%Parameter M was removed from the output of the function 
Screening(1:3) = Output;
%%%%toc
end