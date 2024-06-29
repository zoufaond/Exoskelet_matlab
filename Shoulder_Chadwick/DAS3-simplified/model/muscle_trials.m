% das3_readosim initializes a Matlab struct from an Opensim model file
%
% osimfile: the Opensim model file from which to initialize
% musclepolyfile: a .mat file with polynomials for muscle-tendon lengths
% computed by das3_polynomials.m
% GHpolyfile: a .mat file with polynomials for GH force vectors
% computed by das3_polynomials.m
%
% If this Matlab struct is used to create the Autolev file or the MEX file,
% muscle path  and GH force vector information is not required
% so "musclepolyfile" and "GHpolyfile" are optional inputs
%
% Dimitra Blana, October 2014
% based on gait3d_readosim.m by Ton van den Bogert %
%
% Updated for OpenSim 4.0 by Derek Wolf, October 2019

clear all
% import OpenSim namespace
import org.opensim.modeling.*;

% load the model
osimfile = 'das3_simplified.osim'
Mod = Model(osimfile);

% initialize the system and get the initial state
state = Mod.initSystem();

% equilibrate all muscles
Mod.equilibrateMuscles(state);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BASICS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% basic model parameters
model.name = char(Mod.getName());
model.osimfile = char(Mod.getInputFileName());
file = dir(osimfile);
model.modified = file.date;

% get the set of muscles
MuscleSet = Mod.getMuscles();
model.nMus = MuscleSet.getSize();

% get the metabolic probe, that contains muscle masses
% Probe = Umberger2010MuscleMetabolicsProbe.safeDownCast(Mod.getProbeSet().get('metabolics'));

counter=0;

imus = 4;
    currentMuscle = MuscleSet.get(imus-1);
    
    
    counter=counter+1;
    
    % basic muscle properties
    model.muscles{imus}.name = (char(currentMuscle.getName()));
    model.muscles{imus}.osim_name = char(currentMuscle.getName());
    % model.muscles{imus}.mass = Probe.getMuscleMass(model.muscles{imus}.osim_name);
    model.muscles{imus}.fmax = currentMuscle.getMaxIsometricForce();
    model.muscles{imus}.lceopt = currentMuscle.getOptimalFiberLength();
    model.muscles{imus}.lslack = currentMuscle.getTendonSlackLength();
    model.muscles{imus}.pennopt = currentMuscle.getPennationAngleAtOptimalFiberLength();
    model.muscles{imus}.vmax = currentMuscle.getMaxContractionVelocity();
    
    % act1 = str2double(char(currentMuscle.getPropertyByName('activation1')));
    % act2 = str2double(char(currentMuscle.getPropertyByName('activation2')));
    % 
    % % From Lisa Schutte's dissertation, appendices 2 and 3.
    % % "Using musculoskeletal models to explore strategies for improving
    % % performance in electrical stimulation-induced leg cycle ergometry"
    % % LM Schutte - 1992 - Stanford University
    % % k1 = activation1
    % % k2 = activation2
    % % t_act = 1/(k1+k2)
    % % t_deact = 1/k2
    % model.muscles{imus}.tact = 1/(act1+act2);
    % model.muscles{imus}.tdeact = 1/act2;
    % 
    % model.muscles{imus}.maxact = 1;
    % 
    % currentprop = currentMuscle.getPropertyByName('active_force_length_curve');
    % active_fl_o = currentprop.getValueAsObject;
    % active_fl_spline = SimmSpline.safeDownCast(active_fl_o);
    % num_values = active_fl_spline.getSize;
    % model.muscles{imus}.Xval_active_fl = zeros(num_values,1);
    % model.muscles{imus}.Yval_active_fl = zeros(num_values,1);
    % for i=1:num_values
    %     model.muscles{imus}.Xval_active_fl(i) = active_fl_spline.getX(i-1);
    %     model.muscles{imus}.Yval_active_fl(i) = active_fl_spline.getY(i-1);
    % end
    % 
    % currentprop = currentMuscle.getPropertyByName('passive_force_length_curve');
    % passive_fl_o = currentprop.getValueAsObject;
    % passive_fl_spline = SimmSpline.safeDownCast(passive_fl_o);
    % num_values = passive_fl_spline.getSize;
    % model.muscles{imus}.Xval_passive_fl = zeros(num_values,1);
    % model.muscles{imus}.Yval_passive_fl = zeros(num_values,1);
    % for i=1:num_values
    %     model.muscles{imus}.Xval_passive_fl(i) = passive_fl_spline.getX(i-1);
    %     model.muscles{imus}.Yval_passive_fl(i) = passive_fl_spline.getY(i-1);
    % end
    % model.muscles{imus}.PEEslack = model.muscles{imus}.Xval_passive_fl(4);
    % 
    % get muscle path and DOFs this muscle crosses
    muspath = currentMuscle.getGeometryPath();
    PtSet = muspath.getPathPointSet();
    OrigCoords = muspath.getWrapSet().getSize()
    nPts = PtSet.getSize();
    coord = PtSet.get(0).getLocation(state)
    coord = PtSet.get(1).getLocation(state)

    % 
    origin_segment = PtSet.get(0).getBody();
    % origin_segment_index = find(strcmp(origin_segment,segment_names));
    insertion_segment = PtSet.get(nPts-1).getBody();
    % insertion_segment_index = find(strcmp(insertion_segment,segment_names));