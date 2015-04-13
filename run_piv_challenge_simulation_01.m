function run_piv_challenge_simulation_01;
% This function creates a series of PIV images that closely resemble the
% data that was produced by the PIV Challenge Case E data.

% This is the top level directory to write the simulation data into
top_write_directory='/mnt/current_storage/Projects2/Tomo_PIV/Camera_Simulation_GUI/camera_simulation_package_02/test_directory/';

% This is the directory containing the camera simulation parameters to run
% for each camera
camera_simulation_parameters_read_directory='/mnt/current_storage/Projects2/Tomo_PIV/Camera_Simulation_GUI/camera_simulation_package_02/piv_challenge_simulation_parameters/';

% This is the list of camera parameters to run the simulation over
camera_parameter_list=dir([camera_simulation_parameters_read_directory,'camera*.mat']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Creates Data Write Directories                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This is the name of the particle image top level directory
particle_image_top_directory=[top_write_directory,'camera_images/particle_images/'];
% This creates the top level particle image directory
if not(exist(particle_image_top_directory,'dir'));
    % This creates the particle image directory
    mkdir(particle_image_top_directory);
end;

% This creates the particle image data directories
for camera_index=1:length(camera_parameter_list);
    % This is the current subdirecotry name
    current_subdirectory=[particle_image_top_directory,'camera_',sprintf('%02.0f',camera_index),'/'];
    % This creates the current camera particle image directory
    if not(exist(current_subdirectory,'dir'));
        % This creates the particle image camera directory
        mkdir(current_subdirectory);
    end;
end;

% This is the name of the calibration image top level directory
calibration_image_top_directory=[top_write_directory,'camera_images/calibration_images/'];
% This creates the top level calibration image directory
if not(exist(calibration_image_top_directory,'dir'));
    % This creates the calibration image directory
    mkdir(calibration_image_top_directory);
end;

% This creates the calibration image data directories
for camera_index=1:length(camera_parameter_list);
    % This is the current subdirecotry name
    current_subdirectory=[calibration_image_top_directory,'camera_',sprintf('%02.0f',camera_index),'/'];
    % This creates the current camera calibration image directory
    if not(exist(current_subdirectory,'dir'));
        % This creates the calibration image camera directory
        mkdir(current_subdirectory);
    end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Creates Particle Position Data                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This is a uniform velocity field to simulate
X_Velocity=01.0e3;
Y_Velocity=01.0e3;
Z_Velocity=00.1e3;

% This is the domain of the particles to be simulated
X_Min=-7.5e4;
X_Max=+7.5e4;
Y_Min=-7.5e4;
Y_Max=+7.5e4;
Z_Min=-7.5e3;
Z_Max=+7.5e3;

% This is the number of particles to simulate
total_particle_number=1e7;

% This generates the particle positions
X=(X_Max-X_Min)*rand(total_particle_number,1)+X_Min;
Y=(Y_Max-Y_Min)*rand(total_particle_number,1)+Y_Min;
Z=(Z_Max-Z_Min)*rand(total_particle_number,1)+Z_Min;

% This generates the first frame particle positions
X1=X-X_Velocity/2;
Y1=Y-Y_Velocity/2;
Z1=Z-Z_Velocity/2;

% This generates the second frame particle positions
X2=X+X_Velocity/2;
Y2=Y+Y_Velocity/2;
Z2=Z+Z_Velocity/2;

% This creates the directory to save the particle data positions in
particle_position_data_directory=[top_write_directory,'particle_positions/'];
% This creates the particle position data directory
if not(exist(particle_position_data_directory,'dir'));
    % This creates the particle position directory
    mkdir(particle_position_data_directory);
end;

% This renames the first frame particle position data
X=X1;
Y=Y1;
Z=Z1;
% This writes the first frame particle position data
save([particle_position_data_directory,'particle_data_frame_0001.mat'],'X','Y','Z');

% This renames the second frame particle position data
X=X2;
Y=Y2;
Z=Z2;
% This writes the second frame particle position data
save([particle_position_data_directory,'particle_data_frame_0002.mat'],'X','Y','Z');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Performs Camera Simulation                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This iterates through the different cameras performing the image
% simulation
for camera_index=1:length(camera_parameter_list);
    
    % This displays that the current camera simulation is being ran
    fprintf('\n\n\n\n');
    disp(['Running camera ',num2str(camera_index),' simulation . . . ']);
    
    % This is the current filename of the camera parameter file to load
    parameter_filename_read=[camera_simulation_parameters_read_directory,camera_parameter_list(camera_index).name];
    % This loads the current parameter data
    load(parameter_filename_read);
    
    % This changes the directory containing the particle locations in the
    % parameters structure
    piv_simulation_parameters.particle_field.data_directory=particle_position_data_directory;
    % This changes the vector giving the frames of particle positions to
    % load in
    % the parameters structure (this indexes into the list generated by the
    % command 'dir([data_directory,data_filename_prefix,'*.mat'])')
    piv_simulation_parameters.particle_field.frame_vector=1:2;
    % This changes the number of particles to simulate out of the list of possible
    % particles (if this number is larger than the number of saved particles,
    % an error will be returned)
    piv_simulation_parameters.particle_field.particle_number=2.5e5;
    
    % This changes the directory to save the particle images in parameters
    % structure
    piv_simulation_parameters.output_data.particle_image_directory=[particle_image_top_directory,'camera_',sprintf('%02.0f',camera_index),'/'];
    % This changes the directory to save the calibration grid images in
    % parameters structure
    piv_simulation_parameters.output_data.calibration_grid_image_directory=[calibration_image_top_directory,'camera_',sprintf('%02.0f',camera_index),'/'];
    
    % This runs the camera simulation for the current camera
    run_piv_simulation_02(piv_simulation_parameters);
    
end;
    
    
    
    
    
    