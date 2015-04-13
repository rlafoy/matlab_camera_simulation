function [angle_data,s1_data,s2_data]=mie_scattering_data(n_medium,n_particle,r_particle,lambda,angle_number);
% This function calculates the Mie scattering about a particle within a
% medium with a real refractive index of 'n_medium', with a particle with a
% complex refractive index of 'n_particle', a radius of 'r_particle', with
% light of wavelength 'lambda', and 'angle_number' angles computed between
% 0 degrees and 180 degrees inclusive.  The code assumes that the particle
% is spherical.  The particle radius and light wavelength are in arbitrary
% units.
%
% The output argument 'angle_data' is the vector of length 'angle_number'
% giving the angles between 0 and 180 degrees, the output argument
% 's1_data' gives the magnitude of the scattered light that is
% perpendicular to the scattering plane, and the output argument 's2_data'
% gives the magnitude of the scattered light that is parallel to the
% scattering plane.
%
% The code calls the fortran function 'bhmie_table' to create a table of
% the scattering values; this code uses the fortran function 'bhmie'. These
% two functions must be compiled (in the terminal) using the commands
%   
%   f77 -c -fPIC bhmie.f
%   f77 -o bhmie_table bhmie.o bhmie_table.f
%
% before this matlab function will run correctly.  The code will also 
% correctly compile with the 'f95' and 'gfortran' compilers, but these are 
% less compatible with matlab and may produce unexpected results.
%
% Authors: Rod La Foy
% Created On: 23 August 2013
% Modified On: 19 March 2015
% Notes: This code calls 'bhmie_table' which is based up the code
% 'callbhmie.f' and 'bhmie.f' written by B.T.Draine, Princeton Univ. Obs.

% This extracts the real part of the particle refractive index
n_particle_real=real(n_particle);
% This extracts the imaginary part of the particle refractive index
n_particle_imag=imag(n_particle);

% This converts the medium's refractive index to a string
n_medium_string=num2str(n_medium,6);
% This converts the particle's real refractive index to a string
n_particle_real_string=num2str(n_particle_real,6);
% This converts the particle's imaginary refractive index to a string
n_particle_imag_string=num2str(n_particle_imag,6);
% This converts the particle's radius to a string
r_particle_string=num2str(r_particle,6);
% This converts the indcident light's wavelength to a string
lambda_string=num2str(lambda,6);
% This converts the number of angles to a string
angle_number_string=num2str(angle_number);

% This is the directory path of this m-file
code_filename=mfilename('fullpath');
% This is the index of the last slash in the filename (to remove the actual
% m-file name)
truncate_index=strfind(code_filename,filesep);
truncate_index=truncate_index(end)-1;
% This is the code_directory
code_directory=code_filename(1:truncate_index);

% This is the filename to save the data to (this uses the 'tempname' matlab
% function which creates a (likely) unique temporary file name
table_filename=[tempname,'.out'];


% This checks whether the table data already exists and deletes it if so
if exist(table_filename,'file');
    % This deletes the file containing the Mie scattering data
    delete(table_filename);
end;

% This is the string to execute to create the mie scattering table (the
% change directory command is to ensure that the code is in the current
% system search path - doing this by setting the environment variables
% wasn't working)
shell_command_string=['cd ',code_directory,'; ./bhmie_table ',n_medium_string,' ',n_particle_real_string,' ',n_particle_imag_string,' ',r_particle_string,' ',lambda_string,' ',angle_number_string,' > ',table_filename];

% This runs the system command string to calculate the Mie scattering
system(shell_command_string);

% This initializes a variable stating the number of attempts that have been
% made to load the Mie scattering data
load_attempt_number=0;

% This tries loading the Mie scattering data and if the file doesn't exist,
% this pauses and waits
while true;
    % This checks whether the file exists and if so loads it and if not,
    % the code waits a short period
    if exist(table_filename,'file');
        % This reads in the output file containing the Mie scattering data
        scattering_data=dlmread(table_filename,'',5,0);
        % This breaks the loop since the data was succesfully loaded
        break;
    else;
        % This increments the variable storing the number of load attempts
        load_attempt_number=load_attempt_number+1;
        % If the number of load attempts exceeds 10, this displays an error
        if load_attempt_number>10;
            % This displays an error stating that the maximum number of
            % load attempts has been exceeded
            error('The Mie scattering data has not be created and the maximum number of loading attempts for the data has been exceeded.');
        end;
        % This pauses for a short period
        pause(0.1);
    end;
end;

% This deletes the file containing the Mie scattering data
delete(table_filename);

% This extracts the independent variable giving the angle for the rays
angle_data=scattering_data(:,1);
% This converts the angular data from degrees to radians
angle_data=pi*angle_data/180;

% This extracts the dependent variable giving the unpolarized scattering
% magnitude
s11_data=scattering_data(:,2);
% This extracts the dependent variable giving the quantity of polarization
% (ie this will be zero for unpolarized scattering . . . I think)
pol_data=scattering_data(:,3);

% This computes the dependent variable giving the differences between the
% perpendicular and parallel polarization magnitudes
s12_data=-s11_data.*pol_data;

% This computes the dependent variable giving the scattering magnitude
% perpendicular to the scattering plane
s1_data=(s11_data-s12_data);
% This computes the dependent variable giving the scattering magnitude
% parallel to the scattering plane
s2_data=(s11_data+s12_data);





