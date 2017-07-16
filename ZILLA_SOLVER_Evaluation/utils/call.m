function [] = call( command_line )
%CALL Run a system command from the matlab prompt.
%   Call will spawn a command window and run the dos command line in that
%   window. Return values from the program are not supported, but if an
%   error is detected by the commandline interprater, this call will return
%   an error as well.
%
%   command_line should be a string that the dos commandline interprated
%   can understand. If no input argument is given, CALL will testfor the
%   existence of a command interprater. 
%
%   This mex file works on Unix and Win32.
