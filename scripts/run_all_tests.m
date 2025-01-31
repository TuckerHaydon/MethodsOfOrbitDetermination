clc; clear all; close all;

% Get the directory of this file even if we're not in its directory.
this_file_path = matlab.desktop.editor.getActiveFilename;
[this_file_directory, ~, ~] = fileparts(this_file_path);

% The test directory is ../test
test_directory = fullfile(this_file_directory, "..", "test");

% Run all tests in the test directory
testResults = runtests(test_directory, "IncludeSubfolders", true);