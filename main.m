close all;
clear all;
clc;
%% parameter setting
format long
addpath('utils/');
addpath('metric_utils\');

lambdaL =1.4;
omega =5; 
patchSize = 15;
frame_num = 3;

readPath = '.\test';
savePath = '.\result';
tuneopts.lambdaL = lambdaL;
tuneopts.omega = omega;
tensor.patchSize = patchSize;
tensor.slideStep = patchSize;
frame_length = frame_num;
temporal_step = frame_num;

[All_Num,time] = target_detection(char(readPath), savePath, temporal_step, frame_length, tensor, tuneopts);

