clc
clear all
close all

addpath 'components'

%%% Params
%A = 123; % mode:shortening
%E = 184;
%tc = "tv0";
A = 65; % mode:puncturing
E = 184;
tc = "tv1";
%A = 134; % mode:repetition
%E = 267;
%tc = "tv2";
seed = 0;

polar_enc_ref(A, E, tc, seed)

