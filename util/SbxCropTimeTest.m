% Executing time crop of sbx file
% 
% 
% SbxCropTime(strPathIn, strSbxIn, strPathOut, frameStart, frameEnd)
% 
% Leander de Kraker
% 2024-5-7
% 

clear
clc

strSbxIn = 'Lambrusco_20230404_003';
strPathIn = '\\vs03\VS03-MVP-1\GluA3_VR\Data_Collection\KOinV1SST_2P\Lambrusco\20230404\';
strPathOut = 'D:\2Pdata\Leander\ShiftedLines\';
frameStart = 2;
frameEnd = round(30.9141*60*5); % 5 minutes at 30.9141 Hz

nframes = frameEnd-frameStart+1;

SbxCropTime(strPathIn, strSbxIn, strPathOut, frameStart, frameEnd)