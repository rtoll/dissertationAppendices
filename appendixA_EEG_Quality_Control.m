%%%% INSTRUCTIONS: READ ME %%%%

% 1. CLICK RUN
% 2. CLICK CHANGE FOLDER
% 3. FOLLOW AUDIO PROMPTS TO SELECT CAP SIZE, RAs, AND VHDR FILE
% 4. ALL ELSE IS AUTOMATED
% 5. PROBLEMS OR QUESTIONS TEXT RUSS AT 512.663.3142

%% 00. FIRST TIME SETUP ONLY
 
clear all; close all; clc;
% code version 5.1
% code last updated 11/10/2015
% author: russ, rtoll@stanford.edu, 512.663.3142

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% copy/paste the location of the FastQC on your computer here %%%%
full_path_to_the_FastQC_folder = ''; % Example: C:\Users\User\Desktop\FastQC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 01. Initialize and define thresholds

addpath(genpath(strcat(full_path_to_the_FastQC_folder,'\MatlabFunctions')))

% add/delete people's emails here
qc.report_will_send_to_these_emails = {'email@example1.edu', ...
'email@example2.edu', ... 
'email@example3.edu'};
% to delete someone's email, delete the entire line including the blue dots

if isempty(full_path_to_the_FastQC_folder)
   error('The path to the FastQC folder needs to be set. Type in the path in the FIRST TIME SETUP ONLY part of this program at the top.')
end
% add paths to necessary custom functions
addpath(genpath(full_path_to_the_FastQC_folder));
load('C:/MATLAB/MatlabFunctions/matlab_master_definitions/master_eeg_bv_definitions.mat')
tts('fast-quality-control-check-initiated')
pause(1);
cd(full_path_to_the_FastQC_folder)
qc.version = 5.0;
 
setpref('Internet','E_mail','matlab.russell.toll@gmail.com');
setpref('Internet','SMTP_Server','smtp.gmail.com');
setpref('Internet','SMTP_Username','matlab.russell.toll@gmail.com');
setpref('Internet','SMTP_Password','etkinlab');
props = java.lang.System.getProperties;
props.setProperty('mail.smtp.auth','true');
props.setProperty('mail.smtp.socketFactory.class', 'javax.net.ssl.SSLSocketFactory');
props.setProperty('mail.smtp.socketFactory.port','465');
 
thresh.imp.bracket(1) = 3;
thresh.imp.bracket(2) = 8;
thresh.imp.bracket(3) = 10;
thresh.imp.bracket(4) = 12;
thresh.imp.bracket(5) = 15;
thresh.imp.bracket(6) = 18;
thresh.imp.bracket(7) = 20;
 
thresh.imp.refgndbracket(1) = 1;
thresh.imp.refgndbracket(2) = 3;
thresh.imp.refgndbracket(3) = 5;
thresh.imp.refgndbracket(4) = 7;
thresh.imp.refgndbracket(5) = 8;
thresh.imp.refgndbracket(6) = 9;
thresh.imp.refgndbracket(7) = 10;
 
thresh.flatline.c = 1;
thresh.flatline.f = 2;
 
thresh.numbadchans.rec.c = 5;
thresh.numbadchans.rec.f = 10;
thresh.numbadchans.reo.c = 5;
thresh.numbadchans.reo.f = 12;
 
thresh.sleepy.c = 2.5;
thresh.sleepy.f = 4;
 
thresh.blink.mag = 4;
thresh.blink.gap = 2; % sec
 
thresh.mag.f = 5;
thresh.mag.c = 3;
 
thresh.jitter.f = 10;
thresh.jitter.c = 6;
 
thresh.parox.mag = 5;
thresh.parox.percent.c = 10;
thresh.parox.percent.f = 20;
 
thresh.hb_lo = 30;
thresh.hb_hi = 120;
 
thresh.imp.datamean.c = 15;
thresh.imp.datamean.f = 20;
 
thresh.imp.dataeach.c = 15;
thresh.imp.dataeach.f = 20;
 
thresh.imp.totbad.c = 5;
thresh.imp.totbad.f = 10;
 
thresh.imp.ref.c = 5;
thresh.imp.ref.f = 10;
 
thresh.imp.gnd.c = 5;
thresh.imp.gnd.f = 10;
 
thresh.alpha.mad.c = 1.5;
thresh.alpha.mad.f = 2.0;
 
thresh.alpha.sumpow.c = .15;
thresh.alpha.sumpow.f = .05;
 
thresh.alpha.nopeak.c = 8;
thresh.alpha.nopeak.f = 16;
 
thresh.blinkspermin.standard = 20;
thresh.blinkspermin.perc.rec.c = .15*thresh.blinkspermin.standard;
thresh.blinkspermin.perc.rec.f = .00*thresh.blinkspermin.standard;
thresh.blinkspermin.perc.reo.c = .25*thresh.blinkspermin.standard;
thresh.blinkspermin.perc.reo.f = .00*thresh.blinkspermin.standard;
 
prog.step = 0;
prog.totsteps = 59;
 
clear ans
clc; prog.step=prog.step+1; fprintf('step %g of %g successful \n%3.1f%% complete \n\n\n',prog.step,prog.totsteps,100*prog.step/prog.totsteps)
sound(audioread(strcat(full_path_to_the_FastQC_folder,'\cwb.mp3')),44000); pause(.5);

%% 02. Select the cap size and RAs 
 
tts('select cap size');
cap_size_list = {'54','56','58','60','62','OTHER SIZE'};
cap_size_index = menu('Select Cap Size',cap_size_list);
qc.cap_size = single(str2num(char(cap_size_list(cap_size_index))));
tts('select the are aze running-this session');
list_of_RAs = {'Person1','Person2','ADD A NAME'};
[RA_selection,ok] = listdlg('Name','sel_ra','PromptString','RAs runnning this session: (ctrl+click for multiple select)','SelectionMode','multiple','ListString',list_of_RAs,'ListSize',[300 300]);
addname = '';
if RA_selection(end) == size(list_of_RAs,2);
   addname = inputdlg('Type the name to add');
   RA_selection = RA_selection(1:end-1);
end
RAs = '';
for i=1:length(RA_selection)
   RAs = strcat([RAs,', ',char(list_of_RAs(RA_selection(i)))]); 
end
RAs = RAs(3:end);
if ~isempty(addname)
   RAs = strcat([RAs,', ',char(addname)])
   send_text_message('512-663-3142','Verizon','new RA',sprintf('%s',char(addname)));
end
% text Russ to alert a QC script has been initiated
send_text_message('512-663-3142','Verizon','QC Script Initiated',sprintf('%s are running it',RAs));
qc.RAs = RAs;
 
clear addname cap_size_* i list_of_RAs ok RA*
clc; prog.step=prog.step+1; fprintf('step %g of %g successful \n%3.1f%% complete \n\n\n',prog.step,prog.totsteps,100*prog.step/prog.totsteps)
sound(audioread(strcat(full_path_to_the_FastQC_folder,'\cwb.mp3')),44000); pause(.5);

%% 03. Manually select vhdr file 
 
EEG = custom_pop_fileio();
EEG.chanlocs = master.chanlocs;
tic;

clear nullclear
clc; prog.step=prog.step+1; fprintf('step %g of %g successful \n%3.1f%% complete \n\n\n',prog.step,prog.totsteps,100*prog.step/prog.totsteps)
sound(audioread(strcat(full_path_to_the_FastQC_folder,'\cwb.mp3')),44000); pause(.5);
 
%% 04. Resample and convert all numeric fields to single 
 
% resample to 250 Hz
if EEG.srate > 250 
   EEG = pop_resample(EEG,250);
end
 
% convert all double EEG fields to single
field_names = fieldnames(EEG);
 
for i=1:size(field_names,1)   
   field_contents = getfield(EEG,field_names{i,1});
   if isa(field_contents,'double')      
      EEG = setfield(EEG,field_names{i,1},single(field_contents));
   end
end
 
clear field* i
clc; prog.step=prog.step+1; fprintf('step %g of %g successful \n%3.1f%% complete \n\n\n',prog.step,prog.totsteps,100*prog.step/prog.totsteps)
sound(audioread(strcat(full_path_to_the_FastQC_folder,'\cwb.mp3')),44000); pause(.5);

%% 05. Recording timestamp 
 
fileID = fopen(strcat(EEG.comments(16:end-4),'vmrk'));
EEG.vmrkdata = textscan(fileID,'%s','Whitespace','\b\t');
fclose(fileID);
for i=1:length(EEG.vmrkdata{1})
   if strfind(EEG.vmrkdata{1}{i},'Mk1=New') == 1  
      temp = EEG.vmrkdata{1}{i};
      temp_commas = strfind(temp,',');
      temp_datetimestart = temp_commas(end) + 1;
      temp_datetime = temp(temp_datetimestart:temp_datetimestart+11);
      temp_datestr = datestr(strcat([temp_datetime(1:4),'-',temp_datetime(5:6),'-',temp_datetime(7:8),' ',temp_datetime(9:10),':',temp_datetime(11:12)]));
      EEG.datetimestamp = temp_datestr(1:end-3);
   end 
end
 
clear ans fileID i temp*
clc; prog.step=prog.step+1; fprintf('step %g of %g successful \n%3.1f%% complete \n\n\n',prog.step,prog.totsteps,100*prog.step/prog.totsteps)
sound(audioread(strcat(full_path_to_the_FastQC_folder,'\cwb.mp3')),44000); pause(.5);

%% 06. Impedances 
 
fileID = fopen(EEG.comments(16:end));
EEG.vhdrdata = textscan(fileID,'%s','Whitespace','\b\t');
fclose(fileID);
% thresh.imp = 20;
imp_start_index = [];
for i=1:length(EEG.vhdrdata{1})
   temp = strsplit(EEG.vhdrdata{1}{i});
   if strcmp(char(temp(1)),'Impedance')
      imp_start_index = i;
   end
end
imps = cell(66,2);
for i=1:64
   imps(i,1) = cellstr(num2str(i));
end
imps(65,1) = cellstr('Ref');
imps(66,1) = cellstr('Gnd');
imps(:,2) = cellstr('ERROR'); 
if ~isempty(imp_start_index)
   for j=imp_start_index+1:length(EEG.vhdrdata{1})
      temp = strsplit(EEG.vhdrdata{1}{j},{':','!'});
      for i=1:64
         if ~isempty(find(str2num(char(temp(1))) == i))
            imps(i,1) = temp(1);
            imps(i,2) = strtrim(temp(2));
         end
      end
      if strcmp(char(temp(1)),'Ref')
         imps(65,1) = temp(1);
         imps(65,2) = strtrim(temp(2));
      end
      if strcmp(char(temp(1)),'Gnd')
         imps(66,1) = temp(1);
         imps(66,2) = strtrim(temp(2));
      end
   end
end
temp = imps(1:64,2);
for c=1:64
   if ~isempty(str2num(char(temp(c))))
      imp.e64(c,1) = str2num(char(temp(c)));
   else
      imp.e64(c,1) = NaN;
   end
end
if ~isempty(str2num(char(imps(65,2))))
   imp.ref = str2num(char(imps(65,2)));
else
   imp.ref = NaN;
end
if ~isempty(str2num(char(imps(66,2))))
   imp.gnd = str2num(char(imps(66,2)));
else
   imp.gnd = NaN;
end
 
clear ans c fileID i imp_start_index imps j temp
clc; prog.step=prog.step+1; fprintf('step %g of %g successful \n%3.1f%% complete \n\n\n',prog.step,prog.totsteps,100*prog.step/prog.totsteps)
sound(audioread(strcat(full_path_to_the_FastQC_folder,'\cwb.mp3')),44000); pause(.5);
 
%% 07. Output folder 
 
pid = EEG.comments(strfind(EEG.comments,'231'):strfind(EEG.comments,'231')+6);
if isempty(pid)
pid = EEG.comments(strfind(EEG.comments,'232'):strfind(EEG.comments,'232')+6);
end
if ~exist(strcat(full_path_to_the_FastQC_folder,'\',pid),'dir'), mkdir(strcat(full_path_to_the_FastQC_folder,'\',pid)); end
cd(strcat(full_path_to_the_FastQC_folder,'\',pid))
EEG.pid = pid;
EEG.chans.bad = [];
 
clear pid
clc; prog.step=prog.step+1; fprintf('step %g of %g successful \n%3.1f%% complete \n\n\n',prog.step,prog.totsteps,100*prog.step/prog.totsteps)
sound(audioread(strcat(full_path_to_the_FastQC_folder,'\cwb.mp3')),44000); pause(.5);
 
%% 08. Epoch data 
 
for i=1:size(EEG.event,2)    
   if strcmpi(EEG.event(i).type,'REC') == 1 && strcmp(EEG.event(i+1).value,'Stimulus') == 1        
      recstart = i+1
   end
   if strcmpi(EEG.event(i).type,'REO') == 1 && strcmp(EEG.event(i+1).value,'Stimulus') == 1    
      reostart = i+1
   end
end
REC = pop_epoch(EEG,{},[0 200],'eventindices',recstart);
REO = pop_epoch(EEG,{},[0 200],'eventindices',reostart);    
 
clear EEG i re*
clc; prog.step=prog.step+1; fprintf('step %g of %g successful \n%3.1f%% complete \n\n\n',prog.step,prog.totsteps,100*prog.step/prog.totsteps)
sound(audioread(strcat(full_path_to_the_FastQC_folder,'\cwb.mp3')),44000); pause(.5);
 
%% 09. Process the REC epoch 
 
EEG = REC;
EEG.rawdata = EEG.data;
 
clear REC
clc; prog.step=prog.step+1; fprintf('step %g of %g successful \n%3.1f%% complete \n\n\n',prog.step,prog.totsteps,100*prog.step/prog.totsteps)
sound(audioread(strcat(full_path_to_the_FastQC_folder,'\cwb.mp3')),44000); pause(.5);
 
%% 10. Apply filters 
 
% High pass filter
EEG = pop_eegfiltnew(EEG,1,0);
% Notch 60 Hz
EEG = pop_eegfiltnew(EEG,58.75,61.25,[],1);
% Low pass filter
EEG = pop_eegfiltnew(EEG,0,50);
 
EEG.filtdata = EEG.data;
filtdata = EEG.data;
 
clear nullclear
clc; prog.step=prog.step+1; fprintf('step %g of %g successful \n%3.1f%% complete \n\n\n',prog.step,prog.totsteps,100*prog.step/prog.totsteps)
sound(audioread(strcat(full_path_to_the_FastQC_folder,'\cwb.mp3')),44000); pause(.5);
 
%% 11. Remove flatline channels and robustly rereference 
 
EEG.chans.bad = find(isnan(imp.e64));
 
for c=1:64
   data_mads(c,1) = mad(EEG.data(c,:),1);
   if data_mads(c,1) == 0, data_mads(c,1) = NaN; end
end
data_mads_rob_norm = abs((data_mads - nanmedian(data_mads)) ./ mad(data_mads,1));
EEG.chans.flatline = find(isnan(data_mads_rob_norm));
EEG.chans.bad = sort(unique(cat(1,EEG.chans.bad,EEG.chans.flatline)));
EEG.chans.keep = setdiff([1:64]',EEG.chans.bad);
EEG.data = filtdata;
EEG.data(EEG.chans.bad,:) = NaN;
robref = NaN(1,length(EEG.data));
for t=1:length(EEG.data)
   robref(t) = nanmedian(EEG.data(:,t));
end
for c=1:64
   EEG.data(c,:) = EEG.data(c,:) - robref;
end
 
clear c data* robref t 
clc; prog.step=prog.step+1; fprintf('step %g of %g successful \n%3.1f%% complete \n\n\n',prog.step,prog.totsteps,100*prog.step/prog.totsteps)
sound(audioread(strcat(full_path_to_the_FastQC_folder,'\cwb.mp3')),44000); pause(.5);
 
%% 12. Remove extreme magnitude channels using robust rereferencing 
 
EEG.chans.mag = [];
i=1;
reiterate = true;
while reiterate
   reiterate = false;
   for c=1:64
      data_mads(c,1) = mad(EEG.data(c,:),1);
   end
   data_mads_rob_norm = (data_mads - nanmedian(data_mads)) ./ mad(data_mads,1);
   EEG.magvals(:,i) = data_mads_rob_norm;
   i=i+1;
   mag_chans = find(abs(data_mads_rob_norm) > thresh.mag.f);
   if ~isempty(mag_chans)
      EEG.chans.mag = sort(unique(cat(1,EEG.chans.mag,mag_chans)));
      reiterate = true;
      EEG.chans.bad = sort(unique(cat(1,EEG.chans.bad,EEG.chans.mag)));
      EEG.chans.keep = setdiff([1:64]',EEG.chans.bad);
      EEG.data = filtdata;
      EEG.data(EEG.chans.bad,:) = NaN;
      robref = NaN(1,length(EEG.data));
      for t=1:length(EEG.data)
         robref(t) = nanmedian(EEG.data(:,t));
      end
      for c=1:64
         EEG.data(c,:) = EEG.data(c,:) - robref;
      end
      clear c robref t
   end
   clear c data* mag_chans
end
tempmags = EEG.magvals;
tempinds = find(isnan(tempmags(:,end)))
for i=1:length(tempinds)
   if ~isempty(find(~isnan(tempmags(tempinds(i),:)),1,'last'))
      tempmags(tempinds(i),end) = tempmags(tempinds(i),find(~isnan(tempmags(tempinds(i),:)),1,'last'));   
   else
      tempmags(tempinds(i),end) = NaN;
   end
   
end
EEG.magvals = tempmags(:,end);
 
clear i reiterate temp*
clc; prog.step=prog.step+1; fprintf('step %g of %g successful \n%3.1f%% complete \n\n\n',prog.step,prog.totsteps,100*prog.step/prog.totsteps)
sound(audioread(strcat(full_path_to_the_FastQC_folder,'\cwb.mp3')),44000); pause(.5);
 
%% 13. Remove extreme jitter channels using robust rereferencing 
 
EEG.chans.jitter = [];
i=1;
reiterate = true;
while reiterate
   reiterate = false;
   jitter = diff(EEG.data,1,2);
   for c=1:64
      jitter_mads(c,1) = mad(jitter(c,:),1);
   end
   jitter_mads_rob_norm = (jitter_mads - nanmedian(jitter_mads)) ./ mad(jitter_mads,1);
   EEG.jittvals(:,i) = jitter_mads_rob_norm;
   i=i+1;
   jitter_chans = find(abs(jitter_mads_rob_norm) > thresh.jitter.f);   
   if ~isempty(jitter_chans)
      EEG.chans.jitter = sort(unique(cat(1,EEG.chans.jitter,jitter_chans)));
      reiterate = true;
      EEG.chans.bad = sort(unique(cat(1,EEG.chans.bad,EEG.chans.jitter)));
      EEG.chans.keep = setdiff([1:64]',EEG.chans.bad);
      EEG.data = filtdata;
      EEG.data(EEG.chans.bad,:) = NaN;
      robref = NaN(1,length(EEG.data));
      for t=1:length(EEG.data)
         robref(t) = nanmedian(EEG.data(:,t));
      end
      for c=1:64
         EEG.data(c,:) = EEG.data(c,:) - robref;
      end
      clear c robref t
   end
   clear c jitter*
end
tempjitts = EEG.jittvals;
tempinds = find(isnan(tempjitts(:,end)))
for i=1:length(tempinds)
   if ~isempty(find(~isnan(tempjitts(tempinds(i),:)),1,'last'))
      tempjitts(tempinds(i),end) = tempjitts(tempinds(i),find(~isnan(tempjitts(tempinds(i),:)),1,'last'));   
   else
      tempjitts(tempinds(i),end) = NaN;
   end
   
end
EEG.jittvals = tempjitts(:,end);
 
clear i reiterate temp*
clc; prog.step=prog.step+1; fprintf('step %g of %g successful \n%3.1f%% complete \n\n\n',prog.step,prog.totsteps,100*prog.step/prog.totsteps)
sound(audioread(strcat(full_path_to_the_FastQC_folder,'\cwb.mp3')),44000); pause(.5);
 
%% 14. Restore filtered data, remove and interpolate bad channels, and rereference 
 
EEG.data = filtdata;
EEG = eeg_interp(EEG,EEG.chans.bad);
EEG = pop_reref(EEG,[]);
 
clear filtdata
clc; prog.step=prog.step+1; fprintf('step %g of %g successful \n%3.1f%% complete \n\n\n',prog.step,prog.totsteps,100*prog.step/prog.totsteps)
sound(audioread(strcat(full_path_to_the_FastQC_folder,'\cwb.mp3')),44000); pause(.5);
 
%% 15. Identify paroxysmal segments 
 
parox = false(size(EEG.data));
for c=1:64
   jittnorm(c,:) = convert_raw_to_norm(diff(EEG.data(c,:)));
   parox(c,(find(abs(jittnorm(c,:)) > thresh.parox.mag))+1) = true;
end
twindow = EEG.srate*.2;
for c=1:64
   for t=1:length(EEG.data)/twindow
      if sum(parox(c,[(t-1)*twindow+1:twindow*t])) > 0
         parox(c,[(t-1)*twindow+1:twindow*t]) = 1;
      end
   end
end
paroxq = false(1,length(parox));
paroxq(find(sum(parox) > 16)) = true;
paroxqa=false(size(paroxq));
twindow = EEG.srate*2;
for t=twindow+1:length(paroxq)-twindow
   if paroxq(t)
      paroxqa(t-twindow:t+twindow) = true;
   end
end
if sum(paroxq(1:twindow)) ~= 0
   paroxqa(1:twindow) = true;
end
if sum(paroxq(end-twindow+1:end)) ~= 0
   paroxqa(end-twindow+1:end) = true;
end
EEG.pre_parox_clean_inds = find(paroxqa == false);
 
clear c jittnorm parox* t twindow
clc; prog.step=prog.step+1; fprintf('step %g of %g successful \n%3.1f%% complete \n\n\n',prog.step,prog.totsteps,100*prog.step/prog.totsteps)
sound(audioread(strcat(full_path_to_the_FastQC_folder,'\cwb.mp3')),44000); pause(.5);
 
%% 16. Apply ICA using clean data segments 
 
[EEG.icaweights,EEG.icawinv,EEG.q] = fastQCica(EEG.data(:,EEG.pre_parox_clean_inds),20,0,0,0,0);
EEG.icasphere = eye(64);
EEG.icachansind = [1:64];
EEG.qev = (sum(abs(EEG.q),2) ./ sum(sum(abs(EEG.q)))) .* size(EEG.q,1);
for q=1:size(EEG.icawinv,2)   
   EEG.icawinv_prop(:,q) = EEG.icawinv(:,q) ./ sum(abs(EEG.icawinv(:,q)));
   EEG.icawinv_propev(:,q) =  EEG.icawinv_prop(:,q) * 64;
end
 
clear q
clc; prog.step=prog.step+1; fprintf('step %g of %g successful \n%3.1f%% complete \n\n\n',prog.step,prog.totsteps,100*prog.step/prog.totsteps)
sound(audioread(strcat(full_path_to_the_FastQC_folder,'\cwb.mp3')),44000); pause(.5);
 
%% 17. Calculate the proportional independent components power spectra 
 
[EEG.qspectra,tempfreqs] = spectopo(EEG.q,0,EEG.srate,'plot','off');
EEG.freqs.indices.initial = find(tempfreqs >= 1 & tempfreqs <= 50);
EEG.freqs.values.all = tempfreqs(EEG.freqs.indices.initial);
EEG.freqs.indices.all = 1:length(EEG.freqs.values.all);
EEG.freqs.indices.widealpha = find(EEG.freqs.values.all >= 6 & EEG.freqs.values.all <= 14);
EEG.freqs.values.widealpha = EEG.freqs.values.all(EEG.freqs.indices.widealpha);
EEG.freqs.indices.muscle = find(EEG.freqs.values.all >= 20);
EEG.freqs.values.muscle = EEG.freqs.values.all(EEG.freqs.indices.muscle);
EEG.qspectra = EEG.qspectra(:,EEG.freqs.indices.initial);
EEG.qspectra = EEG.qspectra - min(min(EEG.qspectra));
EEG.qspectra = EEG.qspectra ./ sum(sum(abs(EEG.qspectra)));
 
clear tempfreqs
clc; prog.step=prog.step+1; fprintf('step %g of %g successful \n%3.1f%% complete \n\n\n',prog.step,prog.totsteps,100*prog.step/prog.totsteps)
sound(audioread(strcat(full_path_to_the_FastQC_folder,'\cwb.mp3')),44000); pause(.5);
 
%% 18. Identify muscle artifacts 
 
art.musc = [];
z=1;
for q=1:size(EEG.q,1)
   % is muscle > alpha
   if max(EEG.qspectra(q,EEG.freqs.indices.widealpha)) < max(EEG.qspectra(q,EEG.freqs.indices.muscle))
      % is it spatially localized
      maxchan_mag = find(abs(EEG.icawinv_propev(:,q)) == max(abs(EEG.icawinv_propev(:,q))));
      adj_chans = find(master.adjmat(maxchan_mag,:) ~= 0);
      adjadj_chans = zeros(length(adj_chans),64);
      for i=1:length(adj_chans)
         adjadj_chans(i,:) = master.adjmat(adj_chans(i),:);
      end
      adjadj_chans = setdiff(unique(convert_multi2singleton(adjadj_chans)),cat(2,maxchan_mag,adj_chans,0));
      rem_chans = setdiff(1:64,cat(2,maxchan_mag,adj_chans,adjadj_chans'));
      if abs(EEG.icawinv_propev(maxchan_mag,q)) > mean(abs(EEG.icawinv_propev(adj_chans,q))) && ...
         mean(abs(EEG.icawinv_propev(adj_chans,q))) > mean(abs(EEG.icawinv_propev(adjadj_chans,q))) && ...
         mean(abs(EEG.icawinv_propev(adjadj_chans,q))) > mean(abs(EEG.icawinv_propev(rem_chans,q)))
         art.musc(z,1) = q;
         z=z+1;
      end
   end
end
 
clear adj_chans adjadj_chans maxchan_mag q rem_chans z
clc; prog.step=prog.step+1; fprintf('step %g of %g successful \n%3.1f%% complete \n\n\n',prog.step,prog.totsteps,100*prog.step/prog.totsteps)
sound(audioread(strcat(full_path_to_the_FastQC_folder,'\cwb.mp3')),44000); pause(.5);
 
%% 19. Identify pulse artifacts 
 
art.pulse = [];
rim_adj = repmat(master.chans.rim,[3,1]);
z=1;
for q=1:size(EEG.q,1)
   % max/min in outer ring and are not 1st or 2nd degree adjacent to each other
   minchan = find(EEG.icawinv_propev(:,q) == min(EEG.icawinv_propev(:,q)));
   maxchan = find(EEG.icawinv_propev(:,q) == max(EEG.icawinv_propev(:,q)));
   if ismember(minchan,master.chans.rim) && ismember(maxchan,master.chans.rim)      
      minchan_nn = find(rim_adj == minchan);
      minchan_nn = rim_adj(minchan_nn(2)-3:minchan_nn(2)+3);      
      if ~ismember(maxchan,minchan_nn)
         qtemp = (normdata(EEG.q(q,:)));
         a=find(abs(qtemp)>4);
         b=[];
         j=1;
         for i=1:length(a)
            if mean(qtemp(a)) < 0
               if a(i) ~= 1 && a(i) ~= length(qtemp) && qtemp(a(i)) < 0 && abs(qtemp(a(i)-1)) < abs(qtemp(a(i))) && abs(qtemp(a(i)+1)) < abs(qtemp(a(i)))
                  b(j,1) = a(i);
                  j=j+1;
               end
            else
               if a(i) ~= 1 && a(i) ~= length(qtemp) && qtemp(a(i)) > 0 && abs(qtemp(a(i)-1)) < abs(qtemp(a(i))) && abs(qtemp(a(i)+1)) < abs(qtemp(a(i)))
                  b(j,1) = a(i);
                  j=j+1;
               end
            end
         end         
         if length(b) >= thresh.hb_lo     
            c = diff(b);
            hb = 60*EEG.srate/mean(c(find(c > prctile(c,25) & c < prctile(c,75))));            
            if hb >= thresh.hb_lo && hb <= thresh.hb_hi               
               art.pulse(z,1) = q;
               z=z+1;               
            end
            clear a b c
         end         
      end
   end   
end
 
clear a b c hb i j maxchan minchan* q qtemp rim_adj z
clc; prog.step=prog.step+1; fprintf('step %g of %g successful \n%3.1f%% complete \n\n\n',prog.step,prog.totsteps,100*prog.step/prog.totsteps)
sound(audioread(strcat(full_path_to_the_FastQC_folder,'\cwb.mp3')),44000); pause(.5);
 
%% 20. Identify vertical eye artifacts 
 
art.vert = [];
z=1;
for q=1:size(EEG.q,1)   
   three_max_chans = sortrows([[1:64]' EEG.icawinv_propev(:,q)],2);
   three_max_chans_pos = three_max_chans(end-2:end,1);
   three_max_chans_neg = three_max_chans(1:3,1);
   if (isempty(setdiff(three_max_chans_neg,master.chans.eye.center)) || isempty(setdiff(three_max_chans_pos,master.chans.eye.center))) && ... 
       mean(abs(EEG.icawinv_propev(master.chans.eye.center,q))) > mean(abs(EEG.icawinv_propev(master.chans.eye.centerNN,q))) && ... 
       mean(abs(EEG.icawinv_propev(master.chans.eye.centerNN,q))) > mean(abs(EEG.icawinv_propev(master.chans.eye.centerRem,q)))
      art.vert(z,1) = q;
      z=z+1;
   end
end
 
clear q three* z
clc; prog.step=prog.step+1; fprintf('step %g of %g successful \n%3.1f%% complete \n\n\n',prog.step,prog.totsteps,100*prog.step/prog.totsteps)
sound(audioread(strcat(full_path_to_the_FastQC_folder,'\cwb.mp3')),44000); pause(.5);
 
%% 21. Identify horizontal eye artifacts 
 
art.horiz = [];
z=1;
for q=1:size(EEG.q,1)   
   single_max_chan = sortrows([[1:64]' EEG.icawinv_propev(:,q)],2);
   single_max_chan_pos = single_max_chan(end,1);
   single_max_chan_neg = single_max_chan(1,1);
   left_lat_mean = mean(EEG.icawinv_propev(master.chans.eye.lwNN,q));
   right_lat_mean = mean(EEG.icawinv_propev(master.chans.eye.rwNN,q));
   not_lat_mean = mean(EEG.icawinv_propev(master.chans.eye.notLat,q));
   lat_ratio = left_lat_mean/right_lat_mean;
   six_max_chans = sortrows([[1:64]' EEG.icawinv_propev(:,q)],2);
   six_max_chans_pos = six_max_chans(end-5:end,1);
   six_max_chans_neg = six_max_chans(1:6,1);
   if (((ismember(single_max_chan_neg,master.chans.eye.lwNN) && ismember(single_max_chan_pos,master.chans.eye.rwNN))) || ... 
      ((ismember(single_max_chan_pos,master.chans.eye.lwNN) && ismember(single_max_chan_neg,master.chans.eye.rwNN)))) && ... 
      abs(left_lat_mean) > abs(not_lat_mean) && abs(right_lat_mean) > abs(not_lat_mean) && ... 
      ((left_lat_mean < 0 && length(setdiff(six_max_chans_neg,master.chans.eye.lwNN)) <= 3 && length(setdiff(six_max_chans_pos,master.chans.eye.rwNN)) <= 3) || ... 
      (left_lat_mean > 0 && length(setdiff(six_max_chans_pos,master.chans.eye.lwNN)) <= 3 && length(setdiff(six_max_chans_neg,master.chans.eye.rwNN)) <= 3))      
      art.horiz(z,1) = q;
      z=z+1;
   end
end
 
clear lat_ratio left_lat_mean not_lat_mean q right_lat_mean single* six* z
clc; prog.step=prog.step+1; fprintf('step %g of %g successful \n%3.1f%% complete \n\n\n',prog.step,prog.totsteps,100*prog.step/prog.totsteps)
sound(audioread(strcat(full_path_to_the_FastQC_folder,'\cwb.mp3')),44000); pause(.5);
 
%% 22. Count blink events 
 
verteye = false(1,length(EEG.pre_parox_clean_inds));
horizeye = false(1,length(EEG.pre_parox_clean_inds));
 
if ~isempty(art.vert), verteye = abs(convert_raw_to_norm(EEG.q(art.vert(1),:))); end
if ~isempty(art.horiz), horizeye = abs(convert_raw_to_norm(EEG.q(art.horiz(1),:))); end
 
vertbin = verteye >= thresh.blink.mag;
horizbin = horizeye >= thresh.blink.mag;
EEG.totbin = vertbin | horizbin;
 
indexx = 1:length(EEG.pre_parox_clean_inds);
t=1;
while t <= length(EEG.pre_parox_clean_inds)-1   
   if EEG.totbin(t) && ~EEG.totbin(t+1)      
      nextbin = find(EEG.totbin & indexx > t,1,'first');      
      if nextbin - t < thresh.blink.gap*EEG.srate         
         EEG.totbin(t:nextbin) = true;
         t=nextbin-1;         
      end      
   end   
   t=t+1;   
end
EEG.blinkcount = 0;
for t=1:length(EEG.totbin)-1   
   if EEG.totbin(t) && ~EEG.totbin(t+1)      
      EEG.blinkcount=EEG.blinkcount+1;      
   end   
end
 
clear horizbin horizeye indexx nextbin t vertbin verteye
clc; prog.step=prog.step+1; fprintf('step %g of %g successful \n%3.1f%% complete \n\n\n',prog.step,prog.totsteps,100*prog.step/prog.totsteps)
sound(audioread(strcat(full_path_to_the_FastQC_folder,'\cwb.mp3')),44000); pause(.5);
 
%% 23. Subtract artifactual components, validate projection, calculate alpha peak power, and save epoch 
 
diffin = sum(abs(convert_multi2singleton((diff(EEG.data(:,EEG.pre_parox_clean_inds),1,2)))));
sumin = sum(abs(convert_multi2singleton(((EEG.data(:,EEG.pre_parox_clean_inds))))));
art.all = sort(unique(cat(1,art.musc,art.pulse,art.vert,art.horiz)));
if length(art.all) > 10, art.kill = art.all(1:10); else art.kill = art.all; end
EEG.newdata = EEG.icawinv(:,setdiff(1:size(EEG.icaweights,1),art.kill)) * ...
             (EEG.icaweights(setdiff(1:size(EEG.icaweights,1),art.kill),:) * ...
              EEG.icasphere)*EEG.data(:,EEG.pre_parox_clean_inds);
diffout = sum(abs(convert_multi2singleton((diff(EEG.newdata(:,:),1,2)))));
sumout = sum(abs(convert_multi2singleton(((EEG.newdata(:,:))))));
if sumout < sumin && diffout < diffin, EEG.valid_proj = true; else EEG.valid_proj = false; end
[EEG.cspectra,~] = spectopo(EEG.newdata,0,EEG.srate,'plot','off');
EEG.cspectra = EEG.cspectra(:,EEG.freqs.indices.initial);
EEG.cspectra = EEG.cspectra - min(min(EEG.cspectra))
EEG.cspectra = EEG.cspectra ./ sum(sum(abs(EEG.cspectra)));
EEG.alpha.freqpow = NaN(64,2);
for c=1:64
   i=1;
   for f=2:length(EEG.freqs.indices.widealpha)-1      
      if EEG.cspectra(c,EEG.freqs.indices.widealpha(f)) > EEG.cspectra(c,EEG.freqs.indices.widealpha(f-1)) && ...
         EEG.cspectra(c,EEG.freqs.indices.widealpha(f)) > EEG.cspectra(c,EEG.freqs.indices.widealpha(f+1))
         tempfp(i,1) = EEG.freqs.values.all(EEG.freqs.indices.widealpha(f));
         tempfp(i,2) = EEG.cspectra(c,EEG.freqs.indices.widealpha(f));
         i=i+1;
      end
   end
   if exist('tempfp','var')
      EEG.alpha.freqpow(c,1) = tempfp(find(tempfp(:,2) == max(tempfp(:,2))),1);
      EEG.alpha.freqpow(c,2) = tempfp(find(tempfp(:,2) == max(tempfp(:,2))),2);
      clear tempfp
   end
end
EEG.alpha.sumpow = sum(sum(EEG.cspectra(:,EEG.freqs.indices.widealpha)));
EEG.alpha.mad = mad(EEG.alpha.freqpow(:,1));
EEG.alpha.nopeak = find(isnan(EEG.alpha.freqpow(:,1)));
EEG.art = art;
REC = EEG;
 
clear art c diff* EEG f i propspec sum*
clc; prog.step=prog.step+1; fprintf('step %g of %g successful \n%3.1f%% complete \n\n\n',prog.step,prog.totsteps,100*prog.step/prog.totsteps)
sound(audioread(strcat(full_path_to_the_FastQC_folder,'\cwb.mp3')),44000); pause(.5);
 
%% 24. Process the REO epoch 
 
EEG = REO;
EEG.rawdata = EEG.data;
 
clear REO
clc; prog.step=prog.step+1; fprintf('step %g of %g successful \n%3.1f%% complete \n\n\n',prog.step,prog.totsteps,100*prog.step/prog.totsteps)
sound(audioread(strcat(full_path_to_the_FastQC_folder,'\cwb.mp3')),44000); pause(.5);
 
%% 25. Apply filters 
 
% High pass filter
EEG = pop_eegfiltnew(EEG,1,0);
% Notch 60 Hz
EEG = pop_eegfiltnew(EEG,58.75,61.25,[],1);
% Low pass filter
EEG = pop_eegfiltnew(EEG,0,50);
 
EEG.filtdata = EEG.data;
filtdata = EEG.data;
 
clear nullclear
clc; prog.step=prog.step+1; fprintf('step %g of %g successful \n%3.1f%% complete \n\n\n',prog.step,prog.totsteps,100*prog.step/prog.totsteps)
sound(audioread(strcat(full_path_to_the_FastQC_folder,'\cwb.mp3')),44000); pause(.5);
 
%% 26. Remove flatline channels and robustly rereference 
 
EEG.chans.bad = find(isnan(imp.e64));
 
for c=1:64
   data_mads(c,1) = mad(EEG.data(c,:),1);
   if data_mads(c,1) == 0, data_mads(c,1) = NaN; end
end
data_mads_rob_norm = abs((data_mads - nanmedian(data_mads)) ./ mad(data_mads,1));
EEG.chans.flatline = find(isnan(data_mads_rob_norm));
EEG.chans.bad = sort(unique(cat(1,EEG.chans.bad,EEG.chans.flatline)));
EEG.chans.keep = setdiff([1:64]',EEG.chans.bad);
EEG.data = filtdata;
EEG.data(EEG.chans.bad,:) = NaN;
robref = NaN(1,length(EEG.data));
for t=1:length(EEG.data)
   robref(t) = nanmedian(EEG.data(:,t));
end
for c=1:64
   EEG.data(c,:) = EEG.data(c,:) - robref;
end
 
clear c data* robref t 
clc; prog.step=prog.step+1; fprintf('step %g of %g successful \n%3.1f%% complete \n\n\n',prog.step,prog.totsteps,100*prog.step/prog.totsteps)
sound(audioread(strcat(full_path_to_the_FastQC_folder,'\cwb.mp3')),44000); pause(.5);
 
%% 27. Remove extreme magnitude channels using robust rereferencing 
 
EEG.chans.mag = [];
i=1;
reiterate = true;
while reiterate
   reiterate = false;
   for c=1:64
      data_mads(c,1) = mad(EEG.data(c,:),1);
   end
   data_mads_rob_norm = (data_mads - nanmedian(data_mads)) ./ mad(data_mads,1);
   EEG.magvals(:,i) = data_mads_rob_norm;
   i=i+1;
   mag_chans = find(abs(data_mads_rob_norm) > thresh.mag.f);
   if ~isempty(mag_chans)
      EEG.chans.mag = sort(unique(cat(1,EEG.chans.mag,mag_chans)));
      reiterate = true;
      EEG.chans.bad = sort(unique(cat(1,EEG.chans.bad,EEG.chans.mag)));
      EEG.chans.keep = setdiff([1:64]',EEG.chans.bad);
      EEG.data = filtdata;
      EEG.data(EEG.chans.bad,:) = NaN;
      robref = NaN(1,length(EEG.data));
      for t=1:length(EEG.data)
         robref(t) = nanmedian(EEG.data(:,t));
      end
      for c=1:64
         EEG.data(c,:) = EEG.data(c,:) - robref;
      end
      clear c robref t
   end
   clear c data* mag_chans
end
tempmags = EEG.magvals;
tempinds = find(isnan(tempmags(:,end)))
for i=1:length(tempinds)
   if ~isempty(find(~isnan(tempmags(tempinds(i),:)),1,'last'))
      tempmags(tempinds(i),end) = tempmags(tempinds(i),find(~isnan(tempmags(tempinds(i),:)),1,'last'));   
   else
      tempmags(tempinds(i),end) = NaN;
   end
   
end
EEG.magvals = tempmags(:,end);
 
clear i reiterate temp*
clc; prog.step=prog.step+1; fprintf('step %g of %g successful \n%3.1f%% complete \n\n\n',prog.step,prog.totsteps,100*prog.step/prog.totsteps)
sound(audioread(strcat(full_path_to_the_FastQC_folder,'\cwb.mp3')),44000); pause(.5);
 
%% 28. Remove extreme jitter channels using robust rereferencing 
 
EEG.chans.jitter = [];
i=1;
reiterate = true;
while reiterate
   reiterate = false;
   jitter = diff(EEG.data,1,2);
   for c=1:64
      jitter_mads(c,1) = mad(jitter(c,:),1);
   end
   jitter_mads_rob_norm = (jitter_mads - nanmedian(jitter_mads)) ./ mad(jitter_mads,1);
   EEG.jittvals(:,i) = jitter_mads_rob_norm;
   i=i+1;
   jitter_chans = find(abs(jitter_mads_rob_norm) > thresh.jitter.f);   
   if ~isempty(jitter_chans)
      EEG.chans.jitter = sort(unique(cat(1,EEG.chans.jitter,jitter_chans)));
      reiterate = true;
      EEG.chans.bad = sort(unique(cat(1,EEG.chans.bad,EEG.chans.jitter)));
      EEG.chans.keep = setdiff([1:64]',EEG.chans.bad);
      EEG.data = filtdata;
      EEG.data(EEG.chans.bad,:) = NaN;
      robref = NaN(1,length(EEG.data));
      for t=1:length(EEG.data)
         robref(t) = nanmedian(EEG.data(:,t));
      end
      for c=1:64
         EEG.data(c,:) = EEG.data(c,:) - robref;
      end
      clear c robref t
   end
   clear c jitter*
end
tempjitts = EEG.jittvals;
tempinds = find(isnan(tempjitts(:,end)))
for i=1:length(tempinds)
   if ~isempty(find(~isnan(tempjitts(tempinds(i),:)),1,'last'))
      tempjitts(tempinds(i),end) = tempjitts(tempinds(i),find(~isnan(tempjitts(tempinds(i),:)),1,'last'));   
   else
      tempjitts(tempinds(i),end) = NaN;
   end
   
end
EEG.jittvals = tempjitts(:,end);
 
clear i reiterate temp*
clc; prog.step=prog.step+1; fprintf('step %g of %g successful \n%3.1f%% complete \n\n\n',prog.step,prog.totsteps,100*prog.step/prog.totsteps)
sound(audioread(strcat(full_path_to_the_FastQC_folder,'\cwb.mp3')),44000); pause(.5);
 
%% 29. Restore filtered data, remove and interpolate bad channels, and rereference 
 
EEG.data = filtdata;
EEG = eeg_interp(EEG,EEG.chans.bad);
EEG = pop_reref(EEG,[]);
 
clear filtdata
clc; prog.step=prog.step+1; fprintf('step %g of %g successful \n%3.1f%% complete \n\n\n',prog.step,prog.totsteps,100*prog.step/prog.totsteps)
sound(audioread(strcat(full_path_to_the_FastQC_folder,'\cwb.mp3')),44000); pause(.5);
 
%% 30. Identify paroxysmal segments 
 
parox = false(size(EEG.data));
for c=1:64
   jittnorm(c,:) = convert_raw_to_norm(diff(EEG.data(c,:)));
   parox(c,(find(abs(jittnorm(c,:)) > thresh.parox.mag))+1) = true;
end
twindow = EEG.srate*.2;
for c=1:64
   for t=1:length(EEG.data)/twindow
      if sum(parox(c,[(t-1)*twindow+1:twindow*t])) > 0
         parox(c,[(t-1)*twindow+1:twindow*t]) = 1;
      end
   end
end
paroxq = false(1,length(parox));
paroxq(find(sum(parox) > 16)) = true;
paroxqa=false(size(paroxq));
twindow = EEG.srate*2;
for t=twindow+1:length(paroxq)-twindow
   if paroxq(t)
      paroxqa(t-twindow:t+twindow) = true;
   end
end
if sum(paroxq(1:twindow)) ~= 0
   paroxqa(1:twindow) = true;
end
if sum(paroxq(end-twindow+1:end)) ~= 0
   paroxqa(end-twindow+1:end) = true;
end
EEG.pre_parox_clean_inds = find(paroxqa == false);
 
clear c jittnorm parox* t twindow
clc; prog.step=prog.step+1; fprintf('step %g of %g successful \n%3.1f%% complete \n\n\n',prog.step,prog.totsteps,100*prog.step/prog.totsteps)
sound(audioread(strcat(full_path_to_the_FastQC_folder,'\cwb.mp3')),44000); pause(.5);
 
%% 31. Apply ICA using clean data segments 
 
[EEG.icaweights,EEG.icawinv,EEG.q] = fastQCica(EEG.data(:,EEG.pre_parox_clean_inds),20,0,0,0,0);
EEG.icasphere = eye(64);
EEG.icachansind = [1:64];
EEG.qev = (sum(abs(EEG.q),2) ./ sum(sum(abs(EEG.q)))) .* size(EEG.q,1);
for q=1:size(EEG.icawinv,2)   
   EEG.icawinv_prop(:,q) = EEG.icawinv(:,q) ./ sum(abs(EEG.icawinv(:,q)));
   EEG.icawinv_propev(:,q) =  EEG.icawinv_prop(:,q) * 64;
end
 
clear q
clc; prog.step=prog.step+1; fprintf('step %g of %g successful \n%3.1f%% complete \n\n\n',prog.step,prog.totsteps,100*prog.step/prog.totsteps)
sound(audioread(strcat(full_path_to_the_FastQC_folder,'\cwb.mp3')),44000); pause(.5);
 
%% 32. Calculate the proportional independent components power spectra 
 
[EEG.qspectra,tempfreqs] = spectopo(EEG.q,0,EEG.srate,'plot','off');
EEG.freqs.indices.initial = find(tempfreqs >= 1 & tempfreqs <= 50);
EEG.freqs.values.all = tempfreqs(EEG.freqs.indices.initial);
EEG.freqs.indices.all = 1:length(EEG.freqs.values.all);
EEG.freqs.indices.widealpha = find(EEG.freqs.values.all >= 6 & EEG.freqs.values.all <= 14);
EEG.freqs.values.widealpha = EEG.freqs.values.all(EEG.freqs.indices.widealpha);
EEG.freqs.indices.muscle = find(EEG.freqs.values.all >= 20);
EEG.freqs.values.muscle = EEG.freqs.values.all(EEG.freqs.indices.muscle);
EEG.qspectra = EEG.qspectra(:,EEG.freqs.indices.initial);
EEG.qspectra = EEG.qspectra - min(min(EEG.qspectra));
EEG.qspectra = EEG.qspectra ./ sum(sum(abs(EEG.qspectra)));
 
clear tempfreqs
clc; prog.step=prog.step+1; fprintf('step %g of %g successful \n%3.1f%% complete \n\n\n',prog.step,prog.totsteps,100*prog.step/prog.totsteps)
sound(audioread(strcat(full_path_to_the_FastQC_folder,'\cwb.mp3')),44000); pause(.5);
 
%% 33. Identify muscle artifacts 
 
art.musc = [];
z=1;
for q=1:size(EEG.q,1)
   % is muscle > alpha
   if max(EEG.qspectra(q,EEG.freqs.indices.widealpha)) < max(EEG.qspectra(q,EEG.freqs.indices.muscle))
      % is it spatially localized
      maxchan_mag = find(abs(EEG.icawinv_propev(:,q)) == max(abs(EEG.icawinv_propev(:,q))));
      adj_chans = find(master.adjmat(maxchan_mag,:) ~= 0);
      adjadj_chans = zeros(length(adj_chans),64);
      for i=1:length(adj_chans)
         adjadj_chans(i,:) = master.adjmat(adj_chans(i),:);
      end
      adjadj_chans = setdiff(unique(convert_multi2singleton(adjadj_chans)),cat(2,maxchan_mag,adj_chans,0));
      rem_chans = setdiff(1:64,cat(2,maxchan_mag,adj_chans,adjadj_chans'));
      if abs(EEG.icawinv_propev(maxchan_mag,q)) > mean(abs(EEG.icawinv_propev(adj_chans,q))) && ...
         mean(abs(EEG.icawinv_propev(adj_chans,q))) > mean(abs(EEG.icawinv_propev(adjadj_chans,q))) && ...
         mean(abs(EEG.icawinv_propev(adjadj_chans,q))) > mean(abs(EEG.icawinv_propev(rem_chans,q)))
         art.musc(z,1) = q;
         z=z+1;
      end
   end
end
 
clear adj_chans adjadj_chans i maxchan_mag q rem_chans z
clc; prog.step=prog.step+1; fprintf('step %g of %g successful \n%3.1f%% complete \n\n\n',prog.step,prog.totsteps,100*prog.step/prog.totsteps)
sound(audioread(strcat(full_path_to_the_FastQC_folder,'\cwb.mp3')),44000); pause(.5);
 
%% 34. Identify pulse artifacts 
 
art.pulse = [];
rim_adj = repmat(master.chans.rim,[3,1]);
z=1;
for q=1:size(EEG.q,1)
   % max/min in outer ring and are not 1st or 2nd degree adjacent to each other
   minchan = find(EEG.icawinv_propev(:,q) == min(EEG.icawinv_propev(:,q)));
   maxchan = find(EEG.icawinv_propev(:,q) == max(EEG.icawinv_propev(:,q)));
   if ismember(minchan,master.chans.rim) && ismember(maxchan,master.chans.rim)      
      minchan_nn = find(rim_adj == minchan);
      minchan_nn = rim_adj(minchan_nn(2)-3:minchan_nn(2)+3);      
      if ~ismember(maxchan,minchan_nn)
         qtemp = (normdata(EEG.q(q,:)));
         a=find(abs(qtemp)>4);
         b=[];
         j=1;
         for i=1:length(a)
            if mean(qtemp(a)) < 0
               if a(i) ~= 1 && a(i) ~= length(qtemp) && qtemp(a(i)) < 0 && abs(qtemp(a(i)-1)) < abs(qtemp(a(i))) && abs(qtemp(a(i)+1)) < abs(qtemp(a(i)))
                  b(j,1) = a(i);
                  j=j+1;
               end
            else
               if a(i) ~= 1 && a(i) ~= length(qtemp) && qtemp(a(i)) > 0 && abs(qtemp(a(i)-1)) < abs(qtemp(a(i))) && abs(qtemp(a(i)+1)) < abs(qtemp(a(i)))
                  b(j,1) = a(i);
                  j=j+1;
               end
            end
         end         
         if length(b) >= thresh.hb_lo     
            c = diff(b);
            hb = 60*EEG.srate/mean(c(find(c > prctile(c,25) & c < prctile(c,75))));            
            if hb >= thresh.hb_lo && hb <= thresh.hb_hi               
               art.pulse(z,1) = q;
               z=z+1;               
            end
            clear a b c
         end         
      end
   end   
end
 
clear a b c hb i j maxchan minchan* q qtemp rim_adj z
clc; prog.step=prog.step+1; fprintf('step %g of %g successful \n%3.1f%% complete \n\n\n',prog.step,prog.totsteps,100*prog.step/prog.totsteps)
sound(audioread(strcat(full_path_to_the_FastQC_folder,'\cwb.mp3')),44000); pause(.5);
 
%% 35. Identify vertical eye artifacts 
 
art.vert = [];
z=1;
for q=1:size(EEG.q,1)   
   three_max_chans = sortrows([[1:64]' EEG.icawinv_propev(:,q)],2);
   three_max_chans_pos = three_max_chans(end-2:end,1);
   three_max_chans_neg = three_max_chans(1:3,1);
   if (length(setdiff(three_max_chans_neg,master.chans.eye.center)) <= 1 || length(setdiff(three_max_chans_pos,master.chans.eye.center)) <= 1) && ... 
       mean(abs(EEG.icawinv_propev(master.chans.eye.center,q))) > mean(abs(EEG.icawinv_propev(master.chans.eye.centerNN,q))) && ... 
       mean(abs(EEG.icawinv_propev(master.chans.eye.centerNN,q))) > mean(abs(EEG.icawinv_propev(master.chans.eye.centerRem,q)))
      art.vert(z,1) = q;
      z=z+1;
   end
end
 
clear q three* z
clc; prog.step=prog.step+1; fprintf('step %g of %g successful \n%3.1f%% complete \n\n\n',prog.step,prog.totsteps,100*prog.step/prog.totsteps)
sound(audioread(strcat(full_path_to_the_FastQC_folder,'\cwb.mp3')),44000); pause(.5);
 
%% 36. Identify horizontal eye artifacts 
 
art.horiz = [];
z=1;
for q=1:size(EEG.q,1)   
   single_max_chan = sortrows([[1:64]' EEG.icawinv_propev(:,q)],2);
   single_max_chan_pos = single_max_chan(end,1);
   single_max_chan_neg = single_max_chan(1,1);
   left_lat_mean = mean(EEG.icawinv_propev(master.chans.eye.lwNN,q));
   right_lat_mean = mean(EEG.icawinv_propev(master.chans.eye.rwNN,q));
   not_lat_mean = mean(EEG.icawinv_propev(master.chans.eye.notLat,q));
   lat_ratio = left_lat_mean/right_lat_mean;
   six_max_chans = sortrows([[1:64]' EEG.icawinv_propev(:,q)],2);
   six_max_chans_pos = six_max_chans(end-5:end,1);
   six_max_chans_neg = six_max_chans(1:6,1);
   if (((ismember(single_max_chan_neg,master.chans.eye.lwNN) && ismember(single_max_chan_pos,master.chans.eye.rwNN))) || ... 
      ((ismember(single_max_chan_pos,master.chans.eye.lwNN) && ismember(single_max_chan_neg,master.chans.eye.rwNN)))) && ... 
      abs(left_lat_mean) > abs(not_lat_mean) && abs(right_lat_mean) > abs(not_lat_mean) && ... 
      ((left_lat_mean < 0 && length(setdiff(six_max_chans_neg,master.chans.eye.lwNN)) <= 3 && length(setdiff(six_max_chans_pos,master.chans.eye.rwNN)) <= 3) || ... 
      (left_lat_mean > 0 && length(setdiff(six_max_chans_pos,master.chans.eye.lwNN)) <= 3 && length(setdiff(six_max_chans_neg,master.chans.eye.rwNN)) <= 3))      
      art.horiz(z,1) = q;
      z=z+1;
   end
end
 
clear lat_ratio left_lat_mean not_lat_mean q right_lat_mean single* six* z
clc; prog.step=prog.step+1; fprintf('step %g of %g successful \n%3.1f%% complete \n\n\n',prog.step,prog.totsteps,100*prog.step/prog.totsteps)
sound(audioread(strcat(full_path_to_the_FastQC_folder,'\cwb.mp3')),44000); pause(.5);
 
%% 37. Count blink events 
 
verteye = false(1,length(EEG.pre_parox_clean_inds));
horizeye = false(1,length(EEG.pre_parox_clean_inds));
 
if ~isempty(art.vert), verteye = abs(convert_raw_to_norm(EEG.q(art.vert(1),:))); end
if ~isempty(art.horiz), horizeye = abs(convert_raw_to_norm(EEG.q(art.horiz(1),:))); end
 
vertbin = verteye >= thresh.blink.mag;
horizbin = horizeye >= thresh.blink.mag;
EEG.totbin = vertbin | horizbin;
 
indexx = 1:length(EEG.pre_parox_clean_inds);
t=1;
while t <= length(EEG.pre_parox_clean_inds)-1   
   if EEG.totbin(t) && ~EEG.totbin(t+1)      
      nextbin = find(EEG.totbin & indexx > t,1,'first');      
      if nextbin - t < thresh.blink.gap*EEG.srate         
         EEG.totbin(t:nextbin) = true;
         t=nextbin-1;         
      end      
   end   
   t=t+1;   
end
EEG.blinkcount = 0;
for t=1:length(EEG.totbin)-1   
   if EEG.totbin(t) && ~EEG.totbin(t+1)      
      EEG.blinkcount=EEG.blinkcount+1;      
   end   
end
 
clear horizbin horizeye indexx nextbin t vertbin verteye
clc; prog.step=prog.step+1; fprintf('step %g of %g successful \n%3.1f%% complete \n\n\n',prog.step,prog.totsteps,100*prog.step/prog.totsteps)
sound(audioread(strcat(full_path_to_the_FastQC_folder,'\cwb.mp3')),44000); pause(.5);
 
%% 38. Subtract artifactual components, validate projection, and save epoch 
 
diffin = sum(abs(convert_multi2singleton((diff(EEG.data(:,EEG.pre_parox_clean_inds),1,2)))));
sumin = sum(abs(convert_multi2singleton(((EEG.data(:,EEG.pre_parox_clean_inds))))));
art.all = sort(unique(cat(1,art.musc,art.pulse,art.vert,art.horiz)));
if length(art.all) > 10, art.kill = art.all(1:10); else art.kill = art.all; end
EEG.newdata = EEG.icawinv(:,setdiff(1:size(EEG.icaweights,1),art.kill)) * ...
             (EEG.icaweights(setdiff(1:size(EEG.icaweights,1),art.kill),:) * ...
              EEG.icasphere)*EEG.data(:,EEG.pre_parox_clean_inds);
diffout = sum(abs(convert_multi2singleton((diff(EEG.newdata(:,:),1,2)))));
sumout = sum(abs(convert_multi2singleton(((EEG.newdata(:,:))))));
if sumout < sumin && diffout < diffin, EEG.valid_proj = true; else EEG.valid_proj = false; end
[EEG.cspectra,~] = spectopo(EEG.newdata,0,EEG.srate,'plot','off');
EEG.cspectra = EEG.cspectra(:,EEG.freqs.indices.initial);
EEG.cspectra = EEG.cspectra - min(min(EEG.cspectra))
EEG.cspectra = EEG.cspectra ./ sum(sum(abs(EEG.cspectra)));
EEG.alpha.freqpow = NaN(64,2);
for c=1:64
   i=1;
   for f=2:length(EEG.freqs.indices.widealpha)-1      
      if EEG.cspectra(c,EEG.freqs.indices.widealpha(f)) > EEG.cspectra(c,EEG.freqs.indices.widealpha(f-1)) && ...
         EEG.cspectra(c,EEG.freqs.indices.widealpha(f)) > EEG.cspectra(c,EEG.freqs.indices.widealpha(f+1))
         tempfp(i,1) = EEG.freqs.values.all(EEG.freqs.indices.widealpha(f));
         tempfp(i,2) = EEG.cspectra(c,EEG.freqs.indices.widealpha(f));
         i=i+1;
      end
   end
   if exist('tempfp','var')
      EEG.alpha.freqpow(c,1) = tempfp(find(tempfp(:,2) == max(tempfp(:,2))),1);
      EEG.alpha.freqpow(c,2) = tempfp(find(tempfp(:,2) == max(tempfp(:,2))),2);
      clear tempfp
   end
end
EEG.alpha.sumpow = sum(sum(EEG.cspectra(:,EEG.freqs.indices.widealpha)));
EEG.alpha.mad = mad(EEG.alpha.freqpow(:,1));
EEG.alpha.nopeak = find(isnan(EEG.alpha.freqpow(:,1)));
EEG.art = art;
REO = EEG;
 
clear art c diff* EEG f i propspec sum*
clc; prog.step=prog.step+1; fprintf('step %g of %g successful \n%3.1f%% complete \n\n\n',prog.step,prog.totsteps,100*prog.step/prog.totsteps)
sound(audioread(strcat(full_path_to_the_FastQC_folder,'\cwb.mp3')),44000); pause(.5);
 
%% 39. Calculate QC scores 
 
fail_count = 0;
caution_count = 0;
rec_sleepy = 0;
reo_sleepy = 0;
% flatliners
qc.flatliners = length(unique(cat(1,REC.chans.flatline,REO.chans.flatline)));
if qc.flatliners >= thresh.flatline.c && num_flatliners < thresh.flatline.f
   qc.resultstr.c.flatliners = 'number of flatline channels';
   qc.resultstr.f.flatliners = '';
   caution_count=caution_count+1;
elseif qc.flatliners >= thresh.flatline.f
   qc.resultstr.c.flatliners = '';
   qc.resultstr.f.flatliners = 'too many flatline channels';
   fail_count=fail_count+1;
else
   qc.resultstr.c.flatliners = '';
   qc.resultstr.f.flatliners = '';
end
% mean chan impedance
qc.impmean = nanmean(imp.e64);
if qc.impmean >= thresh.imp.datamean.c && qc.impmean < thresh.imp.datamean.f
   qc.resultstr.c.impdatamean = 'average electrode impedance is high';
   qc.resultstr.f.impdatamean = '';
   caution_count=caution_count+1;
elseif qc.impmean > thresh.imp.datamean.f
   qc.resultstr.c.impdatamean = '';
   qc.resultstr.f.impdatamean = 'average electrode impedance is too high';
   fail_count=fail_count+1;
else
   qc.resultstr.c.impdatamean = '';
   qc.resultstr.f.impdatamean = '';
end
% total bad impedance chans
qc.badimps = sum(isnan(imp.e64)) + length(find(imp.e64 > thresh.imp.dataeach.f)) + .5*length(find(imp.e64 > thresh.imp.dataeach.c));
if qc.badimps >= thresh.imp.totbad.c && qc.badimps < thresh.imp.totbad.f
   qc.resultstr.c.imptotbad = 'several electrodes have high impedance';
   qc.resultstr.f.imptotbad = '';
   caution_count=caution_count+1;
elseif qc.badimps > thresh.imp.totbad.f
   qc.resultstr.c.imptotbad = '';
   qc.resultstr.f.imptotbad = 'too many electrodes have high impedance';
   fail_count=fail_count+1;
else
   qc.resultstr.c.imptotbad = '';
   qc.resultstr.f.imptotbad = '';
end
% ref impedance
qc.impref = imp.ref;
if qc.impref >= thresh.imp.ref.c && qc.impref < thresh.imp.ref.f
   qc.resultstr.c.impref = 'reference impedance is high';
   qc.resultstr.f.impref = '';
   caution_count=caution_count+1;
elseif qc.impref > thresh.imp.ref.f
   qc.resultstr.c.impref = '';
   qc.resultstr.f.impref = 'reference impedance is too high';
   fail_count=fail_count+1;
else
   qc.resultstr.c.impref = '';
   qc.resultstr.f.impref = '';
end
% gnd impedance
qc.impgnd = imp.gnd;
if qc.impgnd >= thresh.imp.gnd.c && qc.impgnd < thresh.imp.gnd.f
   qc.resultstr.c.impgnd = 'ground impedance is high';
   qc.resultstr.f.impgnd = '';
   caution_count=caution_count+1;
elseif qc.impgnd > thresh.imp.gnd.f
   qc.resultstr.c.impgnd = '';
   qc.resultstr.f.impgnd = 'ground impedance is too high';
   fail_count=fail_count+1;
else
   qc.resultstr.c.impgnd = '';
   qc.resultstr.f.impgnd = '';
end
% REC bad chans
qc.recbadchans = length(REC.chans.bad);
if qc.recbadchans >= thresh.numbadchans.rec.c && qc.recbadchans < thresh.numbadchans.rec.f
   qc.resultstr.c.recbadchans = 'number of bad channels during REC is high';
   qc.resultstr.f.recbadchans = '';
   caution_count=caution_count+1;
elseif qc.recbadchans > thresh.numbadchans.rec.f
   qc.resultstr.c.recbadchans = '';
   qc.resultstr.f.recbadchans = 'number of bad channels during REC is too high';
   fail_count=fail_count+1;
else
   qc.resultstr.c.recbadchans = '';
   qc.resultstr.f.recbadchans = '';
end
% REC parox data
qc.recparox = length(REC.newdata)/length(REC.data);
if qc.recparox >= thresh.parox.percent.c && qc.recparox < thresh.parox.percent.f
   qc.resultstr.c.recparox = 'percentage of paroxysmal/contaminated data during REC is high (participant is moving too much)';
   qc.resultstr.f.recparox = '';
   caution_count=caution_count+1;
elseif qc.recparox > thresh.parox.percent.f
   qc.resultstr.c.recparox = '';
   qc.resultstr.f.recparox = 'percentage of paroxysmal/contaminated data during REC is too high (participant is moving too much)';
   fail_count=fail_count+1;
else
   qc.resultstr.c.recparox = '';
   qc.resultstr.f.recparox = '';
end
% REC alpha peak mad
qc.recalphamad = REC.alpha.mad;
if qc.recalphamad >= thresh.alpha.mad.c && qc.recalphamad < thresh.alpha.mad.f
   qc.resultstr.c.recalphamad = 'frequencies where alpha peaks occur during REC are spread out (low coherence), suggests drowsiness';
   qc.resultstr.f.recalphamad = '';
   rec_sleepy = rec_sleepy + .5;
   caution_count=caution_count+1;
elseif qc.recalphamad > thresh.alpha.mad.f
   qc.resultstr.c.recalphamad = '';
   qc.resultstr.f.recalphamad = 'frequencies where alpha peaks occur during REC are too spread out (very low coherence), suggests drowsiness';
   rec_sleepy = rec_sleepy + 1;
   fail_count=fail_count+1;
else
   qc.resultstr.c.recalphamad = '';
   qc.resultstr.f.recalphamad = '';
end
% REC alpha power
qc.recalphasumpow = REC.alpha.sumpow;
if qc.recalphasumpow >= thresh.alpha.sumpow.c && qc.recalphasumpow < thresh.alpha.sumpow.f
   qc.resultstr.c.recalphasumpow = 'percentage of total power in alpha band during REC is low, suggests drowsiness';
   qc.resultstr.f.recalphasumpow = '';
   rec_sleepy = rec_sleepy + .5;
   caution_count=caution_count+1;
elseif qc.recalphasumpow < thresh.alpha.sumpow.f
   qc.resultstr.c.recalphasumpow = '';
   qc.resultstr.f.recalphasumpow = 'percentage of total power in alpha band during REC is too low, suggests drowsiness';
   rec_sleepy = rec_sleepy + 1;
   fail_count=fail_count+1;
else
   qc.resultstr.c.recalphasumpow = '';
   qc.resultstr.f.recalphasumpow = '';
end
% REC alpha no peaks
qc.recalphanopeak = length(REC.alpha.nopeak);
if qc.recalphanopeak >= thresh.alpha.nopeak.c && qc.recalphanopeak < thresh.alpha.nopeak.f
   qc.resultstr.c.recalphanopeak = 'number of peaks in alpha band during REC is low, suggests drowsiness';
   qc.resultstr.f.recalphanopeak = '';
   rec_sleepy = rec_sleepy + .5;
   caution_count=caution_count+1;
elseif qc.recalphanopeak > thresh.alpha.nopeak.f
   qc.resultstr.c.recalphanopeak = '';
   qc.resultstr.f.recalphanopeak = 'number of peaks in alpha band during REC is too low, suggests drowsiness';
   rec_sleepy = rec_sleepy + 1;
   fail_count=fail_count+1;
else
   qc.resultstr.c.recalphanopeak = '';
   qc.resultstr.f.recalphanopeak = '';
end
% REC blinks
qc.recblinks = REC.blinkcount/(length(REC.newdata)/REC.srate/60); % blinks per minute
if qc.recblinks >= thresh.blinkspermin.perc.rec.c && qc.recblinks < thresh.blinkspermin.perc.rec.f
   qc.resultstr.c.recblinks = 'number of eye movements during REC is low, suggests drowsiness';
   qc.resultstr.f.recblinks = '';
   rec_sleepy = rec_sleepy + .5;
   caution_count=caution_count+1;
elseif qc.recblinks < thresh.blinkspermin.perc.rec.f
   qc.resultstr.c.recblinks = '';
   qc.resultstr.f.recblinks = 'number of eye movements during REC is too low, suggests drowsiness';
   rec_sleepy = rec_sleepy + 1;
   fail_count=fail_count+1;
else
   qc.resultstr.c.recblinks = '';
   qc.resultstr.f.recblinks = '';
end
% REC sleepy tally
if (~isempty(REC.art.pulse) && ~isempty(REC.art.vert) && (REC.art.pulse(1) < REC.art.vert(1))) || ... 
   (~isempty(REC.art.pulse) && ~isempty(REC.art.horiz) && (REC.art.pulse(1) < REC.art.horiz(1)))
   rec_sleepy = rec_sleepy +.5;
end
if isempty(REC.art.vert) || isempty(REC.art.horiz)
   rec_sleepy = rec_sleepy + 1;
end
qc.recsleepy = rec_sleepy;
if qc.recsleepy >= thresh.sleepy.c && qc.recsleepy < thresh.sleepy.f
   qc.resultstr.c.recsleepy = 'participant appears drowsy during REC';
   qc.resultstr.f.recsleepy = '';
   caution_count=caution_count+1;
elseif qc.recsleepy > thresh.sleepy.f
   qc.resultstr.c.recsleepy = '';
   qc.resultstr.f.recsleepy = 'participant appears asleep during REC';
   fail_count=fail_count+1;
else
   qc.resultstr.c.recsleepy = '';
   qc.resultstr.f.recsleepy = '';
end
% REO bad chans
qc.reobadchans = length(REO.chans.bad);
if qc.reobadchans >= thresh.numbadchans.reo.c && qc.reobadchans < thresh.numbadchans.reo.f
   qc.resultstr.c.reobadchans = 'number of bad channels during REO is high';
   qc.resultstr.f.reobadchans = '';
   caution_count=caution_count+1;
elseif qc.reobadchans > thresh.numbadchans.reo.f
   qc.resultstr.c.reobadchans = '';
   qc.resultstr.f.reobadchans = 'number of bad channels during REO is too high';
   fail_count=fail_count+1;
else
   qc.resultstr.c.reobadchans = '';
   qc.resultstr.f.reobadchans = '';
end
% REO parox data
qc.reoparox = length(REO.newdata)/length(REO.data);
if qc.reoparox >= thresh.parox.percent.c && qc.reoparox < thresh.parox.percent.f
   qc.resultstr.c.reoparox = 'percentage of paroxysmal/contaminated data during REO is high (participant is moving too much)';
   qc.resultstr.f.reoparox = '';
   caution_count=caution_count+1;
elseif qc.reoparox > thresh.parox.percent.f
   qc.resultstr.c.reoparox = '';
   qc.resultstr.f.reoparox = 'percentage of paroxysmal/contaminated data during REO is too high (participant is moving too much)';
   fail_count=fail_count+1;
else
   qc.resultstr.c.reoparox = '';
   qc.resultstr.f.reoparox = '';
end
% REO alpha peak mad
qc.reoalphamad = REO.alpha.mad;
if qc.reoalphamad >= thresh.alpha.mad.c && qc.reoalphamad < thresh.alpha.mad.f
   qc.resultstr.c.reoalphamad = 'frequencies where alpha peaks occur during REO are spread out (low coherence), suggests drowsiness';
   qc.resultstr.f.reoalphamad = '';
   reo_sleepy = reo_sleepy + .5;
   caution_count=caution_count+1;
elseif qc.reoalphamad > thresh.alpha.mad.f
   qc.resultstr.c.reoalphamad = '';
   qc.resultstr.f.reoalphamad = 'frequencies where alpha peaks occur during REO are too spread out (very low coherence), suggests drowsiness';
   reo_sleepy = reo_sleepy + 1;
   fail_count=fail_count+1;
else
   qc.resultstr.c.reoalphamad = '';
   qc.resultstr.f.reoalphamad = '';
end
% REO alpha power
qc.reoalphasumpow = REO.alpha.sumpow;
if qc.reoalphasumpow >= thresh.alpha.sumpow.c && qc.reoalphasumpow < thresh.alpha.sumpow.f
   qc.resultstr.c.reoalphasumpow = 'percentage of total power in alpha band during REO is low, suggests drowsiness';
   qc.resultstr.f.reoalphasumpow = '';
   reo_sleepy = reo_sleepy + .5;
   caution_count=caution_count+1;
elseif qc.reoalphasumpow < thresh.alpha.sumpow.f
   qc.resultstr.c.reoalphasumpow = '';
   qc.resultstr.f.reoalphasumpow = 'percentage of total power in alpha band during REO is too low, suggests drowsiness';
   reo_sleepy = reo_sleepy + 1;
   fail_count=fail_count+1;
else
   qc.resultstr.c.reoalphasumpow = '';
   qc.resultstr.f.reoalphasumpow = '';
end
% REO alpha no peaks
qc.reoalphanopeak = length(REO.alpha.nopeak);
if qc.reoalphanopeak >= thresh.alpha.nopeak.c && qc.reoalphanopeak < thresh.alpha.nopeak.f
   qc.resultstr.c.reoalphanopeak = 'number of peaks in alpha band during REO is low, suggests drowsiness';
   qc.resultstr.f.reoalphanopeak = '';
   reo_sleepy = reo_sleepy + .5;
   caution_count=caution_count+1;
elseif qc.reoalphanopeak > thresh.alpha.nopeak.f
   qc.resultstr.c.reoalphanopeak = '';
   qc.resultstr.f.reoalphanopeak = 'number of peaks in alpha band during REO is too low, suggests drowsiness';
   reo_sleepy = reo_sleepy + 1;
   fail_count=fail_count+1;
else
   qc.resultstr.c.reoalphanopeak = '';
   qc.resultstr.f.reoalphanopeak = '';
end
% REO blinks
qc.reoblinks = REO.blinkcount/(length(REO.newdata)/REO.srate/60); % blinks per minute
if qc.reoblinks >= thresh.blinkspermin.perc.reo.c && qc.reoblinks < thresh.blinkspermin.perc.reo.f
   qc.resultstr.c.reoblinks = 'number of eye movements during REO is low, suggests drowsiness';
   qc.resultstr.f.reoblinks = '';
   reo_sleepy = reo_sleepy + .5;
   caution_count=caution_count+1;
elseif qc.reoblinks < thresh.blinkspermin.perc.reo.f
   qc.resultstr.c.reoblinks = '';
   qc.resultstr.f.reoblinks = 'number of eye movements during REO is too low, suggests drowsiness';
   reo_sleepy = reo_sleepy + 1;
   fail_count=fail_count+1;
else
   qc.resultstr.c.reoblinks = '';
   qc.resultstr.f.reoblinks = '';
end
% REO sleepy tally
if (~isempty(REO.art.pulse) && ~isempty(REO.art.vert) && (REO.art.pulse(1) < REO.art.vert(1))) || ... 
   (~isempty(REO.art.pulse) && ~isempty(REO.art.horiz) && (REO.art.pulse(1) < REO.art.horiz(1)))
   reo_sleepy = reo_sleepy +.5;
end
if isempty(REO.art.vert) || isempty(REO.art.horiz)
   reo_sleepy = reo_sleepy + 1;
end
qc.reosleepy = reo_sleepy;
if qc.reosleepy >= thresh.sleepy.c && qc.reosleepy < thresh.sleepy.f
   qc.resultstr.c.reosleepy = 'participant appears drowsy during REO';
   qc.resultstr.f.reosleepy = '';
   caution_count=caution_count+1;
elseif qc.reosleepy > thresh.sleepy.f
   qc.resultstr.c.reosleepy = '';
   qc.resultstr.f.reosleepy = 'participant appears asleep during REO';
   fail_count=fail_count+1;
else
   qc.resultstr.c.reosleepy = '';
   qc.resultstr.f.reosleepy = '';
end
 
clear nullclear
clc; prog.step=prog.step+1; fprintf('step %g of %g successful \n%3.1f%% complete \n\n\n',prog.step,prog.totsteps,100*prog.step/prog.totsteps)
sound(audioread(strcat(full_path_to_the_FastQC_folder,'\cwb.mp3')),44000); pause(.5);
 
%% 40. Tally cautions and failures 

% fail aggregate
fail_trim = {''};
fail_cell = struct2cell(qc.resultstr.f);
j=1;
for i=1:length(fail_cell)
   if ~isempty(char(fail_cell(i)))
      fail_trim(j,1) = fail_cell(i,1);
      j=j+1;
   end  
end
% caution aggregate
caution_trim = {''};
caution_cell = struct2cell(qc.resultstr.c);
j=1;
for i=1:length(caution_cell)
   if ~isempty(char(caution_cell(i)))
      caution_trim(j,1) = caution_cell(i,1);
      j=j+1;
   end  
end
badimps = cat(1,find(isnan(imp.e64)),find(imp.e64 > thresh.imp.dataeach.f));
badimpsstr = '';
for i=1:length(badimps)
   badimpsstr = strcat(badimpsstr,',',char(master.chanlocs(badimps(i)).labels)); 
end
badimpsstr = badimpsstr(2:end);
% REC bad chans
REC.chans.badstr = '';
for i=1:length(REC.chans.bad)
   REC.chans.badstr = strcat(REC.chans.badstr,',',char(master.chanlocs(REC.chans.bad(i)).labels)); 
end
REC.chans.badstr = REC.chans.badstr(2:end);
% REO bad chans
REO.chans.badstr = '';
for i=1:length(REO.chans.bad)
   REO.chans.badstr = strcat(REO.chans.badstr,',',char(master.chanlocs(REO.chans.bad(i)).labels)); 
end
REO.chans.badstr = REO.chans.badstr(2:end);
% QC verdict
if fail_count > 0
   qc.verdict = 'FAIL';
elseif caution_count > 0
   qc.verdict = 'PASS WITH CAUTIONS';
else
   qc.verdict = 'PASS';
end
 
clear caution_c* fail_c* i j
clc; prog.step=prog.step+1; fprintf('step %g of %g successful \n%3.1f%% complete \n\n\n',prog.step,prog.totsteps,100*prog.step/prog.totsteps)
sound(audioread(strcat(full_path_to_the_FastQC_folder,'\cwb.mp3')),44000); pause(.5);
 
%% 41. Generate summary info 

linecounter=1;
line1 = sprintf('Figures will be sent in a separate email shortly. \n\nResult: %s \n\n',qc.verdict);
linecounter=linecounter+1;
if ~isempty(char(fail_trim))
   line2 = sprintf('FAILS: \n');
else
   line2 = '';
end
linecounter=linecounter+1;
if length(char(fail_trim)) ~= 0
   for i=1:length(fail_trim)
      assignin('base',(strcat('line',num2str(linecounter))),sprintf('%g. %s \n',i,char(fail_trim(i))));
      if i == length(fail_trim)
         assignin('base',(strcat('line',num2str(linecounter))),sprintf('%g. %s \n\n',i,char(fail_trim(i))));
      end
      linecounter=linecounter+1;
   end
end
if ~isempty(char(caution_trim))
   assignin('base',(strcat('line',num2str(linecounter))),sprintf('CAUTIONS: \n'));
else
   assignin('base',(strcat('line',num2str(linecounter))),'');
end
linecounter=linecounter+1;
if length(char(caution_trim))
   for i=1:length(caution_trim)
      assignin('base',(strcat('line',num2str(linecounter))),sprintf('%g. %s \n',i,char(caution_trim(i))));
      if i == length(caution_trim)
         assignin('base',(strcat('line',num2str(linecounter))),sprintf('%g. %s \n\n',i,char(caution_trim(i))));
      end
      linecounter=linecounter+1;
   end
end
assignin('base',(strcat('line',num2str(linecounter))),sprintf('Ref Imp: %g kOhm \n',imp.ref));
linecounter=linecounter+1;
assignin('base',(strcat('line',num2str(linecounter))),sprintf('Gnd Imp: %g kOhm \n',imp.gnd));
linecounter=linecounter+1;
assignin('base',(strcat('line',num2str(linecounter))),sprintf('Mean Imp: %2.1f kOhm \n',nanmean(imp.e64)));
linecounter=linecounter+1;
if length(badimps) == 0
   assignin('base',(strcat('line',num2str(linecounter))),sprintf('Chans Imp Too High: %s \n','none'));
else
   assignin('base',(strcat('line',num2str(linecounter))),sprintf('Chans Imp Too High: %g (%s) \n',length(badimps),badimpsstr));
end
linecounter=linecounter+1;
assignin('base',(strcat('line',num2str(linecounter))),sprintf('\nResting Eyes Closed: \n'));
linecounter=linecounter+1;
if length(REC.chans.bad) == 0
   assignin('base',(strcat('line',num2str(linecounter))),sprintf('\tBad Chans: %s \n','none'));
else
   assignin('base',(strcat('line',num2str(linecounter))),sprintf('\tBad Chans: %g (%s) \n',length(REC.chans.bad),REC.chans.badstr));
end
linecounter=linecounter+1;
assignin('base',(strcat('line',num2str(linecounter))),sprintf('\tParoxysmal Data: %2.1f%% \n',round(100*(1-length(REC.newdata)/length(REC.data)),1)));
linecounter=linecounter+1;
assignin('base',(strcat('line',num2str(linecounter))),sprintf('\tSleepiness: %g (0 = awake, 5 = asleep) \n',rec_sleepy));
linecounter=linecounter+1;
assignin('base',(strcat('line',num2str(linecounter))),sprintf('\nResting Eyes Open: \n'));
linecounter=linecounter+1;
if length(REO.chans.bad) == 0
   assignin('base',(strcat('line',num2str(linecounter))),sprintf('\tBad Chans: %s \n','none'));
else
   assignin('base',(strcat('line',num2str(linecounter))),sprintf('\tBad Chans: %g (%s) \n',length(REO.chans.bad),REO.chans.badstr));
end
linecounter=linecounter+1;
assignin('base',(strcat('line',num2str(linecounter))),sprintf('\tParoxysmal Data: %2.1f%% \n',round(100*(1-length(REO.newdata)/length(REO.data)),1)));
linecounter=linecounter+1;
assignin('base',(strcat('line',num2str(linecounter))),sprintf('\tSleepiness: %g (0 = awake, 5 = asleep) \n',reo_sleepy));
linecounter=linecounter+1;
assignin('base',(strcat('line',num2str(linecounter))),sprintf('\nCap Size: %g \n',qc.cap_size));
linecounter=linecounter+1;
assignin('base',(strcat('line',num2str(linecounter))),sprintf('The RAs were: %s \n\n',qc.RAs));
linecounter=linecounter+1;
 
clear bad* caution* fail* i rec* reo*
clc; prog.step=prog.step+1; fprintf('step %g of %g successful \n%3.1f%% complete \n\n\n',prog.step,prog.totsteps,100*prog.step/prog.totsteps)
sound(audioread(strcat(full_path_to_the_FastQC_folder,'\cwb.mp3')),44000); pause(.5);
 
%% 42. Generate impedance table 
 
flagC = 'c';
flagF = 'F';
assignin('base',(strcat('line',num2str(linecounter))),sprintf('IMPEDANCES TABLE: \nELEC   IMP   QC\n'));
linecounter=linecounter+1;
for c=1:64   
   if imp.e64(c) >= thresh.imp.dataeach.c && imp.e64(c) < thresh.imp.dataeach.f      
      flagL = flagC;      
   elseif imp.e64(c) >= thresh.imp.dataeach.f      
      flagL = flagF;      
   else       
      flagL = '';      
   end      
   assignin('base',(strcat('line',num2str(linecounter))),sprintf('E%02g   %3g    %c\n',c,imp.e64(c),flagL));
   linecounter=linecounter+1;
end
if imp.ref >= thresh.imp.ref.c && imp.ref < thresh.imp.ref.f      
      flagL = flagC;      
   elseif imp.ref >= thresh.imp.ref.f      
      flagL = flagF;      
   else       
      flagL = '';      
end      
assignin('base',(strcat('line',num2str(linecounter))),sprintf('%s   %3g    %c\n','REF',imp.ref,flagL));   
linecounter=linecounter+1;
if imp.gnd >= thresh.imp.gnd.c && imp.gnd < thresh.imp.gnd.f      
      flagL = flagC;      
   elseif imp.gnd >= thresh.imp.gnd.f      
      flagL = flagF;      
   else       
      flagL = '';      
end      
assignin('base',(strcat('line',num2str(linecounter))),sprintf('%s   %3g    %c\n','GND',imp.gnd,flagL));   
linecounter=linecounter+1;
 
clear c flagL
clc; prog.step=prog.step+1; fprintf('step %g of %g successful \n%3.1f%% complete \n\n\n',prog.step,prog.totsteps,100*prog.step/prog.totsteps)
sound(audioread(strcat(full_path_to_the_FastQC_folder,'\cwb.mp3')),44000); pause(.5);
 
%% 43. Generate flatline summary 
 
tempcatREC = '';
for i=1:length(REC.chans.flatline)
   templine = sprintf('E%02g',REC.chans.flatline(i));
   tempcatREC = strcat(tempcatREC,',',templine);
   if i == length(REC.chans.flatline), tempcatREC = tempcatREC(2:end); end
end
tempcatREO = '';
for i=1:length(REO.chans.flatline)
   templine = sprintf('E%02g',REO.chans.flatline(i));
   tempcatREO = strcat(tempcatREO,',',templine);
   if i == length(REO.chans.flatline), tempcatREO = tempcatREO(2:end); end
end
if isempty(tempcatREC), tempcatREC = 'none'; end
if isempty(tempcatREO), tempcatREO = 'none'; end
assignin('base',(strcat('line',num2str(linecounter))),sprintf('\nFlatline Channels \nREC: %s \nREO: %s\n\n',tempcatREC,tempcatREO));
linecounter=linecounter+1;
 
clear i temp*
clc; prog.step=prog.step+1; fprintf('step %g of %g successful \n%3.1f%% complete \n\n\n',prog.step,prog.totsteps,100*prog.step/prog.totsteps)
sound(audioread(strcat(full_path_to_the_FastQC_folder,'\cwb.mp3')),44000); pause(.5);
 
%% 44. Generate noisy channel table 
 
assignin('base',(strcat('line',num2str(linecounter))),sprintf('NOISINESS  TABLE: \n          MAG Z SCORE          JITTER Z SCORE \nELEC   REC  QC    REO  QC    REC  QC    REO  QC\n'));
linecounter=linecounter+1;
 
for c=1:64   
   if ismember(c,REC.chans.mag)
      flagRECmag = flagF;
   elseif REC.magvals(c) > thresh.mag.c
      flagRECmag = flagC;
   else
      flagRECmag = ' ';
   end
   if ismember(c,REO.chans.mag)
      flagREOmag = flagF;
   elseif REO.magvals(c) > thresh.mag.c
      flagREOmag = flagC;
   else
      flagREOmag = ' ';
   end
   if ismember(c,REC.chans.jitter)
      flagRECjitt = flagF;
   elseif REC.jittvals(c) > thresh.jitter.c
      flagRECjitt = flagC;
   else
      flagRECjitt = ' ';
   end
   if ismember(c,REO.chans.jitter)
      flagREOjitt = flagF;
   elseif REO.jittvals(c) > thresh.jitter.c
      flagREOjitt = flagC;
   else
      flagREOjitt = ' ';
   end
   assignin('base',(strcat('line',num2str(linecounter))),sprintf('E%02g   %4.1f   %c   %4.1f   %c   %4.1f   %c   %4.1f   %c\n',c,REC.magvals(c),flagRECmag,REO.magvals(c),flagREOmag,REC.jittvals(c),flagRECjitt,REO.jittvals(c),flagREOjitt));   
   linecounter=linecounter+1;
end
 
clear c flag*
clc; prog.step=prog.step+1; fprintf('step %g of %g successful \n%3.1f%% complete \n\n\n',prog.step,prog.totsteps,100*prog.step/prog.totsteps)
sound(audioread(strcat(full_path_to_the_FastQC_folder,'\cwb.mp3')),44000); pause(.5);
 
%% 45. Assemble and email the plain text summary report 
 
linescat = '';
for i=1:linecounter-1
   linescat = strcat(linescat,strcat('line',num2str(i)),',');   
end
linescat = strcat('[',linescat);
linescat(end) = ']';                         
ptmsg = eval(eval('linescat'));
tts('now emailing the summary report');
sendmail(qc.report_will_send_to_these_emails,sprintf('Cohen RS EEG Fast QC: %s',REC.pid),ptmsg);
stoptime = toc;
send_text_message('512-663-3142','Verizon','QC Script PT Runtime',sprintf('%02g:%02g',floor(stoptime/60),floor(rem(stoptime,60))));

clear i line* ptmsg stoptime
clc; prog.step=prog.step+1; fprintf('step %g of %g successful \n%3.1f%% complete \n\n\n',prog.step,prog.totsteps,100*prog.step/prog.totsteps)
sound(audioread(strcat(full_path_to_the_FastQC_folder,'\cwb.mp3')),44000); pause(.5);

%% 46. Generate impedances figure 
 
screensize = get(0,'ScreenSize');
imp.color(1,:) = [000 100 000]./255;
imp.color(2,:) = [000 128 000]./255;
imp.color(3,:) = [173 255 047]./255;
imp.color(4,:) = [255 255 000]./255;
imp.color(5,:) = [255 215 000]./255;
imp.color(6,:) = [255 165 000]./255;
imp.color(7,:) = [255 069 000]./255;
imp.color(8,:) = [255 000 000]./255;
figure;
set(gcf,'visible','off','Position',[1 1 min(screensize(3:end)) min(screensize(3:end))],'Color','w')
hold on
axis('square')
for c=1:64
   % determine fill color
   impval = num2str(imp.e64(c));
   numcolor = [192 192 192]./255;
   if isnan(imp.e64(c))
      facecolor = 'k';
      impval = '!!!';
      numcolor = 'w';
   elseif imp.e64(c) < thresh.imp.bracket(1)
      facecolor = imp.color(1,:);
   elseif imp.e64(c) < thresh.imp.bracket(2)
      facecolor = imp.color(2,:);
   elseif imp.e64(c) < thresh.imp.bracket(3)
      facecolor = imp.color(3,:);
   elseif imp.e64(c) < thresh.imp.bracket(4)
      facecolor = imp.color(4,:);
   elseif imp.e64(c) < thresh.imp.bracket(5)
      facecolor = imp.color(5,:);
   elseif imp.e64(c) < thresh.imp.bracket(6)
      facecolor = imp.color(6,:);
   elseif imp.e64(c) < thresh.imp.bracket(7)
      facecolor = imp.color(7,:);
   else
      facecolor = imp.color(8,:);
   end
   subplot('position',master.capplot(c,:))
   rectangle('Curvature',[1 1],'Position',[.1 .1 .9 .9],'FaceColor',facecolor,'EdgeColor','k','LineWidth',2.5)
   text(.55,1,sprintf('E%g',c),'Color',char(master.colors(c)),'Units','normalized','FontSize',12,'FontWeight','b','HorizontalAlignment','center','VerticalAlignment','bottom')
   text(.55,.55,impval,'Color',numcolor,'Units','normalized','FontSize',24,'FontWeight','b','HorizontalAlignment','center','VerticalAlignment','middle')
   axis('square')
   set(gca,'Visible','off')
end
 
% reference
% determine fill color
impval = num2str(imp.ref);
numcolor = [192 192 192]./255;
if isnan(imp.ref)
   facecolor = 'k';
   impval = '!!!';
   numcolor = 'w';
elseif imp.ref < thresh.imp.refgndbracket(1)
   facecolor = imp.color(1,:);
elseif imp.ref < thresh.imp.refgndbracket(2)
   facecolor = imp.color(2,:);
elseif imp.ref < thresh.imp.refgndbracket(3)
   facecolor = imp.color(3,:);
elseif imp.ref < thresh.imp.refgndbracket(4)
   facecolor = imp.color(4,:);
elseif imp.ref < thresh.imp.refgndbracket(5)
   facecolor = imp.color(5,:);
elseif imp.ref < thresh.imp.refgndbracket(6)
   facecolor = imp.color(6,:);
elseif imp.ref < thresh.imp.refgndbracket(7)
   facecolor = imp.color(7,:);
else
   facecolor = imp.color(8,:);
end
subplot('position',master.capplot(66,:))
rectangle('Curvature',[1 1],'Position',[.1 .1 .9 .9],'FaceColor',facecolor,'EdgeColor','c','LineWidth',2.5)
text(.55,1,'REF','Color','k','Units','normalized','FontSize',12,'FontWeight','b','HorizontalAlignment','center','VerticalAlignment','bottom')
text(.55,.55,impval,'Color',numcolor,'Units','normalized','FontSize',24,'FontWeight','b','HorizontalAlignment','center','VerticalAlignment','middle')
axis('square')
set(gca,'Visible','off')
 
% ground
% determine fill color
impval = num2str(imp.gnd);
numcolor = [192 192 192]./255;
if isnan(imp.gnd)
   facecolor = 'k';
   impval = '!!!';
   numcolor = 'w';
elseif imp.gnd < thresh.imp.refgndbracket(1)
   facecolor = imp.color(1,:);
elseif imp.gnd < thresh.imp.refgndbracket(2)
   facecolor = imp.color(2,:);
elseif imp.gnd < thresh.imp.refgndbracket(3)
   facecolor = imp.color(3,:);
elseif imp.gnd < thresh.imp.refgndbracket(4)
   facecolor = imp.color(4,:);
elseif imp.gnd < thresh.imp.refgndbracket(5)
   facecolor = imp.color(5,:);
elseif imp.gnd < thresh.imp.refgndbracket(6)
   facecolor = imp.color(6,:);
elseif imp.gnd < thresh.imp.refgndbracket(7)
   facecolor = imp.color(7,:);
else
   facecolor = imp.color(8,:);
end
subplot('position',master.capplot(65,:))
rectangle('Curvature',[1 1],'Position',[.1 .1 .9 .9],'FaceColor',facecolor,'EdgeColor','c','LineWidth',2.5)
text(.55,1,'GND','Color','k','Units','normalized','FontSize',12,'FontWeight','b','HorizontalAlignment','center','VerticalAlignment','bottom')
text(.55,.55,impval,'Color',numcolor,'Units','normalized','FontSize',24,'FontWeight','b','HorizontalAlignment','center','VerticalAlignment','middle')
axis('square')
set(gca,'Visible','off')
 
export_fig('imp.png');
 
clear c facecolor impval numcolor screensize
clc; prog.step=prog.step+1; fprintf('step %g of %g successful \n%3.1f%% complete \n\n\n',prog.step,prog.totsteps,100*prog.step/prog.totsteps)
sound(audioread(strcat(full_path_to_the_FastQC_folder,'\cwb.mp3')),44000); pause(.5);
 
%% 47. Generate REC before cleaning time series figure 
 
plotdata = NaN(size(REC.data));
plotdata(:,REC.pre_parox_clean_inds) = REC.data(master.chanorderG,REC.pre_parox_clean_inds);
plotdata2 = REC.data;
plotdata2(:,REC.pre_parox_clean_inds) = NaN;
custom_eegplot('noui',plotdata, ... 
        'srate',REC.srate, ... 
        'winlength',round(length(REC.data)/REC.srate), ... 
        'spacing',100, ... 
        'color',master.colorsG, ... 
        'eloc_file',master.chanlocsG, ... 
        'data2',plotdata2, ...
        'position',get(0,'ScreenSize'));
         set(gcf,'visible','off')
clear plotdata*
pause(2)
export_fig('recdirty.png');
 
clear nullclear
clc; prog.step=prog.step+1; fprintf('step %g of %g successful \n%3.1f%% complete \n\n\n',prog.step,prog.totsteps,100*prog.step/prog.totsteps)
sound(audioread(strcat(full_path_to_the_FastQC_folder,'\cwb.mp3')),44000); pause(.5);
 
%% 48. Generate REC after cleaning time series figure 
 
plotdata = NaN(size(REC.data));
plotdata(:,REC.pre_parox_clean_inds) = REC.newdata(master.chanorderG,:);
custom_eegplot('noui',plotdata, ... 
        'srate',REC.srate, ... 
        'winlength',round(length(REC.data)/REC.srate), ... 
        'spacing',100, ... 
        'color',master.colorsG, ... 
        'eloc_file',master.chanlocsG, ... 
        'position',get(0,'ScreenSize'));
         set(gcf,'visible','off')
clear plotdata*
pause(2)
export_fig('recclean.png');
 
clear nullclear
clc; prog.step=prog.step+1; fprintf('step %g of %g successful \n%3.1f%% complete \n\n\n',prog.step,prog.totsteps,100*prog.step/prog.totsteps)
sound(audioread(strcat(full_path_to_the_FastQC_folder,'\cwb.mp3')),44000); pause(.5);
 
%% 49. Generate REO before cleaning time series figure 
 
plotdata = NaN(size(REO.data));
plotdata(:,REO.pre_parox_clean_inds) = REO.data(master.chanorderG,REO.pre_parox_clean_inds);
plotdata2 = REO.data;
plotdata2(:,REO.pre_parox_clean_inds) = NaN;
custom_eegplot('noui',plotdata, ... 
        'srate',REO.srate, ... 
        'winlength',round(length(REO.data)/REO.srate), ... 
        'spacing',100, ... 
        'color',master.colorsG, ... 
        'eloc_file',master.chanlocsG, ... 
        'data2',plotdata2, ...
        'position',get(0,'ScreenSize'));
         set(gcf,'visible','off')
clear plotdata*
pause(2)
export_fig('reodirty.png');
 
clear nullclear
clc; prog.step=prog.step+1; fprintf('step %g of %g successful \n%3.1f%% complete \n\n\n',prog.step,prog.totsteps,100*prog.step/prog.totsteps)
sound(audioread(strcat(full_path_to_the_FastQC_folder,'\cwb.mp3')),44000); pause(.5);
 
%% 50. Generate REO after cleaning time series figure 
 
plotdata = NaN(size(REO.data));
plotdata(:,REO.pre_parox_clean_inds) = REO.newdata(master.chanorderG,:);
custom_eegplot('noui',plotdata, ... 
        'srate',REO.srate, ... 
        'winlength',round(length(REO.data)/REO.srate), ... 
        'spacing',100, ... 
        'color',master.colorsG, ... 
        'eloc_file',master.chanlocsG, ... 
        'position',get(0,'ScreenSize'));
         set(gcf,'visible','off')
clear plotdata*
pause(2)
export_fig('reoclean.png');
 
clear nullclear
clc; prog.step=prog.step+1; fprintf('step %g of %g successful \n%3.1f%% complete \n\n\n',prog.step,prog.totsteps,100*prog.step/prog.totsteps)
sound(audioread(strcat(full_path_to_the_FastQC_folder,'\cwb.mp3')),44000); pause(.5);
 
%% 51. Generate REC spectra figure 
 
figure;
set(gcf,'Position',get(0,'ScreenSize'),'Color','w','visible','off')
hold on
for c=1:64
   plot(REC.freqs.values.all,REC.cspectra(c,:),'Color',char(master.colors(c)),'LineWidth',2)
end
set(gca,'YTick',[])
xlabel('Frequency (Hz)')
title('REC Spectra')
axes('Position',[.5 .6 .4 .3])
hold on
for c=1:64
   plot(REC.freqs.values.widealpha,REC.cspectra(c,REC.freqs.indices.widealpha),'Color',char(master.colors(c)),'LineWidth',2)
end
set(gca,'YTick',[])
xlabel('Frequency (Hz)')
title('Alpha Band')
export_fig('recspectra.png');
 
clear c
clc; prog.step=prog.step+1; fprintf('step %g of %g successful \n%3.1f%% complete \n\n\n',prog.step,prog.totsteps,100*prog.step/prog.totsteps)
sound(audioread(strcat(full_path_to_the_FastQC_folder,'\cwb.mp3')),44000); pause(.5);
 
%% 52. Generate REO spectra figure 
 
figure;
set(gcf,'Position',get(0,'ScreenSize'),'Color','w','visible','off')
hold on
for c=1:64
   plot(REO.freqs.values.all,REO.cspectra(c,:),'Color',char(master.colors(c)),'LineWidth',2)
end
set(gca,'YTick',[])
xlabel('Frequency (Hz)')
title('REO Spectra')
axes('Position',[.5 .6 .4 .3])
hold on
for c=1:64
   plot(REO.freqs.values.widealpha,REO.cspectra(c,REO.freqs.indices.widealpha),'Color',char(master.colors(c)),'LineWidth',2)
end
set(gca,'YTick',[])
xlabel('Frequency (Hz)')
title('Alpha Band')
set(gcf,'Color','w','Visible','off')
export_fig('reospectra.png');
 
clear c master
clc; prog.step=prog.step+1; fprintf('step %g of %g successful \n%3.1f%% complete \n\n\n',prog.step,prog.totsteps,100*prog.step/prog.totsteps)
sound(audioread(strcat(full_path_to_the_FastQC_folder,'\cwb.mp3')),44000); pause(.5);
 
%% 53. Generate ICA topographies figure 
 
figure;
set(gcf,'Position',get(0,'ScreenSize'),'Color','w','Visible','off')
for q=1:20
   subplot(4,10,q)
   set(gcf,'Color','w','Visible','off')
   topoplot(REC.icawinv_propev(:,q),REC.chanlocs(REC.icachansind),'headrad',0,'electrodes','off','maplimits',[-3 3],'style','map', ... 
   'emarker2',{[find(REC.icawinv_propev(:,q) == min(REC.icawinv_propev(:,q)),1) find(REC.icawinv_propev(:,q) == max(REC.icawinv_propev(:,q)),1)],'*','w',10,2});
   if ismember(q,REC.art.musc)
      artstr = 'Muscle';
   elseif ismember(q,REC.art.pulse)
      artstr = 'Pulse';
   elseif ismember(q,REC.art.vert)
      artstr = 'Vert. Eye';
   elseif ismember(q,REC.art.horiz)
      artstr = 'Horiz. Eye';
   else
      artstr = '';
   end
   if isempty(artstr)
      title(strcat('\color{blue}',sprintf('IC %g  EV %0.1f  %s',q,REC.qev(q),artstr)))
   else
      title(strcat('\color{red}',sprintf('IC %g  EV %0.1f  %s',q,REC.qev(q),artstr)))
   end
   set(gcf,'Color','w','Visible','off')
end
for q=1:20
   subplot(4,10,q+20)
   set(gcf,'Color','w','Visible','off')
   topoplot(REO.icawinv_propev(:,q),REO.chanlocs(REO.icachansind),'headrad',0,'electrodes','off','maplimits',[-3 3],'style','map', ... 
   'emarker2',{[find(REO.icawinv_propev(:,q) == min(REO.icawinv_propev(:,q)),1) find(REO.icawinv_propev(:,q) == max(REO.icawinv_propev(:,q)),1)],'*','w',10,2});
   if ismember(q,REO.art.musc)
      artstr = 'Muscle';
   elseif ismember(q,REO.art.pulse)
      artstr = 'Pulse';
   elseif ismember(q,REO.art.vert)
      artstr = 'Vert. Eye';
   elseif ismember(q,REO.art.horiz)
      artstr = 'Horiz. Eye';
   else
      artstr = '';
   end
   if isempty(artstr)
      title(strcat('\color{blue}',sprintf('IC %g  EV %0.1f  %s',q,REO.qev(q),artstr)))
   else
      title(strcat('\color{red}',sprintf('IC %g  EV %0.1f  %s',q,REO.qev(q),artstr)))
   end
   set(gcf,'Color','w','Visible','off')
end
export_fig('topos.png');
 
clear ans artstr q
clc; prog.step=prog.step+1; fprintf('step %g of %g successful \n%3.1f%% complete \n\n\n',prog.step,prog.totsteps,100*prog.step/prog.totsteps)
sound(audioread(strcat(full_path_to_the_FastQC_folder,'\cwb.mp3')),44000); pause(.5);
 
%% 54. Generate alpha topographies figure 
 
figure
set(gcf,'Position',get(0,'ScreenSize'),'Color','w','Visible','off')
subplot(1,2,1)
set(gcf,'Color','w','Visible','off')
topoplot(convert_raw_to_norm((mean(REC.cspectra(:,REC.freqs.indices.widealpha),2))'),REC.chanlocs,'headrad',.5,'electrodes','off','maplimits','minmax'); colorbar
title('REC','FontSize',12,'Color','b')
subplot(1,2,2)
set(gcf,'Color','w','Visible','off')
topoplot(convert_raw_to_norm((mean(REO.cspectra(:,REO.freqs.indices.widealpha),2))'),REO.chanlocs,'headrad',.5,'electrodes','off','maplimits','minmax'); colorbar
title('REO','FontSize',12,'Color','r')
export_fig('alphatopos.png');
 
clear ans
clc; prog.step=prog.step+1; fprintf('step %g of %g successful \n%3.1f%% complete \n\n\n',prog.step,prog.totsteps,100*prog.step/prog.totsteps)
sound(audioread(strcat(full_path_to_the_FastQC_folder,'\cwb.mp3')),44000); pause(.5);
 
%% 55. Attach figures and email report

tts('now emailing the figures')
sendmail(qc.report_will_send_to_these_emails,sprintf('Cohen RS EEG Fast QC Figures: %s',REC.pid),'',{'imp.png','recdirty.png','recclean.png','reodirty.png','reoclean.png','recspectra.png','reospectra.png','topos.png','alphatopos.png'});

clear nullclear
clc; prog.step=prog.step+1; fprintf('step %g of %g successful \n%3.1f%% complete \n\n\n',prog.step,prog.totsteps,100*prog.step/prog.totsteps)
sound(audioread(strcat(full_path_to_the_FastQC_folder,'\cwb.mp3')),44000); pause(.5);

%% 56. Convert figures into help base64f for html 
 
b64.imp = base64file('imp.png');
b64.recdirty = base64file('recdirty.png');
b64.recclean = base64file('recclean.png');
b64.reodirty = base64file('reodirty.png');
b64.reoclean = base64file('reoclean.png');
b64.recspectra = base64file('recspectra.png');
b64.reospectra = base64file('reospectra.png');
b64.topos = base64file('topos.png');
b64.alphatopos = base64file('alphatopos.png');
 
clear nullclear
clc; prog.step=prog.step+1; fprintf('step %g of %g successful \n%3.1f%% complete \n\n\n',prog.step,prog.totsteps,100*prog.step/prog.totsteps)
sound(audioread(strcat(full_path_to_the_FastQC_folder,'\cwb.mp3')),44000); pause(.5);
 
%% 57. Generate, open and email html figures report
 
report_name = strcat(REC.pid,'_REC_REO_FastQC_report.html');
fileID = fopen(report_name,'w');
fprintf(fileID,'%s','<body style="background-color:Pearl">');
fprintf(fileID,'%s',sprintf('<title>EEG RS QC %s</title><font size=5>',REC.pid));
fprintf(fileID,'%s',sprintf('Hello Cohen Team Member, <br><br>These are the Resting State EEG Fast Quality Control check figures for participant <b><font color="cornflowerblue" size=6>%s</font></b>',REC.pid));
% figures
fprintf(fileID,'%s','<br><br> Impedances<br>');
fprintf(fileID,'%s',strcat('<img src="data:image/png;base64,',b64.imp,'">'));
fprintf(fileID,'%s','<br><br> REC Before Cleaning<br>');
fprintf(fileID,'%s',strcat('<img src="data:image/png;base64,',b64.recdirty,'">'));
fprintf(fileID,'%s','<br><br> REC After Cleaning<br>');
fprintf(fileID,'%s',strcat('<img src="data:image/png;base64,',b64.recclean,'">'));
fprintf(fileID,'%s','<br><br> REO Before Cleaning<br>');
fprintf(fileID,'%s',strcat('<img src="data:image/png;base64,',b64.reodirty,'">'));
fprintf(fileID,'%s','<br><br> REO After Cleaning<br>');
fprintf(fileID,'%s',strcat('<img src="data:image/png;base64,',b64.reoclean,'">'));
fprintf(fileID,'%s','<br><br> REC Spectra<br>');
fprintf(fileID,'%s',strcat('<img src="data:image/png;base64,',b64.recspectra,'">'));
fprintf(fileID,'%s','<br><br> REO Spectra<br>');
fprintf(fileID,'%s',strcat('<img src="data:image/png;base64,',b64.reospectra,'">'));
fprintf(fileID,'%s','<br><br> ICA Topographies, REC Top 2 Rows, REO Bottom 2 Rows<br>');
fprintf(fileID,'%s',strcat('<img src="data:image/png;base64,',b64.topos,'">'));
fprintf(fileID,'%s','<br><br> Alpha Band Topographies (Normalized)<br>');
fprintf(fileID,'%s',strcat('<img src="data:image/png;base64,',b64.alphatopos,'">'));
fclose(fileID);
tts('now opening the html figures page in the web browser')
web(report_name,'-browser')
tts('now emailing the html report to Russ')
sendmail('rtoll@stanford.edu',sprintf('Cohen RS EEG Fast QC HTML: %s',REC.pid),'',{report_name});

clear ans b64 fileID report_name
clc; prog.step=prog.step+1; fprintf('step %g of %g successful \n%3.1f%% complete \n\n\n',prog.step,prog.totsteps,100*prog.step/prog.totsteps)
sound(audioread(strcat(full_path_to_the_FastQC_folder,'\cwb.mp3')),44000); pause(.5);
 
%% 58. Save and email data files to Russ

REC.data = [];
REC.rawdata = [];
REC.filtdata = [];
REC.newdata = [];
REC.times = [];
REC.pre_parox_clean_inds = [];
REC.totbin = [];
REO.data = [];
REO.rawdata = [];
REO.filtdata = [];
REO.newdata = [];
REO.times = [];
REO.pre_parox_clean_inds = [];
REO.totbin = [];
save('QC_Data','imp','qc','REC','REO','thresh');
tts('now emailing the data to Russ')
sendmail(qc.report_will_send_to_these_emails,sprintf('Cohen RS EEG Fast QC DATA: %s',REC.pid),'',{'QC_Data.mat'});

clear imp props qc REC REO thresh
clc; prog.step=prog.step+1; fprintf('step %g of %g successful \n%3.1f%% complete \n\n\n',prog.step,prog.totsteps,100*prog.step/prog.totsteps)
sound(audioread(strcat(full_path_to_the_FastQC_folder,'\cwb.mp3')),44000); pause(.5);

%% 59. Terminate and exit

clc
sound(audioread(strcat(full_path_to_the_FastQC_folder,'\cwb.mp3')),44000); pause(.5); 
disp('100% complete')
stoptime = toc;
send_text_message('512-663-3142','Verizon','QC Script Complete',sprintf('Runtime %02g:%02g',floor(stoptime/60),floor(rem(stoptime,60))));
tts('program complete, matlab will now close')

clear full_path_to_the_FastQC_folder prog stoptime
exit
