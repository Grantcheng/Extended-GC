%this data read function is used in read the C++ output results."*.txt"
function Data=Dadaread(fname)
fid=fopen([fname '.txt']);
txt={};
tline=fgetl(fid);
while ischar(tline) % ischar:Determine whether item is character array
    txt{end+1}=tline;
    tline=fgetl(fid);
end
fclose(fid);
Data={};
for i=1:length(txt)
Data{i}=str2num(txt{i});
end