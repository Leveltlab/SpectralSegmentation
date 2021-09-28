function SysMem = getSystemMem()
% get available memory

SysMem = [];
archstr = computer('arch');
%only for fwindows available
if strcmp(archstr, 'win64')
    usr = memory();
    SysMem = usr.MemAvailableAllArrays;
    
elseif strcmp(archstr, 'glnxa64') 
     [~, Memstr] = system("free -b | awk '/^Mem/ {print $4}'");
     SysMem = str2double(Memstr);
    
elseif strcmp(archstr, 'maci64')
    [~, Memstr] = system("sysctl hw.memsize | awk '{print $2}'");
    SysMem = str2double(Memstr);
end