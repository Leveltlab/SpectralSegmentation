function RS = smoothG(RDtemp, se)
%
%RS = smoothG(RDtemp, se)
%se is the standard deviation of the smoothing gaussian
%column wise smoothing!!!
%
%
%Chris van der Togt 
%Dept. of Vision and Cognition 
%Netherlands Institute for Neurosciences
%10/2012


[Segment, Nchan] = size(RDtemp);
r = -3*se:3*se;
G = 1/sqrt(2*pi*se).*exp(-r.^2/(2*se*se));
HS = conv(ones(Segment,1), G,  'same');


 RS = RDtemp;
 for i = 1:Nchan
     RS(:,i) = conv(double(RDtemp(:,i)), G,  'same')./HS;
 end
 
 

