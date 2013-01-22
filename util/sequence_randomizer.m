function [sampling]=sequence_randomizer(length,items)
% generates a random orderring for sequence learning experiments. Enter the
% length of the sequence you need, and the number of response buttons. The
% function will return a pseudorandom order such that each trial type
% occurs and then set of items is selected. This function also eliminates
% the repetitions that would normally occur from one sampling of n
% responses to the next sampling.

repetitions=length/items;

if (mod(length,items)~=0)
    warning('While operationa, be aware that irregular dimensions have been provided. \n Length=%d is not an integer multiple of %d items.',length, items)
    repetitions=ceil(repetitions);
    l=items*repetitions;
end;

% generate random permutations
for i=1:repetitions
    sampling(i,:)=randperm(items);
end;


% check for repeats that will occur if run in order
for i=2:repetitions
    if(sampling(i,1)==sampling(i-1,4))
        %resample if there is a transitional repeat
        sampling(i,:)=randperm(4);
        %reduce i so that any newly generated item is checked
        i=i-1;
    end;
end;

%resize the sampling into column (reshape works down columns)
sampling=reshape(sampling',[repetitions*items,1]);

%trim off any extras due to irregular entry
sampling=sampling(1:length,:);

return;