function r = MC_runs(zeroone)
% r = runs(zeroone)
% r = runs(index)
%    zeroone is a logical matrix or index is a sorted integer matrix (like the output of find)
%    returns a matrix with two columns.  The first is the index of the beginning of
%    each run, and the second is the length of the run.

if isempty(zeroone)
  r = zeros(2,0);;
  return;
end;

if islogical(zeroone)
  index = find(zeroone);
else
  index = zeroone;
end;

d = diff(index);

if size(d,2)==1,
    d=d';
end

not_one = [0 find(d ~= 1) length(d)+1];

for i=1:length(not_one)-1
  r(1,i) = index(not_one(i)+1);
  r(2,i) = index(not_one(i+1))-r(1,i)+1;
end;

if isempty(r)
    r=zeros(2,0);
end

