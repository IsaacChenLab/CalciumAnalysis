function [new_groups,ngroups]=MC_spSetGroups(action,varargin)

if nargin<2
    error('arguments');
end

old_groups=varargin{1};
new_groups=old_groups;
ngr=length(old_groups);

switch action
    case 'add'
        add_groups=varargin{2};
        for k=1:length(add_groups)
            new_groups(ngr+k)=add_groups(k);
        end
        if nargin==4
            instead_group=varargin{3};
            new_groups=MC_spSetGroups('delete',new_groups,instead_group);        
        end
    case 'merge'
        grs=varargin{2};
        for i=2:length(grs)
            if size(new_groups{1},1)==1
                new_groups{grs(1)}=[new_groups{grs(1)} new_groups{grs(i)}];
            else
                new_groups{grs(1)}=[new_groups{grs(1)}; new_groups{grs(i)}];
            end
        end
        if length(grs)>1
            new_groups=MC_spSetGroups('delete',new_groups,grs(2:end));
        end
    case 'delete'
        grs=varargin{2};
        for i=grs
            new_groups{i}=[];
        end
        new_groups=MC_spSetGroups('purge',new_groups);
    case 'purge'
        notempty=[];
        for i=1:length(new_groups),
            if ~isempty(new_groups{i})
                notempty=[notempty i];
            end
        end
        new_groups=new_groups(notempty);
    otherwise
        error('action not recognized');
end

ngroups=length(new_groups);

return;
