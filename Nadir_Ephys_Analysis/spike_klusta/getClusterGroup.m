function myclusters = getClusterGroup(groups,groupnumber)

myclusters = groups(find(groups(:,2)==groupnumber),1);
end
