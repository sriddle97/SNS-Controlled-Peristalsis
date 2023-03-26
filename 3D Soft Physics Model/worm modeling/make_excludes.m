function make_excludes(helix_vecs)
    
    dummy = helix_vecs("left_handed");
    LH_helixes = dummy{1,1};
    dummy = helix_vecs("right_handed");
    RH_helixes = dummy{1,1};
    dummy = helix_vecs("intersection_points");
    intersection_pts = dummy{1,1};
    
    for i = 1:length(LH_helixes)
        for ii = 1:length(RH_helixes)
            [~,b,~] = intersect(LH_helixes{i},RH_helixes{ii},'rows');
            for iii = 1:length(b)
                fprintf('\n<exclude body1="%s" body2="%s" />',append('LH',num2str(i-1),'B_',num2str(b(iii)-2)),append('RH',num2str(ii-1),'B_',num2str(b(iii)-2)))
                fprintf('\n<exclude body1="%s" body2="%s" />',append('LH',num2str(i-1),'B_',num2str(b(iii)-2)),append('RH',num2str(ii-1),'B_',num2str(b(iii)-1)))
                fprintf('\n<exclude body1="%s" body2="%s" />',append('LH',num2str(i-1),'B_',num2str(b(iii)-1)),append('RH',num2str(ii-1),'B_',num2str(b(iii)-2)))
                fprintf('\n<exclude body1="%s" body2="%s" />\n',append('LH',num2str(i-1),'B_',num2str(b(iii)-1)),append('RH',num2str(ii-1),'B_',num2str(b(iii)-1)))            
            end
        end
    end
end