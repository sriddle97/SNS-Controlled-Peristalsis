function make_connects(helix_vecs)
    
    dummy = helix_vecs("left_handed");
    LH_helixes = dummy{1,1};
    dummy = helix_vecs("right_handed");
    RH_helixes = dummy{1,1};
    dummy = helix_vecs("intersection_points");
    intersection_pts = dummy{1,1};

    for i = 1:length(intersection_pts)
        for ii = 1:length(RH_helixes)
            [~,a,b] = intersect(LH_helixes{ii},intersection_pts{i},'rows');
            [~,c,d] = intersect(RH_helixes{ii},intersection_pts{i},'rows');
            for iii = 1:length(a)
                fprintf('\n<connect body1="%s" body2="%s" anchor="0 0 0"/>',append('ball',num2str(i-1),num2str(b(iii)-1)),append('LH',num2str(ii-1),'B_',num2str(a(iii)-2)))
            end
            for iii = 1:length(c)
                fprintf('\n<connect body1="%s" body2="%s" anchor="0 0 0"/>',append('ball',num2str(i-1),num2str(d(iii)-1)),append('RH',num2str(ii-1),'B_',num2str(c(iii)-2)))
            end
        end
    end
    fprintf('\n')
end