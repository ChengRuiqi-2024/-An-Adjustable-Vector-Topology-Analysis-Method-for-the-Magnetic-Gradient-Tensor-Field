function p = collect_points(mask,X,Y,result)
    [p.row,p.col] = find(mask);
    p.x = X(mask);
    p.y = Y(mask);
    p.index = result.index(mask);
    p.score = result.score(mask);
end