function b = beh_analysis(sub)
s = 0;
for ss = sub 
    s = s+1;
cd (['data/sbj' num2str(ss)])
fname = dir('trip_s*');
label = {'color_redblue', 'ori_VH', 'prior_705030', 'bias_no_color_ori_resp' ...
    'expResp', 'tooslow', 'hit', 'resp', 'falsealarm', 'RT_fromtgonset'};
for r = 1:16
    load (fname(r).name);
    if r == 1
        for i = 1:10
            pcat = p;
        end
    else
        for i = 1:10
    pcat.(label{i}) = cat(1, pcat.(label{i}), p.(label{i}));
        end
    end
    
end

for i = 1: size(pcat.hit,1)
    if pcat.bias_no_color_ori_resp(i, 1) == 0 
        pcat.pr_color(i, 1) =2;  pcat.pr_ori(i, 1) =2; pcat.pr_resp(i, 1) =2; 
    elseif pcat.bias_no_color_ori_resp(i, 1) == 1 
        pcat.pr_color(i, 1) = pcat.prior_705030(i,1);  pcat.pr_ori(i, 1) =2; pcat.pr_resp(i, 1) =2; 
    elseif pcat.bias_no_color_ori_resp(i, 1) == 2 
        pcat.pr_ori(i, 1) = pcat.prior_705030(i,1);  pcat.pr_color(i, 1) =2; pcat.pr_resp(i, 1) =2; 
    elseif  pcat.bias_no_color_ori_resp(i, 1) == 3
        pcat.pr_resp(i, 1) = pcat.prior_705030(i,1);  pcat.pr_color(i, 1) =2; pcat.pr_ori(i, 1) =2; 
    end
end


for rb = 1:2 
    for pr = 1:3
b.hit_color(s, rb , pr) = nanmean(pcat.hit(pcat.color_redblue == rb & pcat.pr_color == pr));
b.rt_color(s, rb  , pr) = nanmean(pcat.RT_fromtgonset(pcat.color_redblue == rb & pcat.pr_color == pr & pcat.hit ==1));

b.hit_ori(s, rb , pr) = nanmean(pcat.hit(pcat.ori_VH == rb & pcat.pr_ori == pr));
b.rt_ori(s, rb , pr) = nanmean(pcat.RT_fromtgonset(pcat.ori_VH== rb & pcat.pr_ori == pr & pcat.hit ==1));

b.hit_resp(s, rb , pr) = nanmean(pcat.hit(pcat.expResp == rb & pcat.pr_resp == pr));
b.rt_resp(s, rb , pr) = nanmean(pcat.RT_fromtgonset(pcat.expResp == rb & pcat.pr_resp == pr & pcat.hit ==1));

    end
end


b.hit_color(s, 3 ,:) = nanmean(b.hit_color(s, 1:2 ,:),2); b.rt_color(s, 3 ,:) = nanmean(b.rt_color(s, 1:2 ,:),2);
b.hit_ori(s, 3 ,:) = nanmean(b.hit_ori(s, 1:2 ,:),2); b.rt_ori(s, 3 ,:) = nanmean(b.rt_ori(s, 1:2 ,:),2);
b.hit_resp(s, 3 ,:) = nanmean(b.hit_resp(s, 1:2 ,:),2); b.rt_resp(s, 3 ,:) = nanmean(b.rt_resp(s, 1:2 ,:),2);
cd ../..
end

b.hit_color(:, 3 ,:) = b.hit_color(:, 3, :) - repmat(nanmean(b.hit_color(:, 3, :),3), [1,1,3]) +...
    repmat(nanmean(nanmean(b.hit_color(:, 3, :), 3), 1), [s, 1, 3]);

b.hit_ori(:, 3 ,:) = b.hit_ori(:, 3, :) - repmat(nanmean(b.hit_ori(:, 3, :),3), [1,1,3]) +...
    repmat(nanmean(nanmean(b.hit_ori(:, 3, :), 3), 1), [s, 1, 3]);
b.hit_resp(:, 3 ,:) = b.hit_resp(:, 3, :) - repmat(nanmean(b.hit_resp(:, 3, :),3), [1,1,3]) +...
    repmat(nanmean(nanmean(b.hit_resp(:, 3, :), 3), 1), [s, 1, 3]);

b.rt_color(:, 3 ,:) = b.rt_color(:, 3, :) - repmat(nanmean(b.rt_color(:, 3, :),3), [1,1,3]) +...
    repmat(nanmean(nanmean(b.rt_color(:, 3, :), 3), 1), [s, 1, 3]);
b.rt_ori(:, 3 ,:) = b.rt_ori(:, 3, :) - repmat(nanmean(b.rt_ori(:, 3, :),3), [1,1,3]) +...
    repmat(nanmean(nanmean(b.rt_ori(:, 3, :), 3), 1), [s, 1, 3]);
b.rt_resp(:, 3 ,:) = b.rt_resp(:, 3, :) - repmat(nanmean(b.rt_resp(:, 3, :),3), [1,1,3]) +...
    repmat(nanmean(nanmean(b.rt_resp(:, 3, :), 3), 1), [s, 1, 3]);

