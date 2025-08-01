function [t,p] = p5x2cv(score_1,score_2)

    p_i = score_1 - score_2;
    p_1_1 = p_i(1);
    p_i_bar = (p_i(1,:) + p_i(2,:)) / 2;
    s_i_sqr = (p_i(1,:) - p_i_bar).^2 + (p_i(2,:) - p_i_bar).^2;
    s_sqr = sum(s_i_sqr);
    t = p_1_1 / ((s_sqr / 5)^0.5); 
    p = 2*(1-tcdf(abs(t),5));
end