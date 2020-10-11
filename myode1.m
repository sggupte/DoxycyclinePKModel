function dydt = myode1(t,y,C,t_i,i)

    i = interp1(t_i,i,t);

    dydt = [i-(C(2).*y(1))-(C(1).*y(1)); 
        C(1).*y(1)-C(5).*y(2)-C(4).*y(2)+C(3).*y(3);
        C(4).*y(2)-C(3).*y(3)];
end