function dydt = myode2(t,y,C,t_i,i)

    i = interp1(t_i,i,t);

    dydt = [i-(C(2).*y(1))-(C(1).*y(1)); 
        C(1).*y(1)-C(5).*y(2)-C(4).*y(2)+C(3).*y(3);
        C(4).*y(2)-C(3).*y(3);
        C(6).*y(4)-C(7).*((C(8).*y(3))/(C(8).*y(3)+C(9))).*y(4)-1/C(10).*y(4)];
end