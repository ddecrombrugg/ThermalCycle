function cp_moy = cp_fumes(t_3,t_4,comp_f)
step = 1;
T_3 = t_3-273.15;
T_4 = t_4-273.15;
T_c3 = [300.*ones(floor(300-273.15)+1,1)' 301:step:t_3];
T_c4 = [300.*ones(floor(300-273.15)+1,1)' 301:step:t_4];
c_pa_moy3 = mean(janaf('CO2',T_c3)*comp_f(1) + janaf('H2O',T_c3)*comp_f(2)...
              + janaf('O2',T_c3)*comp_f(3) + janaf('N2',T_c3)*comp_f(4));
c_pa_moy4 = mean(janaf('CO2',T_c4)*comp_f(1) + janaf('H2O',T_c4)*comp_f(2)...
              + janaf('O2',T_c4)*comp_f(3) + janaf('N2',T_c4)*comp_f(4));
          
cp_moy = (c_pa_moy4*T_4-c_pa_moy3*T_3)/(T_4-T_3);
end