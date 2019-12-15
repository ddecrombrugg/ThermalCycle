function [vect vech vecs] = CCGTPlot(data)%use,p_1,p_2,s_1,h_1,h_2,eta_Si
    if data(1) == 1
        vecp=data(2):-0.001:data(3);
        h_2s=zeros(length(vecp),1);
        vecs=zeros(length(vecp),1);
        vect=zeros(length(vecp),1);
        % Isentropic values !
        s_2s=data(4)*ones(length(vecp),1);% same as point 3
        for i=1:length(vecp)
            h_2s(i)=XSteam('H_ps',vecp(i),s_2s(i)); %hs_4
        end
        vech=data(5)-data(7).*(data(5)-h_2s); %h_4
        for i=1:length(vecp)
            vecs(i)=XSteam('s_ph',vecp(i),vech(i)); %s_4
            vect(i)=XSteam('T_hs',vech(i),vecs(i));
        end
    elseif data(1) == 2
        if data(6)>data(5)
            vech=data(1):0.5:data(6);
        else
            vech=data(6):0.5:data(1);
        end
        vecp=data(3)*ones(length(vech),1);
        vect=zeros(length(vecp),1);
        vecs=zeros(length(vecp),1);
        for i=1:length(vecp)
            vect(i)=XSteam('T_ph',vecp(i),vech(i));
            vecs(i)=XSteam('s_ph',vecp(i),vech(i));
        end
    end
end