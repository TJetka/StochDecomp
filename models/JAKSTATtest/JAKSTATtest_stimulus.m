function result = JAKSTATtest_stimulus(t,i)
    result=[];
    stim{1}=5*t-10+sin(t);
    stim{2}=0;
    stim{3}=10*((t>10)&(t<15));
    result=stim{i};
end

