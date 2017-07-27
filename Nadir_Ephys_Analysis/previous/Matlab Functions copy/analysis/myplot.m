function myplot

global SWEEP_DATA TRI
hold off;
plot (SWEEP_DATA(1:3000,end-5:end),'k');
end