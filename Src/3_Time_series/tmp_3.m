%% Calculate DC data
cut_point=(1745-1740)/0.025;
mask=mask(1:cut_point, :);
figure(3)
imagesc(mask);  title('Channel status'); xlabel('Frequency'); ylabel('Time slot');

[fn,tn]=size(mask); dc=zeros(tn,1);
for i=1:tn
    dc(i,1)=sum(mask(i,:))/fn;
end
figure(4)
plot(1:tn,dc); title('Duty cycle');