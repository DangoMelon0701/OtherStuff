function I=trapecio(fu,hi)
    I=hi*(sum(fu)-(fu(1)+fu(length(fu)))/2);
end