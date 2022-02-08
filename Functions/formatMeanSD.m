function converted = formatMeanSD(raw_mean,raw_sd,precision,scaling)
raw_mean = raw_mean*scaling;
raw_sd = raw_sd*scaling;
first_sig_digit = fix(log10(abs(raw_mean)));
if first_sig_digit > 2
    m = fix(raw_mean);
    sd = fix(raw_sd);
    converted = sprintf(['%i' char(177) '%i'],m,sd);
elseif raw_mean == 0
    converted = '0';
else
    round_idx = precision-first_sig_digit-1;
    m = round(raw_mean,precision,'significant');
    sd = round(raw_sd,round_idx,'decimals');
    spec = ['%.' num2str(round_idx) 'f'];
    converted = sprintf([spec ' ' char(177) ' ' spec],m,sd);
end