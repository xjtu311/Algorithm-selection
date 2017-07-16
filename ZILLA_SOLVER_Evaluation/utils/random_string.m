function str = random_string
s = rand('twister');
rn = rand;
rand('twister',s);
str = datestr(now, 'mmmm_dd_yyyy_HH_MM_SS_FFF');
str = strcat(str, '-', num2str(rn));
