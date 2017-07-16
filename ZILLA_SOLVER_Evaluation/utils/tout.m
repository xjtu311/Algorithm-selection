function tout(str)
global log_fid
global traj_fid
global detailed_traj_fid
fprintf(log_fid, str);
fprintf(traj_fid, str);
fprintf(detailed_traj_fid, str);
fprintf(str);