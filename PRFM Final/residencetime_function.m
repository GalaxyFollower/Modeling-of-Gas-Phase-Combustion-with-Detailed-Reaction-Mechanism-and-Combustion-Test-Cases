function func_tr = residencetime_function(h,tr,Vx,tr_prev)

func_tr = tr - tr_prev -h/Vx;
