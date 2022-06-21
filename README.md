# powerLineCalc
power Line wire Calc
//弦截法
extern "C" fn rootsec(x1:f64,x2:f64,eps:f64,f:Callback)->f64
//Aitken法
extern "C" fn root_Aitken(x:f64,interation:usize,eps:f64,f:Callback)->f64
//一维brent法
extern "C" fn zbrent(x1:f64,x2:f64,tol:f64,f:Callback)->f64
//孤立档线长系数计算
extern "C" fn  koljgq(n:i32,span:f32,h:f32,a:f32,em:f32,s_weight_1:f32,s_weight_2:f32,s_length_1:f32,s_length_2:f32,gama_x:f32, a_a:[f32;100], q_q:[f32;100], nq:i32, calmode:i32)->f32  //arr: &mut [[f64;6];3]
//孤立档判断用Fm计算
extern "C" fn  fm_gld(n:i32,span:f32,h:f32,a:f32,em:f32,s_weight_1:f32,s_weight_2:f32,s_length_1:f32,s_length_2:f32,linear_co:f32,gama_m:f32,sigma_m:f32,t_m:f32)->f32
//孤立档安装工况Fm计算
extern "C" fn  fm_gld_e(n:i32,span:f32,h:f32,a:f32,em:f32,s_weight_1:f32,s_weight_2:f32,s_length_1:f32,s_length_2:f32,linear_co:f32,gama_m:f32,sigma_m:f32,t_m:f32)->f32
//孤立档x点弧垂计算
extern "C" fn  sag_x(n:i32,s_weight_1:f32,s_weight_2:f32,s_length_1:f32,s_length_2:f32,span:f32,h:f32,area:f32,sig0:f32, x:f32, gama_x:f32, a_a:[f32;100], q_q:[f32;100], nq:i32, calmode:i32)
-> f32
//孤立档Fm最大值计算
extern "C" fn gld_fm_max(n:i32,span:f32,h:f32,a:f32,em:f32,linear_co:f32,s_ice_co:f32,s_weight_1:f32,s_weight_2:f32,s_length_1:f32,s_length_2:f32,deta_l:f32,az_wj:f32,gama:[f32;5],qxt:[f32;5],sigma_max:f32,sigma_av:f32,calmode:i32,a_a:[f32;100], q_q:[f32;100],nq:i32,kzs:&mut Contr_Str)-> *const c_char  //##calmode==1 竣工 shi aA=[],qQ=[],nq=0
 //孤立档最大弧垂计算
 extern "C" fn gld_sag_max(span:f32,h:f32,area:f32,s_weight_1:f32,s_weight_2:f32,s_length_1:f32,s_length_2:f32,lj:f32,sig0:f32, gama_x:f32,  calmode:i32)
->f32
//ridders改进试位法计算
extern "C" fn ridders(x1:f64,x2:f64,eps:f64,f:Callback)->f64
//不动点迭代法
extern "C" fn rtdeadpoint(x0:f64,tol:f64,f:Callback) -> f64 
//二分法
extern "C" fn rootbis(x1:f64,x2:f64,tol:f64,f:Callback)->f64
//flsp法
extern "C" fn rootflsp(x1:f64,x2:f64,eps:f64,f:Callback)->f64
//黄金分割法
extern "C" fn golden(mut x1:f64,mut x2:f64,tol:f64,f:Callback,res:&mut [f64;2])->bool
golden2( ax:f64, bx:f64, cx:f64,tol:f64,f:Callback,res:&mut [f64;2])->bool
extern "C" fn brent_golden( ax:f64, bx:f64, cx:f64,tol:f64,f:Callback,res:&mut [f64;2])->bool
//brak法
extern "C" fn mnbrak( mut ax:f64, mut bx:f64,glimit:f64,f:Callback,res:&mut [f64;6])->bool
//拉格朗日插值
extern "C" fn lagrange_interp_100( x:f64,n:usize, a:&[f64;100], b:&[f64;100])->f64
//
extern "C" fn lagrange_interp_1000( x:f64,n:usize, a:&[f64;1000], b:&[f64;1000])->f64
extern "C" fn lagrange_interp_10000( x:f64,n:usize, a:&[f64;10000], b:&[f64;10000])->f64
extern "C" fn inlagrn_100( t:f64,n:usize, x:&[f64;100], y:&[f64;100])->f64
extern "C" fn inlagrn_three_100( t:f64,n:usize, x:&[f64;100], y:&[f64;100])->f64
//追赶法求三对角矩阵
extern "C" fn chasing_100( n:usize, a:&[f64;100], b:&[f64;100], c:&[f64;100], d:&[f64;100],x: &mut [f64;100])->bool
extern "C" fn chasing_1000( n:usize, a:&[f64;1000], b:&[f64;1000], c:&[f64;1000], d:&[f64;1000],x: &mut [f64;1000])->bool
extern "C" fn chasing_10000( n:usize, a:&[f64;10000], b:&[f64;10000], c:&[f64;10000], d:&[f64;10000],x: &mut [f64;10000])->bool
//
extern "C" fn m_calc_100( n:usize,ff0:f64,ffn:f64, x_lst:&[f64;100],y_lst:&[f64;100], m: &mut [f64;100])->bool
extern "C" fn m_calc_1000( n:usize,ff0:f64,ffn:f64, x_lst:&[f64;1000],y_lst:&[f64;1000], m: &mut [f64;1000])->bool
extern "C" fn m_calc_10000( n:usize,ff0:f64,ffn:f64, x_lst:&[f64;10000],y_lst:&[f64;10000], m: &mut [f64;10000])->bool
//插值三次样条
extern "C" fn cubic_sp_interpo_10000( n_x_:usize,n_x_lst:usize,m:&[f64;10000], x_:&[f64;10000],y_:&[f64;10000], x_lst: &[f64;10000],y_lst: &mut [f64;10000])->bool
//二维插值
extern "C" fn interpo2d_100(n:usize,m:usize,u:f32,v:f32,x:&[f32;100],y:&[f32;100],z:&[[f32;100];100])->f32
//
extern "C" fn interpo2d(n:usize,m:usize,u:f32,v:f32,x:&[f32;10000],y:&[f32;10000],z:&[[f32;10000];10000])->f32
//最小二乘法拟合多项式，int m = a.size();	//拟合多项式的项数，即拟合多项式的最高次数为m-1
//要求m<=n且m<=20。若m>n或m>20，本函数自动按m=min{n,20}处理。
//a,m-1次多项式系数,dt 0-误差平方和，1-误差绝对值之和，2-误差绝对值的最大值
#[no_mangle]
extern "C" fn fit_curve_least_squares_m(n:usize,mut m:usize,x:&[f64;201],y:&[f64;201],a:&mut [f64;20],dt:&mut [f64;3])->bool

//切比雪夫曲线拟合,最佳拟合多项式，给定点的偏差最大值为最小
	//int n = x.size();	//给定数据点的个数
	//int m = a.size()-1;	//拟合多项式的项数，即拟合多项式的最高次数为m-1
	//要求m<n且m<=20。若m>=n或m>20，本函数自动按m=min{n-1,20}处理。
#[no_mangle]
extern "C" fn fit_curve_chebyshev_m(n:usize,mut m:usize,x:&[f64;201],y:&[f64;201],a:&mut [f64;20])->bool

//最佳一致逼近多项式里米兹法
//[a,b]为区间，n为阶数=p.size()-1,int n = p.size()-1;	//n-1次最佳一致逼近多项式的项数
//要求n<=20; 若n>20，函数自动取n=20
#[no_mangle]
extern "C" fn approximation_remez_m(mut n:usize,eps:f64,a:f64,b:f64,f:Callback,p:&mut [f64;21])->bool
//矩形域的最小二乘曲面拟合
//int n = x.size();			//给定数据点的X坐标个数
//int m = y.size();			//给定数据点的Y坐标个数
//int p = a.GetRowNum();		//拟合多项式中变量x的最高次数加1
//并要求p<=n且p<=20; 否则在本函数中自动取p=min{n,20}处理
//int q = a.GetColNum();		//拟合多项式中变量y的最高次数加1
//并要求q<=m且q<=20; 否则在本函数中自动取q=min{m,20}处理
#[no_mangle]
extern "C" fn fit_surface_least_squares_m( n:usize, m:usize,mut p:usize,mut q:usize,x:&[f64;100],y:&[f64;100],z:&[[f64;100];100],a:&mut [[f64;20];20],dt:&mut [f64;3])->bool
//二维多项式   dCoff为二维多项式系数数组，dX为自变量x，dY为自变量y  stx_nno  //行数   sty_nno//列数
#[no_mangle]
extern "C" fn poly_value_two_dims(stx_nno:usize,sty_nno:usize,dx:f64,dy:f64,dcoff:& [[f64;20];20])->f64
//双线性插值
extern "C" fn bilinear_interpo_m(x:f64,y:f64,x_ary:&[f64;2],y_ary:&[f64;2],z_ary:&mut [[f64;100];100],b:&mut [f64;100])->f64
//线性方程组 高斯消元法
extern "C" fn le_total_choice_gauss_m(n:usize,a:&mut [[f64;100];100],b:&mut [f64;100])->bool
//bla坐标转xyz
extern "C" fn bla2xyz(zone_nnumber:i32,lat:f64,lon:f64,alt:f64,res:&mut [f64;3])->bool
//代表档距计算
extern "C" fn r_span(ary:&mut [u32;100],n:usize)->i32
//有高差的弧垂计算
extern "C" fn sag_h(gama_x:f64,sigma_x:f64,span:f64,h:f64)->f64
//连续档Fm求解
extern "C" fn fm(span:f64,gama_m:f64,sigma_m:f64,t_m:f64,linear_co:f64,elastic_modulus:f64)->f64
//导线参数翻译
extern "C" fn cdrs(json:*const c_char,arr: &mut [[f64;6];3],st:&mut CCdrs)-> bool
//气象参数翻译
extern "C" fn qx(json:*const c_char,arr: &mut [[f64;3];14],rtst:&mut Qx)-> *const c_char
//三次方程解析求解
extern "C" fn root_cubic_simple(b:f64,c:f64,d:f64)->f64
//判断经纬度在中国区内
extern "C" fn out_of_china(lat: f64, lng: f64) -> bool 
//// wgs2gcj convert WGS-84 coordinate(wgs_lat, wgs_lng) to GCJ-02 coordinate.
#[no_mangle]
extern "C" fn wgs2gcj(wgs_lat: f64, wgs_lng: f64) -> (f64,f64) 
// gcj2wgs convert GCJ-02 coordinate(gcj_lat, gcj_lng) to WGS-84 coordinate.
// The output WGS-84 coordinate's accuracy is 1m to 2m. If you want more exactly result, use gcj2wgs_exact.
#[no_mangle]
extern "C" fn gcj2wgs(gcj_lat: f64, gcj_lng: f64) -> (f64,f64)
// gcj2wgs_exact convert GCJ-02 coordinate(gcj_lat, gcj_lng) to WGS-84 coordinate.
// The output WGS-84 coordinate's accuracy is less than 0.5m, but much slower than gcj2wgs.
#[no_mangle]
extern "C" fn gcj2wgs_exact(gcj_lat: f64, gcj_lng: f64,res:&mut [f64;2]) ->bool
// distance calculate the distance between point(lat_a, lng_a) and point(lat_b, lng_b), unit in meter.
#[no_mangle]
extern "C" fn distance(lat_a: f64, lng_a: f64, lat_b: f64, lng_b: f64) -> f64 
 // Converts degrees to radians.
#[no_mangle]
 extern "C" fn deg_to_rad(deg:f64)->f64
 //Converts radians to degrees.
 extern "C" fn rad_to_deg(rad:f64)->f64
 /*
 * arc_length_of_meridian
 *
 * Computes the ellipsoidal distance from the equator to a point at a
 * given latitude.
 *
 * Reference: Hoffmann-Wellenhof, B., Lichtenegger, H., and Collins, J.,
 * GPS: Theory and Practice, 3rd ed.  New York: Springer-Verlag Wien, 1994.
 *
 * Inputs:
 *     phi - Latitude of the point, in radians.
 *
 * Globals:
 *     WGS84_A - Ellipsoid model major axis.
 *     WGS84_B - Ellipsoid model minor axis.
 *
 * Returns:
 *     The ellipsoidal distance of the point from the equator, in meters.
 *
 */
#[no_mangle]
 extern "C" fn arc_length_of_meridian(phi:f64)->f64
 extern "C" fn arc_length_of_meridian_u(a:f64,f:f64,phi:f64)->f64
 extern "C" fn utm_central_meridian(zone:i32)->f64
//foot_point_latitude
extern "C" fn foot_point_latitude(a:f64,f:f64,y:f64)->f64
//map_lat_lon_to_xy
extern "C" fn map_lat_lon_to_xy(a:f64,f:f64,phi:f64,lambda:f64,lambda0:f64,xy:&mut UTMCoor)
//LatLonToUTMXY
extern "C" fn lat_lon_to_utm_xy( lat:f64,  lon:f64,  zone:i32, xy:&mut UTMCoor)
//map_xy_to_lat_lon
extern "C" fn map_xy_to_lat_lon(a:f64,f:f64, x:f64, y:f64, lambda0:f64,philambda:&mut WGS84Corr)
//UTMXYToLatLon
extern "C" fn  utm_xy_to_lat_lon( mut x:f64,  mut y:f64,  zone:i32,  southhemi:bool, latlon:&mut WGS84Corr)

extern "C" fn geo_direct(lat1: f64, lon1: f64, azi1: f64, s12: f64,res:&mut [f64;3])
 

