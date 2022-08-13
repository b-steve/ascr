#include <TMB.hpp>

//DIY functions
template<class Type>
Type bessi0(Type x){
  Type ax;
  Type ans;
  Type y;
  if ((ax = fabs(x)) < 3.75){
    y = x/3.75;
    y *= y;
    ans = 1.0 + y*(3.5156229 + y*(3.0899424 + y*(1.2067492 + y*(0.2659732 + y*(0.360768e-1 + y*0.45813e-2)))));
  } else {
    y = 3.75/ax;
    ans = (exp(ax)/sqrt(ax))*(0.39894228 + y*(0.1328592e-1 + y*(0.225319e-2 + 
      y*(-0.157565e-2 + y*(0.916281e-2 + y*(-0.2057706e-1 + y*(0.2635537e-1 + 
      y*(-0.1647633e-1 + y*0.392377e-2))))))));
  }
  return ans;
}


template<class Type>
Type Gamma(Type z){
  if(z < 0.5){
    z = 1 - z;
    Type y = sqrt(2 * 3.14159265358979323846) * pow(z + 6.5, (z - 0.5)) * exp((-1) * z - 6.5) * 
      (0.99999999999980993 + 676.5203681218851/(z) + (-1259.1392167224028)/(z + 1) + 771.32342877765313/(z + 2) +
      (-176.61502916214059)/(z + 3) + 12.507343278686905/(z + 4) + (-0.13857109526572012)/(z + 5));
    //+ 9.9843695780195716e-6/(z + 6) + 1.5056327351493116e-7/(z + 7));
    
    return 3.14159265358979323846 / (sin(3.14159265358979323846 * (z - 1)) * y);
  } else {
    Type y = sqrt(2 * 3.14159265358979323846) * pow(z + 6.5, (z - 0.5)) * exp((-1) * z - 6.5) * 
      (0.99999999999980993 + 676.5203681218851/(z) + (-1259.1392167224028)/(z + 1) + 771.32342877765313/(z + 2) +
      (-176.61502916214059)/(z + 3) + 12.507343278686905/(z + 4) + (-0.13857109526572012)/(z + 5)); 
    //+ 9.9843695780195716e-6/(z + 6) + 1.5056327351493116e-7/(z + 7));
    
    return y;
  }
}

vector<int> incheck(const vector<int> &a, const vector<int> &b){
  int lengtha = a.size();
  int lengthb = b.size();
  vector<int> ans(lengtha);
  for(int i = 0; i < lengtha; i++){
    ans(i) = 0;
    for(int j = 0; j < lengthb; j++){
      if(a(i) == b(j)){
        ans(i) = 1;
        break;
      }
    }
  }
  return ans;
}

int incheck_scalar(const int &a, const vector<int> &b){
  int lengthb = b.size();
  int ans = 0;
  for(int i = 0; i < lengthb; i++){
    if(a == b(i)){
      ans = 1;
      break;
    }
  }
  return ans;
}

template<class Type>
vector<Type> isNA(const vector<Type> &x){
  int len = x.size();
  vector<Type> ans(len);
  for(int i = 0; i < len; i++){
    if(R_IsNA(asDouble(x(i)))){
      ans(i) = 1;
    } else {
      ans(i) = 0;
    }
  }
  return ans;
}

template<class Type>
void trans(Type *x, const int &link){
  
  if(link == 2){
    *x = exp(*x);
  }
  if(link == 3){
    *x = exp(*x)/(1 + exp(*x));
  }
  
}



//look up series functions------------------------------------------------------------------


//here n_i is n_IDs_for_datafull
//make sure this function is only used to locate the observations in the sessions with detections,
//for example, if there is no detection in session 2, than we should not use this function to
//locate any observations in session 2
int lookup_data_full(const int &is_ani, const int &s, const int &a, const int &i, const int &t,
                     const vector<int> &n_a, const vector<int> &n_i,
                     const vector<int> &n_t, const vector<int> &n_i_each_a){
  int ans = 0;
  if(is_ani == 1){
    //this is index to locate the correct element in n_i_each_a
    int i_a_index = 0;
    if(s > 1){
      for(int s_index = 1; s_index < s; s_index++){
        ans = ans + n_i(s_index - 1) * n_t(s_index - 1);
        i_a_index = i_a_index + n_a(s_index - 1);
      }
    }

    if(a > 1){
      for(int a_index = 1; a_index < a; a_index++){
        ans = ans + n_i_each_a(i_a_index) * n_t(s - 1);
        i_a_index++;
      }
    }

  } else {
    if(s > 1){
      for(int s_index = 1; s_index < s; s_index++){
        ans = ans + n_i(s_index - 1) * n_t(s_index - 1);
      }
    }

  }
  
  ans = ans + (i - 1) * n_t(s - 1) + t - 1;
  return ans;
  
}


int lookup_data_dist_theta(const int &s, const int &t, const int &m,
                           const vector<int> &n_t, const vector<int> &n_m){
  int ans = 0;
  if(s > 1){
    for(int s_index = 1; s_index < s; s_index++){
      ans = ans + n_t(s_index - 1) * n_m(s_index - 1);
    }
  }
  ans = ans + (m - 1) * n_t(s - 1) + t - 1 ;
  
  return ans;
}


int lookup_data_mask(const int &s, const int &m, const vector<int> &n_m){
  int ans = 0;
  if(s > 1){
    for(int s_index = 1; s_index < s; s_index++){
      ans = ans + n_m(s_index - 1);
    }
  }
  ans = ans + m - 1;
  
  return ans;
}

//here n_i is n_IDs

int lookup_data_IDmask(const int &is_ani, const int &s, const int &a, const int &i, const int &m,
                       const vector<int> &n_a, const vector<int> &n_i, 
                       const vector<int> &n_i_each_a, const vector<int> &n_m){
  int ans = 0;
  if(is_ani == 1){
    //this is index to locate the correct element in n_i_each_a
    int i_a_index = 0;

    if(s > 1){
      for(int s_index = 1; s_index < s; s_index++){
        ans = ans + n_i[s_index - 1] * n_m[s_index - 1];
        i_a_index = i_a_index + n_a[s_index - 1];
      }
    }

    if(a > 1){
      for(int a_index = 1; a_index < a; a_index++){
        ans = ans + n_i_each_a[i_a_index] * n_m[s - 1];
        i_a_index++;
      }
    }

    ans = ans + (m - 1) * n_i_each_a[i_a_index] + i - 1;

  } else {

    if(s > 1){
      for(int s_index = 1; s_index < s; s_index++){
        ans = ans + n_i[s_index - 1] * n_m[s_index - 1];
      }
    }
    ans = ans + (i - 1) * n_m[s - 1] + m - 1;
  }
  
  return ans;
}

int lookup_bincapt_uid(const int &s, const int &uid, const int &t, const vector<int> &n_uid_session, 
                       const vector<int> &n_t){
  
  int ans = 0;
  if(s > 1){
    for(int s_index = 1; s_index < s; s_index++){
      ans = ans + n_uid_session(s_index - 1) * n_t(s_index - 1);
    }
  }
  
  ans = ans + (uid - 1) * n_t(s - 1) + t - 1;
  
  return ans;
}

//here n_i is n_IDs as well
int lookup_n_detection(const int &is_ani, const int &s, const int &a, const int &i, const vector<int> &n_a,
                       const vector<int> &n_i, const vector<int> &n_i_each_a){
  int ans = 0;
  if(is_ani == 1){
    //this is index to locate the correct element in n_i_each_a
    int i_a_index = 0;
    if(s > 1){
      for(int s_index = 1; s_index < s; s_index++){
        ans = ans + n_i[s_index - 1];
        i_a_index = i_a_index + n_a[s_index - 1];
      }
    }

    if(a > 1){
      for(int a_index = 1; a_index < a; a_index++){
        ans = ans + n_i_each_a[i_a_index];
        i_a_index++;
      }
    }
    
  } else {
    if(s > 1){
      for(int s_index = 1; s_index < s; s_index++){
        ans = ans + n_i[s_index - 1];
      }
    }
    
  }

  ans = ans + i - 1;
  return ans;
}

//this one can be used to locate the nmuber of detections and the number of ids
//for each uid according to the input vectors respectively
int lookup_n_det_and_id_uid(const int &s, const int &uid, const vector<int> &n_uid_session){
  int ans = 0;
  
  if(s > 1){
    for(int s_index = 1; s_index < s; s_index++){
      ans = ans + n_uid_session(s_index - 1);
    }
  }
  
  ans = ans + uid - 1;
  
  return ans;
}

int look_up_u(const int &s, const int &i, const vector<int> &n_i){
  int ans = 0;
  if(s > 1){
    for(int s_index = 1; s_index < s; s_index++){
      ans += n_i(s_index - 1);
    }
  }
  
  ans += i - 1;
  return ans;
}

//this is the function to look up the ids based on session and uid from "u_id_match"
//this function directly return IDs inside of index in u_id_match
vector<int> look_up_ids_uid(const int &s, const int &uid, const vector<int> &n_i, 
                            const vector<int> n_ids_uid, const vector<int> &n_uid_session, const vector<int> &u_id_match){
  int index = 0;
  
  if(s > 1){
    for(int s_index = 1; s_index < s; s_index++){
      index += n_i(s_index - 1);
    }
  }
  
  int index_n_id_each_uid = lookup_n_det_and_id_uid(s, 1, n_uid_session);
  
  if(uid > 1){
    for(int uid_index = 1; uid_index < uid; uid_index++){
      index += n_ids_uid[index_n_id_each_uid];
      index_n_id_each_uid++;
    }
  }
  
  int n_ids = n_ids_uid[index_n_id_each_uid];
  
  vector<int> ans(n_ids);
  
  for(int i = 0; i < n_ids; i++){
    ans(i) = u_id_match(index);
    index++;
  }
  
  return ans;
}

//this is the function to look up the indices of traps based on session and uid from "index_traps_uid"
vector<int> look_up_traps_uid(const int &s, const int &uid, const int &index_n_det_uid, const int &n_det,
                              const vector<int> n_det_uid, const vector<int> &index_traps_uid){
  
  vector<int> ans(n_det);
  int index = 0;
  
  if(index_n_det_uid > 0){
    for(int j = 0; j < index_n_det_uid; j++){
      index += n_det_uid[j];
    }
  }
  
  for(int k = 0; k < n_det; k++){
    ans(k) = index_traps_uid(index);
    index++;
  }
  
  return ans;
}

//this is the function to lookup the index from capt_bin_uid
int lookup_data_uid(const int &s, const int &uid, const int &t,
                    const vector<int> &n_uid, const vector<int> &n_t){
  int ans = 0;
  
  if(s > 1){
    for(int s_index = 1; s_index < s; s_index++){
      ans = ans + n_uid(s_index - 1) * n_t(s_index - 1);
    }
  }
  ans = ans + (uid - 1) * n_t(s - 1) + t - 1;
  
  
  return ans;
  
}

//this is only used for animal_ID included model to look up the index for n_calls_each_animal
int lookup_n_call_ani(const int &s, const int &a, const vector<int> &n_a){
  int ans = 0;
  if(s > 1){
    for(int s_index = 1; s_index < s; s_index++){
      ans = ans + n_a(s_index - 1);
    }
  }

  ans = ans + a - 1;
  return ans;

}

//end of look up functions-------------------------------------------------------------



//detect functions---------------------------------------------------------------------

//hn
template<class Type>
Type det_hn(const Type &dx, const vector<Type> &param){
  Type ans = 0.0;
  Type g0 = param(0);
  Type sigma = param(1);
  ans = g0 * exp(-1 * pow(dx, 2) / (2 * pow(sigma, 2)));
  return ans;
}

//hhn
template<class Type>
Type det_hhn(const Type &dx, const vector<Type> &param){
  Type ans = 0.0;
  Type lambda0 = param(0);
  Type sigma = param(1);
  Type lambda_d = lambda0 * exp(pow(dx, 2) * (-0.5) / pow(sigma, 2));
  ans = 1 - exp((-1) * lambda_d);
  return ans;
}

//hr
template<class Type>
Type det_hr(const Type &dx, const vector<Type> &param){
  Type ans = 0.0;
  Type g0 = param(0);
  Type sigma = param(1);
  Type z = param(2);
  ans = g0 * (1 - exp((-1) * pow(dx / sigma, (-z))));
  return ans;
}

//lth
template<class Type>
Type det_lth(const Type &dx, const vector<Type> &param){
  Type ans = 0.0;
  Type shape_1 = param(0);
  Type shape_2 = param(1);
  Type scale = param(2);
  ans = 0.5 - 0.5 * erf(shape_1 - exp(shape_2 - scale * dx));
  return ans;
}

//th
template<class Type>
Type det_th(const Type &dx, const vector<Type> &param){
  Type ans = 0.0;
  Type shape = param(0);
  Type scale = param(1);
  ans = 0.5 - 0.5 * erf(dx / scale - shape);
  return ans;
}

//ss
//in ss, we do something different, here we use
//the same input, but output is back-transformed essx
//instead of the probability be detected
template<class Type>
Type essx_ss_identical(const Type &dx, const vector<Type> &param){
  Type b0_ss = param(0);
  Type b1_ss = param(1);
  Type essx = b0_ss - b1_ss * dx;
  return essx;
}

template<class Type>
Type essx_ss_log(const Type &dx, const vector<Type> &param){
  Type b0_ss = param(0);
  Type b1_ss = param(1);
  Type essx = exp(b0_ss - b1_ss * dx);
  return essx;
}

template<class Type>
Type essx_ss_spherical(const Type &dx, const vector<Type> &param){
  Type b0_ss = param(0);
  Type b1_ss = param(1);
  Type essx = 0.0;
  if(dx > 1){
    essx += b0_ss - 20 * log10(dx) - b1_ss * (dx - 1);
  } else {
    essx += b0_ss;
  }
  return essx;
}

//ss_dir
template<class Type>
Type essx_ss_dir_identical(const Type &dx, const vector<Type> &param){
  Type b0_ss = param(0);
  Type b1_ss = param(1);
  Type b2_ss = param(2);
  Type theta = param(3);
  Type essx = b0_ss - (b1_ss - b2_ss * (cos(theta) - 1)) * dx;
  return essx;
}

template<class Type>
Type essx_ss_dir_log(const Type &dx, const vector<Type> &param){
  Type b0_ss = param(0);
  Type b1_ss = param(1);
  Type b2_ss = param(2);
  Type theta = param(3);
  Type essx = exp(b0_ss - (b1_ss - b2_ss * (cos(theta) - 1)) * dx);
  return essx;
}

template<class Type>
Type essx_ss_dir_spherical(const Type &dx, const vector<Type> &param){
  Type b0_ss = param(0);
  Type b1_ss = param(1);
  Type b2_ss = param(2);
  Type theta = param(3);
  Type essx = 0.0;
  if(dx > 1){
    essx += b0_ss;
  } else {
    essx += b0_ss - 20 * log10(dx) - (b1_ss - b2_ss * (cos(theta) - 1)) * (dx - 1);
  }
  return essx;
}

//ss_het
template<class Type>
Type essx_ss_het(const Type &dx, const vector<Type> &param){
  Type b0_ss = param(0);
  Type b1_ss = param(1);
  Type u = param(2);
  Type essx = b0_ss - b1_ss * dx + u;
  return essx;
}

//end of detection functions------------------------------------------------------------------------

template<class Type>
Type objective_function<Type>::operator() ()
{
  Type nll = Type(0.0);
  Type *pointer_nll = &nll;
  Type Inf = std::numeric_limits<double>::infinity();
  Type *pointer_Inf = &Inf;
  Type pi = 3.14159265358979323846;
  DATA_INTEGER(n_sessions);
  DATA_IVECTOR(n_animals);
  DATA_IVECTOR(n_IDs);
  DATA_IVECTOR(n_IDs_for_datafull); 
  DATA_IVECTOR(n_traps);
  DATA_IVECTOR(n_masks);
  DATA_IVECTOR(n_detection);
  DATA_IVECTOR(n_detection_uid);
  DATA_IVECTOR(n_calls_each_animal);
  DATA_IVECTOR(n_uid_session);
  DATA_IVECTOR(n_ids_each_uid);
  DATA_IVECTOR(index_traps_uid);
  
  
  DATA_INTEGER(nrow_data_full);
  DATA_INTEGER(nrow_data_mask);
  DATA_INTEGER(nrow_dx);
  DATA_INTEGER(nrow_id_mask);
  
  DATA_VECTOR(A);
  DATA_VECTOR(survey_length);
  DATA_SCALAR(sound_speed);
  DATA_SCALAR(cue_rates);
  
  DATA_INTEGER(het_method);
  DATA_INTEGER(n_dir_quadpoints);
  DATA_INTEGER(n_het_source_quadpoints);
  DATA_VECTOR(het_source_nodes);
  DATA_VECTOR(het_source_weights);
  DATA_SCALAR(cutoff);
  
  
  DATA_INTEGER(detfn_index);
  DATA_IVECTOR(param_og);
  DATA_IMATRIX(par_n_col);
  DATA_IVECTOR(par_link);
  DATA_INTEGER(ss_link);
  
  
  DATA_INTEGER(is_animalID);
  DATA_INTEGER(is_ss);
  DATA_INTEGER(is_ss_origin);
  DATA_INTEGER(is_ss_het);
  DATA_INTEGER(is_ss_dir);
  DATA_INTEGER(is_bearing);
  DATA_INTEGER(is_toa);
  DATA_INTEGER(is_dist);
  DATA_INTEGER(is_local);
  DATA_INTEGER(is_freqs);
  
  DATA_IVECTOR(u_id_match);
  DATA_VECTOR(capt_bin_uid);
  
  DATA_VECTOR(capt_bin);
  DATA_VECTOR(capt_bearing);
  DATA_VECTOR(capt_dist);
  DATA_VECTOR(capt_ss);
  DATA_VECTOR(capt_toa);
  DATA_VECTOR(dx);
  DATA_VECTOR(theta);
  DATA_VECTOR(index_local);
  DATA_VECTOR(toa_ssq);
  
  
  DATA_MATRIX(g0_DX);
  DATA_MATRIX(sigma_DX);
  DATA_MATRIX(lambda0_DX);
  DATA_MATRIX(z_DX);
  DATA_MATRIX(shape_1_DX);
  DATA_MATRIX(shape_2_DX);
  DATA_MATRIX(shape_DX);
  DATA_MATRIX(scale_DX);
  DATA_MATRIX(b0_ss_DX);
  DATA_MATRIX(b1_ss_DX);
  DATA_MATRIX(b2_ss_DX);
  DATA_MATRIX(sigma_ss_DX);
  DATA_MATRIX(kappa_DX);
  DATA_MATRIX(alpha_DX);
  DATA_MATRIX(sigma_toa_DX);
  //sigma_b0_ss is not extentable
  //DATA_MATRIX(sigma_b0_ss_DX);
  DATA_MATRIX(D_DX);
  DATA_MATRIX(mu_DX);
  
  
  DATA_MATRIX(g0_DX_mask);
  DATA_MATRIX(sigma_DX_mask);
  DATA_MATRIX(lambda0_DX_mask);
  DATA_MATRIX(z_DX_mask);
  DATA_MATRIX(shape_1_DX_mask);
  DATA_MATRIX(shape_2_DX_mask);
  DATA_MATRIX(shape_DX_mask);
  DATA_MATRIX(scale_DX_mask);
  DATA_MATRIX(b0_ss_DX_mask);
  DATA_MATRIX(b1_ss_DX_mask);
  DATA_MATRIX(b2_ss_DX_mask);
  DATA_MATRIX(sigma_ss_DX_mask);
  DATA_MATRIX(kappa_DX_mask);
  DATA_MATRIX(alpha_DX_mask);
  DATA_MATRIX(sigma_toa_DX_mask);
  //sigma_b0_ss is not extentable
  //DATA_MATRIX(sigma_b0_ss_DX_mask);
  DATA_MATRIX(D_DX_mask);
  DATA_MATRIX(mu_DX_mask);
  
  DATA_MATRIX(g0_bound);
  DATA_MATRIX(sigma_bound);
  DATA_MATRIX(lambda0_bound);
  DATA_MATRIX(z_bound);
  DATA_MATRIX(shape_1_bound);
  DATA_MATRIX(shape_2_bound);
  DATA_MATRIX(shape_bound);
  DATA_MATRIX(scale_bound);
  DATA_MATRIX(b0_ss_bound);
  DATA_MATRIX(b1_ss_bound);
  DATA_MATRIX(b2_ss_bound);
  DATA_MATRIX(sigma_ss_bound);
  DATA_MATRIX(kappa_bound);
  DATA_MATRIX(alpha_bound);
  DATA_MATRIX(sigma_toa_bound);
  DATA_MATRIX(sigma_b0_ss_bound);
  DATA_MATRIX(D_bound);
  DATA_MATRIX(mu_bound);
  
  std::cout << "is_ss_origin: " << is_ss_origin << std::endl;
  std::cout << "detfn_index: " << detfn_index << std::endl;
  std::cout << "ss.link: " << ss_link << std::endl;


  //define a pointer variable which point to a detect function
  //based on our input "detfn_index"
  //and also point out the number of parameters for that detfn
  
  Type (*detfn)(const Type &dx, const vector<Type> &param);
  int n_detfn_param;
  if(detfn_index == 1){
    detfn = det_hn;
    n_detfn_param = 2;
  } else if(detfn_index == 2){
    detfn = det_hhn;
    n_detfn_param = 2;
  } else if(detfn_index == 3){
    detfn = det_hr;
    n_detfn_param = 3;
  } else if(detfn_index == 5){
    detfn = det_lth;
    n_detfn_param = 3;
  } else if(detfn_index == 4){
    detfn = det_th;
    n_detfn_param = 2;
  } else if(detfn_index == 6){
    if(ss_link == 1){
      detfn = essx_ss_identical;
    } else if(ss_link == 2){
      detfn = essx_ss_log;
    } else if(ss_link == 4){
      detfn = essx_ss_spherical;
    }
    n_detfn_param = 2;
  } else if(detfn_index == 7){
    if(ss_link == 1){
      detfn = essx_ss_dir_identical;
    } else if(ss_link == 2){
      detfn = essx_ss_dir_log;
    } else if(ss_link == 4){
      detfn = essx_ss_dir_spherical;
    }
    n_detfn_param = 4;
  } else if(detfn_index == 8){
    detfn = essx_ss_het;
    n_detfn_param = 3;
  }
  //define the parameter vector which will be used later
  vector<Type> detfn_param(n_detfn_param);
  
  //declare all parameters and generate its corresponding
  //data vectors based on data_full and data_mask
  //and only report this parameter if it is used in the model
  //to avoid NaN in std, outside of this .cpp, the unused
  //parameters should be fixed in MakeADFun() by "map=list(...)"
  
  
  //in addition, very important point, although all parameters
  //are used here in their original names, such as g0, however,
  //all of them are transformed by their corresponding link function
  //in the likelihood, they need to be backtransformed
  
  int i;
  //g0
  PARAMETER_VECTOR(g0);
  
  for(i = 0; i < (par_n_col(0, 0) + par_n_col(0, 1)); i++){
    if(g0(i) < g0_bound(0, i) || g0(i) > g0_bound(1, i)) *pointer_nll += *pointer_Inf; //infinity is Inf or INFINITY?
  }
  
  vector<Type> g0_full = g0.head(par_n_col(0, 0));
  vector<Type> g0_mask = g0.tail(par_n_col(0, 1));
  vector<Type> g0_vec_full = g0_DX * g0_full;
  vector<Type> g0_vec_mask(nrow_data_mask);
  
  if(par_n_col(0, 1) == 0){
    g0_vec_mask.setZero();
  } else {
    g0_vec_mask = g0_DX_mask * g0_mask;
  }
  if(incheck_scalar(1, param_og) == 1) ADREPORT(g0);
  
  //sigma
  PARAMETER(sigma);
  
  //for(i = 0; i < (par_n_col(1, 0) + par_n_col(1, 1)); i++){
  //  if(sigma(i) < sigma_bound(0, i) || sigma(i) > sigma_bound(1, i)) *pointer_nll += *pointer_Inf;
  //}
  
  //vector<Type> sigma_full = sigma.head(par_n_col(1, 0));
  //vector<Type> sigma_mask = sigma.tail(par_n_col(1, 1));
  //vector<Type> sigma_vec_full = sigma_DX * sigma_full;
  //vector<Type> sigma_vec_mask(nrow_data_mask);
  
  //if(par_n_col(1, 1) == 0){
  //  sigma_vec_mask.setZero();
  //} else {
  //  sigma_vec_mask = sigma_DX_mask * sigma_mask;
  //}
  if(incheck_scalar(2, param_og) == 1) ADREPORT(sigma);
  
  //lambda0
  PARAMETER(lambda0);
  
  //for(i = 0; i < (par_n_col(2, 0) + par_n_col(2, 1)); i++){
  //  if(lambda0(i) < lambda0_bound(0, i) || lambda0(i) > lambda0_bound(1, i)) *pointer_nll += *pointer_Inf;
  //}
  
  //vector<Type> lambda0_full = lambda0.head(par_n_col(2, 0));
  //vector<Type> lambda0_mask = lambda0.tail(par_n_col(2, 1));
  //vector<Type> lambda0_vec_full = lambda0_DX * lambda0_full;
  //vector<Type> lambda0_vec_mask(nrow_data_mask);
  
  //if(par_n_col(2, 1) == 0){
  //  lambda0_vec_mask.setZero();
  //} else {
  //  lambda0_vec_mask = lambda0_DX_mask * lambda0_mask;
  //}
  if(incheck_scalar(3, param_og) == 1) ADREPORT(lambda0);
  
  //z
  PARAMETER_VECTOR(z);
  for(i = 0; i < (par_n_col(3, 0) + par_n_col(3, 1)); i++){
    if(z(i) < z_bound(0, i) || z(i) > z_bound(1, i)) *pointer_nll += *pointer_Inf;
  }
  
  vector<Type> z_full = z.head(par_n_col(3, 0));
  vector<Type> z_mask = z.tail(par_n_col(3, 1));
  vector<Type> z_vec_full = z_DX * z_full;
  vector<Type> z_vec_mask(nrow_data_mask);
  
  if(par_n_col(3, 1) == 0){
    z_vec_mask.setZero();
  } else {
    z_vec_mask = z_DX_mask * z_mask;
  }

  if(incheck_scalar(4, param_og) == 1) ADREPORT(z);
  
  //shape_1
  PARAMETER_VECTOR(shape_1);
  
  for(i = 0; i < (par_n_col(4, 0) + par_n_col(4, 1)); i++){
    if(shape_1(i) < shape_1_bound(0, i) || shape_1(i) > shape_1_bound(1, i)) *pointer_nll += *pointer_Inf;
  }
  
  vector<Type> shape_1_full = shape_1.head(par_n_col(4, 0));
  vector<Type> shape_1_mask = shape_1.tail(par_n_col(4, 1));
  vector<Type> shape_1_vec_full = shape_1_DX * shape_1_full;
  vector<Type> shape_1_vec_mask(nrow_data_mask);
  
  if(par_n_col(4, 1) == 0){
    shape_1_vec_mask.setZero();
  } else {
    shape_1_vec_mask = shape_1_DX_mask * shape_1_mask;
  }
  if(incheck_scalar(5, param_og) == 1) ADREPORT(shape_1);
  
  //shape_2
  PARAMETER_VECTOR(shape_2);
  
  for(i = 0; i < (par_n_col(5, 0) + par_n_col(5, 1)); i++){
    if(shape_2(i) < shape_2_bound(0, i) || shape_2(i) > shape_2_bound(1, i)) *pointer_nll += *pointer_Inf;
  }
  
  vector<Type> shape_2_full = shape_2.head(par_n_col(5, 0));
  vector<Type> shape_2_mask = shape_2.tail(par_n_col(5, 1));
  vector<Type> shape_2_vec_full = shape_2_DX * shape_2_full;
  vector<Type> shape_2_vec_mask(nrow_data_mask);
  
  if(par_n_col(5, 1) == 0){
    shape_2_vec_mask.setZero();
  } else {
    shape_2_vec_mask = shape_2_DX_mask * shape_2_mask;
  }
  if(incheck_scalar(6, param_og) == 1) ADREPORT(shape_2);
  
  //shape
  PARAMETER_VECTOR(shape);
  
  for(i = 0; i < (par_n_col(6, 0) + par_n_col(6, 1)); i++){
    if(shape(i) < shape_bound(0, i) || shape(i) > shape_bound(1, i)) *pointer_nll += *pointer_Inf;
  }
  
  vector<Type> shape_full = shape.head(par_n_col(6, 0));
  vector<Type> shape_mask = shape.tail(par_n_col(6, 1));
  vector<Type> shape_vec_full = shape_DX * shape_full;
  vector<Type> shape_vec_mask(nrow_data_mask);
  
  if(par_n_col(6, 1) == 0){
    shape_vec_mask.setZero();
  } else {
    shape_vec_mask = shape_DX_mask * shape_mask;
  }
  if(incheck_scalar(7, param_og) == 1) ADREPORT(shape);
  
  //scale
  PARAMETER_VECTOR(scale);
  
  for(i = 0; i < (par_n_col(7, 0) + par_n_col(7, 1)); i++){
    if(scale(i) < scale_bound(0, i) || scale(i) > scale_bound(1, i)) *pointer_nll += *pointer_Inf;
  }
  
  vector<Type> scale_full = scale.head(par_n_col(7, 0));
  vector<Type> scale_mask = scale.tail(par_n_col(7, 1));
  vector<Type> scale_vec_full = scale_DX * scale_full;
  vector<Type> scale_vec_mask(nrow_data_mask);
  
  if(par_n_col(7, 1) == 0){
    scale_vec_mask.setZero();
  } else {
    scale_vec_mask = scale_DX_mask * scale_mask;
  }
  if(incheck_scalar(8, param_og) == 1) ADREPORT(scale);
  
  //b0_ss
  PARAMETER_VECTOR(b0_ss);
  
  for(i = 0; i < (par_n_col(8, 0) + par_n_col(8, 1)); i++){
    if(b0_ss(i) < b0_ss_bound(0, i) || b0_ss(i) > b0_ss_bound(1, i)) *pointer_nll += *pointer_Inf;
  }
  
  vector<Type> b0_ss_full = b0_ss.head(par_n_col(8, 0));
  vector<Type> b0_ss_mask = b0_ss.tail(par_n_col(8, 1));
  vector<Type> b0_ss_vec_full = b0_ss_DX * b0_ss_full;
  vector<Type> b0_ss_vec_mask(nrow_data_mask);
  
  if(par_n_col(8, 1) == 0){
    b0_ss_vec_mask.setZero();
  } else {
    b0_ss_vec_mask = b0_ss_DX_mask * b0_ss_mask;
  }
  if(incheck_scalar(9, param_og) == 1) ADREPORT(b0_ss);
  
  //b1_ss
  PARAMETER_VECTOR(b1_ss);
  
  for(i = 0; i < (par_n_col(9, 0) + par_n_col(9, 1)); i++){
    if(b1_ss(i) < b1_ss_bound(0, i) || b1_ss(i) > b1_ss_bound(1, i)) *pointer_nll += *pointer_Inf;
  }
  
  vector<Type> b1_ss_full = b1_ss.head(par_n_col(9, 0));
  vector<Type> b1_ss_mask = b1_ss.tail(par_n_col(9, 1));
  vector<Type> b1_ss_vec_full = b1_ss_DX * b1_ss_full;
  vector<Type> b1_ss_vec_mask(nrow_data_mask);
  
  if(par_n_col(9, 1) == 0){
    b1_ss_vec_mask.setZero();
  } else {
    b1_ss_vec_mask = b1_ss_DX_mask * b1_ss_mask;
  }
  if(incheck_scalar(10, param_og) == 1) ADREPORT(b1_ss);
  
  //b2_ss
  PARAMETER_VECTOR(b2_ss);
  
  for(i = 0; i < (par_n_col(10, 0) + par_n_col(10, 1)); i++){
    if(b2_ss(i) < b2_ss_bound(0, i) || b2_ss(i) > b2_ss_bound(1, i)) *pointer_nll += *pointer_Inf;
  }
  
  vector<Type> b2_ss_full = b2_ss.head(par_n_col(10, 0));
  vector<Type> b2_ss_mask = b2_ss.tail(par_n_col(10, 1));
  vector<Type> b2_ss_vec_full = b2_ss_DX * b2_ss_full;
  vector<Type> b2_ss_vec_mask(nrow_data_mask);
  
  if(par_n_col(10, 1) == 0){
    b2_ss_vec_mask.setZero();
  } else {
    b2_ss_vec_mask = b2_ss_DX_mask * b2_ss_mask;
  }
  if(incheck_scalar(11, param_og) == 1) ADREPORT(b2_ss);
  
  
  //sigma_ss
  PARAMETER_VECTOR(sigma_ss);
  
  for(i = 0; i < (par_n_col(11, 0) + par_n_col(11, 1)); i++){
    if(sigma_ss(i) < sigma_ss_bound(0, i) || sigma_ss(i) > sigma_ss_bound(1, i)) *pointer_nll += *pointer_Inf;
  }
  
  vector<Type> sigma_ss_full = sigma_ss.head(par_n_col(11, 0));
  vector<Type> sigma_ss_mask = sigma_ss.tail(par_n_col(11, 1));
  vector<Type> sigma_ss_vec_full = sigma_ss_DX * sigma_ss_full;
  vector<Type> sigma_ss_vec_mask(nrow_data_mask);
  
  if(par_n_col(11, 1) == 0){
    sigma_ss_vec_mask.setZero();
  } else {
    sigma_ss_vec_mask = sigma_ss_DX_mask * sigma_ss_mask;
  }
  if(incheck_scalar(12, param_og) == 1) ADREPORT(sigma_ss);
  
  //kappa
  PARAMETER_VECTOR(kappa);
  
  for(i = 0; i < (par_n_col(12, 0) + par_n_col(12, 1)); i++){
    if(kappa(i) < kappa_bound(0, i) || kappa(i) > kappa_bound(1, i)) *pointer_nll += *pointer_Inf;
  }
  
  
  vector<Type> kappa_full = kappa.head(par_n_col(12, 0));
  vector<Type> kappa_mask = kappa.tail(par_n_col(12, 1));
  vector<Type> kappa_vec_full = kappa_DX * kappa_full;
  vector<Type> kappa_vec_mask(nrow_data_mask);
  
  if(par_n_col(12, 1) == 0){
    kappa_vec_mask.setZero();
  } else {
    kappa_vec_mask = kappa_DX_mask * kappa_mask;
  }
  if(incheck_scalar(13, param_og) == 1) ADREPORT(kappa);
  
  //alpha
  PARAMETER_VECTOR(alpha);
  
  for(i = 0; i < (par_n_col(13, 0) + par_n_col(13, 1)); i++){
    if(alpha(i) < alpha_bound(0, i) || alpha(i) > alpha_bound(1, i)) *pointer_nll += *pointer_Inf;
  }
  
  
  vector<Type> alpha_full = alpha.head(par_n_col(13, 0));
  vector<Type> alpha_mask = alpha.tail(par_n_col(13, 1));
  vector<Type> alpha_vec_full = alpha_DX * alpha_full;
  vector<Type> alpha_vec_mask(nrow_data_mask);
  
  if(par_n_col(13, 1) == 0){
    alpha_vec_mask.setZero();
  } else {
    alpha_vec_mask = alpha_DX_mask * alpha_mask;
  }
  if(incheck_scalar(14, param_og) == 1) ADREPORT(alpha);
  
  //sigma_toa
  PARAMETER_VECTOR(sigma_toa);
  
  for(i = 0; i < (par_n_col(14, 0) + par_n_col(14, 1)); i++){
    if(sigma_toa(i) < sigma_toa_bound(0, i) || sigma_toa(i) > sigma_toa_bound(1, i)) *pointer_nll += *pointer_Inf;
  }
  
  vector<Type> sigma_toa_full = sigma_toa.head(par_n_col(14, 0));
  vector<Type> sigma_toa_mask = sigma_toa.tail(par_n_col(14, 1));
  vector<Type> sigma_toa_vec_full = sigma_toa_DX * sigma_toa_full;
  vector<Type> sigma_toa_vec_mask(nrow_data_mask);
  
  if(par_n_col(14, 1) == 0){
    sigma_toa_vec_mask.setZero();
  } else {
    sigma_toa_vec_mask = sigma_toa_DX_mask * sigma_toa_mask;
  }
  if(incheck_scalar(15, param_og) == 1) ADREPORT(sigma_toa);
  
  //sigma_b0_ss, this is not extentable, so just declare it as a scalar
  PARAMETER(sigma_b0_ss);
  
  if(sigma_b0_ss < sigma_b0_ss_bound(0, 0) || sigma_b0_ss > sigma_b0_ss_bound(1, 0)) *pointer_nll += *pointer_Inf;
  
  
  if(incheck_scalar(16, param_og) == 1) ADREPORT(sigma_b0_ss);
  
  //settle with "D" as it will be used regardless type of model
  PARAMETER_VECTOR(D);
  
  for(i = 0; i < (par_n_col(16, 0) + par_n_col(16, 1)); i++){
    if(D(i) < D_bound(0, i) || D(i) > D_bound(1, i)) *pointer_nll += *pointer_Inf;
  }
  
  ADREPORT(D);
  
  vector<Type> D_full = D.head(par_n_col(16, 0));
  vector<Type> D_mask = D.tail(par_n_col(16, 1));
  vector<Type> D_vec_full = D_DX * D_full;
  vector<Type> D_vec_mask(nrow_data_mask);
  
  if(par_n_col(16, 1) == 0){
    D_vec_mask.setZero();
  } else {
    D_vec_mask = D_DX_mask * D_mask;
  }
  

  //deal with 'mu', which is the cue.rate for animal_ID included model
  PARAMETER_VECTOR(mu);
  
  for(i = 0; i < (par_n_col(17, 0) + par_n_col(17, 1)); i++){
    if(mu(i) < mu_bound(0, i) || mu(i) > mu_bound(1, i)) *pointer_nll += *pointer_Inf;
  }
  
  
  vector<Type> mu_full = mu.head(par_n_col(17, 0));
  vector<Type> mu_mask = mu.tail(par_n_col(17, 1));
  vector<Type> mu_vec_full = mu_DX * mu_full;
  vector<Type> mu_vec_mask(nrow_data_mask);
  
  if(par_n_col(17, 1) == 0){
    mu_vec_mask.setZero();
  } else {
    mu_vec_mask = mu_DX_mask * mu_mask;
  }
  if(incheck_scalar(18, param_og) == 1) ADREPORT(mu);

  //finally declare the latent variable "u"
  PARAMETER_VECTOR(u);
  
  //declare variables for all backtransformed parameters 
  Type g0_tem;
  Type sigma_tem;
  Type lambda0_tem;
  Type z_tem;
  Type shape_1_tem;
  Type shape_2_tem;
  Type shape_tem;
  Type scale_tem;
  Type b0_ss_tem;
  Type b1_ss_tem;
  Type b2_ss_tem;
  Type sigma_ss_tem;
  Type kappa_tem;
  Type alpha_tem;
  Type sigma_toa_tem;
  Type sigma_b0_ss_tem;
  
  Type *p_g0_tem = &g0_tem;
  Type *p_sigma_tem = &sigma_tem;
  Type *p_lambda0_tem = &lambda0_tem;
  Type *p_z_tem = &z_tem;
  Type *p_shape_1_tem = &shape_1_tem;
  Type *p_shape_2_tem = &shape_2_tem;
  Type *p_shape_tem = &shape_tem;
  Type *p_scale_tem = &scale_tem;
  Type *p_b0_ss_tem = &b0_ss_tem;
  Type *p_b1_ss_tem = &b1_ss_tem;
  Type *p_b2_ss_tem = &b2_ss_tem;
  Type *p_sigma_ss_tem = &sigma_ss_tem;
  Type *p_kappa_tem = &kappa_tem;
  Type *p_alpha_tem = &alpha_tem;
  Type *p_sigma_toa_tem = &sigma_toa_tem;
  Type *p_sigma_b0_ss_tem = &sigma_b0_ss_tem;
  
  //essx = E(ss|x), since it is "session-mask-trap"
  //level data, use data_dist_theta's index
  vector<Type> essx(nrow_dx);
  essx.setZero();
  
  int index_data_mask;
  int index_data_full_D;
  int index_data_full;
  int index_data_dist_theta;
  int index_data_IDmask;
  int index_zi;
  int index_ci;
  int s;
  int m;
  int t;
  int a;
  int n_m;
  int n_i;
  int n_t;
  int n_uid;
  int n_a;
  int ci;
  
  int id;
  int index_data_uid;
  int index_data_uid_inital;
  int index_n_det_and_id_uid;
  int n_det;
  int index_data_full_initial;

  //essentially this is "Type *pointer_g0_vec_full;" below are the same
  Type *p_g0_full;
  Type *p_g0_mask;
  //Type *p_sigma_full;
  //Type *p_sigma_mask;
  //Type *p_lambda0_full;
  //Type *p_lambda0_mask;
  Type *p_z_full;
  Type *p_z_mask;
  Type *p_shape_1_full;
  Type *p_shape_1_mask;
  Type *p_shape_2_full;
  Type *p_shape_2_mask;
  Type *p_shape_full;
  Type *p_shape_mask;
  Type *p_scale_full;
  Type *p_scale_mask;
  Type *p_b0_ss_full;
  Type *p_b0_ss_mask;
  Type *p_b1_ss_full;
  Type *p_b1_ss_mask;
  Type *p_b2_ss_full;
  Type *p_b2_ss_mask;
  Type *p_sigma_ss_full;
  Type *p_sigma_ss_mask;
  Type *p_kappa_full;
  Type *p_kappa_mask;
  Type *p_alpha_full;
  Type *p_alpha_mask;
  Type *p_sigma_toa_full;
  Type *p_sigma_toa_mask;
  Type *p_D_full;
  Type *p_D_mask;
  Type *p_mu_full;
  Type *p_mu_mask;
  
  Type *p_dx;
  Type *p_theta;
  Type *p_essx;
  Type *p_capt_bin;
  Type *p_capt_bearing;
  Type *p_capt_dist;
  Type *p_capt_ss;
  Type *p_capt_toa;
  Type *p_toa_ssq;
  Type *p_index_local;
  
  
  Type area_unit;
  Type servey_len;
  Type lambda_theta;
  Type l_i;
  Type l_w;
  
  Type fy_toa_log;
  Type fy_bear_log;
  Type fy_dist_log;
  Type fy_ss_log;
  
  vector<Type> esa(n_sessions);
  
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //begin of likelihood function

  for(s = 1; s <= n_sessions; s++){
    n_m = n_masks(s - 1);
    n_t = n_traps(s - 1);
    n_a = n_animals(s - 1);
    area_unit = A(s - 1);
    servey_len = survey_length(s - 1);
    n_uid = n_uid_session(s - 1);
	
	//if(s == 105){
	//	std::cout << "session " << s << ", check point1, nll: " << *pointer_nll << std::endl;
	//}
	
    //firstly calculate lambda(theta), the rate of Poisson distribution
    
    //currently we don't have ID - level extension, so p_dot is a vector
    //with length of n_mask and p_k is a matrix with n_mask * n_trap
    //in the future, if we add ID - level extension, both of them should
    //add one dimension
    
    //and incidentally calculate all of the p_k and p_dot
    lambda_theta = Type(0.0);

    //one p_dot for each mask point, the prob of a call be detected by
    //at least one detector
    vector<Type> p_dot(n_m);

    //row index of p_k is for mask, colnum index is for trap/detector
    matrix<Type> p_k(n_m, n_t);

    //only used for animal_ID included model, the prob of an animal
    //has at least one call be detected by at least one detector
    vector<Type> p_m(n_m);
    
    vector<Type> D_tem(n_m);
    Type *p_D_tem = &D_tem[0];
    vector<Type> mu_tem(n_m);
    Type *p_mu_tem = &mu_tem[0];
    
    index_data_mask = lookup_data_mask(s, 1, n_masks);
    
    //"D" is special, not related to trap nor ID, but may relate to session
    //so we set it to the first row of each session

    index_data_full_D = lookup_data_full(is_animalID, s, 1, 1, 1, n_animals, 
                                         n_IDs_for_datafull, n_traps, n_calls_each_animal);
    //wo do not have any "ID-level" parameter extension, so this initial index
    //is the same with the initial index for "D", but we set another variable
    //for this index because this one will loop for the traps, and the one for "D"
    //will not. In the future if we enable the 'animal_ID' level extension, this one
    //will be moved to the loop of animal_ID, or even furthur to the 'ID' if 'ID' level
    //extension is available
    index_data_full = lookup_data_full(is_animalID, s, 1, 1, 1, n_animals, 
                                        n_IDs_for_datafull, n_traps, n_calls_each_animal);
    
    index_data_dist_theta = lookup_data_dist_theta(s, 1, 1, n_traps, n_masks);
    
    //both 'D' and 'mu' are not 'trap-level' extendable, so put them here temporarily
    //in the future, if any 'ID-level' extension available, this may be changed.
    p_D_full = &D_vec_full[index_data_full_D];
    p_mu_full = &mu_vec_full[index_data_full_D];

    p_D_mask = &D_vec_mask[index_data_mask];
    p_g0_mask = &g0_vec_mask[index_data_mask];
    
    //p_sigma_mask = &sigma_vec_mask[index_data_mask];
    
    //p_lambda0_mask = &lambda0_vec_mask[index_data_mask];
    
    p_z_mask = &z_vec_mask[index_data_mask];
    
    p_shape_1_mask = &shape_1_vec_mask[index_data_mask];
    
    p_shape_2_mask = &shape_2_vec_mask[index_data_mask];
    
    p_shape_mask = &shape_vec_mask[index_data_mask];
    
    p_scale_mask = &scale_vec_mask[index_data_mask];
    
    p_b0_ss_mask = &b0_ss_vec_mask[index_data_mask];
    
    p_b1_ss_mask = &b1_ss_vec_mask[index_data_mask];
    
    p_b2_ss_mask = &b2_ss_vec_mask[index_data_mask];
    
    p_sigma_ss_mask = &sigma_ss_vec_mask[index_data_mask];

    p_mu_mask = &mu_vec_mask[index_data_mask];
    
    p_dx = &dx[index_data_dist_theta];
    p_theta = &theta[index_data_dist_theta];
    p_essx = &essx[index_data_dist_theta];
    
    esa(s - 1) = Type(0.0);
    
    for(m = 1; m <= n_m; m++){
      *p_D_tem = *p_D_full + *p_D_mask;
      trans(p_D_tem, par_link(16));
      
      
      p_dot(m - 1) = Type(1.0);
      
      
      p_g0_full = &g0_vec_full[index_data_full];
      
      //p_sigma_full = &sigma_vec_full[index_data_full];
      
      //p_lambda0_full = &lambda0_vec_full[index_data_full];
      
      p_z_full = &z_vec_full[index_data_full];
      
      p_shape_1_full = &shape_1_vec_full[index_data_full];
      
      p_shape_2_full = &shape_2_vec_full[index_data_full];
      
      p_shape_full = &shape_vec_full[index_data_full];
      
      p_scale_full = &scale_vec_full[index_data_full];
      
      p_b0_ss_full = &b0_ss_vec_full[index_data_full];
      
      p_b1_ss_full = &b1_ss_vec_full[index_data_full];
      
      p_b2_ss_full = &b2_ss_vec_full[index_data_full];
      
      p_sigma_ss_full = &sigma_ss_vec_full[index_data_full];
      
      
      for(t = 1; t <= n_t; t++){
        
        if(detfn_index == 1){
          *p_g0_tem = *p_g0_full + *p_g0_mask;
          trans(p_g0_tem, par_link(0));
          *p_sigma_tem = sigma;
          trans(p_sigma_tem, par_link(1));
          
          p_g0_full++;
          //p_sigma_full++;
          
          detfn_param(0) = g0_tem;
          detfn_param(1) = sigma_tem;
          
        } else if(detfn_index == 2){
          *p_lambda0_tem = lambda0;
          trans(p_lambda0_tem, par_link(2));
          *p_sigma_tem = sigma;
          trans(p_sigma_tem, par_link(1));
          
          //p_lambda0_full++;
          //p_sigma_full++;
          
          detfn_param(0) = lambda0_tem;
          detfn_param(1) = sigma_tem;
          
        } else if(detfn_index == 3){
          *p_g0_tem = *p_g0_full + *p_g0_mask;
          trans(p_g0_tem, par_link(0));
          *p_sigma_tem = sigma;
          trans(p_sigma_tem, par_link(1));
          *p_z_tem = *p_z_full + *p_z_mask;
          trans(p_z_tem, par_link(3));
          
          p_g0_full++;
          //p_sigma_full++;
          p_z_full++;
          
          
          detfn_param(0) = g0_tem;
          detfn_param(1) = sigma_tem;	
          detfn_param(2) = z_tem;	
          
        } else if(detfn_index == 5){
          *p_shape_1_tem = *p_shape_1_full + *p_shape_1_mask;
          trans(p_shape_1_tem, par_link(4));
          *p_shape_2_tem = *p_shape_2_full + *p_shape_2_mask;
          trans(p_shape_2_tem, par_link(5));
          *p_scale_tem = *p_scale_full + *p_scale_mask;
          trans(p_scale_tem, par_link(7));
          
          p_shape_1_full++;
          p_shape_2_full++;
          p_scale_full++;
          
          detfn_param(0) = shape_1_tem;
          detfn_param(1) = shape_2_tem;	
          detfn_param(2) = scale_tem;
          
        } else if(detfn_index == 4){
          *p_shape_tem = *p_shape_full + *p_shape_mask;
          trans(p_shape_tem, par_link(6));
          *p_scale_tem = *p_scale_full + *p_scale_mask;
          trans(p_scale_tem, par_link(7));
          
          p_shape_full++;
          p_scale_full++;
          
          
          detfn_param(0) = shape_tem;
          detfn_param(1) = scale_tem;	
        } else if(detfn_index == 6){
          *p_b0_ss_tem = *p_b0_ss_full + *p_b0_ss_mask;
          trans(p_b0_ss_tem, par_link(8));
          *p_b1_ss_tem = *p_b1_ss_full + *p_b1_ss_mask;
          trans(p_b1_ss_tem, par_link(9));
          
          p_b0_ss_full++;
          p_b1_ss_full++;
          
          
          detfn_param(0) = b0_ss_tem;
          detfn_param(1) = b1_ss_tem;
          
        } else if(detfn_index == 7){
          *p_b0_ss_tem = *p_b0_ss_full + *p_b0_ss_mask;
          trans(p_b0_ss_tem, par_link(8));
          *p_b1_ss_tem = *p_b1_ss_full + *p_b1_ss_mask;
          trans(p_b1_ss_tem, par_link(9));
          *p_b2_ss_tem = *p_b2_ss_full + *p_b2_ss_mask;
          trans(p_b2_ss_tem, par_link(10));
          
          p_b0_ss_full++;
          p_b1_ss_full++;
          p_b2_ss_full++;
          
          
          detfn_param(0) = b0_ss_tem;
          detfn_param(1) = b1_ss_tem;
          detfn_param(2) = b2_ss_tem;
          detfn_param(3) = *p_theta;
          
          p_theta++;
        } 
        //het is detfn_index == 8, not sure how to do it yet
        
        
        
        if(is_ss == 0){
          p_k(m - 1, t - 1) = (*detfn)(*p_dx, detfn_param);
        } else if (is_ss_origin == 1){
          *p_sigma_ss_tem = *p_sigma_ss_full + *p_sigma_ss_mask;
          trans(p_sigma_ss_tem, par_link(11));
          
          p_sigma_ss_full++;
          
          *p_essx = (*detfn)(*p_dx, detfn_param);
          p_k(m - 1, t - 1) = 1 - pnorm((cutoff - *p_essx) / sigma_ss_tem);

          p_essx++;
        }



        p_dx++;
        
        p_dot(m - 1) *= 1 - p_k(m - 1, t - 1);
        //end for trap t
      }
      
      p_dot(m - 1) = 1 - p_dot(m - 1);

      if(is_animalID == 0){
        esa(s - 1) += p_dot(m - 1);
        lambda_theta += *p_D_tem * p_dot(m - 1);
      } else {
        *p_mu_tem = *p_mu_full + *p_mu_mask;
        trans(p_mu_tem, par_link(17));
        p_m(m - 1) = 1 - exp((-1) * (*p_mu_tem) * servey_len * p_dot(m - 1));
        esa(s - 1) += p_m(m - 1);
        lambda_theta += *p_D_tem * p_m(m - 1);
      }
      
      
      p_D_mask++;
      p_D_tem++;
      p_g0_mask++;
      //p_sigma_mask++;
      //p_lambda0_mask++;
      p_z_mask++;
      p_shape_1_mask++;
      p_shape_2_mask++;
      p_shape_mask++;
      p_scale_mask++;
      p_b0_ss_mask++;
      p_b1_ss_mask++;
      p_b2_ss_mask++;
      p_sigma_ss_mask++;
      p_mu_mask++;
      p_mu_tem++;
      //end for mask m
    }
    
    //std::cout << "session " << s << ", check point 1, lambda_theta: " << lambda_theta << std::endl;

    esa(s - 1) *= area_unit;
    if(is_animalID == 0){
      lambda_theta *= area_unit * servey_len * cue_rates;
    } else {
      lambda_theta *= area_unit;
    }
    //std::cout << "session " << s << ", check point 2, lambda_theta: " << lambda_theta << std::endl;
    //end of lambda_theta calculation
    
    
    //canceled out original likelihood: nll -= dpois(Type(n_i), lambda_theta, true);
    *pointer_nll += lambda_theta;
	//if(s == 105){
	//	std::cout << "session "<< s <<", check point2, nll: " << *pointer_nll << std::endl;
	//}	
	
    //just like D, sigma_toa is not id nor trap extentable, so just take
    //the first value in this session
    p_sigma_toa_full = &sigma_toa_vec_full[index_data_full_D];
    
    if(is_animalID == 0){
      if(n_uid > 0){
        for(int uid = 1; uid <= n_uid; uid++){
          //ids: a vector contains all IDs under this uid
          vector<int> ids = look_up_ids_uid(s, uid, n_IDs, n_ids_each_uid, n_uid_session, u_id_match);
          
          //calculate and store the capture fw for this unique capture history
          vector<Type> fw(n_m);
          index_data_uid_inital = lookup_data_uid(s, uid, 1, n_uid_session, n_traps);
          if(is_local == 0){
            for(m = 1; m <= n_m; m++){
              index_data_uid = index_data_uid_inital;
              fw(m - 1) = Type(1.0);
              for(t = 1; t <= n_t; t++){
                if(is_ss == 0){
                  fw(m - 1) *= pow(p_k(m - 1, t - 1), capt_bin_uid(index_data_uid)) * 
                    pow((1 - p_k(m - 1, t - 1)), (1 - capt_bin_uid(index_data_uid)));
                } else if(is_ss_het == 0){
                  //pow(p_k(), capt_bin()) term could be cancelled out with fy_ss
                  fw(m - 1) *= pow((1 - p_k(m - 1, t - 1)), (1 -capt_bin_uid(index_data_uid)));
                }
                index_data_uid++;
              }
            }
          } else {
            //all ids under this uid have the same detection history, take the first one is enough
            //to locate the local masks
            id = ids(0);
            //the initial index in data_IDmask for this id
            index_data_IDmask = lookup_data_IDmask(is_animalID, s, 0, id, 1,
                                                  0, n_IDs, 0, n_masks);
            p_index_local = &index_local[index_data_IDmask];
            for(m = 1; m <= n_m; m++){
              if(*p_index_local == 1){
                index_data_uid = index_data_uid_inital;
                fw(m - 1) = Type(1.0);
                for(t = 1; t <= n_t; t++){
                  if(is_ss == 0){
                    fw(m - 1) *= pow(p_k(m - 1, t - 1), capt_bin_uid(index_data_uid)) * 
                      pow((1 - p_k(m - 1, t - 1)), (1 - capt_bin_uid(index_data_uid)));
                  } else if(is_ss_het == 0){
                    //pow(p_k(), capt_bin()) term could be cancelled out with fy_ss
                    fw(m - 1) *= pow((1 - p_k(m - 1, t - 1)), (1 - capt_bin_uid(index_data_uid)));
                  }
                  index_data_uid++;
                }
              }
              p_index_local++;
            }
          }
          
          
          //n_i: the number of ID under this uid
          //n_det: the number of detections under this uid
          //traps: a vector contains all activated detectors under this uid 
          index_n_det_and_id_uid = lookup_n_det_and_id_uid(s, uid, n_uid_session);
          n_i = n_ids_each_uid(index_n_det_and_id_uid);
          n_det = n_detection_uid(index_n_det_and_id_uid);
          
          vector<int> traps = look_up_traps_uid(s, uid, index_n_det_and_id_uid, n_det,
                                                n_detection_uid, index_traps_uid);
          
          //initial index for data_mask
          index_data_mask = lookup_data_mask(s, 1, n_masks);
          
          for(i = 0; i < n_i; i++){
            id = ids(i);
            l_i = Type(0.0);
            
            //the initial index in data_IDmask for this id
            index_data_IDmask = lookup_data_IDmask(is_animalID, s, 0, id, 1,
                                                  0, n_IDs, 0, n_masks);
            
            p_toa_ssq = &toa_ssq[index_data_IDmask];
            p_index_local = &index_local[index_data_IDmask];
            
            p_sigma_toa_mask = &sigma_toa_vec_mask[index_data_mask];
            p_kappa_mask = &kappa_vec_mask[index_data_mask];
            p_alpha_mask = &alpha_vec_mask[index_data_mask];
            p_sigma_ss_mask = &sigma_ss_vec_mask[index_data_mask];
            
            //reset p_D_tem
            p_D_tem = &D_tem[0];
            
            if(is_local == 0){
              for(m = 1; m <= n_m; m++){
                fy_toa_log = Type(0.0);
                fy_bear_log = Type(0.0);
                fy_dist_log = Type(0.0);
                fy_ss_log = Type(0.0);
                
                //toa
                if(is_toa == 1){
                  //"sigma_toa" is not trap extendable nor ID either, so take index_data_full_D as well
                  *p_sigma_toa_tem = *p_sigma_toa_full + *p_sigma_toa_mask;
                  trans(p_sigma_toa_tem, par_link(14));
                  
                  fy_toa_log += (1 - n_det) * log(sigma_toa_tem) + (-0.5) * (*p_toa_ssq) / pow(sigma_toa_tem, 2);
                  
                  p_sigma_toa_mask++;
                  p_toa_ssq++;
                }
                
                for(t = 0; t < n_det; t++){
                  index_data_full = lookup_data_full(0, s, 0, id, traps(t),
                                                    0, n_IDs_for_datafull, n_traps, 0);
                  
                  index_data_dist_theta = lookup_data_dist_theta(s, traps(t), m, n_traps, n_masks);
                  
                  //bearing
                  if(is_bearing == 1){
                    p_capt_bearing = &capt_bearing[index_data_full];
                    p_kappa_full = &kappa_vec_full[index_data_full];
                    p_theta = &theta[index_data_dist_theta];
                    
                    *p_kappa_tem = *p_kappa_full + *p_kappa_mask;
                    trans(p_kappa_tem, par_link(12));
                    
                    fy_bear_log += kappa_tem * cos(*p_capt_bearing - *p_theta) - 
                      log(bessi0(kappa_tem));
                    
                  }
                  
                  //dist
                  if(is_dist == 1){
                    p_capt_dist = &capt_dist[index_data_full];
                    p_alpha_full = &alpha_vec_full[index_data_full];
                    p_dx = &dx[index_data_dist_theta];
                    
                    *p_alpha_tem = *p_alpha_full + *p_alpha_mask;
                    trans(p_alpha_tem, par_link(13));
                    
                    fy_dist_log +=  ((-1) * (alpha_tem * (log(*p_dx) - log(alpha_tem)) + log(Gamma(alpha_tem))) + (alpha_tem - 1) * 
                      log(*p_capt_dist) - alpha_tem * (*p_capt_dist) / *p_dx);
                    
                  }
                  
                  //ss_origin
                  if(is_ss_origin == 1){
                    p_capt_ss = &capt_ss[index_data_full];
                    p_essx = &essx[index_data_dist_theta];
                    p_sigma_ss_full = &sigma_ss_vec_full[index_data_full];
                    
                    *p_sigma_ss_tem = *p_sigma_ss_full + *p_sigma_ss_mask;
                    trans(p_sigma_ss_tem, par_link(11));
                    
                    fy_ss_log += dnorm(*p_capt_ss, *p_essx, sigma_ss_tem, true);
                    
                  }
                  
                  //end of loop of traps
                }
                //for fx=f(x|n;theta)
                //canceled out p_dot & lambda from original likelihood: 
                //fx = D_tem[m-1] * p_dot(m - 1) / lambda_theta = D_tem[m-1];
                //so fx was directly integrated into l_i below
                
                //we sum up likelihood (not log-likelihood) of each mask 
                l_i += (*p_D_tem) * fw(m - 1) * exp(fy_toa_log + fy_bear_log + fy_dist_log + fy_ss_log);
                
                p_D_tem++;
                p_kappa_mask++;
                p_sigma_ss_mask++;
                p_alpha_mask++;
                
                //end of loop of masks
              }
            } else {
              for(m = 1; m <= n_m; m++){
                if(*p_index_local == 1){
                  fy_toa_log = Type(0.0);
                  fy_bear_log = Type(0.0);
                  fy_dist_log = Type(0.0);
                  fy_ss_log = Type(0.0);
                  
                  //toa
                  if(is_toa == 1){
                    //"sigma_toa" is not trap extendable nor ID either, so take index_data_full_D as well
                    *p_sigma_toa_tem = *p_sigma_toa_full + *p_sigma_toa_mask;
                    trans(p_sigma_toa_tem, par_link(14));
                    fy_toa_log += (1 - n_det) * log(sigma_toa_tem) + (-0.5) * (*p_toa_ssq) / pow(sigma_toa_tem, 2);
                  }
                  
                  for(t = 0; t < n_det; t++){
                    index_data_full = lookup_data_full(0, s, 0, id, traps(t),
                                                      0, n_IDs_for_datafull, n_traps, 0);
                    index_data_dist_theta = lookup_data_dist_theta(s, traps(t), m, n_traps, n_masks);
                    //bearing
                    if(is_bearing == 1){
                      p_capt_bearing = &capt_bearing[index_data_full];
                      p_kappa_full = &kappa_vec_full[index_data_full];
                      p_theta = &theta[index_data_dist_theta];
                      
                      *p_kappa_tem = *p_kappa_full + *p_kappa_mask;
                      trans(p_kappa_tem, par_link(12));
                      
                      fy_bear_log += kappa_tem * cos(*p_capt_bearing - *p_theta) - 
                        log(bessi0(kappa_tem));
                    }
                    
                    //dist
                    if(is_dist == 1){
                      p_capt_dist = &capt_dist[index_data_full];
                      p_alpha_full = &alpha_vec_full[index_data_full];
                      p_dx = &dx[index_data_dist_theta];
                      
                      *p_alpha_tem = *p_alpha_full + *p_alpha_mask;
                      trans(p_alpha_tem, par_link(13));
                      
                      fy_dist_log +=  ((-1) * (alpha_tem * (log(*p_dx) - log(alpha_tem)) + log(Gamma(alpha_tem))) + (alpha_tem - 1) * 
                        log(*p_capt_dist) - alpha_tem * (*p_capt_dist) / *p_dx);
                    }
                    
                    //ss_origin
                    if(is_ss_origin == 1){
                      p_capt_ss = &capt_ss[index_data_full];
                      p_essx = &essx[index_data_dist_theta];
                      p_sigma_ss_full = &sigma_ss_vec_full[index_data_full];
                      
                      *p_sigma_ss_tem = *p_sigma_ss_full + *p_sigma_ss_mask;
                      trans(p_sigma_ss_tem, par_link(11));
                      
                      fy_ss_log += dnorm(*p_capt_ss, *p_essx, sigma_ss_tem, true);
                    }
                    //end of loop of traps
                  }
                  //for fx=f(x|n;theta)
                  //canceled out p_dot & lambda from original likelihood: 
                  //fx = D_tem[m-1] * p_dot(m - 1) / lambda_theta = D_tem[m-1];
                  //so fx was directly integrated into l_i below
                  
                  //we sum up likelihood (not log-likelihood) of each mask 
                  l_i += (*p_D_tem) * fw(m - 1) * exp(fy_toa_log + fy_bear_log + fy_dist_log + fy_ss_log);
                }
                
                p_sigma_toa_mask++;
                p_toa_ssq++;
                p_D_tem++;
                p_kappa_mask++;
                p_sigma_ss_mask++;
                p_alpha_mask++;
                p_index_local++;
                //end of loop of masks
              }
            }
            *pointer_nll -= log(l_i * area_unit * servey_len * cue_rates);
            //end of id
          }
          //end of uid
        }
        //end of if(uid > 0)
      }
      //end of if(is_animalID == 0)
    } else {
      //begin of animal_ID model
      if(n_a > 0){

        //initial index for data_mask
        index_data_mask = lookup_data_mask(s, 1, n_masks);
        
        for(a = 1; a <= n_a; a++){
  
		  //likelihood of each animal
		  l_i = Type(0.0);

		  //the number of calls detected by this animal
		  index_ci = lookup_n_call_ani(s, a, n_animals);
		  ci = n_calls_each_animal[index_ci];
		  //for each call, how many detectors detected it
		  vector<int> Zi_vec(ci);
		  for(i = 1; i <= ci; i++){
			index_zi = lookup_n_detection(is_animalID, s, a, i,
			  n_animals, n_IDs, n_calls_each_animal);
			Zi_vec[i - 1] = n_detection[index_zi];
		  }
		  
		  //reset p_D_tem
		  p_D_tem = &D_tem[0];
		  //reset p_mu_tem
		  p_mu_tem = &mu_tem[0];
		  
		  //the initial index in data_full for this animal
		  index_data_full_initial = lookup_data_full(is_animalID, s, a, 1, 1,
			n_animals, n_IDs_for_datafull, n_traps, n_calls_each_animal);

		  //the initial index in data_IDmask for this animal
		  index_data_IDmask = lookup_data_IDmask(is_animalID, s, a, 1, 1,
			n_animals, n_IDs, n_calls_each_animal, n_masks);
		

		  p_toa_ssq = &toa_ssq[index_data_IDmask];
		  p_index_local = &index_local[index_data_IDmask];
		  
		  p_sigma_toa_mask = &sigma_toa_vec_mask[index_data_mask];
		  p_kappa_mask = &kappa_vec_mask[index_data_mask];
		  p_alpha_mask = &alpha_vec_mask[index_data_mask];
		  p_sigma_ss_mask = &sigma_ss_vec_mask[index_data_mask];

		  if(is_local == 0){
			for(m = 1; m <= n_m; m++){
			  index_data_full = index_data_full_initial;
			  //likelihood for capture history
			  l_w = Type(1.0);
			  //std::cout << "session: " << s << ", animal: " << a << ", mask: "<< m << " begins." << std::endl;
			  for(i = 1; i <= ci; i++){
				for(t = 1; t <= n_t; t++){
				  if(is_ss == 0){
					l_w *= pow(p_k(m - 1, t - 1), capt_bin[index_data_full]) *
					  pow((1 - p_k(m - 1, t - 1)), (1 - capt_bin[index_data_full]));
				  } else if(is_ss_het == 0){
					//pow(p_k(), capt_bin()) term could be cancelled out with fy_ss
					l_w *= pow((1 - p_k(m - 1, t - 1)), (1 - capt_bin[index_data_full]));
				  }
				  index_data_full++;
				}
			  }
			  
			  //likelihood for extra information
			  fy_toa_log = Type(0.0);
			  fy_bear_log = Type(0.0);
			  fy_dist_log = Type(0.0);
			  fy_ss_log = Type(0.0);

			  //reset the index_data_full to the first observation of this animal
			  index_data_full = index_data_full_initial;
			  p_capt_bearing = &capt_bearing[index_data_full];
			  p_kappa_full = &kappa_vec_full[index_data_full];
			  p_capt_dist = &capt_dist[index_data_full];
			  p_alpha_full = &alpha_vec_full[index_data_full];
			  p_capt_ss = &capt_ss[index_data_full];
			  p_sigma_ss_full = &sigma_ss_vec_full[index_data_full];


			  index_data_dist_theta = lookup_data_dist_theta(s, 1, m, n_traps, n_masks);
			  for(i = 1; i <= ci; i++){
				//std::cout << "session: " << s << ", animal: " << a << ", mask: "<< m << ", id: " << i << ", loop for extra info begins." << std::endl;
				n_det = Zi_vec[i - 1];
				p_theta = &theta[index_data_dist_theta];
				p_dx = &dx[index_data_dist_theta];
				p_essx = &essx[index_data_dist_theta];
				//std::cout << "check point 5" << std::endl;
				//toa
				if(is_toa == 1){
				  //"sigma_toa" is not trap extendable nor ID either, so take index_data_full_D as well
				  *p_sigma_toa_tem = *p_sigma_toa_full + *p_sigma_toa_mask;
				  trans(p_sigma_toa_tem, par_link(14));

				  fy_toa_log += (1 - n_det) * log(sigma_toa_tem) + (-0.5) * (*p_toa_ssq) / pow(sigma_toa_tem, 2);
				  
				  p_toa_ssq++;
				}

				for(t = 1; t <= n_t; t++){
				  //bearing
				  if(is_bearing == 1){
					if(*p_capt_bearing != 0){
					  *p_kappa_tem = *p_kappa_full + *p_kappa_mask;
					  trans(p_kappa_tem, par_link(12));
					  
					  fy_bear_log += kappa_tem * cos(*p_capt_bearing - *p_theta) - 
						log(bessi0(kappa_tem));
					}

					p_capt_bearing++;
					p_kappa_full++;
					p_theta++;
				  }

				  //dist
				  if(is_dist == 1){
					if(*p_capt_dist != 0){
					  *p_alpha_tem = *p_alpha_full + *p_alpha_mask;
					  trans(p_alpha_tem, par_link(13));

					  fy_dist_log +=  ((-1) * (alpha_tem * (log(*p_dx) - log(alpha_tem)) + log(Gamma(alpha_tem))) + (alpha_tem - 1) * 
						log(*p_capt_dist) - alpha_tem * (*p_capt_dist) / *p_dx);
					}

					p_capt_dist++;
					p_alpha_full++;
					p_dx++;
				  }

				  //ss_origin
				  if(is_ss_origin == 1){
					if(*p_capt_ss != 0){
					  *p_sigma_ss_tem = *p_sigma_ss_full + *p_sigma_ss_mask;
					  trans(p_sigma_ss_tem, par_link(11));
					  
					  //if(s == 105 && a == 1 && i == 1 && m < 50){
						//  std::cout << "mask: " << m  << ", trap" << t << ", ss: " << (*p_capt_ss) << std::endl;
						//  std::cout << "mask: " << m <<  ", trap" << t << ", essx: " << (*p_essx) << std::endl;
						//  std::cout << "mask: " << m <<  ", trap" << t << ", sigma_ss: " << sigma_ss_tem << std::endl;
					  //}
					  

					  fy_ss_log += dnorm(*p_capt_ss, *p_essx, sigma_ss_tem, true);
					}

					p_capt_ss++;
					p_sigma_ss_full++; 
					p_essx++;
				  }

				  //end of trap 't'
				}
				//end of call 'i'
			  }
			  
			  //if(s == 105 && a == 1 && m >900 & m<2000){
				//  std::cout << "until mask: " << m - 1 << ", l_i: " << l_i << std::endl;
				//  std::cout << "mask: " << m << ", D(m): " << (*p_D_tem) << std::endl;
				//  std::cout << "mask: " << m << ", ci: " << ci << std::endl;
				//  std::cout << "mask: " << m << ", mu: " << *p_mu_tem << std::endl;
				//  std::cout << "mask: " << m << ", l_w: " << l_w << std::endl;
				//  std::cout << "mask: " << m << ", p_dot[m - 1]: " << p_dot[m - 1] << std::endl;
				//  std::cout << "mask: " << m << ", fy_ss_log: " << fy_ss_log << std::endl;
			  //}
			  
			  
			  l_i += (*p_D_tem) * dpois(Type(ci), Type((*p_mu_tem) * servey_len)) * l_w * 
				exp((*p_mu_tem) * servey_len * (1 - p_dot[m - 1])) * exp(fy_toa_log + fy_bear_log + fy_dist_log + fy_ss_log);

			  p_D_tem++;
			  p_mu_tem++;
			  p_sigma_toa_mask++;
			  p_kappa_mask++;
			  p_alpha_mask++;
			  p_sigma_ss_mask++;
			  //end of mask 'm'
			}
			//end of if(is_local == 0)
		  } else {
			//begin of 'is_local == 1'

			//the initial index in data_IDmask for this animal
			index_data_IDmask = lookup_data_IDmask(is_animalID, s, a, 1, 1, 
			  n_animals, n_IDs, n_calls_each_animal, n_masks);
			p_index_local = &index_local[index_data_IDmask];

			for(m = 1; m <= n_m; m++){
			  if(*p_index_local == 1){
				index_data_full = index_data_full_initial;
				//likelihood for capture history
				l_w = Type(1.0);
				for(i = 1; i <= ci; i++){
				  for(t = 1; t <= n_t; t++){
					if(is_ss == 0){
					  l_w *= pow(p_k(m - 1, t - 1), capt_bin[index_data_full]) *
						pow((1 - p_k(m - 1, t - 1)), (1 - capt_bin[index_data_full]));
					} else if(is_ss_het == 0){
					  //pow(p_k(), capt_bin()) term could be cancelled out with fy_ss
					  l_w *= pow((1 - p_k(m - 1, t - 1)), (1 - capt_bin[index_data_full]));
					}


					index_data_full++;
				  }
				}

				//likelihood for extra information
				fy_toa_log = Type(0.0);
				fy_bear_log = Type(0.0);
				fy_dist_log = Type(0.0);
				fy_ss_log = Type(0.0);

				//reset the index_data_full to the first observation of this animal
				index_data_full = index_data_full_initial;
				p_capt_bearing = &capt_bearing[index_data_full];
				p_kappa_full = &kappa_vec_full[index_data_full];
				p_capt_dist = &capt_dist[index_data_full];
				p_alpha_full = &alpha_vec_full[index_data_full];
				p_capt_ss = &capt_ss[index_data_full];
				p_sigma_ss_full = &sigma_ss_vec_full[index_data_full];


				index_data_dist_theta = lookup_data_dist_theta(s, 1, m, n_traps, n_masks);

				for(i = 1; i <= ci; i++){
				  n_det = Zi_vec[i - 1];
				  p_theta = &theta[index_data_dist_theta];
				  p_dx = &dx[index_data_dist_theta];
				  p_essx = &essx[index_data_dist_theta];

				  //toa
				  if(is_toa == 1){
					//"sigma_toa" is not trap extendable nor ID either, so take index_data_full_D as well
					*p_sigma_toa_tem = *p_sigma_toa_full + *p_sigma_toa_mask;
					trans(p_sigma_toa_tem, par_link(14));
					
					fy_toa_log += (1 - n_det) * log(sigma_toa_tem) + (-0.5) * (*p_toa_ssq) / pow(sigma_toa_tem, 2);
					
					p_toa_ssq++;
				  }

				  for(t = 1; t <= n_t; t++){
					//bearing
					if(is_bearing == 1){
					  if(*p_capt_bearing != 0){
						*p_kappa_tem = *p_kappa_full + *p_kappa_mask;
						trans(p_kappa_tem, par_link(12));
						
						fy_bear_log += kappa_tem * cos(*p_capt_bearing - *p_theta) - 
						  log(bessi0(kappa_tem));
					  }

					  p_capt_bearing++;
					  p_kappa_full++;
					  p_theta++;
					}

					//dist
					if(is_dist == 1){
					  if(*p_capt_dist != 0){
						*p_alpha_tem = *p_alpha_full + *p_alpha_mask;
						trans(p_alpha_tem, par_link(13));

						fy_dist_log +=  ((-1) * (alpha_tem * (log(*p_dx) - log(alpha_tem)) + log(Gamma(alpha_tem))) + (alpha_tem - 1) * 
						  log(*p_capt_dist) - alpha_tem * (*p_capt_dist) / *p_dx);
					  }

					  p_capt_dist++;
					  p_alpha_full++;
					  p_dx++;
					}

					//ss_origin
					if(is_ss_origin == 1){
					  if(*p_capt_ss != 0){
						*p_sigma_ss_tem = *p_sigma_ss_full + *p_sigma_ss_mask;
						trans(p_sigma_ss_tem, par_link(11));

						fy_ss_log += dnorm(*p_capt_ss, *p_essx, sigma_ss_tem, true);
					  }

					  p_capt_ss++;
					  p_sigma_ss_full++; 
					  p_essx++;
					}

					//end of trap 't'
				  }
				  //end of call 'i'
				}

			  l_i += (*p_D_tem) * dpois(Type(ci), Type((*p_mu_tem) * servey_len)) * l_w * 
				exp((*p_mu_tem) * servey_len * (1 - p_dot[m - 1])) * exp(fy_toa_log + fy_bear_log + fy_dist_log + fy_ss_log);

				//end of if(*p_local_index == 1)

			  } else {
				p_toa_ssq += ci;
			  }

			  p_D_tem++;
			  p_mu_tem++;
			  
			  p_sigma_toa_mask++;
			  p_kappa_mask++;
			  p_alpha_mask++;
			  p_sigma_ss_mask++;

			  //because the sort order of data_IDmask for animal.model is
			  //session - animalID - mask - callID, and the local information is
			  //a session - animalID - mask level data, so each time, we need
			  //to add "ci" which is the number of calls of this animal to locate
			  //the next mask
			  p_index_local += ci;
			  //end of mask 'm'
			}

			//end of 'is_local == 1'
		  }

		  *pointer_nll -= log(l_i * area_unit);

		  //end of animal 'a'
        }
        //end of if(n_a > 0)
      }
      //end of animal_ID model
    }
	
	//if(s == 105){
	//	std::cout << "session "<< s <<", check point3, nll: " << *pointer_nll << std::endl;
	//}
	
	
	//if(s >= 100 && s < 500){
	//	std::cout << "session "<< s <<", check point3, nll: " << *pointer_nll << std::endl;
	//}
    //end for session s
  }
  ADREPORT(esa);
  return nll;
}

