#include "optimize.h"

#include <cassert>
#include <iostream>

using namespace Ipopt;

// constructor
GaitSystemID::GaitSystemID()
{}

//destructor
GaitSystemID::~GaitSystemID()
{}

// returns the size of the problem
bool GaitSystemID::get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
                                Index& nnz_h_lag, IndexStyleEnum& index_style)
{

  // q = 6 # controls
  // p = 12 # sensors
  // There are 9 DoF, 18 states
  // Free model parameters: n x 6 x 12 gains, n x 6 T*'s
  // Let's start with n = 10 (number of gain discretization points)
  // N ~= 36,000 : 100 hz over 6 minutes (about 3/4 of the 8 minute data from
  // each trial)
  // The number of free parameters will then be: 36,000 x 18 + 10 x 6 x 12 + 10
  // x 6 = 648780 (780 model parameter unknowns)
  n = 648780;

  // The constraints are (N - 1) x 18 = 647982
  // The constraints are the equations of motion evaluated for at every time
  // node expect the first one.
  m = 647982;

  // There should be only two per row for the xi and xi-1, then values for the
  // 780 model parameters.
  // num rows = num constraints
  // num cols = num free parameters
  // There are two non zeros per row per state + a nonzero for each free
  // parameter in the dynamic equations (i.e. parameter derivs are zero in the
  // kinematic equations)
  // (2 * 18) * 647982 + 780 * 647982 / 2
  nnz_jac_g = 276040332;

  // We will let IPOPT estimate the Hessian.
  // TODO : I'm not sure if this needs to be set or what.
  nnz_h_lag = 10;

  // use the C style indexing (0-based)
  index_style = TNLP::C_STYLE;

  return true;
}

// returns the variable bounds
bool GaitSystemID::get_bounds_info(Index n, Number* x_l, Number* x_u,
                                   Index m, Number* g_l, Number* g_u)
{
  // here, the n and m we gave IPOPT in get_nlp_info are passed back to us.
  // If desired, we could assert to make sure they are what we think they are.
  assert(n == 648780);
  assert(m == 647982);

  // Could set bounds for the states not to deviate too far from the data. That
  // would mean having a vector of +/- some bound around the measurements

  // the variables have lower bounds of 1
  for (Index i=0; i<4; i++) {
    x_l[i] = 1.0;
  }

  // the variables have upper bounds of 5
  for (Index i=0; i<4; i++) {
    x_u[i] = 5.0;
  }

  // The constraints must all equal zero. (or we could just bound them to be
  // small.
  for (Index i=0; i<m; i++){
      g_l[i] = g_u[i] = 0.0;
  }

  return true;
}

// returns the initial point for the problem
bool GaitSystemID::get_starting_point(Index n, bool init_x, Number* x,
                                      bool init_z, Number* z_L, Number* z_U,
                                      Index m, bool init_lambda,
                                      Number* lambda)
{
  // Here, we assume we only have starting values for x, if you code
  // your own NLP, you can provide starting values for the dual variables
  // if you wish
  assert(init_x == true);
  assert(init_z == false);
  assert(init_lambda == false);

  // We will initialize the free parameters with the measured states and the
  // directly identified gains and T*. This needs to come from loading some
  // data from disk.
  // initialize to the given starting point
  x[0] = 1.0;
  x[1] = 5.0;
  x[2] = 5.0;
  x[3] = 1.0;

  return true;
}

// returns the value of the objective function
bool GaitSystemID::eval_f(Index n, const Number* x, bool new_x, Number& obj_value)
{
  assert(n == 4);

  // This will minimize the sum of the squares of the measured minus the states
  // Extract the states from x and then compare to measured.
  //
  // x is in this order:
  // x11, ..., x1N, x21, ..., x2N, ..., x181, ..., x18N, t*, k
  // TODO : Think of a flattened order for t* and k

  obj_value = x[0] * x[3] * (x[0] + x[1] + x[2]) + x[2];

  return true;
}

// return the gradient of the objective function grad_{x} f(x)
bool HS071_NLP::eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f)
{
  assert(n == 4);

  grad_f[0] = x[0] * x[3] + x[3] * (x[0] + x[1] + x[2]);
  grad_f[1] = x[0] * x[3];
  grad_f[2] = x[0] * x[3] + 1;
  grad_f[3] = x[0] * (x[0] + x[1] + x[2]);

  return true;
}

// return the value of the constraints: g(x)
bool GaitSystemID::eval_g(Index n, const Number* x, bool new_x, Index m, Number* g)
{
  assert(n == 4);
  assert(m == 2);

  g[0] = x[0] * x[1] * x[2] * x[3];
  g[1] = x[0]*x[0] + x[1]*x[1] + x[2]*x[2] + x[3]*x[3];

  return true;
}

// return the structure or values of the jacobian
bool GaitSystemID::eval_jac_g(Index n, const Number* x, bool new_x,
                              Index m, Index nele_jac, Index* iRow, Index *jCol,
                              Number* values)
{
  if (values == NULL) {
    // return the structure of the jacobian

    // this particular jacobian is dense
    iRow[0] = 0;
    jCol[0] = 0;
    iRow[1] = 0;
    jCol[1] = 1;
    iRow[2] = 0;
    jCol[2] = 2;
    iRow[3] = 0;
    jCol[3] = 3;
    iRow[4] = 1;
    jCol[4] = 0;
    iRow[5] = 1;
    jCol[5] = 1;
    iRow[6] = 1;
    jCol[6] = 2;
    iRow[7] = 1;
    jCol[7] = 3;
  }
  else {
    // return the values of the jacobian of the constraints

    values[0] = x[1]*x[2]*x[3]; // 0,0
    values[1] = x[0]*x[2]*x[3]; // 0,1
    values[2] = x[0]*x[1]*x[3]; // 0,2
    values[3] = x[0]*x[1]*x[2]; // 0,3

    values[4] = 2*x[0]; // 1,0
    values[5] = 2*x[1]; // 1,1
    values[6] = 2*x[2]; // 1,2
    values[7] = 2*x[3]; // 1,3
  }

  return true;
}

//return the structure or values of the hessian
bool GaitSystemID::eval_h(Index n, const Number* x, bool new_x,
                          Number obj_factor, Index m, const Number* lambda,
                          bool new_lambda, Index nele_hess, Index* iRow,
                          Index* jCol, Number* values)
{
  return false;
}

void GaitSystemID::finalize_solution(SolverReturn status,
                                     Index n, const Number* x, const Number* z_L, const Number* z_U,
                                     Index m, const Number* g, const Number* lambda,
                                     Number obj_value,
                                     const IpoptData* ip_data,
                                     IpoptCalculatedQuantities* ip_cq)
{
  // here is where we would store the solution to variables, or write to a file, etc
  // so we could use the solution.

  // For this example, we write the solution to the console
  std::cout << std::endl << std::endl << "Solution of the primal variables, x" << std::endl;
  for (Index i=0; i<n; i++) {
     std::cout << "x[" << i << "] = " << x[i] << std::endl;
  }

  std::cout << std::endl << std::endl << "Solution of the bound multipliers, z_L and z_U" << std::endl;
  for (Index i=0; i<n; i++) {
    std::cout << "z_L[" << i << "] = " << z_L[i] << std::endl;
  }
  for (Index i=0; i<n; i++) {
    std::cout << "z_U[" << i << "] = " << z_U[i] << std::endl;
  }

  std::cout << std::endl << std::endl << "Objective value" << std::endl;
  std::cout << "f(x*) = " << obj_value << std::endl;

  std::cout << std::endl << "Final value of the constraints:" << std::endl;
  for (Index i=0; i<m ;i++) {
    std::cout << "g(" << i << ") = " << g[i] << std::endl;
  }
}
