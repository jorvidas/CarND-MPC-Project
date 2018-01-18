#include "MPC.h"
#include <cppad/cppad.hpp>
#include <cppad/ipopt/solve.hpp>
#include "Eigen-3.3/Eigen/Core"

using CppAD::AD;

/*
* TUNABLE
*/

// Set the timestep length and duration
int N = 9;
double dt = 0.10;

// This value assumes the model presented in the classroom is used.
//
// It was obtained by measuring the radius formed by running the vehicle in the
// simulator around in a circle with a constant steering angle and velocity on a
// flat terrain.
//
// Lf was tuned until the the radius formed by the simulating the model
// presented in the classroom matched the previous radius.
//
// This is the length from front to CoG that has a similar radius.
const double Lf = 2.67;

// target velocity
double ref_v = 80.0;

// variable positions for convenience
int x_start = 0;
int y_start = x_start + N;
int psi_start = y_start + N;
int v_start = psi_start + N;
int cte_start = v_start + N;
int epsi_start = cte_start + N;
int delta_start = epsi_start + N;
int a_start = delta_start + N - 1;

class FG_eval {
 public:
  // Fitted polynomial coefficients
  Eigen::VectorXd coeffs;
  FG_eval(Eigen::VectorXd coeffs) {this->coeffs = coeffs;}

  typedef CPPAD_TESTVECTOR(AD<double>) ADvector;
  void operator()(ADvector& fg, const ADvector& vars) {
    /*
    * FG EVAL:
    * - CALCULATE COST
    * - SET MODEL CONSTRAINTS
    */

    /*
    * CALCULATE COST
    */

    /*
    * TUNABLE
    */

    // Cost Coefficients
    // state
    int cte_cc = 2500;
    int epsi_cc = 3000;
    int ref_v_cc = 1;
    // minimize actuators
    int delta_cc = 5;
    int a_cc = 5;
    // delta actuators
    int delta_delta_cc = 200;
    int delta_a_cc = 5;

    // cost update for reference state.
    for (int t = 0; t < N; t++) {
      fg[0] += cte_cc * CppAD::pow(vars[cte_start + t], 2);
      fg[0] += epsi_cc * CppAD::pow(vars[epsi_start + t], 2);
      fg[0] += ref_v_cc * CppAD::pow(vars[v_start + t] - ref_v, 2);
    }

    // cost update for minimizing the use of actuators.
    for (int t = 0; t < N - 1; t++) {
      fg[0] += delta_cc * CppAD::pow(vars[delta_start + t], 2);
      fg[0] += a_cc * CppAD::pow(vars[a_start + t], 2);
    }

    // cost update for minimizing the value gap between sequential actuations.
    for (int t = 0; t < N - 2; t++) {
      fg[0] += delta_delta_cc * CppAD::pow(vars[delta_start + t + 1] -
               vars[delta_start + t], 2);
      fg[0] += delta_a_cc * CppAD::pow(vars[a_start + t + 1] -
               vars[a_start + t], 2);
    }

    /*
    * SET MODEL CONSTRAINTS
    */

    // initial state constraints (given)
    // everything shifted 1 to accomodate cost at fg[0]
    fg[1 + x_start] = vars[x_start];
    fg[1 + y_start] = vars[y_start];
    fg[1 + psi_start] = vars[psi_start];
    fg[1 + v_start] = vars[v_start];
    fg[1 + cte_start] = vars[cte_start];
    fg[1 + epsi_start] = vars[epsi_start];

    // The rest of the constraints (model)
    for (int t = 1; t < N; t++) {
      // PREP
      // The state at time t+1 .
      AD<double> x1 = vars[x_start + t];
      AD<double> y1 = vars[y_start + t];
      AD<double> psi1 = vars[psi_start + t];
      AD<double> v1 = vars[v_start + t];
      AD<double> cte1 = vars[cte_start + t];
      AD<double> epsi1 = vars[epsi_start + t];

      // The state at time t.
      AD<double> x0 = vars[x_start + t - 1];
      AD<double> y0 = vars[y_start + t - 1];
      AD<double> psi0 = vars[psi_start + t - 1];
      AD<double> v0 = vars[v_start + t - 1];
      AD<double> cte0 = vars[cte_start + t - 1];
      AD<double> epsi0 = vars[epsi_start + t - 1];

      // Only consider the actuation at time t.
      AD<double> delta0 = vars[delta_start + t - 1];
      AD<double> a0 = vars[a_start + t - 1];

      AD<double> f0 = coeffs[0] + coeffs[1] * x0 + coeffs[2] * x0 * x0 +
                      coeffs[3] * x0 * x0 * x0;
      AD<double> psides0 = CppAD::atan(3 * coeffs[3] * x0 * x0 +
                                       2 * coeffs[2] * x0 + coeffs[1]);

      // equations for the model:
      // x_[t+1] = x[t] + v[t] * cos(psi[t]) * dt
      // y_[t+1] = y[t] + v[t] * sin(psi[t]) * dt
      // psi_[t+1] = psi[t] + v[t] / Lf * delta[t] * dt
      // v_[t+1] = v[t] + a[t] * dt
      // cte[t+1] = f(x[t]) - y[t] + v[t] * sin(epsi[t]) * dt
      // epsi[t+1] = psi[t] - psides[t] + v[t] * delta[t] / Lf * dt

      // next x
      fg[1 + x_start + t] = x1 - (x0 + v0 * CppAD::cos(psi0) * dt);
      // next y
      fg[1 + y_start + t] = y1 - (y0 + v0 * CppAD::sin(psi0) * dt);
      // next psi
      fg[1 + psi_start + t] = psi1 - (psi0 - v0 * delta0 / Lf * dt);
      // next v
      fg[1 + v_start + t] = v1 - (v0 + a0 * dt);
      // next cte
      fg[1 + cte_start + t] =
          cte1 - ((f0 - y0) + (v0 * CppAD::sin(epsi0) * dt));
      // next epsi
      fg[1 + epsi_start + t] =
          epsi1 - ((psi0 - psides0) - v0 * delta0 / Lf * dt);
    }
  }
};

//
// MPC class definition implementation.
//
MPC::MPC() {}
MPC::~MPC() {}

vector<double> MPC::Solve(Eigen::VectorXd state, Eigen::VectorXd coeffs) {
  bool ok = true;
  typedef CPPAD_TESTVECTOR(double) Dvector;
  /*
  *  MPC components:
  *  - SET UP FOR SOLVER
  *  - SOLVE
  *  - RETURN ACTUATORS
  *
  * Solver will need:
  *  - VARIABLES
  *  - BOUNDING ON VARIABLES
  *  - BOUNDING ON CONSTRAINTS
  *  - COST/CONSTRAINT OBJECT
  *  - OPTIONS
  *  - SOLUTION SPACE
  */

  /*
  * TUNABLE
  */

  double latency = 0.1;                           // also in main.cpp
  int latent_steps = (latency + 0.5 * dt) / dt;

  /*
  * SET UP FOR SOLVER
  */

  // VARIABLES

  // state for convenience
  double x = state[0];
  double y = state[1];
  double psi = state[2];
  double v = state[3];
  double cte = state[4];
  double epsi = state[5];

  // size_of_state_vect * positions + size_actuator_vect * (Position - 1)
  // one less actuator because no actuator for final position
  int n_vars = 6 * N + 2 * (N - 1);
  // states only
  int n_constraints = 6 * N;

  // initial value of the independent variables.
  // future values initialized to 0
  Dvector vars(n_vars);
  for (int i = 0; i < n_vars; i++) {
    vars[i] = 0;
  }

  // current state initialized to current state
  vars[x_start] = x;
  vars[y_start] = y;
  vars[psi_start] = psi;
  vars[v_start] = v;
  vars[cte_start] = cte;
  vars[epsi_start] = epsi;

  // BOUNDING ON VARIABLES
  Dvector vars_lowerbound(n_vars);
  Dvector vars_upperbound(n_vars);

  //  can get to any position if actuators allow it
  for (int i = 0; i < delta_start; i++) {
    vars_lowerbound[i] = -1.0e19;
    vars_upperbound[i] = 1.0e19;
  }

  // steering angle is limited to -25*/25*
  for (int i = delta_start; i < a_start; i++) {
    vars_lowerbound[i] = -0.436332*Lf;
    vars_upperbound[i] = 0.436332*Lf;
  }

  // acceleration is limited to -1/1
  for (int i = a_start; i < n_vars; i++) {
    vars_lowerbound[i] = -1.0;
    vars_upperbound[i] = 1.0;
  }

  // BOUNDING ON CONSTRAINTS
  Dvector constraints_lowerbound(n_constraints);
  Dvector constraints_upperbound(n_constraints);

  // all future states should be constrained to 0 as both halves of model
  // calculations offset
  // Example:
  // future_x = current_x + change_in_x
  // 0 = future_x - (current_x + change_in_x)
  for (int i = 0; i < n_constraints; i++) {
    constraints_lowerbound[i] = 0;
    constraints_upperbound[i] = 0;
  }

  // current state needs to be current state, not calculated by model, is given
  constraints_lowerbound[x_start] = x;
  constraints_lowerbound[y_start] = y;
  constraints_lowerbound[psi_start] = psi;
  constraints_lowerbound[v_start] = v;
  constraints_lowerbound[cte_start] = cte;
  constraints_lowerbound[epsi_start] = epsi;

  constraints_upperbound[x_start] = x;
  constraints_upperbound[y_start] = y;
  constraints_upperbound[psi_start] = psi;
  constraints_upperbound[v_start] = v;
  constraints_upperbound[cte_start] = cte;
  constraints_upperbound[epsi_start] = epsi;

  // COST/CONSTRAINT OBJECT
  FG_eval fg_eval(coeffs);

  // OPTIONS
  // options for IPOPT solver
  std::string options;
  // Uncomment this if you'd like more print information
  options += "Integer print_level  0\n";
  // NOTE: Setting sparse to true allows the solver to take advantage
  // of sparse routines, this makes the computation MUCH FASTER. If you
  // can uncomment 1 of these and see if it makes a difference or not but
  // if you uncomment both the computation time should go up in orders of
  // magnitude.
  options += "Sparse  true        forward\n";
  options += "Sparse  true        reverse\n";
  // NOTE: Currently the solver has a maximum time limit of 0.5 seconds.
  // Change this as you see fit.
  options += "Numeric max_cpu_time          0.5\n";

  // SOLUTION SPACE
  CppAD::ipopt::solve_result<Dvector> solution;

  /*
  * SOLVE
  */

  // uses ipopt for finding best actuator values
  CppAD::ipopt::solve<Dvector, FG_eval>(
      options, vars, vars_lowerbound, vars_upperbound, constraints_lowerbound,
      constraints_upperbound, fg_eval, solution);

  // Check some of the solution values
  ok &= solution.status == CppAD::ipopt::solve_result<Dvector>::success;

  /*
  * RETURN NEXT ACTUATORS
  */

  // for dt < latency
  double delta_to_pass = 0.0;
  double a_to_pass = 0.0;

  // displacement weighted - get back to center
  for (int i = 0; i < latent_steps; i++) {
    double weighting = (latent_steps - i - 0.5);
    delta_to_pass += weighting * solution.x[delta_start + i];
  }
  // normalized
  double normalizer = (latent_steps * latent_steps * 0.5);
  delta_to_pass /= normalizer;

  // end state weighted - adjust speed to necessary speed
  for (int i = 0; i < latent_steps; i++) {
    a_to_pass += solution.x[a_start + i];
  }
  // normalized
  a_to_pass /= latent_steps;

  // return
  vector<double> solved = {delta_to_pass, a_to_pass};
  for (int i = 0; i < N; i++) {
    solved.push_back(solution.x[x_start + i]);
    solved.push_back(solution.x[y_start + i]);
  }

  return solved;
}
