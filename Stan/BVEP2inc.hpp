#include "stan/math/prim/meta/is_rev_matrix.hpp"
#include <stan/math.hpp>
#include <stan/math/rev/core/reverse_pass_callback.hpp>

namespace BVEP3_seeg_noncentered_reparam_model_namespace {

  // namespaces
  using namespace stan;
  using namespace stan::math;
  using namespace Eigen;
  using namespace std;

  /// /// /// /// /// /// /// /// /// debug tools
#define T(e) (e).transpose()
  // https://stackoverflow.com/a/27375675
  template <typename Arg, typename... Args>
    void print(Arg&& arg, Args&&... args)
    {
      std::cout << std::forward<Arg>(arg);
      using expander = int[];
      (void)expander{0, (void(std::cout << std::forward<Args>(args)), 0)...};
      std::cout << std::endl;
    }
#define dbgvar(v) print("[ode_step_rev] " #v " = ", T(v), " adj ", T(v.adj()), " adj_op ", T(v.adj_op()))

  /// /// /// /// /// /// /// /// /// ODE RHS fwd & rev

  // fwd pass
  template <typename Mat, typename Mat2, typename Mat3>
  VectorXd ode_rhs_c(const double time, const Mat3& xz, const Mat2& SC, 
      const double I1, const double tau0, const double K,
      const Mat& eta, ostream* pstream) {
    int nn = xz.size() / 2;
    double rtau0 = 1.0 / tau0;
    VectorXd dxz(xz.size());
    for (int i=0; i<nn; i++)
    {
      double acc = 0.0;
      for (int j=0; j<nn; j++)
        acc += SC(i,j)*xz(j);
      double x = xz(i), z = xz(i+nn);
      dxz(i) = 1.0 - x*x*x - 2*x*x - z + I1;
      dxz(i+nn) = rtau0*(4*(x - eta(i)) - z - K*acc);
    }
    return dxz;
  }

  // reusable reverse mode vjp
  template <typename VarMat, typename VarMat2, typename Var, typename Mat>
  inline auto ode_rhs_rev(
      int nn, VarMat& dxz_a, VarMat2& xz_a, VarMat& eta_a,
      Var& K_a, Mat& SC_a, double rtau0, double I1) {
    VectorXd g = dxz_a.adj_op();
    VectorXd g_dx = g.head(nn), g_dz = g.tail(nn);
    VectorXd x = value_of(xz_a.head(nn));
    eta_a.adj() += -g_dz.tail(nn)*rtau0*4;
    VectorXd gx = SC_a * x;
    double ksum = 0.0;
    for (int i=0;i<nn;i++) ksum += g_dz(i)*gx(i);
    K_a.adj() -= ksum * rtau0;
    VectorXd g_x(xz_a.size());
    g_x.fill(0.0);
    g_x.head(nn).array() += g_dx.array()*(-3*x.array().square() - 4*x.array());
    g_x.head(nn) += g_dz*4*rtau0 - value_of(K_a)*rtau0*(SC_a.transpose()*g_dz);
    g_x.tail(nn) += -g_dx;
    g_x.tail(nn) += -g_dz*rtau0;
    xz_a.adj() += g_x;
  }

  // fwd + reverse when args are var types
  template <typename VarMat, typename Var, typename Mat,
           require_rev_matrix_t<VarMat>* = nullptr >
  VarMat ode_rhs_c(const double time, const VarMat& xz, const Mat& SC,
     const double I1, const double tau0, const Var& K,
     const VarMat& eta, ostream* pstream) { 
   int nn = xz.size() / 2;
   double rtau0 = 1.0 / tau0;
   // call fwd pass impl
   VectorXd dxz_val = ode_rhs_c(
       time, value_of(xz), SC, I1, tau0, 
       value_of(K), value_of(eta), pstream);
   VarMat dxz_(dxz_val);
   // construct bwd pass
   arena_t<VarMat> dxz_a(dxz_);
   arena_t<VarMat> xz_a(xz);
   arena_t<VarMat> eta_a(eta);
   arena_t<Var> K_a(K);
   arena_t<Mat> SC_a(SC);
   reverse_pass_callback([nn, dxz_a, xz_a, eta_a, K_a, SC_a, rtau0, I1]() mutable {
       ode_rhs_rev(nn, dxz_a, xz_a, eta_a, K_a, SC_a, rtau0, I1);
       });
   // return results
   return dxz_;
  }

  /// /// /// /// /// /// /// /// /// Heun step fwd & rev

  // fwd pass
  template <typename Mat, typename Mat2, typename Mat3>
  VectorXd ode_step_c(const double& time, const double& dt, const Mat& xz,
      const Mat2& SC, 
      const double& I1, const double& tau0, const double& K,
      const Mat3& eta, ostream* pstream) {
    VectorXd d1 = ode_rhs_c(time, xz, SC, I1, tau0, K, eta, pstream);
    VectorXd d2 = ode_rhs_c(time, xz + dt * d1, SC, I1, tau0, K, eta, pstream);
    return xz + dt/2*(d1 + d2);
  }

  typedef Matrix<var_value<double>,Dynamic,1> VarVec;

  // reusable reverse mode vjp
  template <typename VVec1, typename VVec2, typename VVec3, typename Var, typename Mat>
  void ode_step_rev(
      int nn, double dt, VVec1& nx, VVec2& xz, VVec3& eta,
      Var& K, Mat& SC, double rtau0, double I1)
  {
    // 1st state of Heun needed for the VJP
    auto d1 = ode_rhs_c(0.0, value_of(xz), SC, I1, 1/rtau0, 
        value_of(K), value_of(eta), nullptr);
    // nx = ode_step_c(), nx.adj() is V of VJP
    // nx = x + dt/2*(d1 + d2)
    arena_t<VarVec> g_d1(nx*dt/2);
    arena_t<VarVec> g_d2(nx*dt/2);
    g_d1.adj() = nx.adj()*dt/2;
    g_d2.adj() = nx.adj()*dt/2;
    xz.adj() += nx.adj();
    // d2 = ode_rhs(x + dt*d1, e, k)
    arena_t<VarVec> xzi(xz + dt*d1);
    ode_rhs_rev(nn, g_d2, xzi, eta, K, SC, rtau0, I1);
    xz.adj() += xzi.adj();
    g_d1.adj() += xzi.adj() * dt;
    // d1 = ode_rhs(x, e, k)
    ode_rhs_rev(nn, g_d1, xz, eta, K, SC, rtau0, I1);
    // xz.adj() += xz.adj()? no need
  }

  // fwd + reverse when args are var types
  template <typename VarMat, typename Var, typename Mat,
           require_rev_matrix_t<VarMat>* = nullptr >
  VarMat ode_step_c(const double time, const double dt, const VarMat& xz,
     const Mat& SC, const double I1, const double tau0, const Var& K,
     const VarMat& eta, ostream* pstream) { 
   int nn = xz.size() / 2;
   double rtau0 = 1.0 / tau0;
   // call fwd pass impl
   VectorXd nx = ode_step_c(
       time, dt, value_of(xz), SC, I1, tau0, 
       value_of(K), value_of(eta), pstream);
   VarMat ret(nx);
   // construct bwd pass
   arena_t<VarMat> nx_a(ret);
   arena_t<VarMat> xz_a(xz);
   arena_t<VarMat> eta_a(eta);
   arena_t<Var> K_a(K);
   arena_t<Mat> SC_a(SC);
   reverse_pass_callback([nn, dt, nx_a, xz_a, eta_a, K_a, SC_a, rtau0, I1]() mutable {
       ode_step_rev(nn, dt, nx_a, xz_a, eta_a, K_a, SC_a, rtau0, I1);
       });
   // return results
   return ret;
  }

  /// /// /// /// /// /// /// /// /// Solve loop fwd & rev
  
  template <typename M1, typename M2, typename M3>
  MatrixXd ode_sol_c(const double& dt, const int& nt, const M1& xz,
    const M2& SC, const double& I1, const double& tau0, const double& K,
    const M3& eta, ostream* pstream)
  {
    MatrixXd sol(xz.size(), nt);
    sol.col(0) = xz;
    for (int t=1; t<nt; t++)
      sol.col(t) = ode_step_c(t*dt, dt, sol.col(t-1), SC, I1, tau0, K, eta, pstream);
    return sol;
  }

  // above, we're a bit loose with template arg VarMat being a vector
  // but below Eigen doesn't like it in `VarMat sol(sol_val)` because
  // it needs to be matrix of (m>1, n>1)
  typedef Matrix<stan::math::var_value<double>, Dynamic, Dynamic> VarMat;
  // TODO revise above functions to use VarVec instead
 
  template <typename VarMat, typename VarVec, typename Var, typename Mat,
            require_rev_matrix_t<VarMat>* = nullptr>
  void ode_sol_rev(const double dt, VarMat& sol, VarVec& xz, Mat& SC,
      Var& K, VarVec& eta, const double I1, const double tau0)
  {
    const int nt = sol.cols(), nn = sol.rows() / 2;
    VarVec nx(xz.size()), x(xz.size());
    for (int t=nt-1; t>0; t--)
    {
      x = sol.col(t-1);
      nx = sol.col(t);
      ode_step_rev(nn,dt,nx,x,eta,K,SC,1/tau0,I1);
    }
    xz.adj() += x.adj();
  }

  template <typename VarVec, typename Var, typename Mat,
           require_rev_matrix_t<VarVec>* = nullptr >
  VarMat ode_sol_c(const double& dt, const int& nt, const VarVec& xz,
    const Mat& SC, const double& I1, const double& tau0, const Var& K,
    const VarVec& eta, ostream* pstream) {
    MatrixXd sol_val = ode_sol_c(dt, nt, value_of(xz), SC, I1, tau0,
        value_of(K), value_of(eta), pstream);
    VarMat sol(sol_val);
    // copy stuff into arena for rev
    arena_t<VarMat> sol_a(sol);
    arena_t<VarVec> xz_a(xz);
    arena_t<Mat> SC_a(SC);
    arena_t<Var> K_a(K);
    arena_t<VarVec> eta_a(eta);
    // push rev onto ad stack
    reverse_pass_callback([dt,sol_a,xz_a,SC_a,K_a,eta_a,I1,tau0]() mutable {
      ode_sol_rev(dt,sol_a,xz_a,SC_a,K_a,eta_a,I1,tau0);
    });
    return sol;
  }

} // *_model_namespace
