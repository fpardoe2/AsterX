#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>

#include "c2p.hxx"
#include "c2p_1DEntropy.hxx"
#include "c2p_1DPalenzuela.hxx"
#include "c2p_1DRePrimAnd.hxx"
#include "c2p_2DNoble.hxx"

#include "c2p_utils.hxx"

#include "setup_eos.hxx"

namespace Con2PrimFactory {

using namespace Arith;
using namespace EOSX;
using namespace AsterUtils;

extern "C" void Con2PrimFactory_Test(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;

  // Get local eos object
  auto eos_3p_ig = global_eos_3p_ig;

  // Set atmo values
  const CCTK_REAL rho_atmo = 1e-10;
  CCTK_REAL eps_atmo = 1e-8;
  const CCTK_REAL Ye_atmo = 0.5;
  const CCTK_REAL press_atmo =
      eos_3p_ig->press_from_rho_eps_ye(rho_atmo, eps_atmo, Ye_atmo);
  const CCTK_REAL temp_atmo =
      eos_3p_ig->temp_from_rho_eps_ye(rho_atmo, eps_atmo, Ye_atmo);
  const CCTK_REAL entropy_atmo =
      eos_3p_ig->kappa_from_rho_eps_ye(rho_atmo, eps_atmo, Ye_atmo);

  // Setting up atmosphere
  const CCTK_REAL rho_atmo_cut = rho_atmo * (1 + 1.0e-3);
  atmosphere atmo(rho_atmo, eps_atmo, Ye_atmo, press_atmo, temp_atmo,
                  entropy_atmo, rho_atmo_cut);

  // Metric

  const CCTK_REAL alp = 1.0;

  const vec<CCTK_REAL, 3> beta{0.0, 0.0, 0.0};

  const smat<CCTK_REAL, 3> g{1.0, 0.0, 0.0,
                             1.0, 0.0, 1.0}; // xx, xy, xz, yy, yz, zz

  // Set BH limiters
  const CCTK_REAL alp_thresh = -1.;
  const CCTK_REAL rho_BH = 1e20;
  const CCTK_REAL eps_BH = 1e20;
  const CCTK_REAL vwlim_BH = 1e20;

  // Mag. limits
  const CCTK_REAL sigma_max = 100.;
  const CCTK_REAL inv_beta_max = 100.;

  // Con2Prim objects
  // (eos_3p_ig, atmo, max_iter, c2p_tol
  //  alp_thresh,
  //  vw_lim, B_lim, rho_BH, eps_BH, vwlim_BH,
  //  Ye_lenient, use_z, use_temperature)
  c2p_2DNoble c2p_Noble(eos_3p_ig, atmo, 100, 1e-8, alp_thresh, 1, 1, rho_BH,
                        eps_BH, vwlim_BH, sigma_max, inv_beta_max, true, false,
                        false, false, false, 1.0);
  c2p_1DPalenzuela c2p_Pal(eos_3p_ig, atmo, 100, 1e-8, alp_thresh, 1, 1, rho_BH,
                           eps_BH, vwlim_BH, sigma_max, inv_beta_max, true,
                           false, false, false, false, 1.0);
  c2p_1DRePrimAnd c2p_RPA(eos_3p_ig, atmo, 100, 1e-8, alp_thresh, 1, 1, rho_BH,
                          eps_BH, vwlim_BH, sigma_max, inv_beta_max, true,
                          false, false, false, false, 1.0);
  c2p_1DEntropy c2p_Ent(eos_3p_ig, atmo, 100, 1e-8, alp_thresh, 1, 1, rho_BH,
                        eps_BH, vwlim_BH, sigma_max, inv_beta_max, true, false,
                        false, false, false, 1.0);

  // Construct error report object:
  c2p_report rep_Noble;
  c2p_report rep_Pal;
  c2p_report rep_RPA;
  c2p_report rep_Ent;

  // Set primitive seeds
  const CCTK_REAL rho_in = 0.125;
  CCTK_REAL eps_in = 0.8;
  const CCTK_REAL Ye_in = 0.5;
  const CCTK_REAL press_in =
      eos_3p_ig->press_from_rho_eps_ye(rho_in, eps_in, Ye_in);
  const CCTK_REAL temp_in =
      eos_3p_ig->temp_from_rho_eps_ye(rho_in, eps_in, Ye_in);
  const CCTK_REAL entropy_in =
      eos_3p_ig->kappa_from_rho_eps_ye(rho_in, eps_in, Ye_in);
  const vec<CCTK_REAL, 3> vup_in = {0.0, 0.0, 0.0};
  const vec<CCTK_REAL, 3> Bup_in = {0.5, -0.5, 0.0};
  const vec<CCTK_REAL, 3> vdown_in = calc_contraction(g, vup_in);
  const CCTK_REAL wlor_in = calc_wlorentz(vdown_in, vup_in);

  prim_vars pv;
  // rho(p.I), eps(p.I), dummy_Ye, press(p.I), entropy, v_up, wlor, Bup
  prim_vars pv_seeds{rho_in,     eps_in, Ye_in,   press_in, temp_in,
                     entropy_in, vup_in, wlor_in, Bup_in};

  // cons_vars cv{dens(p.I), {momx(p.I), momy(p.I), momz(p.I)}, tau(p.I),
  //  dummy_DYe, DEnt, {dBx(p.I), dBy(p.I), dBz(p.I)}};

  cons_vars cv_Noble;
  cons_vars cv_Pal;
  cons_vars cv_Ent;
  cons_vars cv_all;

  // C2P solve may modify cons_vars input
  // Use the following to test C2Ps independently
  cv_Noble.from_prim(pv_seeds, g);
  cv_Pal.from_prim(pv_seeds, g);
  cv_Ent.from_prim(pv_seeds, g);

  // Use the following to test C2Ps in sequence
  // (As done in AsterX, the evolution thorn)
  cv_all.from_prim(pv_seeds, g);

  // Testing C2P Noble
  CCTK_VINFO("Testing C2P Noble...");
  c2p_Noble.solve(eos_3p_ig, pv, pv_seeds, cv_all, alp, beta, g, rep_Noble);

  printf("pv_seeds, pv: \n"
         "rho: %f, %f \n"
         "eps: %f, %f \n"
         "Ye: %f, %f \n"
         "press: %f, %f \n"
         "temperature: %f, %f \n"
         "entropy: %f, %f \n"
         "velx: %f, %f \n"
         "vely: %f, %f \n"
         "velz: %f, %f \n"
         "Bx: %f, %f \n"
         "By: %f, %f \n"
         "Bz: %f, %f \n",
         pv_seeds.rho, pv.rho, pv_seeds.eps, pv.eps, pv_seeds.Ye, pv.Ye,
         pv_seeds.press, pv.press, pv_seeds.temperature, pv.temperature,
         pv_seeds.entropy, pv.entropy, pv_seeds.vel(0), pv.vel(0),
         pv_seeds.vel(1), pv.vel(1), pv_seeds.vel(2), pv.vel(2),
         pv_seeds.Bvec(0), pv.Bvec(0), pv_seeds.Bvec(1), pv.Bvec(1),
         pv_seeds.Bvec(2), pv.Bvec(2));
  printf("cv: \n"
         "dens: %f \n"
         "tau: %f \n"
         "momx: %f \n"
         "momy: %f \n"
         "momz: %f \n"
         "DYe: %f \n"
         "dBx: %f \n"
         "dBy: %f \n"
         "dBz: %f \n"
         "DEnt: %f \n",
         cv_all.dens, cv_all.tau, cv_all.mom(0), cv_all.mom(1), cv_all.mom(2),
         cv_all.DYe, cv_all.dBvec(0), cv_all.dBvec(1), cv_all.dBvec(2),
         cv_all.DEnt);
  /*
    assert(pv.rho == pv_seeds.rho);
    assert(pv.eps == pv_seeds.eps);
    assert(pv.press == pv_seeds.press);
    assert(pv.vel == pv_seeds.vel);
    assert(pv.Bvec == pv_seeds.Bvec);
  */

  rep_Noble.debug_message();

  // Testing C2P Palenzuela
  CCTK_VINFO("Testing C2P Palenzuela...");
  // c2p_Pal.solve(eos_3p_ig, pv, pv_seeds, cv, g, rep_Pal);
  c2p_Pal.solve(eos_3p_ig, pv, cv_all, alp, beta, g, rep_Pal);

  printf("pv_seeds, pv: \n"
         "rho: %f, %f \n"
         "eps: %f, %f \n"
         "Ye: %f, %f \n"
         "press: %f, %f \n"
         "temperature: %f, %f \n"
         "entropy: %f, %f \n"
         "velx: %f, %f \n"
         "vely: %f, %f \n"
         "velz: %f, %f \n"
         "Bx: %f, %f \n"
         "By: %f, %f \n"
         "Bz: %f, %f \n",
         pv_seeds.rho, pv.rho, pv_seeds.eps, pv.eps, pv_seeds.Ye, pv.Ye,
         pv_seeds.press, pv.press, pv_seeds.temperature, pv.temperature,
         pv_seeds.entropy, pv.entropy, pv_seeds.vel(0), pv.vel(0),
         pv_seeds.vel(1), pv.vel(1), pv_seeds.vel(2), pv.vel(2),
         pv_seeds.Bvec(0), pv.Bvec(0), pv_seeds.Bvec(1), pv.Bvec(1),
         pv_seeds.Bvec(2), pv.Bvec(2));
  printf("cv: \n"
         "dens: %f \n"
         "tau: %f \n"
         "momx: %f \n"
         "momy: %f \n"
         "momz: %f \n"
         "DYe: %f \n"
         "dBx: %f \n"
         "dBy: %f \n"
         "dBz: %f \n"
         "DEnt: %f \n",
         cv_all.dens, cv_all.tau, cv_all.mom(0), cv_all.mom(1), cv_all.mom(2),
         cv_all.DYe, cv_all.dBvec(0), cv_all.dBvec(1), cv_all.dBvec(2),
         cv_all.DEnt);
  /*
    assert(pv.rho == pv_seeds.rho);
    assert(pv.eps == pv_seeds.eps);
    assert(pv.press == pv_seeds.press);
    assert(pv.vel == pv_seeds.vel);
    assert(pv.Bvec == pv_seeds.Bvec);
  */
  rep_Pal.debug_message();

  // Testing C2P RePrimAnd
  CCTK_VINFO("Testing C2P RePrimAnd...");
  c2p_RPA.solve(eos_3p_ig, pv, cv_all, alp, beta, g, rep_RPA);

  printf("pv_seeds, pv: \n"
         "rho: %f, %f \n"
         "eps: %f, %f \n"
         "Ye: %f, %f \n"
         "press: %f, %f \n"
         "temperature: %f, %f \n"
         "entropy: %f, %f \n"
         "velx: %f, %f \n"
         "vely: %f, %f \n"
         "velz: %f, %f \n"
         "Bx: %f, %f \n"
         "By: %f, %f \n"
         "Bz: %f, %f \n",
         pv_seeds.rho, pv.rho, pv_seeds.eps, pv.eps, pv_seeds.Ye, pv.Ye,
         pv_seeds.press, pv.press, pv_seeds.temperature, pv.temperature,
         pv_seeds.entropy, pv.entropy, pv_seeds.vel(0), pv.vel(0),
         pv_seeds.vel(1), pv.vel(1), pv_seeds.vel(2), pv.vel(2),
         pv_seeds.Bvec(0), pv.Bvec(0), pv_seeds.Bvec(1), pv.Bvec(1),
         pv_seeds.Bvec(2), pv.Bvec(2));
  printf("cv: \n"
         "dens: %f \n"
         "tau: %f \n"
         "momx: %f \n"
         "momy: %f \n"
         "momz: %f \n"
         "DYe: %f \n"
         "dBx: %f \n"
         "dBy: %f \n"
         "dBz: %f \n"
         "DEnt: %f \n",
         cv_all.dens, cv_all.tau, cv_all.mom(0), cv_all.mom(1), cv_all.mom(2),
         cv_all.DYe, cv_all.dBvec(0), cv_all.dBvec(1), cv_all.dBvec(2),
         cv_all.DEnt);
  rep_RPA.debug_message();

  // Testing C2P Entropy
  CCTK_VINFO("Testing C2P Entropy...");
  c2p_Ent.solve(eos_3p_ig, pv, cv_all, alp, beta, g, rep_Ent);

  printf("pv_seeds, pv: \n"
         "rho: %f, %f \n"
         "eps: %f, %f \n"
         "Ye: %f, %f \n"
         "press: %f, %f \n"
         "temperature: %f, %f \n"
         "entropy: %f, %f \n"
         "velx: %f, %f \n"
         "vely: %f, %f \n"
         "velz: %f, %f \n"
         "Bx: %f, %f \n"
         "By: %f, %f \n"
         "Bz: %f, %f \n",
         pv_seeds.rho, pv.rho, pv_seeds.eps, pv.eps, pv_seeds.Ye, pv.Ye,
         pv_seeds.press, pv.press, pv_seeds.temperature, pv.temperature,
         pv_seeds.entropy, pv.entropy, pv_seeds.vel(0), pv.vel(0),
         pv_seeds.vel(1), pv.vel(1), pv_seeds.vel(2), pv.vel(2),
         pv_seeds.Bvec(0), pv.Bvec(0), pv_seeds.Bvec(1), pv.Bvec(1),
         pv_seeds.Bvec(2), pv.Bvec(2));
  printf("cv: \n"
         "dens: %f \n"
         "tau: %f \n"
         "momx: %f \n"
         "momy: %f \n"
         "momz: %f \n"
         "DYe: %f \n"
         "dBx: %f \n"
         "dBy: %f \n"
         "dBz: %f \n"
         "DEnt: %f \n",
         cv_all.dens, cv_all.tau, cv_all.mom(0), cv_all.mom(1), cv_all.mom(2),
         cv_all.DYe, cv_all.dBvec(0), cv_all.dBvec(1), cv_all.dBvec(2),
         cv_all.DEnt);
  /*
    assert(pv.rho == pv_seeds.rho);
    assert(pv.eps == pv_seeds.eps);
    assert(pv.press == pv_seeds.press);
    assert(pv.vel == pv_seeds.vel);
    assert(pv.Bvec == pv_seeds.Bvec);
  */
  rep_Ent.debug_message();
}

} // namespace Con2PrimFactory
