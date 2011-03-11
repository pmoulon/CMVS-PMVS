
PROJECT(gslLib CXX)

ADD_SUBDIRECTORY(cblas)

INCLUDE_DIRECTORIES(./ ./include)

#defined by hand
ADD_LIBRARY(gsl
#.\config.h
#.\gsl_version.h

./templates_on.h

./version.c
					./block/block.c
					./block/file.c
./block/init.c
./fft/dft.c
./fft/factorize.c
					./fft/fft.c
					./fft/signals.c
					./min/bracketing.c

					./min/brent.c
					./min/convergence.c
					./min/fsolver.c
					./min/golden.c
					./min/quad_golden.c
					./blas/blas.c

					./cdf/beta.c
					./cdf/betainv.c
					./cdf/binomial.c
					./cdf/cauchy.c
					./cdf/cauchyinv.c
					./cdf/chisq.c
					./cdf/chisqinv.c
					./cdf/exponential.c
					./cdf/exponentialinv.c
					./cdf/exppow.c
					./cdf/fdist.c
					./cdf/fdistinv.c
					./cdf/flat.c
					./cdf/flatinv.c
					./cdf/gamma.c
					./cdf/gammainv.c
					./cdf/gauss.c
					./cdf/gaussinv.c
					./cdf/geometric.c
					./cdf/gumbel1.c
					./cdf/gumbel1inv.c
					./cdf/gumbel2.c
					./cdf/gumbel2inv.c
					./cdf/hypergeometric.c
					./cdf/laplace.c
					./cdf/laplaceinv.c
					./cdf/logistic.c
					./cdf/logisticinv.c
					./cdf/lognormal.c
					./cdf/lognormalinv.c
					./cdf/nbinomial.c
					./cdf/pareto.c
					./cdf/paretoinv.c
					./cdf/pascal.c
					./cdf/poisson.c


					./cdf/rayleigh.c



					./cdf/rayleighinv.c



					./cdf/tdist.c



					./cdf/tdistinv.c



					./cdf/weibull.c



					./cdf/weibullinv.c



					./cheb/deriv.c



					./cheb/eval.c



					./cheb/init.c





					./cheb/integ.c





					./combination/combination.c



					./combination/file.c




					./combination/init.c



					./complex/math.c




					./deriv/deriv.c

			

					./dht/dht.c




					./diff/diff.c




					#./eigen/francis.c



					./eigen/gen.c



					./eigen/genherm.c



					./eigen/genhermv.c



					./eigen/gensymm.c



					./eigen/gensymmv.c



					./eigen/genv.c



					./eigen/herm.c



					./eigen/hermv.c



					./eigen/jacobi.c



					./eigen/nonsymm.c



					./eigen/nonsymmv.c



					./eigen/schur.c



					./eigen/sort.c



					./eigen/symm.c



					./eigen/symmv.c




					./err/error.c



					./err/message.c



					./err/stream.c



					./err/strerror.c




					./fit/linear.c



					./histogram/add.c



					./histogram/add2d.c



					./histogram/calloc_range.c



					./histogram/calloc_range2d.c



					./histogram/copy.c



					./histogram/copy2d.c



					./histogram/file.c




					./histogram/file2d.c



					./histogram/get.c



					./histogram/get2d.c



					./histogram/init.c


						



					./histogram/init2d.c



					./histogram/maxval.c



					./histogram/maxval2d.c



					./histogram/oper.c



					./histogram/oper2d.c



					./histogram/params.c



					./histogram/params2d.c



					./histogram/pdf.c



					./histogram/pdf2d.c



					./histogram/reset.c



					./histogram/reset2d.c



					./histogram/stat.c



					./histogram/stat2d.c




					./ieee-utils/endian.c



					./ieee-utils/env.c



					./ieee-utils/fp.c



					./ieee-utils/make_rep.c



					./ieee-utils/print.c



					./ieee-utils/read.c



					./ieee-utils/standardize.c




					./integration/qag.c



					./integration/qagp.c



					./integration/qags.c



					./integration/qawc.c



					./integration/qawf.c



					./integration/qawo.c



					./integration/qaws.c



					./integration/qcheb.c



					./integration/qk.c



					./integration/qk15.c



					./integration/qk21.c



					./integration/qk31.c



					./integration/qk41.c



					./integration/qk51.c



					./integration/qk61.c



					./integration/qmomo.c



					./integration/qmomof.c



					./integration/qng.c



					./integration/workspace.c



					./interpolation/accel.c



					./interpolation/akima.c



					./interpolation/cspline.c



					./interpolation/interp.c



					./interpolation/linear.c




					./interpolation/poly.c



					./interpolation/spline.c


					./linalg/balance.c



					./linalg/balancemat.c



					./linalg/bidiag.c



					./linalg/cholesky.c



					./linalg/choleskyc.c



					./linalg/exponential.c



					./linalg/hermtd.c



					./linalg/hessenberg.c



					./linalg/hesstri.c



					./linalg/hh.c



					./linalg/householder.c



					./linalg/householdercomplex.c



					./linalg/lq.c



					./linalg/lu.c



					./linalg/luc.c



					./linalg/multiply.c



					./linalg/ptlq.c



					./linalg/qr.c



					./linalg/qrpt.c



					./linalg/svd.c



					./linalg/symmtd.c



					./linalg/tridiag.c



					./matrix/copy.c


						



					./matrix/file.c


						




					./matrix/getset.c



					./matrix/init.c


						




					./matrix/matrix.c



					./matrix/minmax.c








						







						







						



					./matrix/oper.c








						







						







						



					./matrix/prop.c








						







						







						



					./matrix/rowcol.c



					./matrix/submatrix.c



					./matrix/swap.c








						







						







						



					./matrix/view.c








						







						







						




				monte


					./monte/miser.c



					./monte/plain.c



					./monte/vegas.c




				multifit


					./multifit/convergence.c








						







						







						



					./multifit/covar.c



					./multifit/fdfsolver.c



					./multifit/fsolver.c








						







						







						



					./multifit/gradient.c



					./multifit/lmder.c



					./multifit/multilinear.c



					./multifit/work.c




				multimin


					./multimin/conjugate_fr.c



					./multimin/conjugate_pr.c



					./multimin/convergence.c


						


					./multimin/diff.c








						







						







						



					./multimin/fdfminimizer.c



					./multimin/fminimizer.c



					./multimin/simplex.c



					./multimin/simplex2.c



					./multimin/steepest_descent.c



					./multimin/vector_bfgs.c



					./multimin/vector_bfgs2.c




				multiroots


					./multiroots/broyden.c



					./multiroots/convergence.c


						




					./multiroots/dnewton.c



					./multiroots/fdfsolver.c








						







						







						



					./multiroots/fdjac.c



					./multiroots/fsolver.c



					./multiroots/gnewton.c



					./multiroots/hybrid.c



					./multiroots/hybridj.c



					./multiroots/newton.c




				ntuple


					./ntuple/ntuple.c




				ode-initval


					./ode-initval/bsimp.c



					./ode-initval/control.c



					./ode-initval/cscal.c








						







						







						



					./ode-initval/cstd.c



					./ode-initval/evolve.c



					./ode-initval/gear1.c



					./ode-initval/gear2.c



					./ode-initval/rk2.c



					./ode-initval/rk2imp.c



					./ode-initval/rk2simp.c



					./ode-initval/rk4.c



					./ode-initval/rk4imp.c



					./ode-initval/rk8pd.c



					./ode-initval/rkck.c



					./ode-initval/rkf45.c



					./ode-initval/step.c




				permutation


					./permutation/canonical.c



					./permutation/file.c


						




					./permutation/init.c


						






					./permutation/permutation.c



					./permutation/permute.c




					./poly/dd.c



					./poly/deriv.c


			

					./poly/eval.c








						







						







						



					./poly/solve_cubic.c



					./poly/solve_quadratic.c



					./poly/zsolve.c



					./poly/zsolve_cubic.c



					./poly/zsolve_init.c



					./poly/zsolve_quadratic.c




					./qrng/halton.c



					./qrng/niederreiter-2.c



					./qrng/qrng.c



					./qrng/reversehalton.c



					./qrng/sobol.c




					./randist/bernoulli.c



					./randist/beta.c




					./randist/bigauss.c



					./randist/binomial.c



					./randist/binomial_tpe.c



					./randist/cauchy.c








						







						







						



					./randist/chisq.c








						







						







						



					./randist/dirichlet.c



					./randist/discrete.c



					./randist/erlang.c



					./randist/exponential.c


						



					./randist/exppow.c








						







						







						



					./randist/fdist.c








						







						







						



					./randist/flat.c








						







						







						



					./randist/gamma.c








						







						







						



					./randist/gauss.c








						







						







						



					./randist/gausstail.c



					./randist/gausszig.c



					./randist/geometric.c



					./randist/gumbel.c



					./randist/hyperg.c



					./randist/landau.c



					./randist/laplace.c








						







						







						



					./randist/levy.c



					./randist/logarithmic.c



					./randist/logistic.c








						







						







						



					./randist/lognormal.c








						







						







						



					./randist/multinomial.c



					./randist/nbinomial.c



					./randist/pareto.c








						







						







						



					./randist/pascal.c



					./randist/poisson.c



					./randist/rayleigh.c






					./randist/shuffle.c



					./randist/sphere.c



					./randist/tdist.c

					



					./randist/weibull.c







					./rng/borosh13.c



					./rng/cmrg.c



					./rng/coveyou.c



					./rng/default.c



					./rng/file.c


						



					./rng/fishman18.c



					./rng/fishman20.c



					./rng/fishman2x.c



					./rng/gfsr4.c



					./rng/knuthran.c



					./rng/knuthran2.c



					./rng/knuthran2002.c



					./rng/lecuyer21.c



					./rng/minstd.c



					./rng/mrg.c



					./rng/mt.c



					./rng/r250.c



					./rng/ran0.c



					./rng/ran1.c



					./rng/ran2.c



					./rng/ran3.c



					./rng/rand.c



					./rng/rand48.c



					./rng/random.c



					./rng/randu.c



					./rng/ranf.c



					./rng/ranlux.c



					./rng/ranlxd.c



					./rng/ranlxs.c



					./rng/ranmar.c



					./rng/rng.c



					./rng/slatec.c



					./rng/taus.c



					./rng/taus113.c



					./rng/transputer.c



					./rng/tt.c



					./rng/types.c



					./rng/uni.c



					./rng/uni32.c



					./rng/vax.c



					./rng/waterman14.c



					./rng/zuf.c




					./roots/bisection.c



					./roots/brent.c








						







						







						



					./roots/convergence.c


						


						



					./roots/falsepos.c



					./roots/fdfsolver.c


						




					./roots/fsolver.c


						





					./roots/newton.c








						







						







						



					./roots/secant.c



					./roots/steffenson.c




			

					./siman/siman.c




					./sort/sort.c








						







						







						



					./sort/sortind.c



					./sort/sortvec.c



					./sort/sortvecind.c



					./sort/subset.c



					./sort/subsetind.c




			

					./specfunc/airy.c



					./specfunc/airy_der.c



					./specfunc/airy_zero.c



					./specfunc/atanint.c



					./specfunc/bessel.c



					./specfunc/bessel_amp_phase.c



					./specfunc/bessel_i.c



					./specfunc/bessel_I0.c



					./specfunc/bessel_I1.c



					./specfunc/bessel_In.c



					./specfunc/bessel_Inu.c



					./specfunc/bessel_j.c



					./specfunc/bessel_J0.c



					./specfunc/bessel_J1.c



					./specfunc/bessel_Jn.c



					./specfunc/bessel_Jnu.c



					./specfunc/bessel_k.c



					./specfunc/bessel_K0.c



					./specfunc/bessel_K1.c



					./specfunc/bessel_Kn.c



					./specfunc/bessel_Knu.c



					./specfunc/bessel_olver.c



					./specfunc/bessel_sequence.c



					./specfunc/bessel_temme.c



					./specfunc/bessel_y.c



					./specfunc/bessel_Y0.c



					./specfunc/bessel_Y1.c



					./specfunc/bessel_Yn.c



					./specfunc/bessel_Ynu.c



					./specfunc/bessel_zero.c



					./specfunc/beta.c


						



					./specfunc/beta_inc.c



					./specfunc/clausen.c



					./specfunc/coulomb.c



					./specfunc/coulomb_bound.c



					./specfunc/coupling.c



					./specfunc/dawson.c



					./specfunc/debye.c



					./specfunc/dilog.c



					./specfunc/elementary.c



					./specfunc/ellint.c



					./specfunc/elljac.c



					./specfunc/erfc.c



					./specfunc/exp.c



					./specfunc/expint.c



					./specfunc/expint3.c



					./specfunc/fermi_dirac.c



					./specfunc/gamma.c


						




					./specfunc/gamma_inc.c



					./specfunc/gegenbauer.c



					./specfunc/hyperg.c








						







						







						



					./specfunc/hyperg_0F1.c



					./specfunc/hyperg_1F1.c



					./specfunc/hyperg_2F0.c



					./specfunc/hyperg_2F1.c



					./specfunc/hyperg_U.c



					./specfunc/laguerre.c



					./specfunc/lambert.c



					./specfunc/legendre_con.c



					./specfunc/legendre_H3d.c



					./specfunc/legendre_poly.c



					./specfunc/legendre_Qn.c



					./specfunc/log.c



					./specfunc/mathieu_angfunc.c



					./specfunc/mathieu_charv.c



					./specfunc/mathieu_coeff.c



					./specfunc/mathieu_radfunc.c



					./specfunc/mathieu_workspace.c



					./specfunc/poch.c



					./specfunc/pow_int.c



					./specfunc/psi.c



					./specfunc/result.c



					./specfunc/shint.c



					./specfunc/sinint.c



					./specfunc/synchrotron.c



					./specfunc/transport.c



					./specfunc/trig.c



					./specfunc/zeta.c




					./statistics/absdev.c



					./statistics/covariance.c



					./statistics/kurtosis.c



					./statistics/lag1.c



					./statistics/mean.c



					./statistics/median.c



					./statistics/minmax.c


						



					./statistics/p_variance.c



					./statistics/quantiles.c



					./statistics/skew.c



					./statistics/ttest.c



					./statistics/variance.c



					./statistics/wabsdev.c



					./statistics/wkurtosis.c



					./statistics/wmean.c



					./statistics/wskew.c



					./statistics/wvariance.c




					./sum/levin_u.c



					./sum/levin_utrunc.c



					./sum/work_u.c



					./sum/work_utrunc.c





					./sys/coerce.c



					./sys/expm1.c



					./sys/fcmp.c



					./sys/fdiv.c



					./sys/hypot.c








						







						







						



					./sys/infnan.c



					./sys/invhyp.c



					./sys/ldfrexp.c



					./sys/log1p.c



					./sys/minmax.c


						



					./sys/pow_int.c








						







						







						



					./sys/prec.c




			

					./test/results.c




					./vector/copy.c



					./vector/file.c



					./vector/init.c



					./vector/minmax.c



					./vector/oper.c



					./vector/prop.c



					./vector/reim.c



					./vector/subvector.c



					./vector/swap.c



					./vector/vector.c



					./vector/view.c





					./wavelet/bspline.c



					./wavelet/daubechies.c



					./wavelet/dwt.c



					./wavelet/haar.c



					./wavelet/wavelet.c




					./bspline/bspline.c
)
