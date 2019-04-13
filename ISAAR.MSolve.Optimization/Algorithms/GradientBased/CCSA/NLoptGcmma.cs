using System;
using System.Collections.Generic;
using System.Text;
using static ISAAR.MSolve.Optimization.Algorithms.GradientBased.CCSA.NLoptApi;
using static ISAAR.MSolve.Optimization.Algorithms.GradientBased.CCSA.NLoptStopping;
using static ISAAR.MSolve.Optimization.Algorithms.GradientBased.CCSA.NLoptOptimize;
using static ISAAR.MSolve.Optimization.Algorithms.GradientBased.CCSA.NLoptOptions;

//TODO: remove all offsets. In 1D arrays they are useless. Fof constraints there are row major 2D arrays. Use C# multidim arrays instead.
namespace ISAAR.MSolve.Optimization.Algorithms.GradientBased.CCSA
{
    /// <summary>
    /// C# port of the implementation of the Globally Convergent Method of Moving Asymptotes (GCMMA) from the NLopt library:
    /// https://nlopt.readthedocs.io/en/latest/NLopt_Algorithms/, 
    /// https://github.com/stevengj/nlopt/blob/master/src/algs/mma/mma.c
    /// That implementation is based on the paper http://citeseer.ist.psu.edu/viewdoc/summary?doi=10.1.1.146.5196:
    /// "A class of globally convergent optimization methods based on conservative convex separable approximations, Svanberg K."
    /// The code in NLopt is under MIT license and was not adapted from professor Svanberg's original Matlab and FORTRAN 
    /// implementations.
    /// Authors: Serafeim Bakalakos
    /// Author of original C code: Steven G. Johnson
    /// </summary>
    internal class NLoptGcmma
    {
        public static int mma_verbose;

        /// <summary>
        /// Magic minimum value for rho in MMA ... the 2002 paper says it should be a "fixed, strictly positive `small' number, 
        /// e.g. 1e-5"... Grrr, I hate these magic numbers, which seem like they should depend on the objective function in some 
        /// way... In particular, note that rho is dimensionful (= dimensions of objective function) 
        /// </summary>
        private const double MMA_RHOMIN = 1e-5;

        /// <summary>
        /// Note that we implement a hidden feature not in the standard nlopt_minimize_constrained interface: whenever the 
        /// constraint function returns NaN, that constraint becomes inactive.
        /// Originally described at
        /// </summary>
        /// <param name="n"></param>
        /// <param name="f"></param>
        /// <param name="f_data"></param>
        /// <param name="m"></param>
        /// <param name="fc"></param>
        /// <param name="lb">Lower bounds. Will not be modified.</param>
        /// <param name="ub">Upper bounds. Will not be modified.</param>
        /// <param name="x">On input: initial guess. On output: minimizer</param>
        /// <param name="minf">
        /// The minimum value of the objective function: <paramref name="minf"/> = <paramref name="f"/>(<paramref name="x"/>).
        /// </param>
        /// <param name="stop"></param>
        /// <param name="dual_opt"></param>
        internal static nlopt_result mma_minimize(int n, nlopt_func f, object[] f_data, int m, nlopt_constraint[] fc, 
            double[] lb, double[] ub, double[] x, ref double minf, nlopt_stopping stop, nlopt_opt dual_opt)
        {
            // Translated from https://github.com/stevengj/nlopt/blob/master/src/algs/mma/mma.c#L145
            nlopt_result ret = nlopt_result.NLOPT_SUCCESS;
            double rho, fcur;
            double[] xcur, sigma, dfdx, dfdx_cur, xprev, xprevprev; 
            double[] dfcdx, dfcdx_cur;
            double[] fcval, fcval_cur, rhoc, gcval, y, dual_lb, dual_ub;
            int i, ifc, j, k = 0;
            dual_data dd;
            bool feasible;
            double infeasibility;
            int mfc;

            m = nlopt_count_constraints(mfc = m, fc);
            if (nlopt_get_dimension(dual_opt) != m)
            {
                nlopt_stop_msg(stop, 
                    string.Format("dual optimizer has wrong dimension %d != %d", nlopt_get_dimension(dual_opt), m));
                return nlopt_result.NLOPT_INVALID_ARGS;
            }
            try
            {
                // The original source malloc'd a large array sigma and set the rest as offsets. 
                sigma = new double[n];
                dfdx = new double[n];
                dfdx_cur = new double[n];
                xcur = new double[n];
                xprev = new double[n];
                xprevprev = new double[n];
                fcval = new double[m];
                fcval_cur = new double[m];
                rhoc = new double[m];
                gcval = new double[m];
                dual_lb = new double[m];
                dual_ub = new double[m];
                y = new double[m];
                dfcdx = new double[m * n];
                dfcdx_cur = new double[m * n];
            }
            catch (OutOfMemoryException)
            {
                return nlopt_result.NLOPT_OUT_OF_MEMORY;
            }

            dd = new dual_data();
            dd.n = n;
            dd.x = x;
            dd.lb = lb;
            dd.ub = ub;
            dd.sigma = sigma;
            dd.dfdx = dfdx;
            dd.dfcdx = dfcdx;
            dd.fcval = fcval;
            dd.rhoc = rhoc;
            dd.xcur = xcur;
            dd.gcval = gcval;

            for (j = 0; j < n; ++j)
            {
                if (nlopt_isinf(ub[j]) || nlopt_isinf(lb[j]))
                    sigma[j] = 1.0; // arbitrary default
                else
                    sigma[j] = 0.5 * (ub[j] - lb[j]);
            }
            rho = 1.0;
            for (i = 0; i < m; ++i)
            {
                rhoc[i] = 1.0;
                dual_lb[i] = y[i] = 0.0;
                dual_ub[i] = double.MaxValue; // Originally dual_ub[i] = HUGE_VAL;
            }

            dd.fval = fcur = minf = f(n, x, 0, dfdx, 0, f_data, 0);
            ++stop.nevals_p;
            Array.Copy(x, xcur, n);
            if (nlopt_stop_forced(stop))
            {
                ret = nlopt_result.NLOPT_FORCED_STOP;
                return ret;
            }

            feasible = true; infeasibility = 0;
            for (i = ifc = 0; ifc < mfc; ++ifc)
            {
                nlopt_eval_constraint(fcval, i, dfcdx, i * n, fc, ifc, n, x, 0);
                i += fc[ifc].m;
                if (nlopt_stop_forced(stop))
                {
                    ret = nlopt_result.NLOPT_FORCED_STOP;
                    return ret;
                }
            }
            for (i = 0; i < m; ++i)
            {
                feasible = feasible && (fcval[i] <= 0 || nlopt_isnan(fcval[i]));
                if (fcval[i] > infeasibility) infeasibility = fcval[i];
            }

            // For non-feasible initial points, set a finite (large) upper-bound on the dual variables.  What this means is that,
            // if no feasible solution is found from the dual problem, it will minimize the dual objective with the unfeasible
            // constraint weighted by 1e40 -- basically, minimizing the unfeasible constraint until it becomes feasible or until 
            // we at least obtain a step towards a feasible point. Svanberg suggested a different approach in his 1987 paper, 
            // basically introducing additional penalty variables for unfeasible constraints, but this is easier to implement 
            // and at least as efficient.

            if (!feasible)
            {
                for (i = 0; i < m; ++i) dual_ub[i] = 1e40;

            }

            nlopt_set_min_objective(dual_opt, dual_func, new object[] { dd }, 0);
            nlopt_set_lower_bounds(dual_opt, dual_lb);
            nlopt_set_upper_bounds(dual_opt, dual_ub);
            nlopt_set_stopval(dual_opt, double.MinValue);
            nlopt_remove_inequality_constraints(dual_opt);
            nlopt_remove_equality_constraints(dual_opt);

            while (true) //outer iterations
            { 
                double fprev = fcur;
                if (nlopt_stop_forced(stop)) ret = nlopt_result.NLOPT_FORCED_STOP;
                else if (nlopt_stop_evals(stop)) ret = nlopt_result.NLOPT_MAXEVAL_REACHED;
                else if (nlopt_stop_time(stop)) ret = nlopt_result.NLOPT_MAXTIME_REACHED;
                else if (feasible && (minf < stop.minf_max)) ret = nlopt_result.NLOPT_MINF_MAX_REACHED;
                if (ret != nlopt_result.NLOPT_SUCCESS) return ret;
                if (++k > 1) Array.Copy(xprev, xprevprev, n);
                Array.Copy(xcur, xprev, n);

                while (true) //inner iterations
                {
                    double min_dual = double.NaN; 
                    double infeasibility_cur;
                    bool feasible_cur, inner_done;
                    int save_verbose;
                    bool new_infeasible_constraint;
                    nlopt_result reti;

                    // Solve dual problem
                    dd.rho = rho; dd.count = 0;
                    save_verbose = mma_verbose;
                    mma_verbose = 0; // no recursive verbosity
                    reti = nlopt_optimize_limited(dual_opt, y, 0, ref min_dual, 0,
                                  stop.maxtime - (stop.watch.Elapsed.TotalSeconds - stop.start));
                    mma_verbose = save_verbose;
                    if (reti < 0 || reti == nlopt_result.NLOPT_MAXTIME_REACHED)
                    {
                        ret = reti;
                        return ret;
                    }

                    dual_func(m, y, 0, null, 0, new object[] { dd }, 0); // evaluate final xcur etc. 
                    if (mma_verbose == 1)
                    {
                        Console.WriteLine("MMA dual converged in %d iterations to g=%g:", dd.count, dd.gval);
                        for (i = 0; i < Math.Min(mma_verbose, m); ++i)
                        {
                            Console.WriteLine("    MMA y[%u]=%g, gc[%u]=%g\n", i, y[i], i, dd.gcval[i]);
                        }
                    }

                    fcur = f(n, xcur, 0, dfdx_cur, 0, f_data, 0);
                    ++stop.nevals_p;
                    if (nlopt_stop_forced(stop))
                    {
                        ret = nlopt_result.NLOPT_FORCED_STOP;
                        return ret;
                    }
                    feasible_cur = true; infeasibility_cur = 0;
                    new_infeasible_constraint = false;
                    inner_done = dd.gval >= fcur;
                    for (i = ifc = 0; ifc < mfc; ++ifc)
                    {
                        nlopt_eval_constraint(fcval_cur, i, dfcdx_cur, i * n, fc, ifc, n, xcur, 0);
                        i += fc[ifc].m;
                        if (nlopt_stop_forced(stop))
                        {
                            ret = nlopt_result.NLOPT_FORCED_STOP;
                            return ret;
                        }
                    }
                    for (i = ifc = 0; ifc < mfc; ++ifc)
                    {
                        int i0 = i, inext = i + fc[ifc].m;
                        for (; i < inext; ++i)
                        {
                            if (!nlopt_isnan(fcval_cur[i]))
                            {
                                feasible_cur = feasible_cur && (fcval_cur[i] <= fc[ifc].tol[i - i0]);
                                if (!nlopt_isnan(fcval[i])) inner_done = inner_done && (dd.gcval[i] >= fcval_cur[i]);
                                else if (fcval_cur[i] > 0) new_infeasible_constraint = true;
                                if (fcval_cur[i] > infeasibility_cur) infeasibility_cur = fcval_cur[i];
                            }
                        }
                    }

                    if ((fcur < minf && (inner_done || feasible_cur || !feasible))
                        || (!feasible && infeasibility_cur < infeasibility))
                    {
                        if (mma_verbose == 1 && !feasible_cur) Console.WriteLine("MMA - using infeasible point?");
                        dd.fval = minf = fcur;
                        infeasibility = infeasibility_cur;
                        Array.Copy(fcval_cur, fcval, m);
                        Array.Copy(xcur, x, n);
                        Array.Copy(dfdx_cur, dfdx, n);
                        Array.Copy(dfcdx_cur, dfcdx, n * m);

                        // once we have reached a feasible solution, the algorithm should never make the solution infeasible
                        // again (if inner_done), although the constraints may be violated slightly by rounding errors etc. so we
                        // must be a little careful about checking feasibility
                        if (infeasibility_cur == 0)
                        {
                            if (!feasible) // reset upper bounds to infin. 
                            {
                                for (i = 0; i < m; ++i) dual_ub[i] = double.MaxValue;
                                nlopt_set_upper_bounds(dual_opt, dual_ub);
                            }
                            feasible = true;
                        }
                        else if (new_infeasible_constraint) feasible = false;
                    }

                    if (nlopt_stop_forced(stop)) ret = nlopt_result.NLOPT_FORCED_STOP;
                    else if (nlopt_stop_evals(stop)) ret = nlopt_result.NLOPT_MAXEVAL_REACHED;
                    else if (nlopt_stop_time(stop)) ret = nlopt_result.NLOPT_MAXTIME_REACHED;
                    else if (feasible && minf < stop.minf_max) ret = nlopt_result.NLOPT_MINF_MAX_REACHED;
                    if (ret != nlopt_result.NLOPT_SUCCESS) return ret;

                    if (inner_done) break;

                    if (fcur > dd.gval) rho = Math.Min(10 * rho, 1.1 * (rho + (fcur - dd.gval) / dd.wval));
                    for (i = 0; i < m; ++i)
                    {
                        if (!nlopt_isnan(fcval_cur[i]) && fcval_cur[i] > dd.gcval[i])
                        {
                            rhoc[i] = Math.Min(10 * rhoc[i],
                                 1.1 * (rhoc[i] + (fcval_cur[i] - dd.gcval[i]) / dd.wval));
                        }
                    }

                    if (mma_verbose == 1) Console.WriteLine("MMA inner iteration: rho -> %g", rho);
                    for (i = 0; i < Math.Min(mma_verbose, m); ++i)
                    {
                        Console.WriteLine("                 MMA rhoc[%u] -> %g", i, rhoc[i]);
                    } 
                }

                if (nlopt_stop_ftol(stop, fcur, fprev)) ret = nlopt_result.NLOPT_FTOL_REACHED;
                if (nlopt_stop_x(stop, xcur, xprev)) ret = nlopt_result.NLOPT_XTOL_REACHED;
                if (ret != nlopt_result.NLOPT_SUCCESS) return ret;

                /* update rho and sigma for iteration k+1 */
                rho = Math.Max(0.1 * rho, MMA_RHOMIN);
                if (mma_verbose == 1) Console.WriteLine("MMA outer iteration: rho -> %g", rho);
                for (i = 0; i < m; ++i) rhoc[i] = Math.Max(0.1 * rhoc[i], MMA_RHOMIN);
                for (i = 0; i < Math.Min(mma_verbose, m); ++i)
                {
                    Console.WriteLine("                 MMA rhoc[%u] -> %g", i, rhoc[i]);
                }
                    
                if (k > 1)
                {
                    for (j = 0; j < n; ++j)
                    {
                        double dx2 = (xcur[j] - xprev[j]) * (xprev[j] - xprevprev[j]);
                        double gam = dx2 < 0 ? 0.7 : (dx2 > 0 ? 1.2 : 1);
                        sigma[j] *= gam;
                        if (!nlopt_isinf(ub[j]) && !nlopt_isinf(lb[j]))
                        {
                            sigma[j] = Math.Min(sigma[j], 10 * (ub[j] - lb[j]));
                            sigma[j] = Math.Max(sigma[j], 0.01 * (ub[j] - lb[j]));
                        }
                    }
                    for (j = 0; j < Math.Min(mma_verbose, n); ++j)
                    {
                        Console.WriteLine("                 MMA sigma[%u] -> %g\n", j, sigma[j]);
                    }
                }
            }

        }

        /// <summary>
        /// Function for MMA's dual solution of the approximate problem.
        /// Originally described at
        /// </summary>
        /// <param name="y">Will not be modified.</param>
        private static double dual_func(int m, double[] y, int offset_y, double[] grad, int grad_offset, 
            object[] d, int offset_d)
        {
            var dd = (dual_data)d[offset_d];
            // Unpack dual_data DTO. TODO: copy comments for these
            int n = dd.n;
            double[] x = dd.x, lb = dd.lb, ub = dd.ub, sigma = dd.sigma, dfdx = dd.dfdx;
            double[] dfcdx = dd.dfcdx;
            double rho = dd.rho, fval = dd.fval;
            double[] rhoc = dd.rhoc, fcval = dd.fcval;
            double[] xcur = dd.xcur;
            double[] gcval = dd.gcval;

            double val;

            dd.count++;

            val = dd.gval = fval;
            dd.wval = 0;
            for (int i = 0; i < m; ++i) val += y[offset_y + i] * (gcval[i] = nlopt_isnan(fcval[i]) ? 0 : fcval[i]);

            for (int j = 0; j < n; ++j)
            {
                double u, v, dx, denominv, c, sigma2, dx2;

                // First, compute xcur[j] for y.  Because this objective is separable, we can minimize over x analytically, 
                // and the minimum dx is given by the solution of a quadratic equation: u dx^2 + 2 v sigma^2 dx + u sigma^2 = 0
                // where u and v are defined by the sums below. Because of the definitions, it is guaranteed that 
                // |u/v| <= sigma, and it follows that the only dx solution with |dx| <= sigma is given by:
                // (v/u) sigma^2 (-1 + sqrt(1 - (u / v sigma)^2)) = (u/v) / (-1 - sqrt(1 - (u / v sigma)^2))
                // (which goes to zero as u -> 0).  The latter expression is less susceptible to roundoff error.

                if (sigma[j] == 0) // special case for lb[i] == ub[i] dims, dx=0 
                { 
                    xcur[j] = x[j];
                    continue;
                }

                u = dfdx[j];
                v = Math.Abs(dfdx[j]) * sigma[j] + 0.5 * rho;
                for (int i = 0; i < m; ++i) 
                {
                    if (!nlopt_isnan(fcval[i]))
                    {
                        u += dfcdx[i * n + j] * y[offset_y + i];
                        v += (Math.Abs(dfcdx[i * n + j]) * sigma[j] + 0.5 * rhoc[i]) * y[offset_y + i];
                    }
                }
                u *= (sigma2 = sigma[j] * sigma[j]);
                dx = (u / v) / (-1 - Math.Sqrt(Math.Abs(1 - Math.Pow(u / (v * sigma[j]), 2))));
                xcur[j] = x[j] + dx;
                if (xcur[j] > ub[j]) xcur[j] = ub[j];
                else if (xcur[j] < lb[j]) xcur[j] = lb[j];
                if (xcur[j] > x[j] + 0.9 * sigma[j]) xcur[j] = x[j] + 0.9 * sigma[j];
                else if (xcur[j] < x[j] - 0.9 * sigma[j]) xcur[j] = x[j] - 0.9 * sigma[j];
                dx = xcur[j] - x[j];

                // Function value:
                dx2 = dx * dx;
                denominv = 1.0 / (sigma2 - dx2);
                val += (u * dx + v * dx2) * denominv;

                // Update gval, wval, gcval (approximant functions)
                c = sigma2 * dx;
                dd.gval += (dfdx[j] * c + (Math.Abs(dfdx[j]) * sigma[j] + 0.5 * rho) * dx2) * denominv;
                dd.wval += 0.5 * dx2 * denominv;
                for (int i = 0; i < m; ++i)
                {
                    if (!nlopt_isnan(fcval[i]))
                    {
                        gcval[i] += (dfcdx[i * n + j] * c + (Math.Abs(dfcdx[i * n + j]) * sigma[j] 
                            + 0.5 * rhoc[i]) * dx2) * denominv;
                    }
                }
                
            }

            // Gradient is easy to compute: since we are at a minimum x (dval/dx=0), we only need the partial derivative 
            // with respect to y, and we negate because we are maximizing: 
            if (grad != null)
            {
                for (int i = 0; i < m; ++i) grad[grad_offset + i] = -gcval[i];
            }
            return -val;
        }

        /// <summary>
        /// Originally described at
        /// </summary>
        private class dual_data //TODO: rename this
        {
            /// <summary>evaluation count, incremented each call</summary>
            internal int count; 

            /// <summary>must be set on input to dimension of x</summary>
            internal int n;

            /// <summary>arrays of length n</summary>
            internal double[] x, lb, ub, sigma, dfdx;  

            /// <summary>m-by-n array (row major) of fc gradients</summary>
            internal double[] dfcdx;  

            /// <summary>must be set on input</summary>
            internal double fval, rho;  

            /// <summary>arrays of length m </summary>
            internal double[] fcval, rhoc; 

            /// <summary>array of length n, output each time</summary>
            internal double[] xcur;  

            /// <summary>output each time</summary>
            internal double gval, wval; 

            /// <summary>output each time (array length m)</summary>
            internal double[] gcval;  
        }
    }
}
