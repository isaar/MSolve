using System;
using System.Collections.Generic;
using System.Text;
using static ISAAR.MSolve.Optimization.Algorithms.GradientBased.CCSA.NLoptApi;
using static ISAAR.MSolve.Optimization.Algorithms.GradientBased.CCSA.NLoptOptions;
using static ISAAR.MSolve.Optimization.Algorithms.GradientBased.CCSA.NLoptStopping;


namespace ISAAR.MSolve.Optimization.Algorithms.GradientBased.CCSA
{
    /// <summary>
    /// Originally described at 
    /// </summary>
    internal class f_max_data
    {
        internal nlopt_func f;
        internal nlopt_precond pre;
        internal object[] f_data;
        internal int f_data_offset;
    }

    internal static class NLoptOptimize
    {
        /// <summary>
        /// Originally described at https://github.com/stevengj/nlopt/blob/master/src/api/optimize.c#L813
        /// </summary>
        internal static nlopt_result nlopt_optimize(nlopt_opt opt, double[] x, int x_offset, ref double opt_f)
        {
            nlopt_func f;
            object[] f_data;
            nlopt_precond pre;
            f_max_data fmd = new f_max_data();
            bool maximize;
            nlopt_result ret;

            nlopt_unset_errmsg(opt);
            if (opt == null /*|| opt_f == null*/ || opt.f == null)
            {
                RETURN_ERR(nlopt_result.NLOPT_INVALID_ARGS, opt, "NULL args to nlopt_optimize");
            }
            f = opt.f;
            f_data = opt.f_data;
            pre = opt.pre;

            // for maximizing, just minimize the f_max wrapper, which flips the sign of everything
            if ((maximize = opt.maximize))
            {
                fmd.f = f;
                fmd.f_data = f_data;
                fmd.pre = pre;
                opt.f = f_max;
                opt.f_data = new object[] { fmd };
                opt.f_data_offset = 0;
                if (opt.pre != null) opt.pre = pre_max;
                opt.stopval = -opt.stopval;
                opt.maximize = false;
            }

            // possibly eliminate lb == ub dimensions for some algorithms
            {
                nlopt_opt elim_opt = opt;
                if (elimdim_wrapcheck(opt))
                {
                    elim_opt = elimdim_create(opt);
                    if (elim_opt == null)
                    {
                        nlopt_set_errmsg(opt, "failure allocating elim_opt");
                        ret = nlopt_result.NLOPT_OUT_OF_MEMORY;
                        goto done;
                    }
                    elimdim_shrink(opt.n, x, x_offset, opt.lb, opt.ub);
                }

                ret = nlopt_optimize_(elim_opt, x, x_offset, ref opt_f);

                if (elim_opt != opt)
                {
                    elimdim_destroy(elim_opt);
                    elimdim_expand(opt.n, x, x_offset, opt.lb, opt.ub);
                }
            }

            done:
            if (maximize)  // restore original signs
            {            
                opt.maximize = maximize;
                opt.stopval = -opt.stopval;
                opt.f = f;
                opt.f_data = f_data;
                opt.pre = pre;
                opt_f = -opt_f;
            }

            return ret;
        }

        /// <summary>
        /// Originally described at https://github.com/stevengj/nlopt/blob/master/src/api/optimize.c#L878
        /// </summary>
        internal static nlopt_result nlopt_optimize_limited(nlopt_opt opt, double[] x, int offset_x, ref double minf, 
            int maxeval, double maxtime)
        {
            int save_maxeval;
            double save_maxtime;
            nlopt_result ret;

            nlopt_unset_errmsg(opt);

            if (opt == null) RETURN_ERR(nlopt_result.NLOPT_INVALID_ARGS, opt, "NULL opt arg");

            save_maxeval = nlopt_get_maxeval(opt);
            save_maxtime = nlopt_get_maxtime(opt);

            /* override opt limits if maxeval and/or maxtime are more stringent */
            if (save_maxeval <= 0 || (maxeval > 0 && maxeval < save_maxeval))
                nlopt_set_maxeval(opt, maxeval);
            if (save_maxtime <= 0 || (maxtime > 0 && maxtime < save_maxtime))
                nlopt_set_maxtime(opt, maxtime);

            ret = nlopt_optimize(opt, x, offset_x, ref minf);

            nlopt_set_maxeval(opt, save_maxeval);
            nlopt_set_maxtime(opt, save_maxtime);

            return ret;
        }

        /// <summary>
        /// Given opt, create a new opt with equal-constraint dimensions eliminated.
        /// Originally described at https://github.com/stevengj/nlopt/blob/master/src/api/optimize.c#L253
        /// </summary>
        private static nlopt_opt elimdim_create(nlopt_opt opt)
        {
            throw new NotImplementedException();
        }

        /// <summary>
        /// Like nlopt_destroy, but also frees elimdim_data.
        /// Originally described at https://github.com/stevengj/nlopt/blob/master/src/api/optimize.c#L316
        /// </summary>
        private static void elimdim_destroy(nlopt_opt opt)
        {
            throw new NotImplementedException();
        }

        /// <summary>
        /// Inverse of elimdim_shrink.
        /// Originally described at https://github.com/stevengj/nlopt/blob/master/src/api/optimize.c#L236
        /// </summary>
        private static void elimdim_expand(int n, double[] v, int offset_v, double[] lb, double[] ub)
        {
            throw new NotImplementedException();
        }

        /// <summary>
        /// Modify v to "shrunk" version, with dimensions for lb[i] == ub[i] elim'ed.
        /// Originally described at https://github.com/stevengj/nlopt/blob/master/src/api/optimize.c#L226
        /// </summary>
        private static void elimdim_shrink(int n, double[] v, int offset_v, double[] lb, double[] ub)
        {
            throw new NotImplementedException();
        }

        /// <summary>
        /// Return whether to use elimdim wrapping.
        /// Originally described at https://github.com/stevengj/nlopt/blob/master/src/api/optimize.c#L340
        /// </summary>
        private static bool elimdim_wrapcheck(nlopt_opt opt)
        {
            throw new NotImplementedException();
        }

        /// <summary>
        /// Wrapper for maximizing: just flip the sign of f and grad. 
        /// Originally described at https://github.com/stevengj/nlopt/blob/master/src/api/optimize.c#L792
        /// </summary>
        private static double f_max(int n, double[] x, int x_offset, double[] grad, int grad_offset, object[] data, 
            int data_offset)
        {
            f_max_data d = (f_max_data)data[data_offset];
            double val = d.f(n, x, x_offset, grad, grad_offset, d.f_data, data_offset);
            if (grad != null)
            {
                for (int i = 0; i<n; ++i) grad[i] = -grad[i];
            }
            return -val;
        }

        /// <summary>
        /// Originally described at https://github.com/stevengj/nlopt/blob/master/src/api/optimize.c#L380
        /// </summary>
        private static nlopt_result nlopt_optimize_(nlopt_opt opt, double[] x, int offset_x, ref double opt_f)
        {
            throw new NotImplementedException();
        }

        /// <summary>
        /// Originally described at https://github.com/stevengj/nlopt/blob/master/src/api/optimize.c#L804
        /// </summary>
        private static void pre_max(int n, double[] x, int x_offset, double[] v, int v_offset, double[] vpre, int vpre_offset,
            object[] data, int data_offset)
        {
            f_max_data d = (f_max_data)data[data_offset];
            d.pre(n, x, x_offset, v, v_offset, vpre, vpre_offset, d.f_data, data_offset);
            for (int i = 0; i < n; ++i) vpre[i] = -vpre[i];
        }
    }
}
