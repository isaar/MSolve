using System;
using System.Collections.Generic;
using System.Text;
using static ISAAR.MSolve.Optimization.Algorithms.GradientBased.CCSA.NLoptStopping;

namespace ISAAR.MSolve.Optimization.Algorithms.GradientBased.CCSA
{
    internal static class NLoptOptions
    {
        /// <summary>
        /// Originally described at https://github.com/stevengj/nlopt/blob/master/src/api/options.c#L792.
        /// </summary>
        internal static string nlopt_set_errmsg(nlopt_opt opt, string msg)
        {
            // This method has been heavily modified, since IO in C# is too different than C.
            //TODO: Is the original intent preserved?
            Console.WriteLine(msg);
            opt.errmsg = msg;
            return opt.errmsg;
        }

        /// <summary>
        /// Originally described at https://github.com/stevengj/nlopt/blob/master/src/api/options.c#L481.
        /// </summary>
        internal static nlopt_result nlopt_remove_equality_constraints(nlopt_opt opt)
        {

            nlopt_unset_errmsg(opt);
            if (opt == null) return nlopt_result.NLOPT_INVALID_ARGS;
            if (opt.munge_on_destroy != null)
            {
                nlopt_munge munge = opt.munge_on_destroy;
                for (int i = 0; i < opt.p; ++i) munge(opt.h[i].f_data, opt.h[i].f_data_offset);
            }
            opt.h = null;
            opt.p = opt.p_alloc = 0;
            return nlopt_result.NLOPT_SUCCESS;
        }

        /// <summary>
        /// Originally described at https://github.com/stevengj/nlopt/blob/master/src/api/options.c#L370.
        /// </summary>
        internal static nlopt_result nlopt_remove_inequality_constraints(nlopt_opt opt)
        {
            nlopt_unset_errmsg(opt);
            if (opt == null) return nlopt_result.NLOPT_INVALID_ARGS;
            if (opt.munge_on_destroy != null)
            {
                nlopt_munge munge = opt.munge_on_destroy;
                for (int i = 0; i < opt.m; ++i) munge(opt.fc[i].f_data, opt.fc[i].f_data_offset);
            }
            opt.fc = null;
            opt.m = opt.m_alloc = 0;
            return nlopt_result.NLOPT_SUCCESS;
        }

        /// <summary>
        /// Originally described at https://github.com/stevengj/nlopt/blob/master/src/api/options.c#L281.
        /// </summary>
        /// <param name="opt"></param>
        /// <param name="lb">Will not be modified</param>
        internal static nlopt_result nlopt_set_lower_bounds(nlopt_opt opt, double[] lb)
        {
            nlopt_unset_errmsg(opt);
            if ((opt != null) && (opt.n == 0 || (lb != null)))
            {
                if (opt.n > 0) Array.Copy(lb, opt.lb, opt.n);

                for (int i = 0; i < opt.n; ++i)
                {
                    if ((opt.lb[i] < opt.ub[i]) && nlopt_istiny(opt.ub[i] - opt.lb[i])) opt.lb[i] = opt.ub[i];
                }
                return nlopt_result.NLOPT_SUCCESS;
            }
            return nlopt_result.NLOPT_INVALID_ARGS;
        }

        /// <summary>
        /// Originally described at https://github.com/stevengj/nlopt/blob/master/src/api/options.c#L252.
        /// </summary>
        internal static nlopt_result nlopt_set_min_objective(nlopt_opt opt, nlopt_func f, object[] f_data, int f_data_offset)
        {
            return nlopt_set_precond_min_objective(opt, f, null, f_data, f_data_offset);
        }

        /// <summary>
        /// Originally described at https://github.com/stevengj/nlopt/blob/master/src/api/options.c#L235.
        /// </summary>
        internal static nlopt_result nlopt_set_precond_min_objective(nlopt_opt opt, nlopt_func f, nlopt_precond pre,
            object[] f_data, int f_data_offset)
        {
            if (opt != null)
            {
                nlopt_unset_errmsg(opt);
                if (opt.munge_on_destroy != null)
                    opt.munge_on_destroy(opt.f_data, opt.f_data_offset);
                opt.f = f;
                opt.f_data = f_data;
                opt.f_data_offset = f_data_offset;
                opt.pre = pre;
                opt.maximize = false;
                if (nlopt_isinf(opt.stopval) && opt.stopval > 0) opt.stopval = double.MinValue; // switch default from max to min
                return nlopt_result.NLOPT_SUCCESS;
            }
            return nlopt_result.NLOPT_INVALID_ARGS;
        }

        /// <summary>
        /// Originally described at https://github.com/stevengj/nlopt/blob/master/src/api/options.c#L321
        /// </summary>
        /// <param name="opt"></param>
        /// <param name="ub">Will not be modified</param>
        internal static nlopt_result nlopt_set_upper_bounds(nlopt_opt opt, double[] ub)
        {
            nlopt_unset_errmsg(opt);
            if ((opt != null) && (opt.n == 0 || (ub != null)))
            {
                if (opt.n > 0) Array.Copy(ub, opt.ub, opt.n);
                for (int i = 0; i < opt.n; ++i)
                    if ((opt.lb[i] < opt.ub[i]) && nlopt_istiny(opt.ub[i] - opt.lb[i])) opt.ub[i] = opt.lb[i];
                return nlopt_result.NLOPT_SUCCESS;
            }
            return nlopt_result.NLOPT_INVALID_ARGS;
        }

        /// <summary>
        /// Originally described at https://github.com/stevengj/nlopt/blob/master/src/api/options.c#L801.
        /// </summary>
        internal static void nlopt_unset_errmsg(nlopt_opt opt)
        {
            if (opt != null) opt.errmsg = null;
        }
    }
}
