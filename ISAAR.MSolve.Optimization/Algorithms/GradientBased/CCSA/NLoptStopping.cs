using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Text;

namespace ISAAR.MSolve.Optimization.Algorithms.GradientBased.CCSA
{
    /// <summary>
    /// Stopping criteria. Originally described at https://github.com/stevengj/nlopt/blob/master/src/util/nlopt-util.h#L75.
    /// </summary>
    internal class nlopt_stopping
    {
        internal int n;
        internal double minf_max;
        internal double ftol_rel;
        internal double ftol_abs;
        internal double xtol_rel;
        internal readonly double[] xtol_abs;
        internal int nevals_p;
        internal int maxeval;
        internal Stopwatch watch; // addition for C#
        internal double maxtime, start;
        internal bool force_stop;

        /// <summary>pointer to msg string to update</summary>
        internal string stop_msg;
    }

    internal static class NLoptStopping
    {
        /// <summary>
        /// Originally described at https://github.com/stevengj/nlopt/blob/master/src/util/stop.c#L108
        /// </summary>
        /// <param name="p"></param>
        /// <param name="c">Will not be modified</param>
        internal static int nlopt_count_constraints(int p, nlopt_constraint[] c)
        {
            int count = 0;
            for (int i = 0; i < p; ++i) count += c[i].m;
            return count;
        }

        /// <summary>
        /// Originally described at https://github.com/stevengj/nlopt/blob/master/src/util/stop.c#L125.
        /// </summary>
        /// <param name="result"></param>
        /// <param name="grad"></param>
        /// <param name="c">Will not be modified</param>
        /// <param name="n"></param>
        /// <param name="x">Will not be modified</param>
        internal static void nlopt_eval_constraint(double[] result, int result_offset, double[] grad, int grad_offset,
            nlopt_constraint[] c, int c_offset, int n, double[] x, int x_offset)
        {
            //offsets should be propagated to methods called here as well.Better to use multidimensional arrays
            nlopt_constraint constraint = c[c_offset];
            if (constraint.f != null) result[result_offset + 0] = constraint.f(n, x, x_offset, grad, grad_offset,
                constraint.f_data, constraint.f_data_offset);
            else constraint.mf(constraint.m, result, result_offset, n, x, x_offset, grad, grad_offset, 
                constraint.f_data, constraint.f_data_offset);
        }

        /// <summary>
        /// Originally described at https://github.com/stevengj/nlopt/blob/master/src/util/stop.c#L166.
        /// Different from the NLopt implementation.
        /// </summary>
        internal static bool nlopt_isinf(double x) => double.IsInfinity(x);

        /// <summary>
        /// Originally described at https://github.com/stevengj/nlopt/blob/master/src/util/stop.c#L204.
        /// Different from the NLopt implementation.
        /// </summary>
        internal static bool nlopt_isnan(double x) => double.IsNaN(x);

        /// <summary>
        /// Originally described at https://github.com/stevengj/nlopt/blob/master/src/util/stop.c#L188.
        /// Somewhat different from the NLopt implementation.
        /// </summary>
        internal static bool nlopt_istiny(double x)
        {
            if (x == 0.0) return true;
            else return Math.Abs(x) < 2.2250738585072014e-308; // assume IEEE 754 double
        }

        /// <summary>
        /// Originally described at https://github.com/stevengj/nlopt/blob/master/src/util/stop.c#L83.
        /// </summary>
        /// <param name="s">Will not be modified</param>
        internal static bool nlopt_stop_evals(nlopt_stopping s)
        {
            return (s.maxeval > 0 && (s.nevals_p >= s.maxeval));
        }

        /// <summary>
        /// Originally described at https://github.com/stevengj/nlopt/blob/master/src/util/stop.c#L39.
        /// </summary>
        /// <param name="s">Will not be modified.</param>
        /// <param name="f"></param>
        /// <param name="oldf"></param>
        internal static bool nlopt_stop_ftol(nlopt_stopping s, double f, double oldf)
        {
            return (relstop(oldf, f, s.ftol_rel, s.ftol_abs));
        }

        /// <summary>
        /// Originally described at https://github.com/stevengj/nlopt/blob/master/src/util/stop.c#L103.
        /// </summary>
        /// <param name="stop">Will not be modified</param>
        internal static bool nlopt_stop_forced(nlopt_stopping stop)
        {
            // This method has been heavily modified, since IO in C# is too different than C.
            //TODO: Is the original intent preserved?
            return stop.force_stop;
        }

        /// <summary>
        /// Originally described at https://github.com/stevengj/nlopt/blob/master/src/util/stop.c#L154
        /// </summary>
        internal static void nlopt_stop_msg(nlopt_stopping s, string msg)
        {
            // This method has been heavily modified, since IO in C# is too different than C.
            //TODO: Is the original intent preserved?
            Console.WriteLine(msg);
            s.stop_msg = msg;
        }

        /// <summary>
        /// Originally described at https://github.com/stevengj/nlopt/blob/master/src/util/stop.c#L93.
        /// </summary>
        /// <param name="s">Will not be modified</param>
        internal static bool nlopt_stop_time(nlopt_stopping s)
        {
            //s.watch.Stop(); // Do not stop this. Other methods will query the stopwatch.Elapsed
            return (s.maxtime > 0 && s.watch.Elapsed.TotalSeconds - s.start >= s.maxtime);
        }

        /// <summary>
        /// Originally described at https://github.com/stevengj/nlopt/blob/master/src/util/stop.c#L49.
        /// </summary>
        /// <param name="s">Will not be modified</param>
        /// <param name="x">Will not be modified</param>
        /// <param name="oldx">Will not be modified</param>
        internal static bool nlopt_stop_x(nlopt_stopping s, double[] x, double[] oldx)
        {
            for (int i = 0; i < s.n; ++i)
                if (!relstop(oldx[i], x[i], s.xtol_rel, s.xtol_abs[i])) return false;
            return true;
        }

        /// <summary>
        /// Originally described at https://github.com/stevengj/nlopt/blob/master/src/util/stop.c#L32.
        /// It is private because it was a static C function
        /// </summary>
        private static bool relstop(double vold, double vnew, double reltol, double abstol) 
        {
            if (nlopt_isinf(vold)) return false;
            return (Math.Abs(vnew - vold) < abstol
                || Math.Abs(vnew - vold) < reltol * (Math.Abs(vnew) + Math.Abs(vold)) * 0.5
                || (reltol > 0 && vnew == vold));        // catch vnew == vold == 0
        }
    }
}
