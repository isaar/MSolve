using System;
using System.Collections.Generic;
using System.Text;
using static ISAAR.MSolve.Optimization.Algorithms.GradientBased.CCSA.NLoptOptions;

namespace ISAAR.MSolve.Optimization.Algorithms.GradientBased.CCSA
{
    internal static class NLoptApi
    {
        /// <summary>
        /// Originally described at https://github.com/stevengj/nlopt/blob/master/src/api/nlopt.h#L206.
        /// </summary>
        internal static int nlopt_get_dimension(nlopt_opt opt) => opt.n; //TODO: could not find the NLopt implementation.

        /// <summary>
        /// Originally described at https://github.com/stevengj/nlopt/blob/master/src/api/nlopt.h#L246.
        /// </summary>
        internal static int nlopt_get_maxeval(nlopt_opt opt) => opt.maxeval; //TODO: could not find the NLopt implementation.

        /// <summary>
        /// Originally described at https://github.com/stevengj/nlopt/blob/master/src/api/nlopt.h#L251
        /// </summary>
        internal static double nlopt_get_maxtime(nlopt_opt opt) => opt.maxtime; //TODO: could not find the NLopt implementation.

        /// <summary>
        /// Originally described at https://github.com/stevengj/nlopt/blob/master/src/api/nlopt.h#L245.
        /// </summary>
        internal static nlopt_result nlopt_set_maxeval(nlopt_opt opt, int maxeval) //TODO: could not find the NLopt implementation.
        {
            if (opt != null)
            {
                opt.maxeval = maxeval;
                return nlopt_result.NLOPT_SUCCESS;
            }
            return nlopt_result.NLOPT_INVALID_ARGS;
        }

        /// <summary>
        /// Originally described at https://github.com/stevengj/nlopt/blob/master/src/api/nlopt.h#L250/
        /// </summary>
        internal static nlopt_result nlopt_set_maxtime(nlopt_opt opt, double maxtime) //TODO: could not find the NLopt implementation.
        {
            if (opt != null)
            {
                opt.maxtime = maxtime;
                return nlopt_result.NLOPT_SUCCESS;
            }
            return nlopt_result.NLOPT_INVALID_ARGS;
        }

        /// <summary>
        /// Originally described at https://github.com/stevengj/nlopt/blob/master/src/api/nlopt.h#L231.
        /// </summary>
        internal static nlopt_result nlopt_set_stopval(nlopt_opt opt, double stopval) //TODO: could not find the NLopt implementation.
        {
            if (opt != null)
            {
                opt.stopval = stopval;
                return nlopt_result.NLOPT_SUCCESS;
            }
            return nlopt_result.NLOPT_INVALID_ARGS;
        }

        /// <summary>
        /// Originally described at https://github.com/stevengj/nlopt/blob/master/src/api/nlopt-internal.h#L94.
        /// </summary>
        internal static nlopt_result RETURN_ERR(nlopt_result err, nlopt_opt opt, string msg)
        {
            nlopt_set_errmsg(opt, msg);
            return err;
        }
    }
}
