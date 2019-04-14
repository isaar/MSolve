using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Text;

namespace ISAAR.MSolve.Optimization.Algorithms.GradientBased.CCSA
{
    /// <summary>
    /// Originally described at https://github.com/stevengj/nlopt/blob/master/src/api/nlopt.h#L158
    /// </summary>
    internal enum nlopt_result
    {
        /// <summary>generic failure code</summary>
        NLOPT_FAILURE = -1,
        NLOPT_INVALID_ARGS = -2,
        NLOPT_OUT_OF_MEMORY = -3,
        NLOPT_ROUNDOFF_LIMITED = -4,
        NLOPT_FORCED_STOP = -5,


        /// <summary>generic success code</summary>
        NLOPT_SUCCESS = 1,
        NLOPT_STOPVAL_REACHED = 2,
        NLOPT_FTOL_REACHED = 3,
        NLOPT_XTOL_REACHED = 4,
        NLOPT_MAXEVAL_REACHED = 5,
        NLOPT_MAXTIME_REACHED = 6,

        NLOPT_MINF_MAX_REACHED = NLOPT_STOPVAL_REACHED
    }

    /// <summary>
    /// Originally described at https://github.com/stevengj/nlopt/blob/master/src/api/nlopt.h#L60
    /// </summary>
    /// <param name = "x" >Will not be modified.</param>
    /// <param name="gradient">Null if not needed.</param>
    /// <param name="func_data">This was originally a void pointer so cast away and pray.</param>
    internal delegate double nlopt_func(int n, double[] x, int x_offset, double[] gradient, int gradient_offset, 
        object[] func_data, int func_data_offset);

    /// <summary>
    /// ginally described at https://github.com/stevengj/nlopt/blob/master/src/api/nlopt.h#L64
    /// </summary>
    /// <param name="x" >Will not be modified.</param>
    /// <param name="gradient">Null if not needed.</param>
    /// <param name="func_data">This was originally a void pointer so cast away and pray.</param>
    internal delegate void nlopt_mfunc(int m, double[] result, int result_offset, int n, double[] x, int x_offset, 
        double[] gradient, int gradient_offset, object[] func_data, int func_data_offset);

    /// <summary>
    /// Originally described at https://github.com/stevengj/nlopt/blob/master/src/api/nlopt.h#L276
    /// </summary>
    internal delegate object nlopt_munge(object[] p, int p_offset);

    /// <summary>
    /// A preconditioner, which preconditions v at x to return vpre. (The meaning of "preconditioning" is algorithm-dependent.)
    /// Originally described at https://github.com/stevengj/nlopt/blob/master/src/api/nlopt.h#L70
    /// </summary>
    /// <param name="x">Will not be modified.</param>
    /// <param name="v">Will not be modified.</param>
    /// <param name="data">This was originally a void pointer so cast away and pray.</param>
    internal delegate void nlopt_precond(int n, double[] x, int x_offset, double[] v, int v_offset, double[] vpre, 
        int vpre_offset, object[] data, int data_offset);

    /// <summary>
    /// Originally described at https://github.com/stevengj/nlopt/blob/master/src/util/nlopt-util.h#L114
    /// </summary>
    internal class nlopt_constraint
    {
        /// <summary>dimension of constraint: mf maps R^n -> R^m</summary>
        internal int m;

        /// <summary>one-dimensional constraint, requires m == 1</summary>
        internal nlopt_func f;
        internal nlopt_mfunc mf;

        /// <summary>preconditioner for f (NULL if none or if mf)</summary>
        internal nlopt_precond pre;  
        internal object[] f_data;
        internal int f_data_offset;
        internal double[] tol;
    }

    /// <summary>
    /// Originally described at https://github.com/stevengj/nlopt/blob/master/src/api/nlopt-internal.h#L35, 
    /// https://github.com/stevengj/nlopt/blob/master/src/api/nlopt.h#L188
    /// </summary>
    internal class nlopt_opt
    {
        /// <summary>
        /// The optimization algorithm (immutable)
        /// </summary>
        //nlopt_algorithm algorithm;

        /// <summary>the dimension of the problem (immutable)</summary>
        internal int n;

        internal nlopt_func f;

        /// <summary>objective function to minimize</summary>
        internal object[] f_data;

        internal int f_data_offset;

        /// <summary>optional preconditioner for f (NULL if none)</summary>
        internal nlopt_precond pre;

        /// <summary>nonzero if we are maximizing, not minimizing</summary>
        internal bool maximize;

        /// <summary>lower and upper bounds (length n) </summary>
        internal double[] lb, ub;

        /// <summary>number of inequality constraints</summary>
        internal int m;

        /// <summary>number of inequality constraints allocated</summary>
        internal int m_alloc;

        /// <summary>inequality constraints, length m_alloc</summary>
        internal nlopt_constraint[] fc;

        /// <summary>number of equality constraints</summary>
        internal int p;

        /// <summary>number of equality constraints allocated</summary>
        internal int p_alloc;

        /// <summary>equality constraints, length p_alloc</summary>
        internal nlopt_constraint[] h;

        /// <summary>hack for wrappers</summary>
        internal nlopt_munge munge_on_destroy, munge_on_copy;

        #region stopping criteria
        /// <summary>stop when f reaches stopval or better</summary>
        internal double stopval;

        /// <summary>relative/absolute f tolerances</summary>
        internal double ftol_rel, ftol_abs;

        /// <summary>rel x tolerance</summary>
        internal double xtol_rel;

        /// <summary>abs x tolerances</summary>
        internal double[] xtol_abs;

        /// <summary>max # evaluations</summary>
        internal int maxeval;

        /// <summary>number of evaluations</summary>
        internal int numevals;

        /// <summary>max time (seconds)</summary>
        internal double maxtime;

        /// <summary>if nonzero, force a halt the next time we try to evaluate the objective during optimization</summary>
        internal int force_stop;

        /// <summary>
        /// when local optimization is used, we need a force_stop in the parent object to force a stop in child optimizations
        /// </summary>
        internal nlopt_opt force_stop_child;
        #endregion

        #region algorithm-specific parameters 
        /// <summary>local optimizer</summary>
        internal nlopt_opt local_opt;

        /// <summary>population size for stochastic algs </summary>
        internal int stochastic_population;

        /// <summary>initial step sizes (length n) for nonderivative algs</summary>
        internal double[] dx;

        /// <summary>max subspace dimension (0 for default)</summary>
        internal int vector_storage;

        /// <summary>algorithm-specific workspace during optimization</summary>
        internal object work;

        internal int work_offset;

        /// <summary>description of most recent error</summary>
        internal string errmsg; 
        #endregion
    }
}
